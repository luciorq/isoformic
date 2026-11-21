# Workers ----------------------------------------------------------------------

#' The number of workers for `n_tasks` tasks
#'
#' `workers` is a single whole number of at least 1 or `"auto"`, which is
#' `parallelly::availableCores()` (it follows Slurm, the CPU affinity mask and
#' the cgroup CPU quota; without parallelly, 1 with a message). The result is
#' capped at the number of tasks. `option_hint` adds to the error that the
#' default is the `isoformic.workers` option.
#' @keywords internal
#' @noRd
resolve_workers <- function(
  workers,
  n_tasks,
  arg = "workers",
  option_hint = identical(arg, "workers")
) {
  if (identical(workers, "auto")) {
    if (rlang::is_installed("parallelly")) {
      workers <- as.integer(parallelly::availableCores())
    } else {
      cli::cli_inform(
        message = c(
          i = paste(
            "{.code {arg} = \"auto\"} needs the {.pkg parallelly} package",
            "to count the CPUs; running with one worker."
          )
        ),
        class = "isoformic_workers_fallback"
      )
      workers <- 1L
    }
  }
  if (
    !is.numeric(workers) ||
      length(workers) != 1L ||
      is.na(workers) ||
      !is.finite(workers) ||
      workers > .Machine$integer.max ||
      workers < 1 ||
      workers != trunc(workers)
  ) {
    cli::cli_abort(
      message = c(
        x = "{.arg {arg}} must be a whole number of at least 1 or {.val auto}.",
        i = "Use {.val auto} for all the CPUs available to the session.",
        i = if (isTRUE(option_hint)) {
          "The default is the {.code isoformic.workers} option, or 1."
        }
      ),
      class = "isoformic_invalid_workers"
    )
  }
  as.integer(max(1L, min(workers, n_tasks)))
}

#' A function that a worker can run without isoformic
#'
#' Its environment becomes the base environment and its source reference is
#' dropped, so that it serializes as a few bytes and a worker does not load
#' isoformic (which may not be installed, under `devtools::load_all()`). It
#' may only call base R and `pkg::fun()`.
#' @keywords internal
#' @noRd
worker_function <- function(fun) {
  environment(fun) <- baseenv()
  attr(fun, "srcref") <- NULL
  fun
}

#' Run one task with its own seed
#'
#' The session's RNG kind is set together with the seed, so that a task gives
#' the same result in the session and in a worker, whose default kind
#' differs. The session's library paths `lib` come first while the task runs
#' (daemons that the user started may lack them; paths that do not exist on
#' the worker are dropped by `.libPaths()`). Warnings are collected and
#' muffled, because a daemon does not send them back: the result is
#' `list(value = , warnings = )`, and `run_tasks()` raises the warnings again
#' in the session.
#' @keywords internal
#' @noRd
seeded_task <- worker_function(function(x, fun, args, kind, lib) {
  old_lib <- .libPaths()
  on.exit(.libPaths(old_lib), add = TRUE)
  .libPaths(c(lib, old_lib))
  set.seed(
    x[["seed"]],
    kind = kind[[1L]],
    normal.kind = kind[[2L]],
    sample.kind = kind[[3L]]
  )
  seen <- new.env(parent = emptyenv())
  seen$warnings <- character()
  value <- withCallingHandlers(
    do.call(fun, c(x[["task"]], args)),
    warning = function(cnd) {
      seen$warnings <- c(seen$warnings, conditionMessage(cnd))
      invokeRestart("muffleWarning")
    }
  )
  list(value = value, warnings = seen$warnings)
})

#' Raise the warnings that seeded tasks collected, naming the task
#'
#' Each distinct message of a task is raised once, in task order, with the
#' class `isoformic_worker_warning`.
#' @keywords internal
#' @noRd
relay_task_warnings <- function(res, labels) {
  for (i in seq_along(res)) {
    for (msg in unique(res[[i]][["warnings"]])) {
      rlang::warn(
        message = c(
          msg,
          i = paste0("Raised by the task on ", labels[[i]], ".")
        ),
        class = "isoformic_worker_warning"
      )
    }
  }
  invisible(res)
}

#' Run `fun` on every task, serially or on mirai workers
#'
#' Each task is a named list of arguments of `fun`; `args` are the arguments
#' common to all tasks. `fun` must be self-contained (see
#' `worker_function()`).
#'
#' `workers` is the requested number of workers, not capped at the number of
#' tasks: the tasks run on workers when both are above one. With one worker
#' (or one task) the tasks run in the session, in order, as calls of `fun`:
#' they use the session's random number stream as before this helper
#' existed. Otherwise one seed per task is drawn in the session before
#' dispatch and every task sets it, so the result is the same for any number
#' of workers above one, and whether or not the daemons were already
#' running, but not the same as with one worker. `seeded = TRUE` runs the
#' seeded tasks in the session (used by the tests). Warnings of seeded tasks
#' are raised again in the session after all tasks return, once per distinct
#' message and task, naming the task (`labels`).
#'
#' The workers are the mirai daemons of the compute profile `"isoformic"`
#' (see `ensure_workers()`), which stay running after the call. Tasks that
#' have not finished when the call exits (an interrupt) are cancelled. The
#' session's library paths are put first in every task, for daemons the
#' user started. An error in a task is raised again naming the task
#' (`labels`), after the warnings of the tasks that succeeded.
#' @keywords internal
#' @noRd
run_tasks <- function(
  tasks,
  fun,
  args = list(),
  workers = 1L,
  labels = NULL,
  seeded = NULL
) {
  requested <- resolve_workers(workers, Inf)
  used <- min(requested, length(tasks))
  seeded <- seeded %||% (used > 1L)
  if (!seeded) {
    return(lapply(tasks, function(task) do.call(fun, c(task, args))))
  }
  labels <- labels %||% paste("task", seq_along(tasks))
  seeds <- sample.int(.Machine$integer.max, length(tasks))
  items <- Map(
    function(task, seed) list(task = task, seed = seed),
    tasks,
    seeds
  )
  names(items) <- NULL
  task_args <- list(
    fun = worker_function(fun),
    args = args,
    kind = RNGkind(),
    lib = .libPaths()
  )
  if (used <= 1L) {
    res <- lapply(items, function(item) {
      withr::with_preserve_seed(do.call(seeded_task, c(list(item), task_args)))
    })
    relay_task_warnings(res, labels)
    return(lapply(res, function(r) r[["value"]]))
  }
  check_mirai()
  # mirai draws from the session's stream when it starts daemons: the
  # stream moves by the seeds only
  withr::with_preserve_seed(ensure_workers(requested, used))
  res <- withr::with_preserve_seed(collect_tasks(items, task_args))
  failed <- vapply(res, mirai::is_error_value, logical(1L))
  relay_task_warnings(res[!failed], labels[!failed])
  failed <- which(failed)
  if (length(failed)) {
    i <- failed[[1L]]
    message <- if (mirai::is_mirai_error(res[[i]])) {
      attr(res[[i]], "message", exact = TRUE) %||% as.character(res[[i]])
    } else {
      paste("the worker returned error value", as.character(res[[i]]))
    }
    cli::cli_abort(
      message = c(
        x = "A worker failed on {labels[[i]]}.",
        i = "{message}",
        if (length(failed) > 1L) {
          c(i = "{length(failed) - 1L} more task{?s} failed.")
        }
      ),
      class = "isoformic_worker_error"
    )
  }
  lapply(res, function(r) r[["value"]])
}

#' Send the seeded tasks to the daemons and wait for all of them
#'
#' When the wait is interrupted, the tasks still queued or running are
#' cancelled (`mirai::stop_mirai()` also interrupts them on the daemons), so
#' that kept daemons do not go on with work nobody collects and the next
#' call does not queue behind it.
#' @keywords internal
#' @noRd
collect_tasks <- function(items, task_args) {
  m <- mirai::mirai_map(
    items,
    seeded_task,
    .args = task_args,
    .compute = workers_profile
  )
  done <- FALSE
  on.exit(if (!done) mirai::stop_mirai(m), add = TRUE)
  res <- wait_tasks(m)
  done <- TRUE
  res
}

#' Wait for all the tasks of a `mirai_map()` and return their values
#'
#' Fails when tasks are still waiting and no daemon has been connected for
#' `seconds` (all of them died during the call, for example killed for
#' their memory): the profile would otherwise wait for a daemon forever.
#' The caller cancels the remaining tasks, and the next call starts
#' isoformic's daemons again.
#' @keywords internal
#' @noRd
wait_tasks <- function(m, seconds = the_workers$grace) {
  lost_since <- NULL
  checked <- Sys.time()
  while (mirai::unresolved(m)) {
    Sys.sleep(0.01)
    # The dispatcher is asked for the connections once a second
    if (as.numeric(difftime(Sys.time(), checked, units = "secs")) < 1) {
      next
    }
    checked <- Sys.time()
    if (profile_connections() > 0L) {
      lost_since <- NULL
    } else if (is.null(lost_since)) {
      lost_since <- checked
    } else if (
      as.numeric(difftime(checked, lost_since, units = "secs")) > seconds
    ) {
      cli::cli_abort(
        message = c(
          x = "The workers stopped before their tasks were done.",
          i = paste(
            "No daemon of the mirai compute profile",
            "{.val {workers_profile}} has been connected for {seconds}",
            "seconds (were they killed for their memory?)."
          )
        ),
        class = "isoformic_worker_error"
      )
    }
  }
  m[]
}

#' Check that mirai can run the workers
#' @keywords internal
#' @noRd
check_mirai <- function() {
  rlang::check_installed(
    "mirai",
    version = "2.3.0",
    reason = "to run with more than one worker."
  )
}

# Daemons ----------------------------------------------------------------------

#' The mirai compute profile of the workers
#' @keywords internal
#' @noRd
workers_profile <- "isoformic"

#' The daemons that isoformic started
#'
#' `url` is the address that `mirai::status()` reports for them (`NULL` when
#' isoformic has not started any), which tells them apart from daemons the
#' user started on the profile with `mirai::daemons()`; `n` is their number;
#' `auto` is `TRUE` when `run_tasks()` started them (it resizes them) and
#' `FALSE` when `isoformic_workers()` did (they are used as they are).
#' `grace` is how many seconds isoformic waits for daemons to connect, and
#' for a daemon to come back during a call, before it fails (the tests
#' lower it).
#' @keywords internal
#' @noRd
the_workers <- new.env(parent = emptyenv())
the_workers$url <- NULL
the_workers$n <- 0L
the_workers$auto <- FALSE
the_workers$grace <- 30

#' The address of the daemons of the profile, or `NULL` when it has none
#' @keywords internal
#' @noRd
profile_url <- function() {
  if (!mirai::daemons_set(.compute = workers_profile)) {
    return(NULL)
  }
  mirai::status(.compute = workers_profile)[["daemons"]]
}

#' Make the profile's daemons ready for `requested` workers and `used` tasks
#'
#' - Daemons that the user started with `mirai::daemons()` are used as they
#'   are; when none of them is connected the call fails instead of waiting
#'   for one forever (see `check_connected()`).
#' - Daemons started with `isoformic_workers(n)` are used as they are, but
#'   the ones that died (an OOM kill, a crash) are launched again, so that
#'   `n` stay connected.
#' - Otherwise `used` daemons are started for the session (the session's
#'   library paths, one BLAS / OpenMP thread each). Later calls reuse them:
#'   more are launched when a call uses more, and they are started again
#'   with `used` when a call requests fewer than are running, so that a call
#'   never runs more tasks at once than its `workers`.
#'
#' The daemons isoformic counts are the connected ones
#' (`mirai::status()$connections`), not the ones it launched: daemons that
#' died between calls are replaced, and when none is left they are all
#' started again, because the profile would otherwise queue the tasks for
#' no daemon.
#' @keywords internal
#' @noRd
ensure_workers <- function(requested, used) {
  url <- profile_url()
  if (!is.null(url) && !identical(url, the_workers$url)) {
    check_connected()
    return(invisible(FALSE))
  }
  auto <- is.null(url) || the_workers$auto
  # isoformic_workers(n) keeps its count, run_tasks() follows the call
  want <- if (auto) used else the_workers$n
  connected <- if (is.null(url)) 0L else profile_connections()
  if (connected == 0L || (auto && connected > requested)) {
    start_workers(want, auto = auto)
  } else if (connected < want) {
    with_worker_env(
      mirai::launch_local(want - connected, .compute = workers_profile)
    )
    the_workers$n <- as.integer(want)
    wait_connections(want)
    check_connected(seconds = 0)
  } else if (auto) {
    the_workers$n <- connected
  }
  invisible(TRUE)
}

#' The number of daemons connected to the profile
#' @keywords internal
#' @noRd
profile_connections <- function() {
  as.integer(mirai::status(.compute = workers_profile)[["connections"]])
}

#' Fail when no daemon of the profile is connected
#'
#' Waits up to `seconds` for one to connect first (remote daemons, or local
#' ones still starting).
#' @keywords internal
#' @noRd
check_connected <- function(seconds = the_workers$grace) {
  if (wait_connections(1L, seconds = seconds)) {
    return(invisible(TRUE))
  }
  cli::cli_abort(
    message = c(
      x = paste(
        "No daemon of the mirai compute profile {.val {workers_profile}}",
        "is connected."
      ),
      i = paste(
        "isoformic does not restart daemons that it did not start: stop",
        "them with {.code isoformic_workers(0)} (or",
        "{.code mirai::daemons(0, .compute = \"isoformic\")}) and call",
        "again."
      )
    ),
    class = "isoformic_worker_error"
  )
}

#' Wait until `n` daemons are connected to the profile, for at most
#' `seconds`
#'
#' `mirai::launch_local()` and `mirai::daemons()` return before the daemons
#' connect: waiting for them gives every task of the call a daemon, lets
#' the next call count the connected daemons, and a daemon stopped before
#' it connects prints an error.
#' @keywords internal
#' @noRd
wait_connections <- function(n, seconds = the_workers$grace) {
  deadline <- Sys.time() + seconds
  repeat {
    if (profile_connections() >= n) {
      return(invisible(TRUE))
    }
    if (Sys.time() >= deadline) {
      return(invisible(FALSE))
    }
    Sys.sleep(0.05)
  }
}

#' Start `n` daemons on the profile and record them as isoformic's
#'
#' Daemons are fresh `Rscript` processes: they get the session's library
#' paths through `R_LIBS` (they do not inherit `.libPaths()`), and one BLAS,
#' OpenMP and data.table thread each, so that `n` workers use `n` cores.
#' Daemons already on the profile are stopped first. Returns when the `n`
#' daemons are connected (or after `the_workers$grace` seconds), and fails
#' when none connects.
#' @keywords internal
#' @noRd
start_workers <- function(n, auto = FALSE) {
  stop_workers(own_only = FALSE)
  with_worker_env(mirai::daemons(n, .compute = workers_profile))
  the_workers$url <- profile_url()
  the_workers$n <- as.integer(n)
  the_workers$auto <- auto
  wait_connections(n)
  check_connected(seconds = 0)
  invisible(n)
}

#' Evaluate `code` with the environment that daemons start with
#' @keywords internal
#' @noRd
with_worker_env <- function(code) {
  withr::with_envvar(
    c(
      R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep),
      OMP_NUM_THREADS = "1",
      OPENBLAS_NUM_THREADS = "1",
      MKL_NUM_THREADS = "1",
      VECLIB_MAXIMUM_THREADS = "1",
      R_DATATABLE_NUM_THREADS = "1"
    ),
    code
  )
}

#' Stop the daemons of the profile
#'
#' With `own_only = TRUE`, only daemons that isoformic started. Never fails:
#' it also runs when the session ends and when the package is unloaded.
#' Returns whether daemons were stopped.
#' @keywords internal
#' @noRd
stop_workers <- function(own_only = TRUE) {
  stopped <- FALSE
  if (isNamespaceLoaded("mirai")) {
    stopped <- tryCatch(
      {
        url <- profile_url()
        if (!is.null(url) && (!own_only || identical(url, the_workers$url))) {
          mirai::daemons(0L, .compute = workers_profile)
          TRUE
        } else {
          FALSE
        }
      },
      error = function(cnd) FALSE
    )
  }
  if (stopped || !own_only) {
    the_workers$url <- NULL
    the_workers$n <- 0L
    the_workers$auto <- FALSE
  }
  stopped
}

#' Start or stop the enrichment workers
#'
#' @description
#' `isoformic_workers()` starts `n` `mirai` daemons ahead for the enrichment
#' functions ([run_enrichment()] and the other functions of that topic),
#' stopping those of the `"isoformic"` compute profile that were running.
#' Every enrichment call then uses them as they are (replacing those that
#' died), whatever its `workers`, as long as `workers` is above 1. The daemons
#' get the session's library paths and one BLAS and OpenMP thread each, and
#' starting them leaves the session's random number stream where it was.
#' `isoformic_workers(0)` stops the daemons of the `"isoformic"` profile,
#' whoever started them; they are also stopped when the session ends or
#' isoformic is unloaded.
#'
#' The section *Workers* of [run_enrichment()] describes how the calls use
#' the daemons, their memory and the reproducibility of their p-values.
#'
#' @param n Number of daemons to start: a whole number of at least 1,
#'   `"auto"` for `parallelly::availableCores()`, or 0 to stop the daemons.
#'
#' @return The number of daemons now running, invisibly.
#'
#' @seealso [run_enrichment()]; the `isoformic.workers` option in
#'   [isoformic-options].
#'
#' @examples
#' \dontrun{
#' isoformic_workers(4)
#' # ... calls of run_enrichment() with workers = 4 ...
#' isoformic_workers(0)
#' }
#' @export
isoformic_workers <- function(n) {
  if (!identical(n, 0) && !identical(n, 0L)) {
    n <- resolve_workers(n, Inf, arg = "n")
  }
  if (n == 0L) {
    stop_workers(own_only = FALSE)
    return(invisible(0L))
  }
  check_mirai()
  withr::with_preserve_seed(start_workers(n, auto = FALSE))
  invisible(n)
}
