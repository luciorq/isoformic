# Stops the daemons of the "isoformic" profile now and when the test ends,
# so that every test starts without kept daemons
local_no_workers <- function(envir = parent.frame()) {
  stop_workers(own_only = FALSE)
  withr::defer(stop_workers(own_only = FALSE), envir = envir)
}

testthat::test_that("resolve_workers validates, caps and resolves auto", {
  testthat::expect_identical(resolve_workers(1L, 5L), 1L)
  testthat::expect_identical(resolve_workers(4, 5L), 4L)
  testthat::expect_identical(resolve_workers(8L, 3L), 3L)
  testthat::expect_identical(resolve_workers(8L, 0L), 1L)
  bad_values <- list(
    0L,
    -1,
    1.5,
    NA_integer_,
    Inf,
    1e10,
    c(1L, 2L),
    "two",
    TRUE,
    NULL
  )
  for (bad in bad_values) {
    testthat::expect_error(
      resolve_workers(bad, 5L),
      class = "isoformic_invalid_workers"
    )
  }
  testthat::expect_error(
    resolve_workers(Inf, Inf),
    class = "isoformic_invalid_workers"
  )
  testthat::local_mocked_bindings(
    is_installed = function(...) FALSE,
    .package = "rlang"
  )
  testthat::expect_message(
    n <- resolve_workers("auto", 5L),
    class = "isoformic_workers_fallback"
  )
  testthat::expect_identical(n, 1L)
})

testthat::test_that("workers = 'auto' follows the scheduler allocation", {
  testthat::skip_if_not_installed("parallelly")
  withr::local_envvar(SLURM_CPUS_PER_TASK = "2")
  withr::local_options(mc.cores = NULL)
  testthat::expect_identical(resolve_workers("auto", 10L), 2L)
  testthat::expect_identical(resolve_workers("auto", 1L), 1L)
})

testthat::test_that("one worker runs the tasks in the session's random stream", {
  tasks <- list(list(n = 3L), list(n = 2L))
  set.seed(42)
  res <- run_tasks(tasks, function(n) stats::runif(n))
  set.seed(42)
  ref <- list(stats::runif(3L), stats::runif(2L))
  testthat::expect_identical(res, ref)
})

testthat::test_that("seeded tasks give the same result in the session and on workers", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("mirai")
  local_no_workers()
  tasks <- lapply(1:4, function(n) list(n = n))
  draw <- function(n) stats::runif(n)
  set.seed(7)
  in_session <- run_tasks(tasks, draw, seeded = TRUE)
  after_session <- stats::runif(1L)
  set.seed(7)
  on_two <- run_tasks(tasks, draw, workers = 2L)
  after_two <- stats::runif(1L)
  # The daemons started by the first call are kept and reused
  testthat::expect_true(mirai::daemons_set(.compute = "isoformic"))
  set.seed(7)
  on_two_kept <- run_tasks(tasks, draw, workers = 2L)
  after_two_kept <- stats::runif(1L)
  set.seed(7)
  on_three <- run_tasks(tasks, draw, workers = 3L)
  testthat::expect_identical(on_two, in_session)
  testthat::expect_identical(on_two_kept, in_session)
  testthat::expect_identical(on_three, in_session)
  # The session's stream moves by the seeds only, whatever runs the tasks
  # and whether the daemons were started or already running
  testthat::expect_identical(after_two, after_session)
  testthat::expect_identical(after_two_kept, after_session)
  # Started ahead with isoformic_workers(), after set.seed()
  set.seed(7)
  isoformic_workers(2L)
  testthat::expect_identical(run_tasks(tasks, draw, workers = 4L), in_session)
  testthat::expect_identical(stats::runif(1L), after_session)
  testthat::expect_identical(isoformic_workers(0), 0L)
  testthat::expect_false(mirai::daemons_set(.compute = "isoformic"))
})

testthat::test_that("kept daemons are reused, grown and shrunk to the request", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("mirai")
  local_no_workers()
  pids <- function(n_tasks, workers) {
    tasks <- lapply(seq_len(n_tasks), function(i) list(i = i))
    unlist(run_tasks(
      tasks,
      function(i) {
        Sys.sleep(0.2)
        Sys.getpid()
      },
      workers = workers
    ))
  }
  first <- pids(2L, 4L)
  # Started for the tasks of the call, not for all the requested workers
  testthat::expect_identical(the_workers$n, 2L)
  testthat::expect_true(the_workers$auto)
  url <- the_workers$url
  # Reused by the next call
  testthat::expect_true(all(pids(2L, 4L) %in% first))
  testthat::expect_identical(the_workers$url, url)
  # Grown, not restarted, when a call uses more
  grown <- pids(6L, 4L)
  testthat::expect_identical(the_workers$n, 4L)
  testthat::expect_identical(the_workers$url, url)
  # The launched daemons are connected when the call returns
  testthat::expect_identical(
    mirai::status(.compute = "isoformic")[["connections"]],
    4L
  )
  before <- unique(c(grown, pids(8L, 4L)))
  # Started again when a call requests fewer than are running
  shrunk <- pids(6L, 2L)
  testthat::expect_identical(the_workers$n, 2L)
  testthat::expect_lte(length(unique(shrunk)), 2L)
  testthat::expect_false(any(shrunk %in% before))
  # Daemons started ahead are used as they are, whatever `workers` is
  isoformic_workers(3L)
  testthat::expect_false(the_workers$auto)
  testthat::expect_lte(length(unique(pids(6L, 2L))), 3L)
  testthat::expect_identical(the_workers$n, 3L)
  testthat::expect_identical(
    mirai::status(.compute = "isoformic")[["connections"]],
    3L
  )
  testthat::expect_error(
    isoformic_workers(-1),
    class = "isoformic_invalid_workers"
  )
})

testthat::test_that("an error on a worker is raised again naming the task", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("mirai")
  local_no_workers()
  tasks <- lapply(1:3, function(n) list(n = n))
  fail_on_two <- function(n) {
    if (n == 2L) {
      stop("no draw for two", call. = FALSE)
    }
    n
  }
  err <- testthat::expect_error(
    run_tasks(
      tasks,
      fail_on_two,
      workers = 2L,
      labels = c("first", "second", "third")
    ),
    class = "isoformic_worker_error"
  )
  testthat::expect_match(conditionMessage(err), "second", fixed = TRUE)
  testthat::expect_match(conditionMessage(err), "no draw for two", fixed = TRUE)
  # The daemons are kept and still work
  testthat::expect_true(mirai::daemons_set(.compute = "isoformic"))
  testthat::expect_identical(
    run_tasks(tasks, function(n) n, workers = 2L),
    as.list(1:3)
  )
})

testthat::test_that("daemons that the user started are used and left running", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("mirai")
  local_no_workers()
  mirai::daemons(1L, .compute = "isoformic")
  url <- profile_url()
  tasks <- lapply(1:3, function(n) list(n = n))
  draw <- function(n) stats::runif(n)
  set.seed(3)
  ref <- run_tasks(tasks, draw, seeded = TRUE)
  set.seed(3)
  res <- run_tasks(tasks, draw, workers = 2L)
  testthat::expect_identical(res, ref)
  # Not resized, not recorded as isoformic's, and not stopped on unload
  testthat::expect_identical(profile_url(), url)
  testthat::expect_null(the_workers$url)
  testthat::expect_false(stop_workers())
  testthat::expect_true(mirai::daemons_set(.compute = "isoformic"))
})

testthat::test_that("enrichment results do not depend on the number of workers", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("fgsea")
  testthat::skip_if_not_installed("mirai")
  local_no_workers()
  iso <- example_isoformic()
  sets <- read_gmt(system.file(
    "extdata",
    "reactome_subset.gmt",
    package = "isoformic"
  ))
  enrich <- function(...) {
    set.seed(11)
    suppressWarnings(run_enrichment(
      iso,
      list(first = sets, second = rev(sets)),
      pval_cutoff = NULL,
      min_size = 3,
      ...
    ))
  }
  serial <- enrich()
  # workers = 1 is the serial run of earlier versions
  testthat::expect_identical(enrich(workers = 1L), serial)
  withr::with_options(
    list(isoformic.workers = 1L),
    testthat::expect_identical(enrich(), serial)
  )
  on_two <- enrich(workers = 2L)
  testthat::expect_identical(enrich(workers = 3L), on_two)
  withr::with_options(
    list(isoformic.workers = 3L),
    testthat::expect_identical(enrich(), on_two)
  )
  # Each call has its own seed: the sampled p-values differ from a serial
  # run, the deterministic columns do not
  testthat::expect_identical(on_two$pathway, serial$pathway)
  testthat::expect_identical(on_two$experiment, serial$experiment)
  testthat::expect_identical(on_two$collection, serial$collection)
  testthat::expect_identical(on_two$ES, serial$ES)
  testthat::expect_identical(on_two$size, serial$size)
  testthat::expect_identical(on_two$leadingEdge, serial$leadingEdge)

  usage <- apa_usage(iso, group = "condition")
  three <- list(a = sets, b = rev(sets), c = sets[1:5])
  apa <- function(workers) {
    set.seed(5)
    suppressWarnings(run_apa_enrichment(
      usage,
      three,
      min_size = 2,
      workers = workers
    ))
  }
  testthat::expect_identical(apa(3L), apa(2L))
  testthat::expect_identical(apa(1L)$ES, apa(2L)$ES)
})

testthat::test_that("enrichment errors on workers name the collection and category", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("fgsea")
  testthat::skip_if_not_installed("mirai")
  local_no_workers()
  det <- tibble::tibble(
    transcript_name = paste0("t", 1:6),
    transcript_type = rep(c("protein_coding", "retained_intron"), each = 3L),
    log2FC = c(2, -1, 0.5, -3, 1, 0.2)
  )
  map <- tibble::tibble(
    transcript_name = paste0("t", 1:6),
    gene_name = paste0("g", 1:6)
  )
  sets <- list(a = c("g1", "g2"), b = c("g4", "g5"))
  # An argument fgsea does not have fails inside every call
  testthat::expect_error(
    run_enrichment(
      det,
      list(first = sets, second = sets),
      tx_to_gene = map,
      tx_types = list(
        protein_coding = "protein_coding",
        retained_intron = "retained_intron"
      ),
      min_size = 1,
      pval_cutoff = NULL,
      workers = 2L,
      not_an_fgsea_argument = TRUE
    ),
    regexp = "collection \"first\", transcript category \"protein_coding\"",
    fixed = TRUE,
    class = "isoformic_worker_error"
  )
})

testthat::test_that("warnings of seeded tasks reach the session naming the task", {
  tasks <- lapply(1:3, function(n) list(n = n))
  warn_on_odd <- function(n) {
    if (n %% 2L == 1L) {
      warning("odd draw ", n)
      warning("odd draw ", n)
    }
    n
  }
  collect <- function(...) {
    seen <- new.env(parent = emptyenv())
    seen$msgs <- character()
    res <- withCallingHandlers(
      run_tasks(tasks, warn_on_odd, labels = c("one", "two", "three"), ...),
      isoformic_worker_warning = function(cnd) {
        seen$msgs <- c(seen$msgs, conditionMessage(cnd))
        invokeRestart("muffleWarning")
      }
    )
    list(res = res, msgs = seen$msgs)
  }
  in_session <- collect(seeded = TRUE)
  testthat::expect_identical(in_session$res, as.list(1:3))
  # One warning per distinct message and task, in task order
  testthat::expect_length(in_session$msgs, 2L)
  testthat::expect_match(in_session$msgs[[1L]], "odd draw 1", fixed = TRUE)
  testthat::expect_match(in_session$msgs[[1L]], "one", fixed = TRUE)
  testthat::expect_match(in_session$msgs[[2L]], "three", fixed = TRUE)

  testthat::skip_on_cran()
  testthat::skip_if_not_installed("mirai")
  local_no_workers()
  on_workers <- collect(workers = 2L)
  testthat::expect_identical(on_workers, in_session)
})

testthat::test_that("tasks on daemons the user started get the session's library paths", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("mirai")
  local_no_workers()
  # Started without R_LIBS: the daemons only know the default library
  mirai::daemons(1L, .compute = "isoformic")
  tasks <- list(list(n = 1L), list(n = 2L))
  libs <- run_tasks(tasks, function(n) .libPaths(), workers = 2L)
  session <- normalizePath(.libPaths(), mustWork = FALSE)
  for (lib in libs) {
    testthat::expect_identical(
      normalizePath(lib, mustWork = FALSE)[seq_along(session)],
      session
    )
  }
})

testthat::test_that("gene sets cut down to the ranking keep names, order and empty sets", {
  sets <- list(a = c("t1", "t9", "t2"), b = "t8", c = c("t3", "t1"))
  testthat::expect_identical(
    sets_in_ranking(sets, c("t1", "t2", "t3")),
    list(a = c("t1", "t2"), b = character(), c = c("t3", "t1"))
  )
  testthat::expect_identical(sets_in_ranking(list(), "t1"), list())
})

testthat::test_that("cutting the gene sets to the ranking does not change enrichment results", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("fgsea")
  testthat::skip_if_not_installed("mirai")
  local_no_workers()
  iso <- example_isoformic()
  sets <- read_gmt(system.file(
    "extdata",
    "reactome_subset.gmt",
    package = "isoformic"
  ))
  collections <- list(first = sets, second = rev(sets))
  usage <- apa_usage(iso, group = "condition")
  enrich <- function() {
    set.seed(13)
    tx <- suppressWarnings(run_enrichment(
      iso,
      collections,
      pval_cutoff = NULL,
      min_size = 3,
      workers = 2L
    ))
    set.seed(5)
    gene <- suppressWarnings(run_apa_enrichment(
      usage,
      collections,
      min_size = 2,
      workers = 2L
    ))
    list(tx = tx, gene = gene)
  }
  cut <- enrich()
  calls <- new.env(parent = emptyenv())
  calls$n <- 0L
  testthat::local_mocked_bindings(sets_in_ranking = function(sets, ids) {
    calls$n <- calls$n + 1L
    sets
  })
  testthat::expect_identical(cut, enrich())
  testthat::expect_gt(calls$n, 0L)
})

testthat::test_that("fgsea warnings are the same with one and more workers", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("fgsea")
  testthat::skip_if_not_installed("mirai")
  local_no_workers()
  det <- tibble::tibble(
    transcript_name = paste0("t", 1:20),
    transcript_type = rep(c("protein_coding", "retained_intron"), each = 10L),
    log2FC = rep(c(1, -1), 10L)
  )
  map <- tibble::tibble(
    transcript_name = paste0("t", 1:20),
    gene_name = paste0("g", 1:20)
  )
  sets <- list(a = paste0("g", 1:6), b = paste0("g", 11:16))
  warnings_of <- function(workers) {
    seen <- new.env(parent = emptyenv())
    seen$msgs <- character()
    set.seed(1)
    withCallingHandlers(
      run_enrichment(
        det,
        list(first = sets, second = sets),
        tx_to_gene = map,
        tx_types = list(
          protein_coding = "protein_coding",
          retained_intron = "retained_intron"
        ),
        min_size = 2,
        pval_cutoff = NULL,
        workers = workers
      ),
      warning = function(cnd) {
        seen$msgs <- c(seen$msgs, conditionMessage(cnd))
        invokeRestart("muffleWarning")
      }
    )
    seen$msgs
  }
  serial <- warnings_of(1L)
  testthat::expect_true(any(grepl("ties", serial, fixed = TRUE)))
  parallel <- warnings_of(2L)
  # The same messages in the same order (one per call here), each followed
  # by the call that raised it
  testthat::expect_length(parallel, length(serial))
  testthat::expect_identical(substr(parallel, 1L, nchar(serial)), serial)
  testthat::expect_match(
    parallel[[1L]],
    "collection \"first\", transcript category \"protein_coding\"",
    fixed = TRUE
  )
})

# Runs `mode` of a child R session that loads isoformic and starts kept
# daemons; returns the paths of its files once `wait_for` exists in them
child_session <- function(mode, envir = parent.frame()) {
  out <- withr::local_tempdir(.local_envir = envir)
  ns_path <- getNamespaceInfo("isoformic", "path")
  src <- if (file.exists(file.path(ns_path, "R", "parallel.R"))) ns_path else ""
  if (nzchar(src)) {
    testthat::skip_if_not_installed("pkgload")
  }
  script <- file.path(out, "child.R")
  writeLines(
    c(
      "args <- commandArgs(trailingOnly = TRUE)",
      "src <- args[[1L]]",
      "out <- args[[2L]]",
      "mode <- args[[3L]]",
      "if (nzchar(src)) {",
      "  pkgload::load_all(src, quiet = TRUE, helpers = FALSE)",
      "} else {",
      "  suppressPackageStartupMessages(library(isoformic))",
      "}",
      "ns <- asNamespace('isoformic')",
      "put <- function(x, name) {",
      "  writeLines(as.character(x), file.path(out, paste0(name, '.tmp')))",
      "  invisible(file.rename(",
      "    file.path(out, paste0(name, '.tmp')),",
      "    file.path(out, name)",
      "  ))",
      "}",
      "put(Sys.getpid(), 'session')",
      "tasks <- lapply(1:2, function(i) list(i = i))",
      "pid <- function(i) {",
      "  Sys.sleep(0.5)",
      "  Sys.getpid()",
      "}",
      "daemons <- unique(unlist(ns$run_tasks(tasks, pid, workers = 2L)))",
      "put(daemons, 'daemons')",
      "if (mode == 'error') {",
      "  try(ns$run_tasks(tasks, function(i) stop('boom'), workers = 2L))",
      "  stop('the session ends on an error')",
      "}",
      "if (mode == 'kill') {",
      "  put('', 'ready')",
      "  Sys.sleep(120)",
      "}",
      "if (mode == 'dead') {",
      "  workers <- ns$the_workers",
      "  workers$grace <- 5",
      "  draws <- lapply(1:4, function(n) list(n = n))",
      "  draw <- function(n) stats::runif(n)",
      "  set.seed(1)",
      "  ref <- ns$run_tasks(draws, draw, seeded = TRUE)",
      "  connected_below <- function(n) {",
      "    while (ns$profile_connections() >= n) Sys.sleep(0.05)",
      "  }",
      "  # All the kept daemons die between two calls (an OOM kill)",
      "  for (p in daemons) tools::pskill(p, tools::SIGKILL)",
      "  connected_below(1L)",
      "  set.seed(1)",
      "  put(identical(ns$run_tasks(draws, draw, workers = 2L), ref), 'all')",
      "  # One of the two dies",
      "  pids <- unique(unlist(ns$run_tasks(tasks, pid, workers = 2L)))",
      "  tools::pskill(pids[[1L]], tools::SIGKILL)",
      "  connected_below(2L)",
      "  pids <- unique(unlist(ns$run_tasks(tasks, pid, workers = 2L)))",
      "  put(length(pids), 'one')",
      "  # All of them die during a call, with tasks still queued",
      "  die <- function(n) tools::pskill(Sys.getpid(), tools::SIGKILL)",
      "  err <- tryCatch(",
      "    ns$run_tasks(draws, die, workers = 2L),",
      "    isoformic_worker_error = function(cnd) 'error'",
      "  )",
      "  put(err, 'during')",
      "  set.seed(1)",
      "  res <- ns$run_tasks(draws, draw, workers = 2L)",
      "  put(identical(res, ref), 'after')",
      "  pids <- unique(unlist(ns$run_tasks(tasks, pid, workers = 2L)))",
      "  put(pids, 'last')",
      "  put('', 'done')",
      "}"
    ),
    script
  )
  system2(
    file.path(R.home("bin"), "Rscript"),
    c(shQuote(script), shQuote(src), shQuote(out), mode),
    stdout = file.path(out, "log"),
    stderr = file.path(out, "log"),
    wait = FALSE,
    env = paste0(
      "R_LIBS=",
      shQuote(paste(.libPaths(), collapse = .Platform$path.sep))
    )
  )
  out
}

wait_for_file <- function(path, seconds = 120) {
  for (i in seq_len(seconds * 10L)) {
    if (file.exists(path)) {
      return(TRUE)
    }
    Sys.sleep(0.1)
  }
  FALSE
}

pid_alive <- function(pid) isTRUE(tools::pskill(pid, 0L))

# Waits up to `seconds` for every process of `pids` to end
processes_end <- function(pids, seconds = 20) {
  for (i in seq_len(seconds * 10L)) {
    if (!any(vapply(pids, pid_alive, logical(1L)))) {
      return(TRUE)
    }
    Sys.sleep(0.1)
  }
  FALSE
}

testthat::test_that("an interrupt cancels the tasks and keeps the daemons", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("mirai")
  local_no_workers()
  tasks <- lapply(1:2, function(i) list(i = i))
  pid <- function(i) {
    Sys.sleep(0.3)
    Sys.getpid()
  }
  daemons <- unique(unlist(run_tasks(tasks, pid, workers = 2L)))
  long <- function(i) {
    Sys.sleep(60)
    i
  }
  # The wait is interrupted (a SIGINT reaches it only on the thread that
  # waits, so the test raises the interrupt itself)
  testthat::local_mocked_bindings(wait_tasks = function(m) {
    Sys.sleep(1)
    rlang::interrupt()
  })
  res <- tryCatch(
    run_tasks(lapply(1:4, function(i) list(i = i)), long, workers = 2L),
    interrupt = function(cnd) "interrupted"
  )
  testthat::expect_identical(res, "interrupted")
  testthat::local_mocked_bindings(wait_tasks = function(m) m[])
  # The next call does not wait for the cancelled 60-second tasks, and runs
  # on the same daemons
  start <- Sys.time()
  again <- unique(unlist(run_tasks(tasks, pid, workers = 2L)))
  testthat::expect_lt(
    as.numeric(difftime(Sys.time(), start, units = "secs")),
    30
  )
  testthat::expect_true(all(again %in% daemons))
  testthat::expect_identical(
    unname(mirai::info(.compute = "isoformic")[c("awaiting", "executing")]),
    c(0L, 0L)
  )
})

testthat::test_that("the daemons end with a session that ends normally", {
  testthat::skip_on_cran()
  testthat::skip_on_os("windows")
  testthat::skip_if_not_installed("mirai")
  out <- child_session("normal")
  testthat::expect_true(wait_for_file(file.path(out, "daemons")))
  session <- as.integer(readLines(file.path(out, "session")))
  daemons <- as.integer(readLines(file.path(out, "daemons")))
  testthat::expect_length(daemons, 2L)
  testthat::expect_true(processes_end(c(session, daemons), seconds = 60))
})

testthat::test_that("the daemons end with a session that ends on an error", {
  testthat::skip_on_cran()
  testthat::skip_on_os("windows")
  testthat::skip_if_not_installed("mirai")
  out <- child_session("error")
  testthat::expect_true(wait_for_file(file.path(out, "daemons")))
  session <- as.integer(readLines(file.path(out, "session")))
  daemons <- as.integer(readLines(file.path(out, "daemons")))
  testthat::expect_length(daemons, 2L)
  testthat::expect_true(processes_end(c(session, daemons), seconds = 60))
  testthat::expect_match(
    paste(readLines(file.path(out, "log")), collapse = "\n"),
    "the session ends on an error",
    fixed = TRUE
  )
})

testthat::test_that("the daemons end with a session that is killed", {
  testthat::skip_on_cran()
  testthat::skip_on_os("windows")
  testthat::skip_if_not_installed("mirai")
  out <- child_session("kill")
  testthat::expect_true(wait_for_file(file.path(out, "ready")))
  session <- as.integer(readLines(file.path(out, "session")))
  daemons <- as.integer(readLines(file.path(out, "daemons")))
  testthat::expect_true(all(vapply(daemons, pid_alive, logical(1L))))
  tools::pskill(session, tools::SIGKILL)
  testthat::expect_true(processes_end(c(session, daemons)))
})

testthat::test_that("kept daemons that died are started again, without a hang", {
  testthat::skip_on_cran()
  testthat::skip_on_os("windows")
  testthat::skip_if_not_installed("mirai")
  # A child session, so that a hang fails the test instead of blocking it
  out <- child_session("dead")
  testthat::expect_true(wait_for_file(file.path(out, "session")))
  session <- as.integer(readLines(file.path(out, "session")))
  # Until the child is done, or ends early on an error
  for (i in seq_len(1800L)) {
    if (file.exists(file.path(out, "done")) || !pid_alive(session)) {
      break
    }
    Sys.sleep(0.1)
  }
  done <- file.exists(file.path(out, "done"))
  if (!done) {
    tools::pskill(session, tools::SIGKILL)
  }
  testthat::expect_true(
    done,
    label = paste(readLines(file.path(out, "log")), collapse = "\n")
  )
  read <- function(name) readLines(file.path(out, name))
  # Same result as before the daemons died
  testthat::expect_identical(read("all"), "TRUE")
  # The dead daemon is replaced: the call runs on two again
  testthat::expect_identical(read("one"), "2")
  # A call whose daemons all died fails instead of waiting forever, and the
  # next call starts them again
  testthat::expect_identical(read("during"), "error")
  testthat::expect_identical(read("after"), "TRUE")
  # No daemon outlives the session
  last <- as.integer(read("last"))
  testthat::expect_true(processes_end(c(session, last), seconds = 60))
})

testthat::test_that("daemons the user started that are not connected fail the call", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("mirai")
  local_no_workers()
  old <- the_workers$grace
  the_workers$grace <- 1
  withr::defer(assign("grace", old, envir = the_workers))
  # A profile that listens, with no daemon connected to it
  mirai::daemons(url = mirai::local_url(), .compute = "isoformic")
  tasks <- lapply(1:2, function(n) list(n = n))
  err <- testthat::expect_error(
    run_tasks(tasks, function(n) n, workers = 2L),
    class = "isoformic_worker_error"
  )
  testthat::expect_match(
    conditionMessage(err),
    "isoformic_workers(0)",
    fixed = TRUE
  )
  # Left running: isoformic does not stop daemons it did not start
  testthat::expect_true(mirai::daemons_set(.compute = "isoformic"))
})
