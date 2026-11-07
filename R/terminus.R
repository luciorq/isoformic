#' Read the transcript groups that Terminus builds
#'
#' Terminus groups the transcripts whose abundances the reads cannot estimate
#' separately and reports one consensus grouping for a set of samples. Run it
#' on Salmon output made with `--numGibbsSamples` (or `--numBootstraps`),
#' `--dumpEq` and `--dumpEqWeights`:
#'
#' ```
#' terminus group -d quant/sample1 -o terminus/out
#' terminus collapse -d quant/sample1 quant/sample2 -o terminus/out
#' ```
#'
#' `collapse` writes `clusters.txt` in the output folder of every sample: one
#' line per group, the name of the group followed by its transcripts, separated
#' by commas. `read_terminus_groups()` reads that file into a long table that
#' [add_uncertainty()] accepts in its `groups` argument, so that the groups of
#' Terminus replace the groups built from the correlation of the replicates.
#'
#' @param path The `clusters.txt` file, or a folder that holds it (the output
#'   folder of one sample, or the folder given to `-o` when it holds one sample
#'   folder).
#' @param ignore_after_bar Logical. Keep only the part of transcript names
#'   before the first `|`.
#'
#' @returns A tibble with `group_id`, `transcript_id` and `group_size`.
#'
#' @seealso [add_uncertainty()], [isoform_ambiguity()]
#'
#' @examples
#' clusters <- tempfile(fileext = ".txt")
#' writeLines(
#'   c(
#'     "NewTr1,ENST00000831217.1,ENST00000831165.1",
#'     "NewTr2,ENST00000831423.1,ENST00000831157.1,ENST00000831246.1"
#'   ),
#'   clusters
#' )
#' read_terminus_groups(clusters)
#'
#' @export
read_terminus_groups <- function(path, ignore_after_bar = TRUE) {
  rlang::check_required(path)
  path <- as.character(path)
  file <- terminus_clusters_file(path)
  lines <- readLines(file, warn = FALSE)
  lines <- lines[nzchar(trimws(lines))]
  if (!length(lines)) {
    cli::cli_abort(
      message = c(x = "{.file {file}} has no groups."),
      class = "isoformic_invalid_terminus_groups"
    )
  }
  fields <- strsplit(lines, ",", fixed = TRUE)
  sizes <- lengths(fields) - 1L
  if (any(sizes < 2L)) {
    cli::cli_abort(
      message = c(
        x = "{.file {file}} is not a Terminus {.file clusters.txt}.",
        i = "Each line needs a group name and at least two transcripts."
      ),
      class = "isoformic_invalid_terminus_groups"
    )
  }
  transcript_id <- unlist(lapply(fields, `[`, -1L), use.names = FALSE)
  if (isTRUE(ignore_after_bar)) {
    transcript_id <- sub("\\|.*$", "", transcript_id)
  }
  res <- tibble::tibble(
    group_id = rep(vapply(fields, `[`, character(1L), 1L), times = sizes),
    transcript_id = transcript_id,
    group_size = rep(sizes, times = sizes)
  )
  attr(res, "source") <- "terminus"
  res
}

#' Locate the `clusters.txt` of a Terminus run
#' @keywords internal
#' @noRd
terminus_clusters_file <- function(path) {
  if (isTRUE(fs::is_file(path))) {
    return(path)
  }
  assert_dir_exists(path, arg = "path")
  direct <- fs::path(path, "clusters.txt")
  if (isTRUE(fs::is_file(direct))) {
    return(as.character(direct))
  }
  found <- fs::dir_ls(path, recurse = 1L, glob = "*clusters.txt", type = "file")
  if (!length(found)) {
    cli::cli_abort(
      message = c(
        x = "No {.file clusters.txt} found in {.path {path}}.",
        i = "Run {.code terminus collapse} to write it."
      ),
      class = "isoformic_invalid_terminus_groups"
    )
  }
  as.character(found[[1L]])
}

#' Transcript to group table for the `groups` argument of add_uncertainty()
#' @keywords internal
#' @noRd
resolve_groups <- function(groups) {
  if (is.null(groups)) {
    return(NULL)
  }
  if (is.character(groups) && length(groups) == 1L) {
    groups <- read_terminus_groups(groups)
  }
  if (
    !is.data.frame(groups) ||
      !all(c("transcript_id", "group_id") %in% names(groups))
  ) {
    cli::cli_abort(
      message = c(
        x = paste0(
          "{.arg groups} needs a table with {.field transcript_id} and ",
          "{.field group_id}, or the path to Terminus output."
        )
      ),
      class = "isoformic_invalid_groups"
    )
  }
  groups <- groups[!is.na(groups$group_id) & !is.na(groups$transcript_id), ]
  if (anyDuplicated(groups$transcript_id)) {
    cli::cli_abort(
      message = c(x = "A transcript is in more than one group."),
      class = "isoformic_invalid_groups"
    )
  }
  res <- tibble::tibble(
    transcript_id = as.character(groups$transcript_id),
    group_id = as.character(groups$group_id)
  )
  attr(res, "source") <- attr(groups, "source")
  attr(res, "cutoff") <- attr(groups, "cutoff")
  res
}

#' Run Terminus on Salmon output from R
#'
#' `run_terminus()` builds the consensus groups of Terminus for a set of
#' samples without leaving R: it runs `terminus group` on every Salmon output
#' directory, then `terminus collapse` on all of them, and returns the groups
#' read by [read_terminus_groups()], ready for `add_uncertainty(groups = )`.
#' Terminus is installed from Bioconda into a conda environment managed by
#' `condathis` the first time it is needed.
#'
#' Every directory needs inferential replicates (`--numGibbsSamples` or
#' `--numBootstraps`) and equivalence classes with their weights (`--dumpEq
#' --dumpEqWeights`); [salmon_quant()] writes both. Without the weights
#' Terminus stops with an unhelpful error, so the directories are checked
#' first. Terminus names the folder of each sample after the basename of its
#' Salmon directory, which must therefore be unique.
#'
#' `terminus group` writes a log of several gigabytes per human sample
#' (`delta.log`) that `collapse` does not use; it is deleted as soon as each
#' sample is grouped, and the other logs after the consensus, unless
#' `keep_logs = TRUE`.
#'
#' @param paths Salmon output directories (or paths to their `quant.sf`).
#' @param output_dir Folder for the Terminus output; one subfolder per sample
#'   holds `groups.txt` and, after the consensus, `clusters.txt` and the
#'   collapsed quantification.
#' @param consensus_thresh Share of samples in which two transcripts must be
#'   grouped for the consensus to group them (`--consensus-thresh`).
#' @param min_spread Minimum spread of a transcript for its edges to be
#'   collapse candidates (`--min-spread`).
#' @param threads Threads for writing the collapsed quantifications.
#' @param workers Number of samples grouped at the same time (in parallel
#'   processes on Unix-like systems). Each one needs several gigabytes of
#'   memory for a human transcriptome.
#' @param keep_logs Logical, keep the logs of `terminus group`.
#' @param env_name Name of the conda environment with Terminus.
#' @param verbose Passed to `condathis::run()`.
#'
#' @returns The consensus groups, as returned by [read_terminus_groups()],
#'   with the output folder in the `output_dir` attribute.
#'
#' @seealso [read_terminus_groups()], [add_uncertainty()], [salmon_quant()]
#'
#' @examples
#' \dontrun{
#' groups <- run_terminus(quant_dirs, output_dir = "terminus")
#' iso <- add_uncertainty(iso, paths = quant_dirs, groups = groups)
#' tx_uncertainty(iso)
#' }
#'
#' @export
run_terminus <- function(
  paths,
  output_dir,
  consensus_thresh = 0.5,
  min_spread = 0.1,
  threads = 4L,
  workers = 1L,
  keep_logs = FALSE,
  env_name = "terminus-env",
  verbose = "output"
) {
  rlang::check_required(paths)
  rlang::check_required(output_dir)
  rlang::check_installed("condathis", reason = "to run Terminus.")
  dirs <- resolve_quant_dirs(paths)
  if (!identical(unique(dirs$kind), "salmon")) {
    cli::cli_abort(
      message = c(
        x = "Terminus works on Salmon output, not on Oarfish output.",
        i = "Pass Salmon output directories in {.arg paths}."
      ),
      class = "isoformic_terminus_input"
    )
  }
  folders <- as.character(fs::path_file(dirs$dir))
  if (anyDuplicated(folders)) {
    cli::cli_abort(
      message = c(
        x = paste0(
          "Terminus names each sample after its directory, and these repeat: ",
          "{.val {unique(folders[duplicated(folders)])}}."
        ),
        i = "Give every Salmon output directory a different name."
      ),
      class = "isoformic_terminus_input"
    )
  }
  for (dir in dirs$dir) {
    check_terminus_input(dir)
  }
  check_terminus_env(env_name)
  output_dir <- as.character(fs::path_abs(output_dir))
  fs::dir_create(output_dir)

  group_one <- function(i) {
    condathis::run(
      "terminus",
      "group",
      "--dir",
      dirs$dir[[i]],
      "--out",
      output_dir,
      "--min-spread",
      min_spread,
      env_name = env_name,
      verbose = verbose
    )
    if (!isTRUE(keep_logs)) {
      unlink(fs::path(output_dir, folders[[i]], "delta.log"))
    }
    invisible(folders[[i]])
  }
  if (workers > 1L && .Platform$OS.type == "unix") {
    # A failed sample comes back as a try-error, reported by the error below;
    # mclapply's generic warning about it adds nothing.
    done <- suppressWarnings(parallel::mclapply(
      seq_len(nrow(dirs)),
      group_one,
      mc.cores = workers
    ))
    failed <- vapply(done, inherits, logical(1L), what = "try-error")
    if (any(failed)) {
      cli::cli_abort(
        message = c(
          x = "{.code terminus group} failed for {.val {folders[failed]}}.",
          i = as.character(done[failed][[1L]])
        ),
        class = "isoformic_terminus_failed"
      )
    }
  } else {
    lapply(seq_len(nrow(dirs)), group_one)
  }
  condathis::run(
    "terminus",
    "collapse",
    "--dirs",
    dirs$dir,
    "--out",
    output_dir,
    "--consensus-thresh",
    consensus_thresh,
    "--threads",
    threads,
    env_name = env_name,
    verbose = verbose
  )
  if (!isTRUE(keep_logs)) {
    logs <- c("delta.log", "die_roll.log", "infrv.log")
    unlink(fs::path(
      output_dir,
      rep(folders, times = length(logs)),
      rep(logs, each = length(folders))
    ))
  }
  groups <- read_terminus_groups(fs::path(output_dir, folders[[1L]]))
  cli::cli_inform(c(
    v = paste0(
      "Terminus grouped {length(unique(groups$transcript_id))} transcripts ",
      "into {length(unique(groups$group_id))} groups ({nrow(dirs)} sample{?s})."
    )
  ))
  attr(groups, "output_dir") <- output_dir
  groups
}

#' Create the conda environment that holds Terminus
#' @keywords internal
#' @noRd
create_terminus_env <- function(
  env_name = "terminus-env",
  version = "0.1.0",
  force = FALSE
) {
  rlang::check_installed("condathis")
  condathis::create_env(
    packages = paste0("bioconda::terminus==", version),
    env_name = env_name,
    additional_channels = "bioconda",
    verbose = "output",
    overwrite = force
  )
}

#' Create the Terminus environment when it is missing
#' @keywords internal
#' @noRd
check_terminus_env <- function(env_name = "terminus-env") {
  if (isFALSE(condathis::env_exists(env_name))) {
    rlang::try_fetch(
      create_terminus_env(env_name = env_name),
      error = function(e) {
        cli::cli_abort(
          message = c(
            x = "The conda environment {.val {env_name}} could not be created.",
            i = "Original error message: {conditionMessage(e)}"
          ),
          class = "isoformic_terminus_env_create_error"
        )
      }
    )
  }
  invisible(TRUE)
}

#' Check that a Salmon directory has what Terminus needs
#'
#' Replicates, equivalence classes, and the weights of the classes: a class
#' line holds the class size k, k transcript indices, k weights when they were
#' written, and the count, so 2k + 2 fields with weights and k + 2 without.
#' `cmd_info.json` does not record whether the weights were requested.
#' @keywords internal
#' @noRd
check_terminus_input <- function(dir) {
  boot <- fs::path(dir, "aux_info", "bootstrap", "bootstraps.gz")
  if (!isTRUE(fs::file_exists(boot))) {
    cli::cli_abort(
      message = c(
        x = "No inferential replicates in {.path {dir}}.",
        i = "Quantify with {.code --numGibbsSamples 30} (the default of {.fn salmon_quant})."
      ),
      class = "isoformic_no_inferential_replicates"
    )
  }
  eq <- fs::path(dir, "aux_info", "eq_classes.txt.gz")
  if (!isTRUE(fs::file_exists(eq))) {
    eq <- fs::path(dir, "aux_info", "eq_classes.txt")
  }
  if (!isTRUE(fs::file_exists(eq))) {
    cli::cli_abort(
      message = c(
        x = "No equivalence classes in {.path {dir}}.",
        i = "Quantify with {.code --dumpEq --dumpEqWeights}."
      ),
      class = "isoformic_no_equivalence_classes"
    )
  }
  con <- gzfile(as.character(eq), open = "r")
  on.exit(close(con), add = TRUE)
  n_tx <- as.integer(readLines(con, n = 1L))
  readLines(con, n = 1L + n_tx)
  first <- readLines(con, n = 1L)
  fields <- strsplit(first, "\t", fixed = TRUE)[[1L]]
  size <- as.integer(fields[[1L]])
  if (length(fields) != 2L * size + 2L) {
    cli::cli_abort(
      message = c(
        x = "The equivalence classes of {.path {dir}} have no weights.",
        i = "Terminus needs them: quantify with {.code --dumpEq --dumpEqWeights} ({.fn salmon_quant} does)."
      ),
      class = "isoformic_no_equivalence_weights"
    )
  }
  invisible(TRUE)
}
