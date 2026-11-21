# Shared helpers of the uncertainty evidence scripts (sourced, not run).
#
# Purpose: one place for the paths, the frozen code snapshot, the provenance
# columns and the Parquet writers that every script of
# data-raw/uncertainty_evidence/ uses, so that all tables carry the same
# provenance and are read from the same frozen package code.
#
# Env vars:
#   ISOFORMIC_UE_OUT      output root of the evaluation (required: every
#                         script stops with a message when it is unset)
#   ISOFORMIC_V50_DIR     work directory with the GENCODE 50 references and
#                         the Salmon runs (required by the scripts that read
#                         quantifications; layout in README.md)
#   ISOFORMIC_UE_RLIB     optional R library prepended to .libPaths() by the
#                         WP-D scripts, default $ISOFORMIC_V50_DIR/rlib when
#                         that directory exists
#   ISOFORMIC_UE_THREADS  threads for DuckDB and data.table, default 4
#   ISOFORMIC_UE_PKG      path of the package snapshot to load; default: the
#                         snapshot recorded in $ISOFORMIC_UE_OUT/prereg.yaml
#                         (written by 01_preregister.R). Another path (the
#                         repository, for the `productivity` stage of
#                         30_audit_prevalence.R, amendment A18) stamps the
#                         tables with "<HEAD sha>+working-tree" and no
#                         snapshot sha256
#
# Use from a script run at the repository root:
#   source("data-raw/uncertainty_evidence/ue_common.R")
#   ue <- ue_setup("10_placenta_prepare")   # loads the snapshot
#
# The live tree under R/ is loaded only where ISOFORMIC_UE_PKG says so
# (amendment A18): another workflow edits it, and the evidence must come from
# a fixed commit. 01_preregister.R creates the snapshot with `git archive`
# (no checkout).

# Value of a required environment variable; stops when it is unset
ue_env <- function(name, what) {
  value <- Sys.getenv(name)
  if (!nzchar(value)) {
    stop(
      "Set the environment variable ",
      name,
      " (",
      what,
      "); see ",
      "data-raw/uncertainty_evidence/README.md",
      call. = FALSE
    )
  }
  value
}

ue_out_root <- function() {
  ue_env("ISOFORMIC_UE_OUT", "output root of the uncertainty evaluation")
}

ue_v50_dir <- function() {
  ue_env(
    "ISOFORMIC_V50_DIR",
    "work directory with the GENCODE 50 references and Salmon runs"
  )
}

# Prepend the optional extra R library (call before ue_setup())
ue_use_rlib <- function() {
  rlib <- Sys.getenv("ISOFORMIC_UE_RLIB")
  v50 <- Sys.getenv("ISOFORMIC_V50_DIR")
  if (!nzchar(rlib) && nzchar(v50)) {
    rlib <- file.path(v50, "rlib")
  }
  if (nzchar(rlib) && dir.exists(rlib)) {
    .libPaths(c(rlib, .libPaths()))
  }
  invisible(rlib)
}

ue_threads <- function() {
  as.integer(Sys.getenv("ISOFORMIC_UE_THREADS", "4"))
}

ue_prereg <- function() {
  path <- file.path(ue_out_root(), "prereg.yaml")
  if (!file.exists(path)) {
    stop(
      "Run data-raw/uncertainty_evidence/01_preregister.R first: ",
      path,
      call. = FALSE
    )
  }
  yaml::read_yaml(path)
}

# Load the frozen package code and return the provenance of the run
ue_setup <- function(step) {
  pre <- ue_prereg()
  pkg <- Sys.getenv("ISOFORMIC_UE_PKG", pre$snapshot$path)
  if (!dir.exists(pkg)) {
    stop("Package snapshot not found: ", pkg, call. = FALSE)
  }
  suppressMessages(devtools::load_all(pkg, quiet = TRUE, export_all = TRUE))
  data.table::setDTthreads(ue_threads())
  # Code other than the snapshot (amendment A18): the tables say so, with the
  # commit of that tree and no snapshot checksum
  is_snapshot <- identical(
    normalizePath(pkg, mustWork = FALSE),
    normalizePath(pre$snapshot$path, mustWork = FALSE)
  )
  git_sha <- pre$snapshot$git_sha
  snapshot_sha256 <- pre$snapshot$tar_sha256
  if (!is_snapshot) {
    head_sha <- tryCatch(
      system2("git", c("-C", pkg, "rev-parse", "HEAD"), stdout = TRUE),
      error = function(e) NA_character_,
      warning = function(w) NA_character_
    )
    git_sha <- paste0(head_sha[1L], "+working-tree")
    snapshot_sha256 <- NA_character_
  }
  list(
    step = step,
    pkg_path = pkg,
    git_sha = git_sha,
    snapshot_sha256 = snapshot_sha256,
    seed_base = pre$seeds$base,
    started = Sys.time(),
    prereg = pre
  )
}

# Provenance columns added to every table written
ue_stamp <- function(x, ue, seed = NA_integer_) {
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  x$git_sha <- ue$git_sha
  x$snapshot_sha256 <- ue$snapshot_sha256
  x$seed <- as.integer(seed)
  x$script <- ue$step
  x
}

# Parquet through the package writer (DuckDB)
ue_write <- function(x, path, ue, seed = NA_integer_) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write_parquet(ue_stamp(x, ue, seed), path)
  invisible(path)
}

# Parquet through DuckDB with numeric columns cast to FLOAT (single
# precision), for the large draw tables; write_parquet() keeps doubles.
ue_write_float <- function(x, path, float_cols, ue) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  x <- ue_stamp(x, ue)
  con <- DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  DBI::dbExecute(con, sprintf("SET threads = %d", ue_threads()))
  duckdb::duckdb_register(con, "ue_tbl", x)
  cols <- vapply(
    colnames(x),
    function(col) {
      q <- DBI::dbQuoteIdentifier(con, col)
      if (col %in% float_cols) sprintf("CAST(%s AS FLOAT) AS %s", q, q) else q
    },
    character(1)
  )
  DBI::dbExecute(
    con,
    sprintf(
      "COPY (SELECT %s FROM ue_tbl) TO %s (FORMAT parquet, COMPRESSION zstd)",
      paste(cols, collapse = ", "),
      DBI::dbQuoteString(con, normalizePath(path, mustWork = FALSE))
    )
  )
  invisible(path)
}

# Session information as a one-column table, written next to each output set
ue_write_session <- function(dir, ue) {
  si <- utils::capture.output(utils::sessionInfo())
  ue_write(data.frame(line = si), file.path(dir, "session_info.parquet"), ue)
}

# Peak resident memory of this R process (Linux), in GB
ue_peak_rss_gb <- function() {
  st <- readLines("/proc/self/status")
  hwm <- grep("^VmHWM", st, value = TRUE)
  as.numeric(gsub("[^0-9]", "", hwm)) / 1024^2
}

# Overdispersion (Baldoni et al. 2024, edgeR::catchSalmon) of every row from
# the per sample moments of the draws: the package estimator, frozen in the
# snapshot. `mean`, `var` are transcript x sample matrices, `n` the draws.
ue_overdispersion <- function(mean, var, n, prior_df = 3) {
  overdispersion_from_moments(
    list(mean = mean, var = var, n = n),
    prior_df = prior_df
  )
}

# InfRV (fishpond form) per transcript and sample from the moments
ue_infrv <- function(mean, var) {
  inferential_rv(list(mean = mean, var = var))
}
