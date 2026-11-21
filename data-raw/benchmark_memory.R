# Memory benchmark of the inferential replicate (draws) and annotation code
# paths
#
# Purpose: measure how much memory the steps that carry the Salmon Gibbs
# draws, and the genome-wide annotation steps, take on a genome-wide
# quantification, so that changes can be compared with a recorded baseline.
# For each step the script reports the peak resident set size (RSS) above the
# RSS at the start of the step, the time it took and, for the draws steps,
# that peak as a multiple of the size of one copy of the draws (transcripts
# kept x samples x draws, 8 bytes per count). The draws steps:
#
#   read_inferential_replicates  read the draws of the expressed transcripts
#   add_uncertainty              add_uncertainty(paths = ) on the object
#   productive_fraction          productive_fraction(group = ) with the draws
#   subset_samples               x[, j] to a contrast of 3 + 3 samples (the
#                                draws of the subset are a view of those of
#                                the object, not a copy)
#   write_isoformic              write_isoformic() of the object, which
#                                writes the draws to inf_reps/
#   read_isoformic               read_isoformic() of that directory, which
#                                reads the draws back
#   inferential_overdispersion   inferential_overdispersion() of the object
#   subset_overdispersion        inferential_overdispersion() of the 3 + 3
#                                subset, which reads the view one draw at a
#                                time
#   subset_write                 write_isoformic() of the 3 + 3 subset, which
#                                writes the selected draws
#
# The annotation steps (no draws; their multiple column is not meaningful):
#
#   prepare_annotation           prepare_isoformic_annotation() of the GFF3
#                                into a directory under the output directory
#                                (deleted at the end)
#   transcript_features          transcript_features() of the annotation
#                                directory, genome-wide
#   splice_site_strength         splice_site_strength() of the annotation
#                                directory, genome-wide, by transcript
#   splice_site_strength_gene    splice_site_strength(gene_name = "RBPJ") of
#                                the object
#   poison_exons                 poison_exons() of the object, all genes
#   transcript_sequences         transcript_sequences() of the transcripts of
#                                six genes
#
# Each step runs in its own R process (this script calls itself through
# Rscript), so the peaks of different steps do not mix. Inside the process,
# the setup (load_all(), reading the quantifications and, for the steps after
# add_uncertainty, attaching the draws) runs first; then the kernel peak RSS
# counter (VmHWM) is reset through /proc/self/clear_refs (Linux) and the step
# runs. Where the counter cannot be reset, the peak since the start of the
# process is used and the `peak_reset` column is FALSE.
#
# Usage (from the package root; `just benchmark-memory` runs the same):
#   ISOFORMIC_FULL=true \
#     ISOFORMIC_BENCHMARK_QUANTS=/path/to/salmon/quants \
#     ISOFORMIC_BENCHMARK_OUT=/path/to/output_dir \
#     Rscript data-raw/benchmark_memory.R
# Environment variables:
#   ISOFORMIC_FULL               must be "true": the script downloads the
#                                GENCODE annotation when it is not cached
#   ISOFORMIC_BENCHMARK_QUANTS   directory with one Salmon output directory
#                                per sample (named by sample ID), quantified
#                                with --numGibbsSamples (required)
#   ISOFORMIC_BENCHMARK_OUT      output directory (required); receives
#                                benchmark_memory.parquet, one Parquet file
#                                per step and the samples used
#   ISOFORMIC_BENCHMARK_SAMPLES  number of samples, split evenly between the
#                                two conditions of the sample sheet
#                                (default 12)
#   ISOFORMIC_BENCHMARK_DRAWS    number of draws kept (default 25)
#   ISOFORMIC_BENCHMARK_SHEET    tab-separated sample sheet with columns
#                                `run` and `condition` (default
#                                data-raw/placenta/placenta_samples.tsv, the
#                                PRJNA622166 placentas)
#   ISOFORMIC_BENCHMARK_GENCODE  GENCODE release of the annotation (default:
#                                the package default)
#   ISOFORMIC_BENCHMARK_ANNOT    an isoformic annotation directory to use
#                                instead of the GENCODE one
#                                (the ORF cache of an annotation directory is
#                                tied to its path and file times: in a fresh
#                                copy, the first step that classifies
#                                productivity (productive_fraction) also
#                                recomputes the ORFs genome-wide (about 4 more
#                                minutes and 4 to 6 GB more at GENCODE 50).
#                                Run the benchmark twice, or warm the copy
#                                first, to compare that step).
#   ISOFORMIC_BENCHMARK_GFF      GFF3 file of the prepare_annotation step
#                                (default: download_reference() of the
#                                GENCODE release, from the cache)
#   ISOFORMIC_BENCHMARK_STEPS    comma-separated names of the steps to run
#                                (default: all)
#   ISOFORMIC_BENCHMARK_PKG      package source loaded with load_all(), for
#                                example an older copy of the package to
#                                record a baseline (default: the current
#                                directory)
#
# Expected memory: one copy of the draws of the expressed transcripts of
# GENCODE 50 for 12 placenta samples and 25 draws is about 0.44 GB (about
# 200,000 transcripts with a mean count of at least 1). Each process also
# holds the annotation and the assays (about 1 to 2 GB), and the steps peak
# at up to about 10 copies of the draws (about 6 GB in all). The annotation
# steps peak at up to about 6 GB (prepare_isoformic_annotation() of the
# GENCODE 50 GFF3 before its parse was split), so a 20 GB cap is ample. Only
# one step process runs at a time. On a shared machine run it under a memory
# cap.

devtools::load_all(
  Sys.getenv("ISOFORMIC_BENCHMARK_PKG", unset = "."),
  quiet = TRUE
)

env_or <- function(name, default = NULL) {
  value <- Sys.getenv(name, unset = "")
  if (nzchar(value)) value else default
}

# Memory counters of this process ---------------------------------------------

proc_status_kb <- function(field) {
  lines <- readLines("/proc/self/status", warn = FALSE)
  line <- grep(paste0("^", field, ":"), lines, value = TRUE)
  if (!length(line)) {
    return(NA_real_)
  }
  as.numeric(gsub("[^0-9]", "", line[[1L]]))
}

reset_peak_rss <- function() {
  tryCatch(
    {
      con <- file("/proc/self/clear_refs", open = "w")
      on.exit(close(con))
      writeLines("5", con)
      TRUE
    },
    error = function(cnd) FALSE,
    warning = function(cnd) FALSE
  )
}

draws_bytes <- function(reps) {
  if (is.null(reps) || !length(reps)) {
    return(NA_real_)
  }
  per <- vapply(
    reps,
    function(m) {
      length(m) * if (is.integer(m)) 4 else 8
    },
    numeric(1L)
  )
  sum(per)
}

# Inputs -----------------------------------------------------------------------

out_dir <- env_or("ISOFORMIC_BENCHMARK_OUT")
quants_dir <- env_or("ISOFORMIC_BENCHMARK_QUANTS")
if (is.null(out_dir) || is.null(quants_dir)) {
  stop(
    "Set ISOFORMIC_BENCHMARK_QUANTS and ISOFORMIC_BENCHMARK_OUT.",
    call. = FALSE
  )
}
n_draws <- as.integer(env_or("ISOFORMIC_BENCHMARK_DRAWS", "25"))
samples_file <- fs::path(out_dir, "samples.parquet")
written_dir <- fs::path(out_dir, "written_experiment")
written_subset_dir <- fs::path(out_dir, "written_subset")
prepared_dir <- fs::path(out_dir, "prepared_annotation")
benchmark_genes <- c("TP53", "BRCA1", "EGFR", "RBPJ", "GAPDH", "MYC")

load_samples <- function() {
  sel <- read_parquet(samples_file)
  dirs <- stats::setNames(sel$quant_dir, sel$sample_id)
  list(col_data = sel[, c("sample_id", "condition")], dirs = dirs)
}

base_object <- function(s) {
  read_salmon(
    s$dirs,
    annot_path = Sys.getenv("ISOFORMIC_BENCHMARK_ANNOT"),
    col_data = s$col_data
  )
}

with_draws <- function(s) {
  add_uncertainty(
    base_object(s),
    paths = s$dirs,
    group = "condition",
    n_draws = n_draws
  )
}

# The first 3 samples of each condition, as x[, j] of the object with draws:
# its draws are a view of those of the object (a `selection` attribute)
contrast_subset <- function(s) {
  iso <- with_draws(s)
  cd <- as.data.frame(iso@col_data)
  j <- unlist(lapply(split(cd$sample_id, cd$condition), utils::head, n = 3L))
  sub <- iso[, unname(j)]
  rm(iso)
  sub
}

# Steps: `setup` builds the input, `run` is measured --------------------------

steps <- list(
  read_inferential_replicates = list(
    setup = function(s) {
      counts <- tx_assay(base_object(s), "counts")
      list(s = s, keep = rownames(counts)[rowMeans(counts) >= 1])
    },
    run = function(x) {
      read_inferential_replicates(
        x$s$dirs,
        transcripts = x$keep,
        n_draws = n_draws
      )
    },
    draws = function(x, res) res
  ),
  add_uncertainty = list(
    setup = function(s) list(s = s, iso = base_object(s)),
    run = function(x) {
      add_uncertainty(
        x$iso,
        paths = x$s$dirs,
        group = "condition",
        n_draws = n_draws
      )
    },
    draws = function(x, res) res@metadata$inf_reps
  ),
  productive_fraction = list(
    setup = function(s) with_draws(s),
    run = function(x) productive_fraction(x, group = "condition"),
    draws = function(x, res) x@metadata$inf_reps
  ),
  subset_samples = list(
    setup = function(s) {
      iso <- with_draws(s)
      cd <- as.data.frame(iso@col_data)
      j <- unlist(lapply(
        split(cd$sample_id, cd$condition),
        utils::head,
        n = 3L
      ))
      list(iso = iso, j = unname(j))
    },
    run = function(x) x$iso[, x$j],
    draws = function(x, res) x$iso@metadata$inf_reps
  ),
  write_isoformic = list(
    setup = function(s) with_draws(s),
    run = function(x) {
      write_isoformic(x, path = written_dir, overwrite = TRUE)
    },
    draws = function(x, res) x@metadata$inf_reps,
    disk = written_dir
  ),
  read_isoformic = list(
    setup = function(s) {
      if (!fs::dir_exists(written_dir)) {
        stop("Run the write_isoformic step first.", call. = FALSE)
      }
      written_dir
    },
    run = function(x) read_isoformic(x),
    draws = function(x, res) res@metadata$inf_reps
  ),
  inferential_overdispersion = list(
    setup = function(s) with_draws(s),
    run = function(x) inferential_overdispersion(x),
    draws = function(x, res) x@metadata$inf_reps
  ),
  subset_overdispersion = list(
    setup = function(s) contrast_subset(s),
    run = function(x) inferential_overdispersion(x),
    draws = function(x, res) x@metadata$inf_reps
  ),
  subset_write = list(
    setup = function(s) contrast_subset(s),
    run = function(x) {
      write_isoformic(x, path = written_subset_dir, overwrite = TRUE)
    },
    draws = function(x, res) x@metadata$inf_reps,
    disk = written_subset_dir
  ),
  prepare_annotation = list(
    setup = function(s) Sys.getenv("ISOFORMIC_BENCHMARK_GFF"),
    run = function(x) {
      fs::dir_create(prepared_dir)
      prepare_isoformic_annotation(
        x,
        output_path = prepared_dir,
        annot_name = "annot",
        overwrite = TRUE
      )
    },
    draws = function(x, res) NULL
  ),
  transcript_features = list(
    setup = function(s) Sys.getenv("ISOFORMIC_BENCHMARK_ANNOT"),
    run = function(x) transcript_features(x),
    draws = function(x, res) NULL
  ),
  splice_site_strength = list(
    setup = function(s) Sys.getenv("ISOFORMIC_BENCHMARK_ANNOT"),
    run = function(x) splice_site_strength(x),
    draws = function(x, res) NULL
  ),
  splice_site_strength_gene = list(
    setup = function(s) base_object(s),
    run = function(x) splice_site_strength(x, gene_name = "RBPJ"),
    draws = function(x, res) NULL
  ),
  poison_exons = list(
    setup = function(s) base_object(s),
    run = function(x) poison_exons(x),
    draws = function(x, res) NULL
  ),
  transcript_sequences = list(
    setup = function(s) {
      annot <- Sys.getenv("ISOFORMIC_BENCHMARK_ANNOT")
      genes <- read_parquet(
        fs::path(annot, "annot_data_genes.parquet"),
        columns = c("gene_id", "gene_name")
      )
      tx <- read_parquet(
        fs::path(annot, "annot_data_transcripts.parquet"),
        columns = c("transcript_id", "gene_id")
      )
      gene_ids <- genes$gene_id[genes$gene_name %in% benchmark_genes]
      list(annot = annot, ids = tx$transcript_id[tx$gene_id %in% gene_ids])
    },
    run = function(x) transcript_sequences(x$annot, transcript_id = x$ids),
    draws = function(x, res) NULL
  )
)

run_step <- function(name) {
  step <- steps[[name]]
  s <- load_samples()
  x <- step$setup(s)
  invisible(gc(full = TRUE))
  start_kb <- proc_status_kb("VmRSS")
  peak_reset <- reset_peak_rss()
  t0 <- proc.time()[["elapsed"]]
  res <- step$run(x)
  elapsed <- proc.time()[["elapsed"]] - t0
  peak_kb <- proc_status_kb("VmHWM")
  reps <- step$draws(x, res)
  # Size on disk of the draws a step wrote (computed outside tibble(), whose
  # `step` column would mask the step)
  disk_gb <- NA_real_
  if (!is.null(step$disk)) {
    inf_dir <- fs::path(step$disk, "inf_reps")
    disk_gb <- if (fs::dir_exists(inf_dir)) {
      sum(as.numeric(fs::dir_info(inf_dir, recurse = TRUE)$size)) / 1024^3
    } else {
      0
    }
  }
  row <- tibble::tibble(
    step = name,
    n_result_rows = if (is.data.frame(res)) nrow(res) else NA_integer_,
    n_samples = nrow(s$col_data),
    n_draws = if (is.null(reps)) NA_integer_ else length(reps),
    n_transcripts_draws = if (is.null(reps)) {
      NA_integer_
    } else {
      nrow(reps[[1L]])
    },
    draws_gb = draws_bytes(reps) / 1024^3,
    start_rss_gb = start_kb / 1024^2,
    peak_rss_gb = peak_kb / 1024^2,
    peak_above_start_gb = (peak_kb - start_kb) / 1024^2,
    elapsed_s = elapsed,
    peak_reset = peak_reset,
    draws_is_view = !is.null(attr(reps, "selection", exact = TRUE)),
    result_draws_is_view = S7::S7_inherits(res) &&
      !is.null(attr(res@metadata$inf_reps, "selection", exact = TRUE)),
    inf_reps_disk_gb = disk_gb
  )
  write_parquet(row, fs::path(out_dir, paste0("step_", name, ".parquet")))
  invisible(row)
}

# Driver -----------------------------------------------------------------------

step_name <- env_or("ISOFORMIC_BENCHMARK_STEP")
if (!is.null(step_name)) {
  run_step(step_name)
  quit(save = "no", status = 0L)
}

if (!identical(Sys.getenv("ISOFORMIC_FULL"), "true")) {
  stop(
    "This benchmark reads a genome-wide annotation; set ISOFORMIC_FULL=true.",
    call. = FALSE
  )
}
fs::dir_create(out_dir)

step_names <- names(steps)
selected <- env_or("ISOFORMIC_BENCHMARK_STEPS")
if (!is.null(selected)) {
  step_names <- trimws(strsplit(selected, ",", fixed = TRUE)[[1L]])
  unknown <- setdiff(step_names, names(steps))
  if (length(unknown)) {
    stop("Unknown steps: ", paste(unknown, collapse = ", "), call. = FALSE)
  }
}

annot_dir <- env_or("ISOFORMIC_BENCHMARK_ANNOT")
gff <- env_or("ISOFORMIC_BENCHMARK_GFF")
if (
  is.null(gff) && (is.null(annot_dir) || "prepare_annotation" %in% step_names)
) {
  release <- env_or("ISOFORMIC_BENCHMARK_GENCODE")
  gff <- download_reference(version = release, file_type = "gff")
}
if (is.null(annot_dir)) {
  annot_dir <- prepare_isoformic_annotation(gff)
}

n_samples <- as.integer(env_or("ISOFORMIC_BENCHMARK_SAMPLES", "12"))
sheet <- utils::read.delim(
  env_or("ISOFORMIC_BENCHMARK_SHEET", "data-raw/placenta/placenta_samples.tsv"),
  colClasses = "character"
)
sheet <- sheet[fs::dir_exists(fs::path(quants_dir, sheet$run)), ]
conditions <- unique(sheet$condition)
if (length(conditions) != 2L) {
  stop("The sample sheet needs exactly two conditions.", call. = FALSE)
}
per_group <- n_samples %/% 2L
picked <- unlist(lapply(conditions, function(cond) {
  utils::head(which(sheet$condition == cond), per_group)
}))
if (length(picked) < 2L * per_group) {
  stop("Not enough quantified samples in each condition.", call. = FALSE)
}
sel <- tibble::tibble(
  sample_id = sheet$run[picked],
  condition = sheet$condition[picked],
  quant_dir = as.character(fs::path_abs(fs::path(
    quants_dir,
    sheet$run[picked]
  )))
)
write_parquet(sel, samples_file)

by_condition <- table(sel$condition)
by_condition <- paste(by_condition, names(by_condition), collapse = ", ")
cli::cli_inform(c(
  i = "Annotation: {.path {annot_dir}}",
  i = "{nrow(sel)} samples ({by_condition}), {n_draws} draws"
))

script <- normalizePath(sub(
  "^--file=",
  "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[[1L]]
))
rscript <- file.path(R.home("bin"), "Rscript")
status <- integer()
for (name in step_names) {
  cli::cli_inform(c(i = "Step {.val {name}}"))
  status[[name]] <- system2(
    rscript,
    shQuote(script),
    env = c(
      paste0("ISOFORMIC_BENCHMARK_STEP=", name),
      paste0("ISOFORMIC_BENCHMARK_ANNOT=", shQuote(annot_dir)),
      paste0("ISOFORMIC_BENCHMARK_GFF=", shQuote(if (is.null(gff)) "" else gff))
    )
  )
}

rows <- lapply(step_names, function(name) {
  f <- fs::path(out_dir, paste0("step_", name, ".parquet"))
  if (status[[name]] == 0L && fs::file_exists(f)) {
    read_parquet(f)
  } else {
    tibble::tibble(step = name)
  }
})
res <- dplyr::bind_rows(rows)
res$exit_status <- unname(status[res$step])
one_copy <- res$draws_gb[res$step == "read_inferential_replicates"]
one_copy <- if (length(one_copy) && !is.na(one_copy)) {
  one_copy
} else {
  stats::median(res$draws_gb, na.rm = TRUE)
}
res$draws_copy_gb <- one_copy
res$peak_multiple_of_draws <- res$peak_above_start_gb / one_copy
res$git_head <- tryCatch(
  system2("git", c("rev-parse", "--short", "HEAD"), stdout = TRUE)[[1L]],
  error = function(cnd) NA_character_
)
write_parquet(res, fs::path(out_dir, "benchmark_memory.parquet"))
for (d in c(written_dir, written_subset_dir, prepared_dir)) {
  if (fs::dir_exists(d)) {
    fs::dir_delete(d)
  }
}

print(
  as.data.frame(res[, c(
    "step",
    "start_rss_gb",
    "peak_above_start_gb",
    "peak_multiple_of_draws",
    "elapsed_s",
    "exit_status"
  )]),
  digits = 3L,
  row.names = FALSE
)
n_tx <- res$n_transcripts_draws[[1L]]
cli::cli_inform(c(
  i = paste0(
    "One copy of the draws: {round(one_copy, 3)} GB ({n_tx} transcripts x ",
    "{nrow(sel)} samples x {n_draws} draws)."
  ),
  v = "Results in {.path {fs::path(out_dir, 'benchmark_memory.parquet')}}."
))
