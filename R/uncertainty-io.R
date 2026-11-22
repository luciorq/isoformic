#' Read the inferential uncertainty that Salmon writes
#'
#' @description
#' Salmon quantifies how well each transcript abundance is determined by the
#' reads, and writes that information next to `quant.sf`:
#'
#' * `read_inferential_replicates()` reads the Gibbs (`--numGibbsSamples`) or
#'   bootstrap (`--numBootstraps`) replicates from `aux_info/bootstrap/`
#'   (the format is identical in Salmon 1.x and in the Rust rewrite, Salmon
#'   2). It returns a list with one transcript x sample matrix of counts per
#'   replicate, the layout used by [productive_fraction()] and
#'   [add_uncertainty()]. The samples are read one at a time, so memory holds
#'   the result and the replicates of a single sample. When the samples have
#'   different numbers of replicates, the smallest number is kept. The mean
#'   and variance over the kept replicates of every transcript, including
#'   those left out by `transcripts`, are recorded in the `moments`
#'   attribute: [inferential_overdispersion()] estimates its prior from them,
#'   so it does not depend on which transcripts are kept.
#' * `equivalence_support()` reads the equivalence classes
#'   (`--dumpEq`, `aux_info/eq_classes.txt.gz`) and reports, for every
#'   transcript, how many of its fragments are compatible with no other
#'   transcript.
#' * `salmon_run_qc()` reads `aux_info/meta_info.json` and reports the
#'   mapping rate, the library type, whether the optimiser converged and how
#'   many replicates were drawn.
#' * `inferential_replicates_table()` and `as_inferential_replicates()`
#'   convert the list of replicates to and from a long table (one row per
#'   transcript, sample and replicate), which can be stored as Parquet with
#'   [write_parquet()] and read back with [read_parquet()], for example to
#'   exchange them with other tools. `inferential_moments_table()` converts
#'   the `moments` attribute to a table in the same way (one row per
#'   transcript and sample); pass it to the `moments` argument of
#'   `as_inferential_replicates()` to restore it. The long table takes
#'   several times the memory of the replicates: to save the replicates of
#'   an object, use [write_isoformic()], which writes one wide file per draw.
#'   Both apply the `selection` of the replicates of a subset `x[i, j]`
#'   (see [IsoformicExperiment]): the replicate table holds the kept
#'   transcripts and samples only, the moment table the kept samples (the
#'   moments keep every transcript). The selection is an attribute of the
#'   unmodified `x@metadata$inf_reps`: a list operation such as [lapply()]
#'   or `[` drops it, and the list then holds the draws of the parent. Pass
#'   the object as `x` (alone, or with such a list as `reps`) to read the
#'   draws of its transcripts and samples in every case.
#'
#' The functions also read Oarfish output (see [read_oarfish()]): its
#' bootstrap replicates (`.infreps.pq`), its unique and ambiguous read counts
#' (`.ambig_info.tsv`, `n_classes` is then `NA`) and its run metadata
#' (`.meta_info.json`). One call reads the output of one tool.
#'
#' @param paths Character vector of Salmon output directories (containing
#'   `quant.sf`) or paths to `quant.sf` files, or Oarfish output prefixes.
#'   Names, when present, are the sample IDs; otherwise the directory name (or
#'   the prefix) is used.
#' @param transcripts Optional transcript IDs to keep. Reading every
#'   transcript of a human annotation for many samples takes gigabytes of
#'   memory, so [add_uncertainty()] only reads the expressed ones.
#' @param n_draws Number of replicates to keep (the first ones). `NULL` keeps
#'   all of them; a number above the replicates available keeps all of them,
#'   with a message.
#' @param ignore_after_bar Logical. Keep only the part of transcript names
#'   before the first `|` (GENCODE FASTA headers used as Salmon references).
#' @param reps A list of transcript x sample count matrices, one per replicate.
#'   In `inferential_replicates_table()` and `inferential_moments_table()`,
#'   it defaults to `x@metadata$inf_reps` when `x` is given.
#' @param x Optional [IsoformicExperiment]. When given, the replicates are
#'   read as those of its transcripts and samples, matched by ID, even when
#'   `reps` lost the `selection` of a subset.
#' @param table A long table with `transcript_id`, `sample_id`, `draw` and
#'   `count` columns, as returned by `inferential_replicates_table()`.
#' @param samp_type Replicate type recorded on the result, `"gibbs"` or
#'   `"bootstrap"`.
#' @param moments Optional table of `inferential_moments_table()` (columns
#'   `transcript_id`, `sample_id`, `n_draws`, `mean` and `var`) to record as
#'   the `moments` attribute of the result.
#' @param min_percent_mapped Samples below this mapping rate (percent) get the
#'   `low_mapping_rate` flag in `salmon_run_qc()`.
#'
#' @returns
#' * `read_inferential_replicates()`: a list of numeric matrices (transcripts
#'   by samples, counts) with the replicate type (`"gibbs"` or `"bootstrap"`)
#'   in the `samp_type` attribute and the moments of every transcript in the
#'   `moments` attribute (a list with the `mean` and `var` matrices,
#'   transcripts by samples, and the number of replicates `n`).
#' * `inferential_replicates_table()`: a tibble with `transcript_id`,
#'   `sample_id`, `draw` and `count`. `as_inferential_replicates()` is its
#'   inverse. `inferential_moments_table()`: a tibble with `transcript_id`,
#'   `sample_id`, `n_draws`, `mean` and `var`.
#' * `equivalence_support()`: a tibble with one row per transcript:
#'   `transcript_id`, `n_classes` (mean number of equivalence classes that
#'   contain the transcript, per sample), `n_fragments` and
#'   `n_unique_fragments` (summed over samples), `unique_fraction` (share of
#'   the fragments of the transcript that map to it alone) and `n_samples`
#'   (samples with equivalence classes).
#' * `salmon_run_qc()`: a tibble with one row per sample: `sample_id`,
#'   `salmon_version`, `library_type`, `detected_library_type`,
#'   `num_processed`, `num_mapped`, `percent_mapped`, `replicate_type`,
#'   `num_replicates`, `optimiser`, `em_converged`, `seq_bias_correct`,
#'   `gc_bias_correct`, `pos_bias_correct`, `frag_length_mean`,
#'   `frag_length_sd` and `qc_flag` (`low_mapping_rate`, `no_replicates`
#'   and `em_not_converged`, separated by `;`, or `NA`). For Oarfish output
#'   the columns are `sample_id`, `num_processed`, `num_mapped`,
#'   `percent_mapped`, `replicate_type`, `num_replicates`, `optimiser`,
#'   `alignment_source`, `prob_model` and `qc_flag`, and the `tool`
#'   attribute is `"oarfish"`.
#'
#' @examples
#' quant_dirs <- list.files(
#'   system.file("extdata", "mini_quants", package = "isoformic"),
#'   full.names = TRUE
#' )
#' names(quant_dirs) <- sub("_quant$", "", basename(quant_dirs))
#' # The bundled quantifications were made without replicates
#' salmon_run_qc(quant_dirs)[, c("sample_id", "percent_mapped", "qc_flag")]
#'
#' @name inferential-io
#' @aliases read_inferential_replicates equivalence_support salmon_run_qc
NULL

#' @rdname inferential-io
#' @export
read_inferential_replicates <- function(
  paths,
  transcripts = NULL,
  n_draws = NULL,
  ignore_after_bar = TRUE
) {
  rlang::check_required(paths)
  dirs <- resolve_quant_dirs(paths)
  n_sample <- nrow(dirs)
  ids <- NULL
  keep <- NULL
  reps <- NULL
  samp_type <- NULL
  available <- NULL
  mom_mean <- NULL
  mom_var <- NULL
  mom_n <- integer(n_sample)
  collect <- gc_pacer()
  # One sample is in memory at a time: its kept transcripts and replicates go
  # into the matrices of the replicates, which are allocated once
  for (i in seq_len(n_sample)) {
    x <- read_sample_replicates(dirs, i)
    sample_ids <- replicate_transcript_ids(x$counts, ignore_after_bar)
    if (i == 1L) {
      ids <- sample_ids
      samp_type <- x$samp_type
      keep <- kept_replicate_rows(ids, transcripts)
      rows <- keep
      all_rows <- NULL
      available <- ncol(x$counts)
      n_kept <- if (is.null(n_draws)) available else min(available, n_draws)
      # Every replicate shares the row and column names
      shared_names <- list(ids[keep], dirs$sample_id)
      empty <- vector(typeof(x$counts), 1L)
      reps <- lapply(seq_len(n_kept), function(k) {
        matrix(
          empty,
          nrow = length(keep),
          ncol = n_sample,
          dimnames = shared_names
        )
      })
      # Two allocations: `mom_var <- mom_mean` would share the matrix, and
      # the first assignment into one of them would copy it
      mom_dimnames <- list(ids, dirs$sample_id)
      mom_mean <- matrix(
        0,
        nrow = length(ids),
        ncol = n_sample,
        dimnames = mom_dimnames
      )
      mom_var <- matrix(
        0,
        nrow = length(ids),
        ncol = n_sample,
        dimnames = mom_dimnames
      )
    } else if (identical(sample_ids, ids)) {
      # The usual case, the same index for every sample: no matching needed
      rows <- keep
      all_rows <- NULL
      available <- min(available, ncol(x$counts))
      if (available < length(reps)) {
        reps <- reps[seq_len(available)]
      }
    } else {
      if (!setequal(sample_ids, ids)) {
        cli::cli_abort(
          message = c(
            x = paste0(
              "Sample {.val {dirs$sample_id[i]}} was quantified against a ",
              "different set of transcripts than {.val {dirs$sample_id[1L]}}."
            ),
            i = "Quantify every sample with the same Salmon index."
          ),
          class = "isoformic_inconsistent_transcripts"
        )
      }
      all_rows <- match(ids, sample_ids)
      rows <- all_rows[keep]
      available <- min(available, ncol(x$counts))
      if (available < length(reps)) {
        reps <- reps[seq_len(available)]
      }
    }
    for (k in seq_along(reps)) {
      reps[[k]][, i] <- x$counts[rows, k]
    }
    # Moments of every transcript, kept or not, for the overdispersion prior
    mom <- replicate_moments(x$counts, length(reps))
    if (is.null(all_rows)) {
      mom_mean[, i] <- mom$mean
      mom_var[, i] <- mom$var
    } else {
      mom_mean[, i] <- mom$mean[all_rows]
      mom_var[, i] <- mom$var[all_rows]
    }
    mom_n[[i]] <- length(reps)
    cells <- length(x$counts)
    rm(x, mom)
    # The replicates of a sample (hundreds of MB for a human annotation) are
    # garbage now: collect them before the next sample is read, so that the
    # garbage of several samples does not pile up
    collect(cells)
  }
  # A later sample with fewer replicates: recompute the earlier samples on
  # the replicates that are kept
  for (i in which(mom_n > length(reps))) {
    x <- read_sample_replicates(dirs, i)
    sample_ids <- replicate_transcript_ids(x$counts, ignore_after_bar)
    mom <- replicate_moments(x$counts, length(reps))
    all_rows <- match(ids, sample_ids)
    mom_mean[, i] <- mom$mean[all_rows]
    mom_var[, i] <- mom$var[all_rows]
    cells <- length(x$counts)
    rm(x, mom)
    collect(cells)
  }
  if (!is.null(n_draws) && n_draws > available) {
    cli::cli_inform(c(
      i = "Only {available} replicate{?s} available per sample; using all of them."
    ))
  }
  attr(reps, "samp_type") <- samp_type
  attr(reps, "moments") <- list(
    mean = mom_mean,
    var = mom_var,
    n = length(reps)
  )
  reps
}

#' Replicates of sample `i` of the table of [resolve_quant_dirs()]
#' @keywords internal
#' @noRd
read_sample_replicates <- function(dirs, i) {
  if (identical(dirs$kind[[i]], "oarfish")) {
    read_oarfish_replicates(dirs$dir[[i]])
  } else {
    read_bootstrap_dir(dirs$dir[[i]])
  }
}

#' Transcript IDs of a matrix of replicates, cut at the first `|` if asked
#' @keywords internal
#' @noRd
replicate_transcript_ids <- function(counts, ignore_after_bar) {
  ids <- rownames(counts)
  if (isTRUE(ignore_after_bar)) {
    ids <- sub("\\|.*$", "", ids)
  }
  ids
}

#' Rows of `ids` to keep: all of them, or those of `transcripts`
#' @keywords internal
#' @noRd
kept_replicate_rows <- function(ids, transcripts) {
  if (is.null(transcripts)) {
    return(seq_along(ids))
  }
  keep <- match(unique(transcripts), ids)
  keep <- keep[!is.na(keep)]
  if (!length(keep)) {
    cli::cli_abort(
      message = c(
        x = "None of the requested transcripts are present in the Salmon replicates."
      ),
      class = "isoformic_no_annotated_transcripts"
    )
  }
  keep
}

#' Mean and variance of each row over the first `n` replicates of a sample
#'
#' Same sums as `slice_moments()`, so the values match those computed from
#' the kept replicates. `.rowSums()` sums the first `n` columns in place (the
#' matrix is stored column by column), with the arithmetic of `rowSums()`,
#' so the only copy is the vector of squares; `^ 2` squares a temporary in
#' place and gives the same values as `x * x`.
#' @keywords internal
#' @noRd
replicate_moments <- function(counts, n) {
  n_row <- nrow(counts)
  s1 <- .rowSums(counts, n_row, n)
  squares <- if (n < ncol(counts)) {
    counts[seq_len(as.double(n_row) * n)]^2
  } else {
    counts^2
  }
  s2 <- .rowSums(squares, n_row, n)
  rm(squares)
  mean_draw <- s1 / n
  var_draw <- if (n > 1L) {
    pmax(s2 - n * mean_draw^2, 0) / (n - 1L)
  } else {
    mean_draw * 0
  }
  list(mean = unname(mean_draw), var = unname(var_draw))
}

#' @rdname inferential-io
#' @export
equivalence_support <- function(
  paths,
  transcripts = NULL,
  ignore_after_bar = TRUE
) {
  rlang::check_required(paths)
  dirs <- resolve_quant_dirs(paths)
  total <- NULL
  unique_total <- NULL
  n_classes <- NULL
  n_samples <- 0L
  ids <- NULL
  oarfish <- identical(unique(dirs$kind), "oarfish")
  if (oarfish) {
    return(oarfish_support(dirs$dir, transcripts, ignore_after_bar))
  }
  # One connection for all samples; only per-transcript sums reach R
  con <- duckdb_connect()
  for (dir in dirs$dir) {
    sup <- eq_sample_support(con, dir)
    if (is.null(sup)) {
      next
    }
    tx <- sup$transcripts
    if (isTRUE(ignore_after_bar)) {
      tx <- sub("\\|.*$", "", tx)
    }
    if (is.null(ids)) {
      ids <- tx
      total <- numeric(length(ids))
      unique_total <- numeric(length(ids))
      n_classes <- numeric(length(ids))
    }
    tid <- match(tx, ids)[sup$sums$tid]
    keep <- !is.na(tid)
    sums <- sup$sums[keep, , drop = FALSE]
    tid <- tid[keep]
    total <- total + sum_by_index(sums$total, tid, length(ids))
    unique_total <- unique_total +
      sum_by_index(sums$unique_total, tid, length(ids))
    n_classes <- n_classes + sum_by_index(sums$n_classes, tid, length(ids))
    n_samples <- n_samples + 1L
  }
  if (n_samples == 0L) {
    cli::cli_abort(
      message = c(
        x = "No equivalence classes found in the Salmon output.",
        i = "Run {.code salmon quant --dumpEq} to write {.file aux_info/eq_classes.txt.gz}."
      ),
      class = "isoformic_no_equivalence_classes"
    )
  }
  res <- tibble::tibble(
    transcript_id = ids,
    n_classes = n_classes / n_samples,
    n_fragments = total,
    n_unique_fragments = unique_total,
    unique_fraction = ifelse(total > 0, unique_total / total, NA_real_),
    n_samples = n_samples
  )
  if (!is.null(transcripts)) {
    res <- res[res$transcript_id %in% transcripts, ]
  }
  res
}

#' @rdname inferential-io
#' @export
salmon_run_qc <- function(paths, min_percent_mapped = 50) {
  rlang::check_required(paths)
  dirs <- resolve_quant_dirs(paths)
  if (identical(unique(dirs$kind), "oarfish")) {
    return(oarfish_run_qc(dirs, min_percent_mapped))
  }
  rows <- lapply(dirs$dir, function(dir) {
    info <- read_meta_info(dir)
    tibble::tibble(
      salmon_version = info$salmon_version,
      library_type = info$library_type,
      detected_library_type = info$detected_library_type,
      num_processed = info$num_processed,
      num_mapped = info$num_mapped,
      percent_mapped = info$percent_mapped,
      replicate_type = info$samp_type,
      num_replicates = info$num_bootstraps,
      optimiser = info$opt_type,
      em_converged = info$em_converged,
      seq_bias_correct = info$seq_bias_correct,
      gc_bias_correct = info$gc_bias_correct,
      pos_bias_correct = info$pos_bias_correct,
      frag_length_mean = info$frag_length_mean,
      frag_length_sd = info$frag_length_sd
    )
  })
  res <- dplyr::bind_rows(rows)
  res <- dplyr::bind_cols(tibble::tibble(sample_id = dirs$sample_id), res)
  res$qc_flag <- vapply(
    seq_len(nrow(res)),
    function(i) {
      flags <- character()
      if (isTRUE(res$percent_mapped[i] < min_percent_mapped)) {
        flags <- c(flags, "low_mapping_rate")
      }
      if (isTRUE(res$num_replicates[i] < 1)) {
        flags <- c(flags, "no_replicates")
      }
      if (isFALSE(res$em_converged[i])) {
        flags <- c(flags, "em_not_converged")
      }
      if (length(flags)) paste(flags, collapse = ";") else NA_character_
    },
    character(1L)
  )
  res
}

# Helpers ----------------------------------------------------------------------

#' Resolve Salmon output directories and sample IDs
#' @keywords internal
#' @noRd
resolve_quant_dirs <- function(paths, arg = rlang::caller_arg(paths)) {
  kinds <- vapply(paths, quant_kind, character(1L))
  if (length(unique(kinds)) > 1L) {
    cli::cli_abort(
      message = c(
        x = "{.arg {arg}} mixes quantifiers: {.val {unique(kinds)}}.",
        i = "Read the output of one tool at a time."
      ),
      class = "isoformic_mixed_quantifiers"
    )
  }
  dirs <- vapply(
    seq_along(paths),
    function(i) {
      p <- as.character(paths[[i]])
      if (identical(kinds[[i]], "oarfish")) {
        return(oarfish_prefix(p))
      }
      if (isTRUE(fs::is_file(p))) {
        p <- as.character(fs::path_dir(p))
      }
      assert_dir_exists(p, arg = arg)
      as.character(fs::path_abs(p))
    },
    character(1L)
  )
  ids <- names(paths)
  if (is.null(ids) || !all(nzchar(ids))) {
    ids <- as.character(fs::path_file(dirs))
  }
  if (anyDuplicated(ids)) {
    cli::cli_abort(
      message = c(x = "Sample IDs derived from {.arg {arg}} are not unique."),
      class = "isoformic_duplicated_sample_ids"
    )
  }
  data.frame(
    sample_id = ids,
    dir = unname(dirs),
    kind = unname(kinds),
    stringsAsFactors = FALSE
  )
}

#' Which tool wrote a quantification: a Salmon directory or Oarfish files
#'
#' Oarfish writes `<prefix>.quant` next to `<prefix>.infreps.pq`; a path is
#' Oarfish output when it is that prefix or the `.quant` file itself.
#' @keywords internal
#' @noRd
quant_kind <- function(path) {
  path <- as.character(path)
  if (grepl("\\.quant$", path) && isTRUE(fs::is_file(path))) {
    return("oarfish")
  }
  if (
    !isTRUE(fs::is_dir(path)) && isTRUE(fs::is_file(paste0(path, ".quant")))
  ) {
    return("oarfish")
  }
  if (
    isTRUE(fs::is_dir(path)) &&
      !isTRUE(fs::is_file(fs::path(path, "quant.sf"))) &&
      length(fs::dir_ls(path, glob = "*.quant", type = "file")) == 1L
  ) {
    return("oarfish")
  }
  "salmon"
}

#' Extract one field of a (flat) JSON text without a JSON parser
#' @keywords internal
#' @noRd
json_field <- function(text, key) {
  pattern <- paste0('"', key, '"\\s*:\\s*(\\[[^\\]]*\\]|"[^"]*"|[^,}\\s]+)')
  hit <- regmatches(text, regexec(pattern, text, perl = TRUE))[[1L]]
  if (length(hit) < 2L) {
    return(NA_character_)
  }
  value <- hit[2L]
  if (identical(value, "null")) {
    return(NA_character_)
  }
  value <- gsub('^"|"$', "", value)
  if (startsWith(value, "[")) {
    value <- gsub('[][" \n\r\t]', "", value)
  }
  value
}

#' Read the fields of `aux_info/meta_info.json` used by the package
#' @keywords internal
#' @noRd
read_meta_info <- function(dir) {
  file <- fs::path(dir, "aux_info", "meta_info.json")
  if (!isTRUE(fs::file_exists(file))) {
    cli::cli_abort(
      message = c(
        x = "{.file {file}} not found.",
        i = "Is {.path {dir}} a Salmon output directory?"
      ),
      class = "isoformic_no_meta_info"
    )
  }
  text <- paste(readLines(file, warn = FALSE), collapse = "\n")
  num <- function(key) suppressWarnings(as.numeric(json_field(text, key)))
  field_chr <- function(key) json_field(text, key)
  field_lgl <- function(key) {
    value <- json_field(text, key)
    if (is.na(value)) NA else identical(value, "true")
  }
  list(
    salmon_version = field_chr("salmon_version"),
    samp_type = field_chr("samp_type"),
    opt_type = field_chr("opt_type"),
    library_type = field_chr("library_types"),
    detected_library_type = field_chr("detected_library_type"),
    num_processed = num("num_processed"),
    num_mapped = num("num_mapped"),
    percent_mapped = num("percent_mapped"),
    num_bootstraps = num("num_bootstraps"),
    num_valid_targets = num("num_valid_targets"),
    em_converged = field_lgl("em_converged"),
    seq_bias_correct = field_lgl("seq_bias_corrected"),
    gc_bias_correct = field_lgl("gc_bias_corrected"),
    pos_bias_correct = field_lgl("pos_bias_corrected"),
    frag_length_mean = num("frag_length_mean"),
    frag_length_sd = num("frag_length_sd")
  )
}

#' Read the replicates of one Salmon output directory
#'
#' `bootstraps.gz` holds raw little-endian doubles: `num_bootstraps` replicates
#' written one after the other, each with the abundance of every transcript in
#' the order of `names.tsv.gz`.
#' @keywords internal
#' @noRd
read_bootstrap_dir <- function(dir) {
  boot_dir <- fs::path(dir, "aux_info", "bootstrap")
  boot_file <- fs::path(boot_dir, "bootstraps.gz")
  names_file <- fs::path(boot_dir, "names.tsv.gz")
  if (!all(fs::file_exists(c(boot_file, names_file)))) {
    cli::cli_abort(
      message = c(
        x = "No inferential replicates found in {.path {dir}}.",
        i = paste0(
          "Quantify with {.code --numGibbsSamples 30} (or ",
          "{.code --numBootstraps}); {.fn salmon_quant} does it by default."
        )
      ),
      class = "isoformic_no_inferential_replicates"
    )
  }
  info <- read_meta_info(dir)
  tx <- readLines(gzfile(names_file), warn = FALSE)
  tx <- unlist(strsplit(paste(tx, collapse = "\t"), "\t", fixed = TRUE))
  tx <- tx[nzchar(tx)]
  n_tx <- length(tx)
  n_boot <- info$num_bootstraps
  if (is.na(n_boot) || n_boot < 1) {
    cli::cli_abort(
      message = c(x = "{.path {dir}} records no inferential replicates."),
      class = "isoformic_no_inferential_replicates"
    )
  }
  con <- gzfile(boot_file, "rb")
  on.exit(close(con), add = TRUE)
  values <- readBin(con, what = "double", n = n_tx * n_boot, endian = "little")
  if (length(values) != n_tx * n_boot) {
    cli::cli_abort(
      message = c(
        x = "{.file {boot_file}} is truncated: expected {n_tx * n_boot} values, read {length(values)}."
      ),
      class = "isoformic_invalid_inferential_replicates"
    )
  }
  # In place: `matrix()` would hold a second copy of the replicates
  dim(values) <- c(n_tx, n_boot)
  dimnames(values) <- list(tx, NULL)
  list(counts = values, samp_type = info$samp_type)
}

#' The equivalence class file of one Salmon output directory
#' @returns The path, or `NULL` when the directory has none.
#' @keywords internal
#' @noRd
eq_classes_file <- function(dir) {
  candidates <- fs::path(
    dir,
    "aux_info",
    c("eq_classes.txt.gz", "eq_classes.txt")
  )
  file <- candidates[fs::file_exists(candidates)]
  if (!length(file)) {
    return(NULL)
  }
  as.character(file[[1L]])
}

#' Open an equivalence class file and read its header
#'
#' The text format is the number of transcripts, the number of classes, the
#' transcript names and one line per class: the class size, the (0-based)
#' transcript indices, optionally their weights (`--dumpEqWeights`), and
#' the fragment count.
#' @returns A list with the open connection `con` (positioned at the first
#'   class line; the caller closes it), `transcripts` and `n_classes`.
#' @keywords internal
#' @noRd
open_eq_classes <- function(file) {
  con <- if (grepl("\\.gz$", file)) gzfile(file, "rt") else file(file, "rt")
  n_tx <- as.integer(readLines(con, n = 1L))
  n_classes <- as.integer(readLines(con, n = 1L))
  tx <- readLines(con, n = n_tx)
  list(con = con, transcripts = tx, n_classes = n_classes)
}

#' Per-transcript fragment sums of the equivalence classes of one sample
#'
#' The class lines are parsed and summed in DuckDB, so only one row per
#' transcript that appears in a class reaches R. The last field of a line is
#' the count whether or not the file has weights.
#'
#' @param con A DuckDB connection.
#' @returns `NULL` without equivalence classes, or a list with `transcripts`
#'   (the names of the header) and `sums`, a data frame with the 1-based
#'   transcript index `tid`, `total` (fragments of all its classes),
#'   `unique_total` (fragments of its classes of size 1) and `n_classes`.
#' @keywords internal
#' @noRd
eq_sample_support <- function(con, dir) {
  file <- eq_classes_file(dir)
  if (is.null(file)) {
    return(NULL)
  }
  header <- open_eq_classes(file)
  close(header$con)
  n_tx <- length(header$transcripts)
  sums <- data.frame(
    tid = integer(),
    total = numeric(),
    unique_total = numeric(),
    n_classes = numeric()
  )
  if (header$n_classes > 0L) {
    # One VARCHAR column per line: the delimiter never occurs in the class
    # lines, which hold numbers and tabs only
    sql <- glue::glue_sql(
      "WITH lines AS (
         SELECT string_split(rtrim(line, chr(9)), chr(9)) AS f
         FROM read_csv(
           {as.character(fs::path_abs(file))},
           header = false, auto_detect = false,
           columns = {{'line': 'VARCHAR'}},
           delim = '|', quote = '', escape = '',
           skip = {n_tx + 2L}, max_line_size = 8388608
         )
         LIMIT {header$n_classes}
       ),
       parsed AS (
         SELECT TRY_CAST(f[1] AS INTEGER) AS size, f
         FROM lines
       ),
       classes AS (
         SELECT size,
           CASE WHEN size IS NULL OR size < 0 OR len(f) < size + 2
             THEN error('isoformic: class line with fewer fields than its size needs')
             ELSE CAST(f[len(f)] AS DOUBLE)
           END AS count, f
         FROM parsed
       ),
       members AS (
         SELECT size, count,
           CAST(unnest(list_slice(f, 2, 1 + size)) AS INTEGER) + 1 AS tid
         FROM classes
       )
       SELECT tid,
         sum(count) AS total,
         coalesce(sum(count) FILTER (WHERE size = 1), 0) AS unique_total,
         CAST(count(*) AS DOUBLE) AS n_classes
       FROM members
       GROUP BY tid
       ORDER BY tid",
      .con = con
    )
    sums <- tryCatch(
      DBI::dbGetQuery(con, sql),
      error = function(e) {
        cli::cli_abort(
          message = c(
            x = "Could not read the class lines of {.file {file}}.",
            i = "A class line is malformed or has fewer fields than its size needs."
          ),
          class = "isoformic_invalid_equivalence_classes",
          parent = e
        )
      }
    )
    sums$tid <- as.integer(sums$tid)
  }
  list(transcripts = header$transcripts, sums = sums)
}

#' Read the members of the equivalence classes of one Salmon output directory
#'
#' The class lines are parsed as numbers with [scan()]; the number of fields
#' of each line gives the bounds of the classes, so files with and without
#' weights (and a mix of both) are read alike.
#'
#' @returns `NULL` without equivalence classes, or a list with `transcripts`,
#'   `members` (a data frame with the 1-based transcript index `tid` and the
#'   class index `class` of every class member) and `classes` (a data frame
#'   with the `size` and fragment `count` of every class).
#' @keywords internal
#' @noRd
read_eq_class_members <- function(dir) {
  file <- eq_classes_file(dir)
  if (is.null(file)) {
    return(NULL)
  }
  header <- open_eq_classes(file)
  on.exit(close(header$con), add = TRUE)
  n_classes <- header$n_classes
  n_tx <- length(header$transcripts)
  members <- data.frame(tid = integer(), class = integer())
  classes <- data.frame(size = integer(), count = numeric())
  if (n_classes > 0L) {
    # Fields per line, read in chunks so that the lines of the whole file
    # are never held at once
    n_fields <- integer(n_classes)
    done <- 0L
    while (done < n_classes) {
      lines <- readLines(header$con, n = min(100000L, n_classes - done))
      if (!length(lines)) {
        break
      }
      # A trailing tab adds no field: scan() returns no token for it
      lines <- sub("\t+$", "", lines)
      n_fields[done + seq_along(lines)] <- nchar(lines, type = "bytes") -
        nchar(gsub("\t", "", lines, fixed = TRUE), type = "bytes") +
        1L
      done <- done + length(lines)
    }
    rm(lines)
    close(header$con)
    on.exit()
    n_values <- sum(as.numeric(n_fields))
    con <- if (grepl("\\.gz$", file)) gzfile(file, "rt") else file(file, "rt")
    on.exit(close(con), add = TRUE)
    values <- if (done == n_classes) {
      scan(con, what = double(), n = n_values, skip = n_tx + 2L, quiet = TRUE)
    } else {
      numeric()
    }
    if (done < n_classes || length(values) != n_values) {
      cli::cli_abort(
        message = c(
          x = "{.file {file}} is truncated: expected {n_classes} class line{?s} and {n_values} value{?s}."
        ),
        class = "isoformic_invalid_equivalence_classes"
      )
    }
    starts <- cumsum(c(1, n_fields))[seq_len(n_classes)]
    size <- as.integer(values[starts])
    short <- which(is.na(size) | size < 0L | n_fields < size + 2L)
    if (length(short)) {
      cli::cli_abort(
        message = c(
          x = "Class line {short[[1L]]} of {.file {file}} has fewer fields than its size needs."
        ),
        class = "isoformic_invalid_equivalence_classes"
      )
    }
    count <- values[starts + n_fields - 1]
    # The members follow the size field of each class
    tid <- as.integer(values[sequence(size, from = starts + 1)]) + 1L
    rm(values)
    cls <- rep.int(seq_len(n_classes), size)
    members <- data.frame(tid = tid, class = cls)
    classes <- data.frame(size = size, count = count)
  }
  list(transcripts = header$transcripts, members = members, classes = classes)
}

#' Read the equivalence classes of one Salmon output directory
#'
#' One row per class member, with the size and count of its class repeated.
#' The package itself reads `read_eq_class_members()`, which does not
#' repeat them.
#' @keywords internal
#' @noRd
read_eq_classes <- function(dir) {
  ec <- read_eq_class_members(dir)
  if (is.null(ec)) {
    return(NULL)
  }
  cls <- ec$members$class
  list(
    transcripts = ec$transcripts,
    members = data.frame(
      tid = ec$members$tid,
      size = ec$classes$size[cls],
      count = ec$classes$count[cls],
      class = cls
    )
  )
}

#' Sum values by a positive integer index of known length
#' @keywords internal
#' @noRd
sum_by_index <- function(x, index, n) {
  out <- numeric(n)
  if (length(index)) {
    sums <- rowsum(x, index, reorder = TRUE)
    out[as.integer(rownames(sums))] <- sums[, 1L]
  }
  out
}

#' @rdname inferential-io
#' @export
inferential_replicates_table <- function(reps = NULL, x = NULL) {
  reps <- table_inf_reps(reps, x)
  if (!is.list(reps) || !length(reps)) {
    cli::cli_abort(
      message = c(x = "{.arg reps} must be a non-empty list of matrices."),
      class = "isoformic_invalid_inferential_replicates"
    )
  }
  # The selection of the replicates of a subset is applied to each draw
  read_draw <- inf_reps_reader(reps)
  first <- read_draw(1L)
  counts <- vector("list", length(reps))
  same <- logical(length(reps))
  for (k in seq_along(reps)) {
    m <- if (k == 1L) first else read_draw(k)
    same[[k]] <- identical(dimnames(m), dimnames(first))
    counts[[k]] <- as.vector(m)
  }
  if (!all(same) || is.null(rownames(first)) || is.null(colnames(first))) {
    cli::cli_abort(
      message = c(
        x = "All replicates need the same named transcripts and samples."
      ),
      class = "isoformic_invalid_inferential_replicates"
    )
  }
  n_tx <- nrow(first)
  n_sample <- ncol(first)
  tibble::tibble(
    transcript_id = rep(rownames(first), times = n_sample * length(reps)),
    sample_id = rep(
      rep(colnames(first), each = n_tx),
      times = length(reps)
    ),
    draw = rep(seq_along(reps), each = n_tx * n_sample),
    count = unlist(counts, use.names = FALSE)
  )
}

#' @rdname inferential-io
#' @export
inferential_moments_table <- function(reps = NULL, x = NULL) {
  mom <- inf_reps_moments(table_inf_reps(reps, x))
  if (is.null(mom)) {
    cli::cli_abort(
      message = c(
        x = "{.arg reps} has no moments of the replicates.",
        i = "They are recorded by {.fn read_inferential_replicates}."
      ),
      class = "isoformic_no_replicate_moments"
    )
  }
  n_tx <- nrow(mom$mean)
  n_sample <- ncol(mom$mean)
  tibble::tibble(
    transcript_id = rep(rownames(mom$mean), times = n_sample),
    sample_id = rep(colnames(mom$mean), each = n_tx),
    n_draws = rep(as.integer(mom$n), n_tx * n_sample),
    mean = as.vector(mom$mean),
    var = as.vector(mom$var)
  )
}

#' The replicates given to the table helpers, read as those of `x`
#'
#' Without `x`, `reps` as it is. With `x` (an IsoformicExperiment), `reps`
#' defaults to its stored replicates, and the selection of its transcripts
#' and samples is rebuilt when a list operation dropped it.
#' @keywords internal
#' @noRd
table_inf_reps <- function(reps, x, call = rlang::caller_env()) {
  if (is.null(x)) {
    return(reps)
  }
  if (!S7::S7_inherits(x, IsoformicExperiment)) {
    cli::cli_abort(
      message = c(
        x = "{.arg x} must be an {.cls IsoformicExperiment} or {.code NULL}."
      ),
      class = "isoformic_invalid_inferential_replicates",
      call = call
    )
  }
  object_inf_reps(x, reps %||% x@metadata$inf_reps)
}

#' Moments of the replicates from a table of `inferential_moments_table()`
#' @keywords internal
#' @noRd
moments_from_table <- function(moments) {
  required <- c("transcript_id", "sample_id", "n_draws", "mean", "var")
  if (!is.data.frame(moments) || !all(required %in% colnames(moments))) {
    cli::cli_abort(
      message = c(
        x = "{.arg moments} needs the columns {.field {required}}."
      ),
      class = "isoformic_invalid_inferential_replicates"
    )
  }
  n_draws <- unique(moments$n_draws)
  ids <- unique(as.character(moments$transcript_id))
  samples <- unique(as.character(moments$sample_id))
  if (
    length(n_draws) != 1L ||
      nrow(moments) != length(ids) * length(samples) ||
      anyNA(moments$mean) ||
      anyNA(moments$var)
  ) {
    cli::cli_abort(
      message = c(
        x = paste0(
          "{.arg moments} must have one mean and variance for every ",
          "transcript and sample, over one number of replicates."
        )
      ),
      class = "isoformic_invalid_inferential_replicates"
    )
  }
  slot <- match(moments$transcript_id, ids) +
    length(ids) * (match(moments$sample_id, samples) - 1)
  # With exactly one row per cell expected, a repeated cell means another
  # cell is missing
  dup <- anyDuplicated(slot)
  if (dup > 0L) {
    abort_duplicated_cell(
      "moments",
      moments$transcript_id[[dup]],
      moments$sample_id[[dup]]
    )
  }
  mean_draw <- matrix(0, length(ids), length(samples))
  var_draw <- mean_draw
  mean_draw[slot] <- moments$mean
  var_draw[slot] <- moments$var
  dimnames(mean_draw) <- dimnames(var_draw) <- list(ids, samples)
  list(mean = mean_draw, var = var_draw, n = as.integer(n_draws))
}

#' The layout of a replicate table in the order of
#' `inferential_replicates_table()`
#'
#' Transcripts fastest, then samples, then draws in increasing order. The
#' transcripts, samples and draws are read from the first rows of each block
#' and the order is checked one draw at a time, so no index or hash over the
#' whole table is built.
#'
#' @returns `NULL` when the table is not in that order, or a list with the
#'   `ids`, `samples` and `draws`, which are then those of [unique()] and
#'   [sort()] over the whole columns.
#' @keywords internal
#' @noRd
canonical_replicates_layout <- function(table) {
  n <- nrow(table)
  if (!n) {
    return(NULL)
  }
  tx_col <- table$transcript_id
  sample_col <- table$sample_id
  draw_col <- table$draw
  per_draw <- first_change(draw_col, 1, n) - 1
  if (per_draw < 1) {
    return(NULL)
  }
  n_tx <- first_change(sample_col, 1, per_draw) - 1
  if (n_tx < 1 || n %% per_draw != 0 || per_draw %% n_tx != 0) {
    return(NULL)
  }
  ids <- as.character(tx_col[seq_len(n_tx)])
  samples <- as.character(sample_col[seq(1, per_draw, by = n_tx)])
  draws <- draw_col[seq(1, n, by = per_draw)]
  distinct <- !anyNA(ids) &&
    !anyNA(samples) &&
    !anyNA(draws) &&
    anyDuplicated(ids) == 0L &&
    anyDuplicated(samples) == 0L &&
    !is.unsorted(draws, strictly = TRUE)
  if (!distinct || !draw_blocks_in_order(table, ids, samples, draws)) {
    return(NULL)
  }
  list(ids = ids, samples = samples, draws = draws)
}

#' Does every draw block repeat the transcripts and samples of the first?
#' @keywords internal
#' @noRd
draw_blocks_in_order <- function(table, ids, samples, draws) {
  per_draw <- length(ids) * length(samples)
  tx_order <- rep(ids, times = length(samples))
  sample_order <- rep(samples, each = length(ids))
  for (k in seq_along(draws)) {
    rows <- (k - 1) * per_draw + seq_len(per_draw)
    same <- isTRUE(all(table$draw[rows] == draws[[k]])) &&
      isTRUE(all(as.character(table$transcript_id[rows]) == tx_order)) &&
      isTRUE(all(as.character(table$sample_id[rows]) == sample_order))
    if (!same) {
      return(FALSE)
    }
  }
  TRUE
}

#' Position of the first element of `x[from:to]` that differs from
#' `x[from]`, or `to + 1` when none does; scanned in chunks
#' @keywords internal
#' @noRd
first_change <- function(x, from, to, chunk = 1e6) {
  value <- x[[from]]
  start <- from
  while (start <= to) {
    end <- min(to, start + chunk - 1)
    differs <- which(x[start:end] != value | is.na(x[start:end]))
    if (length(differs)) {
      return(start + differs[[1L]] - 1)
    }
    start <- end + 1
  }
  to + 1
}

#' Abort on a transcript and sample (and draw) given twice
#' @keywords internal
#' @noRd
abort_duplicated_cell <- function(arg, transcript, sample, draw = NULL) {
  cli::cli_abort(
    message = c(
      x = "{.arg {arg}} must have one row for every transcript, sample{if (is.null(draw)) '' else ' and draw'}.",
      i = paste0(
        "Transcript {.val {transcript}}, sample {.val {sample}}",
        if (!is.null(draw)) ", draw {.val {draw}}",
        " has more than one row, so another is missing."
      )
    ),
    class = "isoformic_invalid_inferential_replicates"
  )
}

#' @rdname inferential-io
#' @export
as_inferential_replicates <- function(
  table,
  samp_type = "gibbs",
  moments = NULL
) {
  required <- c("transcript_id", "sample_id", "draw", "count")
  if (!all(required %in% colnames(table))) {
    cli::cli_abort(
      message = c(
        x = "{.arg table} needs the columns {.field {required}}."
      ),
      class = "isoformic_invalid_inferential_replicates"
    )
  }
  layout <- canonical_replicates_layout(table)
  canonical <- !is.null(layout)
  if (!canonical) {
    layout <- list(
      ids = unique(as.character(table$transcript_id)),
      samples = unique(as.character(table$sample_id)),
      draws = sort(unique(table$draw))
    )
  }
  ids <- layout$ids
  samples <- layout$samples
  draws <- layout$draws
  n_tx <- length(ids)
  per_draw <- as.numeric(n_tx) * length(samples)
  n_cells <- per_draw * length(draws)
  if (nrow(table) != n_cells || anyNA(table$count)) {
    cli::cli_abort(
      message = c(
        x = "{.arg table} must have one count for every transcript, sample and draw."
      ),
      class = "isoformic_invalid_inferential_replicates"
    )
  }
  dims <- list(ids, samples)
  reps <- if (canonical) {
    # Written by inferential_replicates_table(): each draw is a block of
    # rows in matrix order, so no index over the whole table is needed
    lapply(seq_along(draws), function(k) {
      rows <- (k - 1) * per_draw + seq_len(per_draw)
      matrix(as.double(table$count[rows]), nrow = n_tx, dimnames = dims)
    })
  } else {
    # Integer arithmetic while the cells fit, as before
    n_tx_slot <- if (n_cells > .Machine$integer.max) as.numeric(n_tx) else n_tx
    slot <- match(table$transcript_id, ids) +
      n_tx_slot *
        ((match(table$sample_id, samples) - 1L) +
          length(samples) * (match(table$draw, draws) - 1L))
    if (anyNA(slot)) {
      cli::cli_abort(
        message = c(
          x = "{.arg table} must have no missing transcript, sample or draw."
        ),
        class = "isoformic_invalid_inferential_replicates"
      )
    }
    # With exactly one row per cell expected, a repeated cell means another
    # cell is missing
    dup <- if (n_cells <= .Machine$integer.max) {
      repeated <- which(tabulate(slot, n_cells) > 1L)
      if (length(repeated)) match(repeated[[1L]], slot) else 0L
    } else {
      anyDuplicated(slot)
    }
    if (dup > 0L) {
      abort_duplicated_cell(
        "table",
        table$transcript_id[[dup]],
        table$sample_id[[dup]],
        table$draw[[dup]]
      )
    }
    # The row of every cell, in matrix order: slot is a permutation
    row_of_cell <- integer(length(slot))
    row_of_cell[slot] <- seq_along(slot)
    rm(slot)
    lapply(seq_along(draws), function(k) {
      rows <- row_of_cell[(k - 1) * per_draw + seq_len(per_draw)]
      matrix(as.double(table$count[rows]), nrow = n_tx, dimnames = dims)
    })
  }
  attr(reps, "samp_type") <- samp_type
  if (!is.null(moments)) {
    mom <- moments_from_table(moments)
    if (
      mom$n != length(reps) ||
        !setequal(colnames(mom$mean), samples) ||
        !all(ids %in% rownames(mom$mean))
    ) {
      cli::cli_abort(
        message = c(
          x = paste0(
            "{.arg moments} must cover the transcripts, samples and number ",
            "of replicates of {.arg table}."
          )
        ),
        class = "isoformic_invalid_inferential_replicates"
      )
    }
    attr(reps, "moments") <- mom
  }
  reps
}
