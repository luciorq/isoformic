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
#'   [add_uncertainty()].
#' * `equivalence_support()` reads the equivalence classes
#'   (`--dumpEq`, `aux_info/eq_classes.txt.gz`) and reports, for every
#'   transcript, how many of its fragments are compatible with no other
#'   transcript.
#' * `salmon_run_qc()` reads `aux_info/meta_info.json` and reports the
#'   mapping rate, the library type, whether the optimiser converged and how
#'   many replicates were drawn.
#' * `inferential_replicates_table()` and `as_inferential_replicates()`
#'   convert the list of replicates to and from a long table (one row per
#'   transcript, sample and replicate), which can be stored as Parquet.
#'
#' @param paths Character vector of Salmon output directories (containing
#'   `quant.sf`) or paths to `quant.sf` files. Names, when present, are the
#'   sample IDs; otherwise the directory name is used.
#' @param transcripts Optional transcript IDs to keep. Reading every
#'   transcript of a human annotation for many samples takes gigabytes of
#'   memory, so [add_uncertainty()] only reads the expressed ones.
#' @param n_draws Number of replicates to keep (the first ones). `NULL` keeps
#'   all of them.
#' @param ignore_after_bar Logical. Keep only the part of transcript names
#'   before the first `|` (GENCODE FASTA headers used as Salmon references).
#' @param reps A list of transcript x sample count matrices, one per replicate.
#' @param table A long table with `transcript_id`, `sample_id`, `draw` and
#'   `count` columns, as returned by `inferential_replicates_table()`.
#' @param samp_type Replicate type recorded on the result, `"gibbs"` or
#'   `"bootstrap"`.
#' @param min_percent_mapped Samples below this mapping rate (percent) get the
#'   `low_mapping_rate` flag in `salmon_run_qc()`.
#'
#' @returns
#' * `read_inferential_replicates()`: a list of numeric matrices (transcripts
#'   by samples, counts) with the replicate type (`"gibbs"` or `"bootstrap"`)
#'   in the `samp_type` attribute.
#' * `inferential_replicates_table()`: a tibble with `transcript_id`,
#'   `sample_id`, `draw` and `count`. `as_inferential_replicates()` is its
#'   inverse.
#' * `equivalence_support()`: a tibble with one row per transcript:
#'   `transcript_id`, `n_classes` (mean number of equivalence classes per
#'   sample), `n_fragments` and `n_unique_fragments` (summed over samples),
#'   `unique_fraction` (share of the fragments of the transcript that map to
#'   it alone) and `n_samples`.
#' * `salmon_run_qc()`: a tibble with one row per sample.
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
  # One sample is in memory at a time: its kept transcripts and replicates go
  # into the matrices of the replicates, which are allocated once
  for (i in seq_len(n_sample)) {
    x <- if (identical(dirs$kind[[i]], "oarfish")) {
      read_oarfish_replicates(dirs$dir[[i]])
    } else {
      read_bootstrap_dir(dirs$dir[[i]])
    }
    sample_ids <- rownames(x$counts)
    if (isTRUE(ignore_after_bar)) {
      sample_ids <- sub("\\|.*$", "", sample_ids)
    }
    if (i == 1L) {
      ids <- sample_ids
      samp_type <- x$samp_type
      keep <- seq_along(ids)
      if (!is.null(transcripts)) {
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
      }
      rows <- keep
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
      rows <- match(ids, sample_ids)[keep]
      available <- min(available, ncol(x$counts))
      if (available < length(reps)) {
        reps <- reps[seq_len(available)]
      }
    }
    for (k in seq_along(reps)) {
      reps[[k]][, i] <- x$counts[rows, k]
    }
    rm(x)
  }
  if (!is.null(n_draws) && n_draws > available) {
    cli::cli_inform(c(
      i = "Only {available} replicate{?s} available per sample; using all of them."
    ))
  }
  attr(reps, "samp_type") <- samp_type
  reps
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
  for (dir in dirs$dir) {
    ec <- read_eq_classes(dir)
    if (is.null(ec)) {
      next
    }
    if (isTRUE(ignore_after_bar)) {
      ec$transcripts <- sub("\\|.*$", "", ec$transcripts)
    }
    if (is.null(ids)) {
      ids <- ec$transcripts
      total <- numeric(length(ids))
      unique_total <- numeric(length(ids))
      n_classes <- numeric(length(ids))
    }
    tid <- match(ec$transcripts, ids)[ec$members$tid]
    keep <- !is.na(tid)
    count <- ec$members$count[keep]
    is_unique <- ec$members$size[keep] == 1L
    tid <- tid[keep]
    total <- total + sum_by_index(count, tid, length(ids))
    unique_total <- unique_total +
      sum_by_index(count * is_unique, tid, length(ids))
    n_classes <- n_classes + tabulate(tid, length(ids))
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

#' Read the equivalence classes of one Salmon output directory
#'
#' The text format is the number of transcripts, the number of classes, the
#' transcript names and one line per class: the class size, the (0-based)
#' transcript indices, optionally their weights, and the fragment count.
#' @keywords internal
#' @noRd
read_eq_classes <- function(dir) {
  candidates <- fs::path(
    dir,
    "aux_info",
    c("eq_classes.txt.gz", "eq_classes.txt")
  )
  file <- candidates[fs::file_exists(candidates)]
  if (!length(file)) {
    return(NULL)
  }
  file <- file[[1L]]
  con <- if (grepl("\\.gz$", file)) gzfile(file, "rt") else file(file, "rt")
  on.exit(close(con), add = TRUE)
  n_tx <- as.integer(readLines(con, n = 1L))
  n_classes <- as.integer(readLines(con, n = 1L))
  tx <- readLines(con, n = n_tx)
  lines <- if (n_classes > 0L) readLines(con, n = n_classes) else character()
  members <- data.frame(
    tid = integer(),
    size = integer(),
    count = numeric(),
    class = integer()
  )
  if (length(lines)) {
    parts <- strsplit(lines, "\t", fixed = TRUE)
    lens <- lengths(parts)
    flat <- unlist(parts, use.names = FALSE)
    starts <- cumsum(c(1, lens))[seq_along(lens)]
    size <- as.integer(flat[starts])
    count <- as.numeric(flat[starts + lens - 1])
    cls <- rep.int(seq_along(size), size)
    pos <- starts[cls] + sequence(size)
    members <- data.frame(
      tid = as.integer(flat[pos]) + 1L,
      size = size[cls],
      count = count[cls],
      class = cls
    )
  }
  list(transcripts = tx, members = members)
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
inferential_replicates_table <- function(reps) {
  if (!is.list(reps) || !length(reps)) {
    cli::cli_abort(
      message = c(x = "{.arg reps} must be a non-empty list of matrices."),
      class = "isoformic_invalid_inferential_replicates"
    )
  }
  first <- reps[[1L]]
  same <- vapply(
    reps,
    function(m) identical(dimnames(m), dimnames(first)),
    logical(1L)
  )
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
    count = unlist(lapply(reps, as.vector), use.names = FALSE)
  )
}

#' @rdname inferential-io
#' @export
as_inferential_replicates <- function(table, samp_type = "gibbs") {
  required <- c("transcript_id", "sample_id", "draw", "count")
  if (!all(required %in% colnames(table))) {
    cli::cli_abort(
      message = c(
        x = "{.arg table} needs the columns {.field {required}}."
      ),
      class = "isoformic_invalid_inferential_replicates"
    )
  }
  ids <- unique(as.character(table$transcript_id))
  samples <- unique(as.character(table$sample_id))
  draws <- sort(unique(table$draw))
  n_cells <- length(ids) * length(samples) * length(draws)
  if (nrow(table) != n_cells || anyNA(table$count)) {
    cli::cli_abort(
      message = c(
        x = "{.arg table} must have one count for every transcript, sample and draw."
      ),
      class = "isoformic_invalid_inferential_replicates"
    )
  }
  slot <- match(table$transcript_id, ids) +
    length(ids) *
      ((match(table$sample_id, samples) - 1L) +
        length(samples) * (match(table$draw, draws) - 1L))
  values <- numeric(n_cells)
  values[slot] <- table$count
  per_draw <- length(ids) * length(samples)
  reps <- lapply(seq_along(draws), function(k) {
    matrix(
      values[(k - 1L) * per_draw + seq_len(per_draw)],
      nrow = length(ids),
      dimnames = list(ids, samples)
    )
  })
  attr(reps, "samp_type") <- samp_type
  reps
}
