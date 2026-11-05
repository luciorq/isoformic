#' Read Oarfish Quantifications Into an IsoformicExperiment
#'
#' Oarfish quantifies long reads (Oxford Nanopore and PacBio) against a
#' transcriptome. For every sample it writes, next to a common output prefix,
#' `<prefix>.quant` (transcript name, length and estimated read count),
#' `<prefix>.ambig_info.tsv` (reads that map to the transcript alone and reads
#' shared with others), `<prefix>.meta_info.json` and, with
#' `--num-bootstraps`, `<prefix>.infreps.pq` (the bootstrap replicates as a
#' Parquet table). `read_oarfish()` reads them like [read_salmon()] reads
#' Salmon output: an [IsoformicExperiment] with `counts`, `length` and
#' `abundance` assays, and, with `inferential = TRUE`, the uncertainty tables
#' of [add_uncertainty()].
#'
#' Oarfish does not write TPM; the `abundance` assay is computed from the
#' read counts and the transcript lengths (reads per base, scaled to one
#' million). The ambiguity counts of `.ambig_info.tsv` play the role of the
#' equivalence classes of Salmon in [equivalence_support()]: `unique_fraction`
#' is the share of reads mapped to the transcript alone, while `n_classes` is
#' not available. The functions [read_inferential_replicates()],
#' [equivalence_support()], [salmon_run_qc()] and [add_uncertainty()] accept
#' Oarfish prefixes as `paths`, but not together with Salmon directories.
#'
#' @param paths Output prefixes of Oarfish (the value of `-o`), paths to
#'   the `.quant` files, or directories holding a single `.quant` file. Names,
#'   when present, are the sample IDs; otherwise the prefix is used.
#' @inheritParams read_salmon
#' @param inferential Logical. Read the bootstrap replicates and the ambiguity
#'   counts and compute the uncertainty tables with [add_uncertainty()]. The
#'   samples must have been quantified with `--num-bootstraps`.
#'
#' @returns An [IsoformicExperiment].
#'
#' @examples
#' prefixes <- file.path(
#'   system.file("extdata", "oarfish", package = "isoformic"),
#'   c("control_1", "control_2", "treated_1")
#' )
#' names(prefixes) <- basename(prefixes)
#' annot <- system.file(
#'   "extdata", "example_experiment", "annotation",
#'   package = "isoformic"
#' )
#' iso <- read_oarfish(prefixes, annot_path = annot, inferential = TRUE)
#' iso
#' head(tx_uncertainty(iso))
#'
#' @export
read_oarfish <- function(
  paths,
  annot_path,
  col_data = NULL,
  experiment_name = NULL,
  ignore_after_bar = TRUE,
  drop_unannotated = FALSE,
  inferential = FALSE
) {
  rlang::check_required(paths)
  rlang::check_required(annot_path)
  prefixes <- vapply(paths, oarfish_prefix, character(1L))
  sample_ids <- names(paths)
  if (is.null(sample_ids) || !all(nzchar(sample_ids))) {
    sample_ids <- as.character(fs::path_file(prefixes))
  }
  if (anyDuplicated(sample_ids)) {
    cli::cli_abort(
      message = c(x = "Sample IDs derived from {.arg paths} are not unique."),
      class = "isoformic_duplicated_sample_ids"
    )
  }
  quants <- lapply(prefixes, read_oarfish_quant)
  ids <- quants[[1L]]$transcript_id
  for (i in seq_along(quants)[-1L]) {
    if (!setequal(quants[[i]]$transcript_id, ids)) {
      cli::cli_abort(
        message = c(
          x = paste0(
            "Sample {.val {sample_ids[i]}} was quantified against a ",
            "different set of transcripts than {.val {sample_ids[1L]}}."
          ),
          i = "Quantify every sample against the same transcriptome."
        ),
        class = "isoformic_inconsistent_transcripts"
      )
    }
  }
  column <- function(name) {
    mat <- vapply(
      quants,
      function(q) q[[name]][match(ids, q$transcript_id)],
      numeric(length(ids))
    )
    dimnames(mat) <- list(ids, sample_ids)
    mat
  }
  counts <- column("num_reads")
  lengths <- column("len")
  per_base <- counts / lengths
  assays <- list(
    abundance = sweep(per_base, 2L, colSums(per_base), "/") * 1e6,
    counts = counts,
    length = lengths
  )
  assays <- lapply(assays, function(m) {
    m[is.nan(m)] <- 0
    m
  })
  if (isTRUE(ignore_after_bar)) {
    assays <- lapply(assays, function(m) {
      rownames(m) <- sub("\\|.*$", "", rownames(m))
      m
    })
  }
  if (is.null(col_data)) {
    col_data <- tibble::tibble(sample_id = sample_ids)
  }
  annot_dir <- resolve_annot_path(annot_path)
  if (isTRUE(drop_unannotated)) {
    assays <- drop_unannotated_rows(assays, annot_dir)
  }
  iso <- IsoformicExperiment(
    experiment_name = experiment_name,
    annot_path = annot_dir,
    assays = assays,
    col_data = col_data
  )
  if (isTRUE(inferential)) {
    iso <- add_uncertainty(iso, paths = stats::setNames(prefixes, sample_ids))
  }
  iso
}

#' Output prefix of an Oarfish sample
#' @keywords internal
#' @noRd
oarfish_prefix <- function(path, arg = "paths") {
  path <- as.character(path)
  if (isTRUE(fs::is_dir(path))) {
    found <- fs::dir_ls(path, glob = "*.quant", type = "file")
    if (length(found) != 1L) {
      cli::cli_abort(
        message = c(
          x = "Expected one {.file .quant} file in {.path {path}}, found {length(found)}.",
          i = "Pass the output prefix of Oarfish instead."
        ),
        class = "isoformic_no_oarfish_quant"
      )
    }
    path <- as.character(found)
  }
  prefix <- sub("\\.quant$", "", path)
  assert_file_exists(paste0(prefix, ".quant"), arg = arg)
  as.character(fs::path_abs(prefix))
}

#' Read the `.quant` table of one Oarfish sample
#' @keywords internal
#' @noRd
read_oarfish_quant <- function(prefix) {
  quant <- read_delimited(paste0(prefix, ".quant"), null_strings = "NA")
  tibble::tibble(
    transcript_id = as.character(quant$tname),
    len = as.numeric(quant$len),
    num_reads = as.numeric(quant$num_reads)
  )
}

#' Read the bootstrap replicates of one Oarfish sample
#'
#' The Parquet table has one column per replicate (`bootstrap.0`, ...) and one
#' row per transcript, in the order of the `.quant` file.
#' @keywords internal
#' @noRd
read_oarfish_replicates <- function(prefix) {
  file <- paste0(prefix, ".infreps.pq")
  if (!isTRUE(fs::file_exists(file))) {
    cli::cli_abort(
      message = c(
        x = "No inferential replicates found for {.path {prefix}}.",
        i = "Run {.code oarfish --num-bootstraps 30} to write {.file {fs::path_file(file)}}."
      ),
      class = "isoformic_no_inferential_replicates"
    )
  }
  reps <- read_parquet(file)
  quant <- read_oarfish_quant(prefix)
  if (nrow(reps) != nrow(quant)) {
    cli::cli_abort(
      message = c(
        x = "{.file {fs::path_file(file)}} has {nrow(reps)} rows, {.file .quant} has {nrow(quant)}."
      ),
      class = "isoformic_invalid_inferential_replicates"
    )
  }
  order <- order(as.integer(sub("^.*\\.", "", colnames(reps))))
  counts <- as.matrix(reps[, order, drop = FALSE])
  storage.mode(counts) <- "double"
  dimnames(counts) <- list(quant$transcript_id, NULL)
  list(counts = counts, samp_type = "bootstrap")
}

#' Read support from the unique and ambiguous read counts of Oarfish
#' @keywords internal
#' @noRd
oarfish_support <- function(prefixes, transcripts, ignore_after_bar) {
  ids <- NULL
  total <- unique_total <- NULL
  for (prefix in prefixes) {
    file <- paste0(prefix, ".ambig_info.tsv")
    if (!isTRUE(fs::file_exists(file))) {
      cli::cli_abort(
        message = c(
          x = "No ambiguity counts found for {.path {prefix}}.",
          i = "Oarfish writes {.file {fs::path_file(file)}} next to the {.file .quant} file."
        ),
        class = "isoformic_no_equivalence_classes"
      )
    }
    quant <- read_oarfish_quant(prefix)
    ambig <- read_delimited(file, null_strings = "NA")
    if (is.null(ids)) {
      ids <- quant$transcript_id
      total <- unique_total <- numeric(length(ids))
    }
    total <- total + as.numeric(ambig$total_reads)
    unique_total <- unique_total + as.numeric(ambig$unique_reads)
  }
  if (isTRUE(ignore_after_bar)) {
    ids <- sub("\\|.*$", "", ids)
  }
  res <- tibble::tibble(
    transcript_id = ids,
    n_classes = NA_real_,
    n_fragments = total,
    n_unique_fragments = unique_total,
    unique_fraction = ifelse(total > 0, unique_total / total, NA_real_),
    n_samples = length(prefixes)
  )
  if (!is.null(transcripts)) {
    res <- res[res$transcript_id %in% transcripts, ]
  }
  res
}

#' Run metadata of Oarfish samples, in the layout of `salmon_run_qc()`
#' @keywords internal
#' @noRd
oarfish_run_qc <- function(dirs, min_percent_mapped) {
  rows <- lapply(dirs$dir, function(prefix) {
    file <- paste0(prefix, ".meta_info.json")
    if (!isTRUE(fs::file_exists(file))) {
      cli::cli_abort(
        message = c(
          x = "{.file {fs::path_file(file)}} not found.",
          i = "Is {.path {prefix}} an Oarfish output prefix?"
        ),
        class = "isoformic_no_meta_info"
      )
    }
    text <- paste(readLines(file, warn = FALSE), collapse = "\n")
    num <- function(key) suppressWarnings(as.numeric(json_field(text, key)))
    mapped <- num("valid_best_aln")
    processed <- sum(
      mapped,
      num("no_mapping"),
      num("no_valid_aln"),
      na.rm = TRUE
    )
    tibble::tibble(
      num_processed = processed,
      num_mapped = mapped,
      percent_mapped = if (processed > 0) {
        100 * mapped / processed
      } else {
        NA_real_
      },
      replicate_type = "bootstrap",
      num_replicates = num("num_bootstraps"),
      optimiser = "em",
      alignment_source = json_field(text, "alignment_source"),
      prob_model = json_field(text, "prob_model")
    )
  })
  res <- dplyr::bind_cols(
    tibble::tibble(sample_id = dirs$sample_id),
    dplyr::bind_rows(rows)
  )
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
      if (length(flags)) paste(flags, collapse = ";") else NA_character_
    },
    character(1L)
  )
  attr(res, "tool") <- "oarfish"
  res
}

#' Read Oarfish single-cell quantifications
#'
#' With `--single-cell`, Oarfish quantifies a transcriptome-aligned BAM whose
#' reads carry a cell barcode (`CB:Z` tag) and writes, next to its output
#' prefix, `<prefix>.count.mtx` (a Matrix Market table of cells by
#' transcripts), `<prefix>.barcodes.txt`, `<prefix>.features.txt` and
#' `<prefix>.meta_info.json`. It draws no replicates in this mode.
#'
#' * `read_oarfish_sc_counts()` returns the non-zero counts as a long table,
#'   one row per cell and transcript.
#' * `read_oarfish_sc()` sums the cells of each sample (a pseudobulk) and
#'   returns an [IsoformicExperiment] with `counts`, `length` (from the
#'   annotation) and `abundance` (TPM computed from both) assays, so the
#'   isoform analyses of the package apply to single-cell data at the level of
#'   samples or cell types. `cells` says which sample each barcode belongs to;
#'   barcodes absent from it are dropped. To compare cell types, give each
#'   combination of sample and cell type its own `sample_id` and describe it in
#'   `col_data`.
#'
#' @param path The output prefix of Oarfish (the value of `-o`), the
#'   `.count.mtx` file, or a directory holding a single `.count.mtx` file.
#' @param cells A data frame with a `barcode` column and a `sample_id` column.
#' @param col_data Optional sample metadata with a `sample_id` column; the
#'   number of cells and reads of each sample are added as `n_cells` and
#'   `n_reads`.
#' @param min_cells Samples with fewer cells are dropped.
#' @inheritParams read_salmon
#'
#' @returns `read_oarfish_sc_counts()`: a tibble with `barcode`,
#'   `transcript_id` and `count`. `read_oarfish_sc()`: an [IsoformicExperiment].
#'
#' @examples
#' dir <- system.file("extdata", "oarfish_sc", package = "isoformic")
#' counts <- read_oarfish_sc_counts(file.path(dir, "out"))
#' head(counts)
#' cells <- read.delim(file.path(dir, "cells.tsv"))
#' iso <- read_oarfish_sc(
#'   file.path(dir, "out"),
#'   cells = cells,
#'   annot_path = system.file(
#'     "extdata", "example_experiment", "annotation",
#'     package = "isoformic"
#'   )
#' )
#' col_data(iso)
#'
#' @export
read_oarfish_sc_counts <- function(path, ignore_after_bar = TRUE) {
  rlang::check_required(path)
  prefix <- oarfish_sc_prefix(as.character(path))
  barcodes <- readLines(paste0(prefix, ".barcodes.txt"), warn = FALSE)
  features <- readLines(paste0(prefix, ".features.txt"), warn = FALSE)
  entries <- read_matrix_market(paste0(prefix, ".count.mtx"))
  if (
    !identical(entries$dims, c(length(barcodes), length(features))) &&
      !identical(entries$dims, c(length(features), length(barcodes)))
  ) {
    cli::cli_abort(
      message = c(
        x = paste0(
          "{.file {fs::path_file(prefix)}.count.mtx} is {entries$dims[1]} by ",
          "{entries$dims[2]}, but there are {length(barcodes)} barcodes and ",
          "{length(features)} features."
        )
      ),
      class = "isoformic_invalid_oarfish_sc"
    )
  }
  # Oarfish writes cells as rows; accept the transposed layout as well
  cells_as_rows <- identical(entries$dims[[1L]], length(barcodes))
  cell <- if (cells_as_rows) entries$table$i else entries$table$j
  feature <- if (cells_as_rows) entries$table$j else entries$table$i
  if (isTRUE(ignore_after_bar)) {
    features <- sub("\\|.*$", "", features)
  }
  res <- tibble::tibble(
    barcode = barcodes[cell],
    transcript_id = features[feature],
    count = entries$table$value
  )
  res[res$count != 0, ]
}

#' @rdname read_oarfish_sc_counts
#' @export
read_oarfish_sc <- function(
  path,
  cells,
  annot_path,
  col_data = NULL,
  experiment_name = NULL,
  min_cells = 1L,
  ignore_after_bar = TRUE,
  drop_unannotated = FALSE
) {
  rlang::check_required(path)
  rlang::check_required(cells)
  rlang::check_required(annot_path)
  if (
    !is.data.frame(cells) || !all(c("barcode", "sample_id") %in% names(cells))
  ) {
    cli::cli_abort(
      message = c(
        x = "{.arg cells} needs {.field barcode} and {.field sample_id} columns."
      ),
      class = "isoformic_invalid_oarfish_sc"
    )
  }
  if (anyDuplicated(cells$barcode)) {
    cli::cli_abort(
      message = c(x = "A barcode is listed more than once in {.arg cells}."),
      class = "isoformic_invalid_oarfish_sc"
    )
  }
  prefix <- oarfish_sc_prefix(as.character(path))
  counts_long <- read_oarfish_sc_counts(
    prefix,
    ignore_after_bar = ignore_after_bar
  )
  features <- readLines(paste0(prefix, ".features.txt"), warn = FALSE)
  if (isTRUE(ignore_after_bar)) {
    features <- sub("\\|.*$", "", features)
  }
  barcodes <- readLines(paste0(prefix, ".barcodes.txt"), warn = FALSE)
  cells <- cells[as.character(cells$barcode) %in% barcodes, , drop = FALSE]
  unknown <- setdiff(barcodes, as.character(cells$barcode))
  if (length(unknown)) {
    cli::cli_inform(c(
      i = "Dropping {length(unknown)} barcode{?s} not listed in {.arg cells}."
    ))
  }
  n_cells <- table(as.character(cells$sample_id))
  samples <- names(n_cells)[n_cells >= min_cells]
  if (!length(samples)) {
    cli::cli_abort(
      message = c(x = "No sample has at least {min_cells} cell{?s}."),
      class = "isoformic_invalid_oarfish_sc"
    )
  }
  sample_of <- as.character(cells$sample_id)[match(
    counts_long$barcode,
    as.character(cells$barcode)
  )]
  keep <- sample_of %in% samples
  summed <- tapply(
    counts_long$count[keep],
    list(
      factor(counts_long$transcript_id[keep], levels = features),
      factor(sample_of[keep], levels = samples)
    ),
    sum
  )
  summed[is.na(summed)] <- 0
  counts <- matrix(
    summed,
    nrow = length(features),
    dimnames = list(features, samples)
  )
  annot_dir <- resolve_annot_path(annot_path)
  tx <- read_annotation_table(
    annot_dir,
    "transcripts",
    columns = c("transcript_id", "tx_length")
  )
  len <- as.numeric(tx$tx_length[match(features, tx$transcript_id)])
  lengths <- matrix(
    len,
    nrow = length(features),
    ncol = length(samples),
    dimnames = dimnames(counts)
  )
  per_base <- counts / lengths
  per_base[!is.finite(per_base)] <- 0
  assays <- list(
    abundance = sweep(per_base, 2L, colSums(per_base), "/") * 1e6,
    counts = counts,
    length = lengths
  )
  assays$abundance[!is.finite(assays$abundance)] <- 0
  if (isTRUE(drop_unannotated)) {
    assays <- drop_unannotated_rows(assays, annot_dir)
  }
  summary <- tibble::tibble(
    sample_id = samples,
    n_cells = as.integer(n_cells[samples]),
    n_reads = unname(colSums(counts))
  )
  col_data <- if (is.null(col_data)) {
    summary
  } else {
    dplyr::left_join(
      tibble::as_tibble(col_data)[
        as.character(col_data$sample_id) %in% samples,
      ],
      summary,
      by = "sample_id"
    )
  }
  IsoformicExperiment(
    experiment_name = experiment_name,
    annot_path = annot_dir,
    assays = assays,
    col_data = col_data
  )
}

#' Output prefix of an Oarfish single-cell run
#' @keywords internal
#' @noRd
oarfish_sc_prefix <- function(path) {
  if (isTRUE(fs::is_dir(path))) {
    found <- fs::dir_ls(path, glob = "*.count.mtx", type = "file")
    if (length(found) != 1L) {
      cli::cli_abort(
        message = c(
          x = "Expected one {.file .count.mtx} file in {.path {path}}, found {length(found)}.",
          i = "Pass the output prefix of Oarfish instead."
        ),
        class = "isoformic_invalid_oarfish_sc"
      )
    }
    path <- as.character(found)
  }
  prefix <- sub("\\.count\\.mtx$", "", path)
  for (ext in c(".count.mtx", ".barcodes.txt", ".features.txt")) {
    assert_file_exists(paste0(prefix, ext), arg = "path")
  }
  as.character(fs::path_abs(prefix))
}

#' Read a coordinate Matrix Market file with DuckDB
#' @keywords internal
#' @noRd
read_matrix_market <- function(path, envir = parent.frame()) {
  con_file <- file(path, open = "r")
  on.exit(close(con_file), add = TRUE)
  n_skip <- 0L
  repeat {
    line <- readLines(con_file, n = 1L)
    if (!length(line)) {
      break
    }
    n_skip <- n_skip + 1L
    if (!startsWith(line, "%")) {
      break
    }
  }
  if (!length(line) || startsWith(line, "%")) {
    cli::cli_abort(
      message = c(x = "{.file {path}} is not a Matrix Market file."),
      class = "isoformic_invalid_oarfish_sc"
    )
  }
  dims <- as.integer(strsplit(trimws(line), "\\s+")[[1L]])
  if (length(dims) != 3L || anyNA(dims)) {
    cli::cli_abort(
      message = c(
        x = "{.file {path}} has no size line of a coordinate Matrix Market file."
      ),
      class = "isoformic_invalid_oarfish_sc"
    )
  }
  if (dims[[3L]] == 0L) {
    empty <- data.frame(i = integer(), j = integer(), value = numeric())
    return(list(dims = dims[1:2], table = empty))
  }
  con <- duckdb_connect(envir = envir)
  path <- as.character(fs::path_abs(path))
  sql <- glue::glue_sql(
    "SELECT * FROM read_csv({path}, delim = ' ', header = false, skip = {n_skip},
       columns = {{'i': 'INTEGER', 'j': 'INTEGER', 'value': 'DOUBLE'}})",
    .con = con
  )
  table <- DBI::dbGetQuery(con, sql)
  if (nrow(table) != dims[[3L]]) {
    cli::cli_abort(
      message = c(
        x = "{.file {path}} declares {dims[3]} entries but has {nrow(table)}."
      ),
      class = "isoformic_invalid_oarfish_sc"
    )
  }
  list(dims = dims[1:2], table = table)
}
