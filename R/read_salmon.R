#' Read Salmon Quantifications Into an IsoformicExperiment
#'
#' Reads `quant.sf` files produced by Salmon (`salmon quant`) directly with
#' DuckDB, without going through `tximport`, and returns an
#' [IsoformicExperiment] with transcript-level `abundance` (TPM), `counts`
#' (`NumReads`) and `length` (`EffectiveLength`) assays. Gene-level assays
#' are derived lazily by the object (see [gene_assays()]).
#'
#' All files are read in one DuckDB query and pivoted to wide matrices, so
#' reading tens of samples with 200,000 transcripts each takes seconds.
#' Inferential replicates (Gibbs/bootstrap samples) are read only with
#' `inferential = TRUE`, which also computes the uncertainty tables of
#' [add_uncertainty()]; for `swish` they are imported with
#' `tximeta`/`fishpond` and converted with [as_isoformic()].
#'
#' @param paths Character vector of Salmon output directories (containing
#'   `quant.sf`) or paths to `quant.sf` files. Names, when present, are used
#'   as sample IDs; otherwise the directory name is used.
#' @param annot_path Annotation directory or GFF3/GTF file (see
#'   [prepare_isoformic_annotation()]).
#' @param col_data Optional data frame of sample metadata whose first column
#'   (or a `sample_id` column) matches the sample IDs. When `NULL`, a
#'   minimal table with just `sample_id` is created.
#' @param experiment_name Optional experiment name.
#' @param ignore_after_bar Logical. Keep only the part of transcript names
#'   before the first `|` (GENCODE FASTA headers used as Salmon references).
#' @param drop_unannotated Logical. Drop transcripts absent from the
#'   annotation (with a message) instead of failing.
#' @param inferential Logical. Read the inferential replicates, equivalence
#'   classes and run metadata of the Salmon directories and compute the
#'   uncertainty tables with [add_uncertainty()] (default settings; call
#'   that function for other cutoffs or a contrast). The directories must
#'   have been quantified with `--numGibbsSamples` or `--numBootstraps`.
#'
#' @returns An [IsoformicExperiment].
#'
#' @examples
#' quant_dirs <- list.files(
#'   system.file("extdata", "mini_quants", package = "isoformic"),
#'   full.names = TRUE
#' )
#' names(quant_dirs) <- sub("_quant$", "", basename(quant_dirs))
#' annot <- system.file("extdata", "mini_annotation", package = "isoformic")
#' iso <- read_salmon(quant_dirs, annot_path = annot)
#' iso
#'
#' @export
read_salmon <- function(
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
  quant_files <- vapply(
    paths,
    function(p) {
      if (isTRUE(fs::is_dir(p))) {
        p <- fs::path(p, "quant.sf")
      }
      assert_file_exists(p, arg = "paths")
      as.character(fs::path_abs(p))
    },
    character(1L)
  )
  sample_ids <- names(paths)
  if (is.null(sample_ids) || !all(nzchar(sample_ids))) {
    sample_ids <- fs::path_file(fs::path_dir(quant_files))
  }
  if (anyDuplicated(sample_ids)) {
    cli::cli_abort(
      message = c(x = "Sample IDs derived from {.arg paths} are not unique."),
      class = "isoformic_duplicated_sample_ids"
    )
  }
  file_map <- data.frame(
    filename = unname(quant_files),
    sample_id = unname(sample_ids),
    stringsAsFactors = FALSE
  )

  con <- duckdb_connect()
  duckdb::duckdb_register(con, "file_map", file_map)
  files_sql <- DBI::dbQuoteString(con, unname(quant_files))
  DBI::dbExecute(
    con,
    glue::glue_sql(
      "CREATE TEMP TABLE quants AS
       SELECT m.sample_id, q.Name AS transcript_id,
              q.TPM, q.NumReads, q.EffectiveLength
       FROM read_csv(
         [{files_sql*}], delim = '\t', header = true, filename = true,
         columns = {{'Name': 'VARCHAR', 'Length': 'INTEGER',
                     'EffectiveLength': 'DOUBLE', 'TPM': 'DOUBLE',
                     'NumReads': 'DOUBLE'}}
       ) q
       JOIN file_map m ON q.filename = m.filename",
      files_sql = files_sql,
      .con = con
    )
  )
  # The `IN` list makes the PIVOT a single static SELECT; without it DuckDB
  # plans several statements and duckdb >= 2.0 warns while counting the rows.
  pivot <- function(value_col) {
    sql <- glue::glue_sql(
      "PIVOT (SELECT sample_id, transcript_id, {DBI::SQL(value_col)} AS v FROM quants)
       ON sample_id IN ({unname(sample_ids)*}) USING first(v)
       ORDER BY transcript_id",
      .con = con
    )
    df <- DBI::dbGetQuery(con, sql)
    mat <- as.matrix(df[, sample_ids, drop = FALSE])
    rownames(mat) <- df$transcript_id
    storage.mode(mat) <- "double"
    mat
  }
  assays <- list(
    abundance = pivot("TPM"),
    counts = pivot("NumReads"),
    length = pivot("EffectiveLength")
  )
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
    iso <- add_uncertainty(
      iso,
      paths = stats::setNames(dirname(quant_files), sample_ids)
    )
  }
  iso
}

#' Drop assay rows whose transcript is absent from the annotation
#' @keywords internal
#' @noRd
drop_unannotated_rows <- function(assays, annot_dir) {
  annotated <- read_annotation_table(
    annot_dir,
    "transcripts",
    columns = "transcript_id"
  )$transcript_id
  keep <- rownames(assays[[1]]) %in% annotated
  if (!any(keep)) {
    cli::cli_abort(
      message = c(
        x = "None of the {length(keep)} transcript IDs are present in the annotation.",
        i = "Check that the annotation matches the reference used for quantification."
      ),
      class = "isoformic_no_annotated_transcripts"
    )
  }
  n_drop <- sum(!keep)
  if (n_drop > 0L) {
    cli::cli_inform(c(
      i = "Dropping {n_drop} transcript{?s} not present in the annotation."
    ))
    assays <- lapply(assays, function(m) m[keep, , drop = FALSE])
  }
  assays
}
