#' Create a Transcript-to-Gene Mapping Table
#'
#' Builds a transcript-to-gene table from a GENCODE transcriptome FASTA
#' file (parsing the `|`-separated headers) or from a GFF3/GTF annotation
#' (through [prepare_isoformic_annotation()]). Parsing is done by DuckDB,
#' so a full GENCODE FASTA (about 250,000 headers, gzip-compressed) takes a
#' couple of seconds.
#'
#' With an [IsoformicExperiment] object, use [tx_to_gene()] or [tx_annot()],
#' which read the same information from the annotation directory.
#'
#' @param file_path Path to a FASTA (`.fa`, `.fasta`, optionally `.gz`),
#'   GFF3 or GTF file, or to an annotation directory.
#' @param file_type `"auto"` (default) detects the format from the file
#'   extension; otherwise one of `"fasta"`, `"gff"`, `"gtf"`.
#' @param output_path Passed to [prepare_isoformic_annotation()] for
#'   GFF3/GTF inputs.
#'
#' @returns A tibble with columns `transcript_id`, `gene_id`,
#'   `havana_gene_id`, `havana_transcript_id`, `transcript_name`,
#'   `gene_name`, `tx_length` and `transcript_type`. The Havana columns are
#'   `NA` for GFF3/GTF inputs.
#'
#' @examples
#' gff <- system.file("extdata", "example.gff3", package = "isoformic")
#' tx2gene <- make_tx_to_gene(gff, output_path = tempdir())
#' head(tx2gene)
#'
#' @export
make_tx_to_gene <- function(
  file_path,
  file_type = c("auto", "fasta", "gff", "gtf"),
  output_path = NULL
) {
  file_type <- rlang::arg_match(file_type)
  if (is_annotation_dir(file_path)) {
    return(tx_to_gene_from_annotation(file_path))
  }
  assert_file_exists(file_path)
  if (identical(file_type, "auto")) {
    ext <- tolower(fs::path_ext(sub("\\.(gz|bgz|zst|bz2)$", "", file_path)))
    file_type <- if (ext %in% c("fa", "fasta", "fna")) {
      "fasta"
    } else {
      detect_annotation_type(file_path)
    }
  }
  if (identical(file_type, "fasta")) {
    return(tx_to_gene_from_fasta(file_path))
  }
  annot_dir <- prepare_isoformic_annotation(
    input_path = file_path,
    output_path = output_path,
    file_type = file_type
  )
  tx_to_gene_from_annotation(annot_dir)
}

tx_to_gene_cols <- function() {
  c(
    "transcript_id",
    "gene_id",
    "havana_gene_id",
    "havana_transcript_id",
    "transcript_name",
    "gene_name",
    "tx_length",
    "transcript_type"
  )
}

#' Parse GENCODE FASTA headers with DuckDB
#' @keywords internal
#' @noRd
tx_to_gene_from_fasta <- function(file_path) {
  file_path <- as.character(fs::path_abs(file_path))
  con <- duckdb_connect()
  sql <- glue::glue_sql(
    "WITH headers AS (
      SELECT substr(line, 2) AS header
      FROM read_csv(
        {file_path},
        delim = '\t', header = false, auto_detect = false,
        quote = '', escape = '', columns = {{'line': 'VARCHAR'}}
      )
      WHERE line LIKE '>%'
    )
    SELECT
      split_part(header, '|', 1) AS transcript_id,
      nullif(split_part(header, '|', 2), '') AS gene_id,
      nullif(split_part(header, '|', 3), '') AS havana_gene_id,
      nullif(split_part(header, '|', 4), '') AS havana_transcript_id,
      nullif(split_part(header, '|', 5), '') AS transcript_name,
      nullif(split_part(header, '|', 6), '') AS gene_name,
      TRY_CAST(split_part(header, '|', 7) AS INTEGER) AS tx_length,
      nullif(split_part(header, '|', 8), '') AS transcript_type
    FROM headers",
    .con = con
  )
  res <- tibble::as_tibble(DBI::dbGetQuery(con, sql))
  if (!nrow(res)) {
    cli::cli_abort(
      message = c(x = "No FASTA headers found in {.path {file_path}}."),
      class = "isoformic_fasta_no_headers"
    )
  }
  res
}

#' Transcript-to-gene table from an annotation directory
#' @keywords internal
#' @noRd
tx_to_gene_from_annotation <- function(annot_dir) {
  transcripts <- read_annotation_table(annot_dir, "transcripts")
  genes <- read_annotation_table(annot_dir, "genes")
  transcripts |>
    dplyr::left_join(
      genes[, c("gene_id", "gene_name")],
      by = "gene_id"
    ) |>
    dplyr::mutate(
      havana_gene_id = NA_character_,
      havana_transcript_id = NA_character_
    ) |>
    dplyr::select(dplyr::all_of(tx_to_gene_cols()))
}
