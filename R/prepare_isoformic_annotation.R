#' Write Feature Annotation to Parquet Files
#'
#' This function reads an annotation file and parse feature annotation to
#' Parquet files each level of required
#' feature (i.e. gene, transcript, and exon).
#'
#' @param input_path Character string specifying the path to the input GFF file.
#' @param output_path Character string specifying the path to the output
#' directory where Parquet files are written. If `NULL` or an empty string,
#' the cache directory will be used.
#' @param file_type Character string specifying the type of the input file.
#' Currently, only "gff" is supported (default is "gff").
#'
#' @return Invisible path to the created Parquet file.
#'
#' @export
prepare_isoformic_annotation <- function(
  input_path,
  output_path = NULL,
  file_type = c("gff")
) {
  file_type <- stringr::str_to_lower(file_type)
  file_type <- rlang::arg_match(file_type)

  if (!isTRUE(fs::file_exists(input_path))) {
    cli::cli_abort(
      message = c(
        x = "{.path {input_path}} do {.strong not} exist."
      ),
      class = "isoformic_annot_file_dont_exist"
    )
  }
  if (isTRUE(rlang::is_null(output_path) || identical(output_path, ""))) {
    parquet_path <- get_isoformic_cache()
  } else {
    parquet_path <- output_path
  }

  annot_name <- input_path |>
    fs::path_ext_remove() |>
    fs::path_file() |>
    stringr::str_to_lower() |>
    stringr::str_replace_all(
      pattern = stringr::fixed("."),
      replacement = "_"
    )

  parquet_dir <- fs::path(parquet_path, annot_name)

  if (!isTRUE(fs::dir_exists(parquet_dir))) {
    fs::dir_create(parquet_dir, recurse = TRUE)
  }

  parquet_gene_file_path <- fs::path(parquet_dir, "row_data_genes.parquet")
  parquet_transcript_file_path <- fs::path(
    parquet_dir,
    "row_data_transcripts.parquet"
  )
  parquet_exon_file_path <- fs::path(parquet_dir, "row_data_exons.parquet")

  sql_str <- r"---(SET enable_progress_bar = false;
CREATE TEMPORARY TABLE parsed_annotations AS
  SELECT
    type,
    seqid,
    start_pos,
    end_pos,
    strand,
    -- regexp_extract(attributes, 'ID=([^;]+)', 1) AS id,
    regexp_extract(attributes, 'gene_id=([^;]+)', 1) AS gene_id,
    regexp_extract(attributes, 'gene_name=([^;]+)', 1) AS gene_name,
    regexp_extract(attributes, 'gene_type=([^;]+)', 1) AS gene_type,
    regexp_extract(attributes, 'transcript_id=([^;]+)', 1) AS transcript_id,
    regexp_extract(attributes, 'transcript_name=([^;]+)', 1) AS transcript_name,
    regexp_extract(attributes, 'transcript_type=([^;]+)', 1) AS transcript_type,
    regexp_extract(attributes, 'exon_id=([^;]+)', 1) AS exon_id,
    regexp_extract(attributes, 'exon_number=([^;]+)', 1) AS exon_number
  FROM read_csv(
    {`input_path`},
    delim = '\t',
    comment = '#',
    encoding = 'utf-8',
    header = false,
    names = [
      'seqid', 'source', 'type', 'start_pos', 'end_pos',
      'score', 'strand', 'phase', 'attributes'
    ]
  )
WHERE type IN ('gene', 'transcript', 'exon');
COPY (
  SELECT
    gene_id,
    gene_name,
    gene_type,
    seqid,
    start_pos,
    end_pos,
    strand
  FROM parsed_annotations
  WHERE type = 'gene'
)
TO {`parquet_gene_file_path`} (FORMAT parquet, COMPRESSION zstd);
COPY (
  SELECT
    transcript_id,
    transcript_name,
    transcript_type,
    gene_id,
    seqid,
    start_pos,
    end_pos,
    strand
  FROM parsed_annotations
  WHERE type = 'transcript'
)
TO {`parquet_transcript_file_path`} (FORMAT parquet, COMPRESSION zstd);
COPY (
  SELECT
    exon_id,
    transcript_id,
    exon_number,
    seqid,
    start_pos,
    end_pos,
    strand
  FROM parsed_annotations
  WHERE type = 'exon'
)
TO {`parquet_exon_file_path`} (FORMAT parquet, COMPRESSION zstd);
)---"

  duckdb_run(
    sql_string = sql_str,
    db_type = "duckdb_memory",
    db_file_path = NULL,
    read_only = TRUE
  )
  return(invisible(parquet_dir))
}
