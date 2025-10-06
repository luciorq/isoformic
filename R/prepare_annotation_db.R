#' Write Parquet File from GFF
#'
#' This function reads a GFF file and writes its contents to a Parquet file
#' using DuckDB.
#'
#' @param input_path Character string specifying the path to the input GFF file.
#' @param output_path Character string specifying the path to the output Parquet file.
#' If `NULL` or an empty string, a temporary file will be created.
#' @param file_type Character string specifying the type of the input file.
#' Currently, only "gff" is supported (default is "gff").
#'
#' @return Invisible path to the created Parquet file.
#'
#' @export
prepare_annotation_db <- function(
    input_path,
    output_path = NULL,
    file_type = c("gff")) {
  file_type <- stringr::str_to_lower(file_type)
  file_type <- rlang::arg_match(file_type)
  if (isTRUE(rlang::is_null(output_path) || identical(output_path, ""))) {
    parquet_file_path <- "none"
  } else {
    parquet_file_path <- fs::path(
      fs::path_ext_remove(output_path),
      ext = "parquet"
    )
  }

  if (!isTRUE(fs::file_exists(input_path))) {
    cli::cli_abort(
      message = c(
        x = "{.path {input_path}} do {.strong not} exist."
      ),
      class = "isoformic_annot_file_dont_exist"
    )
  }

  sql_str <- r"---(SET enable_progress_bar = false;
COPY (
SELECT t.seqid, t.source, t.type, t.start_pos,
  t.end_pos, t.score, t.strand, t.phase,
  split_part(regexp_split_to_table(t.attributes, ';'), '=', 1) AS key,
  split_part(regexp_split_to_table(t.attributes, ';'), '=', 2) AS value
FROM (
  SELECT *
  FROM read_csv(
    {`input_path`},
    delim = '\t', comment = '#',
    encoding='utf-8',
    header = false,
    names = [
      'seqid', 'source', 'type', 'start_pos', 'end_pos',
      'score', 'strand', 'phase', 'attributes'
    ]
  )
  WHERE type IN ('gene', 'transcript', 'exon')
) t
) TO {`parquet_file_path`} (FORMAT parquet, COMPRESSION zstd);
)---"
  duckdb_run(
    sql_string = sql_str,
    db_type = "duckdb_memory",
    db_file_path = NULL,
    read_only = TRUE
  )
  return(invisible(parquet_file_path))
}
