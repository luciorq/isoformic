#' Create Transcript-to-Gene Relationship Table
#'
#' Extracts a transcript-to-gene mapping table from GENCODE annotation files, such as the transcriptome FASTA file.
#' Currently, only FASTA files are supported.
#'
#' @param file_path A character string specifying the path to the reference file (e.g., GENCODE FASTA file).
#' @param file_type A character string specifying the type of the reference file. Currently, only `"fasta"` is supported.
#'   Default is `"fasta"`.
#'
#' @returns A `tibble` containing the transcript-to-gene mapping information, including transcript IDs, gene IDs,
#'   transcript names, gene names, and transcript types.
#'
#' @details
#' The function reads the headers of the FASTA file and extracts relevant information to create a mapping table.
#' For GTF or GFF3 files, support is not yet implemented.
#'
#' @examples
#' \dontrun{
#' # Assuming you have downloaded the GENCODE transcriptome FASTA file:
#' fasta_file <- download_reference(
#'   version = "43",
#'   organism = "human",
#'   file_type = "fasta",
#'   output_path = "data-raw"
#' )
#'
#' # Create the transcript-to-gene mapping table
#' tx_to_gene <- make_tx_to_gene(file_path = fasta_file, file_type = "fasta")
#'
#' # View the first few rows
#' utils::head(tx_to_gene)
#' }
#'
#' @export
make_tx_to_gene <- function(
  file_path,
  file_type = c("fasta", "gff", "gtf")
) {
  .data <- rlang::.data
  file_type <- rlang::arg_match(file_type)
  if (!isTRUE(fs::file_exists(file_path))) {
    cli::cli_abort(
      c(`x` = "{.path {file_path}} do not exist.")
    )
  }
  if (isTRUE(file_type == "fasta")) {
    fasta_lines <- readr::read_lines(file_path)
    vector_detect <- stringr::str_detect(fasta_lines, "^>.")
    header_table <- fasta_lines[vector_detect]
    header_table <- header_table |>
      tibble::as_tibble() |>
      tidyr::separate(
        .data[["value"]],
        into = paste0("col_", 1:9),
        sep = "\\|"
      ) |>
      dplyr::select(-c("col_9"))
    header_table$col_1 <- stringr::str_remove(header_table$col_1, ">")
  } else if (isTRUE(file_type == "gff")) {
    # TODO: @luciorq Actually implement parsing GFF3 and GTF files
    cli::cli_abort(
      c(
        `x` = "{.var GTF} and {.var GFF3} files are not supported yet.",
        `!` = "Use the GENCODE {.var FASTA} file."
      )
    )
  }
  names(header_table) <- c(
    "transcript_id",
    "gene_id",
    "havanna_gene_id",
    "havanna_transcript_id",
    "transcript_name",
    "gene_name",
    "tx_length",
    "transcript_type"
  )
  return(header_table)
}
