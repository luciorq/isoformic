#' Prepare Annotation
#'
#' Prepare annotation to be imported as `rowRanges` and `rowData` for both
#'   Genes, Transcripts and Exons based Position Annotation Table.
#'   From a GTF or GFF3 annotation file.
#'
#' @param file_path Path to annotation file.
#'
#' @param file_type Character indicating the type of file to download.
#'   One of `"gtf"` or `"gff"`. Defaults to `"gtf"`.
#'
#' @export
prepare_annotation <- function(file_path,
                               file_type = c("gtf", "gff")) {
  file_type <- stringr::str_to_lower(file_type)
  file_type <- rlang::arg_match(file_type)
  .data <- rlang::.data
  if (!isTRUE(fs::file_exists(file_path))) {
    cli::cli_abort(
      message = c(
        x = "{.path {file_path}} do {.strong not} exist."
      ),
      class = "isoformic_annot_file_dont_exist"
    )
  }
  annot_df <- vroom::vroom(
    file = file_path,
    delim = "\t",
    col_names = FALSE,
    col_types = list(
      X1 = "c",
      X2 = "_",
      X3 = "c",
      X4 = "i",
      X5 = "i",
      X6 = "_",
      X7 = "c",
      X8 = "_",
      X9 = "c"
    ),
    comment = "#",
    # n_max = 10,
    progress = FALSE
  )

  colnames(annot_df) <- c(
    "chr", "type", "start", "stop", "strand", "attributes"
  )

  annot_list <- list()
  for (i in c("gene", "transcript", "exon")) {
    attr_df <- annot_df |>
      dplyr::select("type", "attributes") |>
      dplyr::filter(.data[["type"]] %in% i) |>
      dplyr::pull(.data[["attributes"]]) |>
      parse_annot_attributes(feature_type = i, file_type = file_type)
    annot_list[[i]] <- annot_df |>
      dplyr::filter(.data[["type"]] %in% i) |>
      dplyr::select("chr", "start", "stop", "strand") |>
      dplyr::bind_cols(attr_df)
  }
  return(annot_list)
}
