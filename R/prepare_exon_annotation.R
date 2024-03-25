#' Prepare Exon based Position Annotation Table
#' @param gene_names String or vector of gene names to extract.
#' @param file_path Path to annotation file.
#' @inheritParams download_reference
#' @export
prepare_exon_annotation <- function(gene_name,
                                    file_path,
                                    file_type = "gff") {
  if (isTRUE(file_type == "gff")) {
    gff_file_path <- fs::path(file_path)
    annot_df <- readr::read_table(
      file = gff_file_path,
      col_names = FALSE,
      col_types = readr::cols(
        X1 = readr::col_skip(),
        X2 = readr::col_skip(),
        X3 = "c",
        X4 = "c",
        X5 = "c",
        X6 = readr::col_skip(),
        X7 = "c",
        X8 = readr::col_skip(),
        X9 = "c"
      ),
      comment = "#",
      progress = FALSE
    )
  }
  gene_annot_df <- tibble::tibble()
  for (i in gene_name) {
    gene_annot_df_temp <- annot_df |>
      dplyr::filter(
        stringr::str_detect(X9, paste0("gene_name=", i, ";"))
      )
    gene_annot_df <- dplyr::bind_rows(gene_annot_df, gene_annot_df_temp)
  }


  tx_id_vector <- gene_annot_df |>
    dplyr::filter(X3 %in% "transcript") |>
    dplyr::select(X9) |>
    dplyr::mutate(
      tx_id = stringr::str_extract(X9, "transcript_id=.*?;") |>
        stringr::str_remove("^transcript_id=") |>
        stringr::str_remove(";")
    ) |>
    dplyr::select(tx_id) |>
    dplyr::arrange() |>
    dplyr::distinct() |>
    dplyr::pull(tx_id)

  # parent_id <- tx_id_vector[1]
  tx_exon_table <- tx_id_vector |>
    purrr::map_dfr(
      .f = \(parent_id) {
        gene_annot_df |>
          dplyr::filter(
            X3 %in% "exon",
            stringr::str_detect(X9, glue::glue("Parent={parent_id};"))
          ) |>
          dplyr::mutate(
            tx_name = stringr::str_extract(X9, "gene_name=.*?;") |>
              stringr::str_remove("^gene_name=") |>
              stringr::str_remove(";")
          ) |>
          dplyr::select(c(X4, X5, X7, tx_name)) |>
          dplyr::mutate(tx_id = parent_id) |>
          dplyr::relocate(tx_name, .after = tx_id)
      }
    ) #|>
  # dplyr::mutate(tx_name = gene_name)
  tx_exon_table <- tx_exon_table |>
    dplyr::rename(
      exon_left = X4,
      exon_right = X5,
      strand = X7
    )
  return(tx_exon_table)
}
