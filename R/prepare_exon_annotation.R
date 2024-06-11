#' Prepare Exon based Position Annotation Table
#' @param gene_names String or vector of gene names to extract.
#' @param file_path Path to annotation file.
#' @inheritParams download_reference
#' @export
prepare_exon_annotation <- function(gene_name,
                                    file_path,
                                    file_type = c("gff")) {
  file_type <- stringr::str_to_lower(file_type)
  file_type <- rlang::arg_match(file_type)
  if (isTRUE(file_type %in% c("gff"))) {
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

  colnames(annot_df) <- c("chr", "type", "start", "stop", "strand", "attributes")

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

parse_annot_attributes <- function(column_vector,
                                   feature_type = c("gene", "transcript", "exon"),
                                   file_type = c("gtf", "gff")) {
  feature_type <- rlang::arg_match(feature_type)
  file_type <- rlang::arg_match(file_type)

  if (isTRUE(feature_type == "gene")) {
    gene_id_vector <- extract_attribute(column_vector, "gene_id", file_type)
    gene_name_vector <- extract_attribute(column_vector, "gene_name", file_type)
    gene_type_vector <- extract_attribute(column_vector, "gene_type", file_type)
    attr_df <- tibble::tibble(
      gene_id = gene_id_vector,
      gene_name = gene_name_vector,
      gene_type = gene_type_vector
    )
  } else if (isTRUE(feature_type == "exon")) {
    exon_id_vector <- extract_attribute(column_vector, "exon_id", file_type)
    exon_number_vector <- extract_attribute(column_vector, "exon_number", file_type)
    transcript_id_vector <- extract_attribute(column_vector, "transcript_id", file_type)
    attr_df <- tibble::tibble(
      exon_id = exon_id_vector,
      exon_number = exon_number_vector,
      transcript_id = transcript_id_vector
    )
  } else if (isTRUE(feature_type == "transcript")) {
    transcript_id_vector <- extract_attribute(column_vector, "transcript_id", file_type)
    transcript_name_vector <- extract_attribute(column_vector, "transcript_name", file_type)
    transcript_type_vector <- extract_attribute(column_vector, "transcript_type", file_type)
    gene_id_vector <- extract_attribute(column_vector, "gene_id", file_type)
    attr_df <- tibble::tibble(
      transcript_id = transcript_id_vector,
      transcript_name = transcript_name_vector,
      transcript_type = transcript_type_vector,
      gene_id = gene_id_vector
    )
  }
  return(attr_df)
}

extract_attribute <- function(column_vector, attr_name, file_type) {
  if (isTRUE(file_type == "gtf")) {
    attr_vector <- column_vector |>
      stringr::str_extract(pattern = paste0(attr_name, ".*?;")) |>
      stringr::str_remove(paste0("^", attr_name, " \"")) |>
      stringr::str_remove("\";$") |>
      stringr::str_squish()
  } else if (isTRUE(file_type == "gff")) {
    attr_vector <- column_vector |>
      stringr::str_extract(pattern = paste0(attr_name, "=.*?;")) |>
      stringr::str_remove(paste0("^", attr_name, "=")) |>
      stringr::str_remove(";$") |>
      stringr::str_squish()
  }
  return(attr_vector)
}
