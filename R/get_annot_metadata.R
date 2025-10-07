get_annot_metadata <- function(gff_file) {
  .data <- rlang::.data
  metadata_list <- list()

  line_schema <- arrow::schema(
    arrow::field("linestr", arrow::string(), nullable = TRUE)
  )

  metadata_df <- arrow::read_delim_arrow(
    gff_file,
    schema = line_schema,
    delim = "°",
    as_data_frame = FALSE
  ) |>
    dplyr::filter(stringr::str_detect(.data$linestr, "^#")) |>
    dplyr::collect()

  # Comments
  comment_df <- metadata_df |>
    dplyr::filter(stringr::str_detect(.data$linestr, "^#[^#]")) |>
    dplyr::mutate(
      linestr = stringr::str_remove(.data$linestr, "^#"),
      key = stringr::str_extract(.data$linestr, ".+?(?=:)"),
      value = stringr::str_remove(.data$linestr, ".+?:(\\s+)?")
    ) |>
    dplyr::select(-c("linestr"))

  metadata_list[["comments"]] <- comment_df

  # Pragmas
  pragma_df <- metadata_df |>
    dplyr::filter(stringr::str_detect(.data$linestr, "^##[^#]")) |>
    dplyr::mutate(
      linestr = stringr::str_remove(.data$linestr, "^##"),
      key = stringr::str_extract(.data$linestr, "^[^\\s]+"),
      value = stringr::str_remove(.data$linestr, "^[^\\s]+\\s*")
    ) |>
    dplyr::select(-c("linestr"))

  metadata_list[["pragmas"]] <- pragma_df

  metadata_list[["gff_version"]] <- pragma_df |>
    get_gff_version_metadata()
  metadata_list[["seq_region"]] <- pragma_df |>
    get_seq_region_metadata()

  return(metadata_list)
}

# Utils
get_seq_region_metadata <- function(pragma_df) {
  .data <- rlang::.data
  pragma_df |>
    dplyr::filter(.data$key == "sequence-region") |>
    tidyr::separate_wider_delim(
      cols = "value",
      delim = stringr::regex("\\s+"),
      names = c("seqid", "start_pos", "end_pos"),
      too_few = "align_start",
      too_many = "merge"
    ) |>
    dplyr::mutate(
      start_pos = as.integer(.data$start_pos),
      end_pos = as.integer(.data$end_pos)
    ) |>
    dplyr::select("seqid", "start_pos", "end_pos")
}

get_gff_version_metadata <- function(pragma_df) {
  .data <- rlang::.data
  pragma_df |>
    dplyr::filter(.data$key == "gff-version") |>
    dplyr::slice_head(n = 1) |>
    dplyr::pull("value")
}
