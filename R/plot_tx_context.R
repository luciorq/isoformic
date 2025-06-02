#' Plot Transcript Genomic Context
#'
#' This function plots the genomic context of all transcripts of given genes.
#'
#' @param exon_table a tibble with exon information.
#'   Must contain columns `tx_id`, `exon_left`, and `exon_right`.
#'
#' @param custom_colors a vector of colors to use for each transcript. If not
#'    provided, the function will use the default colors. Actually, this
#'    argument is ***NOT implemented** yet.
#'
#' @export
plot_tx_context <- function(
    exon_table,
    custom_colors = NULL) {
  if (isFALSE(tibble::is_tibble(exon_table))) {
    exon_table <- tibble::as_tibble(exon_table)
  }
  .data <- rlang::.data
  .env <- rlang::.env
  exon_table <- exon_table |>
    dplyr::mutate(
      exon_left = as.integer(.data$exon_left),
      exon_right = as.integer(.data$exon_right)
    )

  exon_table <- exon_table |>
    dplyr::group_by(tx_id) |>
    dplyr::arrange(.data$tx_id, .data$exon_right) |>
    dplyr::ungroup()

  # pos_left <- min(exon_table$exon_left, exon_table$exon_right)
  segment_end_vector <- exon_table$exon_left
  segment_end_vector <- segment_end_vector[2:length(segment_end_vector)]
  segment_end_vector <- c(segment_end_vector, NA_integer_)
  plot_data <- exon_table |>
    dplyr::mutate(
      # new_left = segment_left - pos_left,
      # new_right = segment_right - pos_left,
      segment_start = .data$exon_right,
      segment_end = .env$segment_end_vector
    ) |>
    dplyr::mutate(
      segment_middle = as.integer(
        floor((.data$segment_end + .data$segment_start) / 2L)
      )
    ) |>
    dplyr::mutate(tpm = 5L)

  plot_data[nrow(plot_data), "segment_start"] <- NA_integer_

  tx_id_vector <- plot_data$tx_id

  # TODO: @luciorq pass DEG_DET_table as an argument, or remove if not needed
  # DEG_DET_table_nogene <- DEG_DET_table |>
  #   dplyr::filter(transcript_type != "gene")
  # DEG_DET_table_nogene$transcript_type <- as.factor(DEG_DET_table_nogene$transcript_type)

  # tx_id = tx_id_vector[1]
  for (tx_id in tx_id_vector) {
    tx_id_data <- plot_data[plot_data$tx_id %in% tx_id, ]
    exon_right_max <- base::max(tx_id_data$exon_right, na.rm = TRUE)
    plot_data[plot_data$tx_id %in% tx_id & plot_data$exon_right == exon_right_max, "segment_start"] <- NA_integer_
    plot_data[plot_data$tx_id %in% tx_id & plot_data$exon_right == exon_right_max, "segment_middle"] <- NA_integer_
    plot_data[plot_data$tx_id %in% tx_id & plot_data$exon_right == exon_right_max, "segment_end"] <- NA_integer_
  }

  # TODO: @luciorq tx_to_gene should be passed as an argument?
  tx_id_to_name <- tx_to_gene |>
    dplyr::select(transcript_id, transcript_name, transcript_type)

  plot_data |>
    dplyr::arrange(tx_id) |>
    dplyr::left_join(tx_id_to_name, by = c("tx_id" = "transcript_id")) |>
    ggplot2::ggplot() +
    ggplot2::geom_rect(
      mapping = ggplot2::aes(
        xmin = .data$exon_left,
        xmax = .data$exon_right,
        ymin = 0L,
        ymax = .data$tpm,
        fill = .data$transcript_type
      ) # ,
      # color = NA
    ) +
    ggplot2::scale_fill_manual(values = tx_type_palette()) +
    ggplot2::scale_color_manual(values = tx_type_palette()) +
    ggplot2::geom_segment(
      mapping = ggplot2::aes(
        x = .data$segment_start,
        xend = .data$segment_middle,
        y = 0L,
        yend = .data$tpm,
        col = .data$transcript_type,
        alpha = 0.3
      ),
      na.rm = TRUE
    ) +
    ggplot2::geom_segment(
      mapping = ggplot2::aes(
        x = .data$segment_middle,
        xend = .data$segment_end,
        y = .data$tpm,
        yend = 0L,
        col = .data$transcript_type,
        alpha = 0.3
      ),
      na.rm = TRUE
    ) +
    ggplot2::facet_wrap(
      ggplot2::vars(.data$transcript_name),
      ncol = 1
    ) +
    ggplot2::labs(
      x = "Genomic Coordinate",
      y = "TPM"
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.title.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.line.y = ggplot2::element_blank()
    )
}
