#' Plot Transcript Genomic Context
#' @export
plot_tx_context <- function(exon_table, custom_colors = NULL) {
  exon_table <- exon_table |>
    dplyr::mutate(
      exon_left = as.numeric(exon_left),
      exon_right = as.numeric(exon_right)
    )

  exon_table <- exon_table |>
    dplyr::group_by(tx_id) |>
    dplyr::arrange(tx_id, exon_right) |>
    dplyr::ungroup()

  # pos_left <- min(exon_table$exon_left, exon_table$exon_right)
  segment_end_vector <- exon_table$exon_left
  segment_end_vector <- segment_end_vector[2:length(segment_end_vector)]
  segment_end_vector <- c(segment_end_vector, NA)
  plot_data <- exon_table |>
    dplyr::mutate(
      # new_left = segment_left - pos_left,
      # new_right = segment_right - pos_left,
      segment_start = exon_right,
      segment_end = segment_end_vector
    ) |>
    dplyr::mutate(
      segment_middle = floor((segment_end + segment_start) / 2)
    ) |>
    dplyr::mutate(tpm = 5)

  plot_data[nrow(plot_data),"segment_start"] <- NA_real_

  tx_id_vector <- plot_data$tx_id

  DEG_DET_table_nogene <-DEG_DET_table |> filter(transcript_type != "gene")
  DEG_DET_table_nogene$transcript_type <- as.factor(DEG_DET_table_nogene$transcript_type)
  
  for (tx_id in tx_id_vector) {
    tx_id_data <- plot_data[plot_data$tx_id %in% tx_id, ]
    exon_right_max <- max(tx_id_data$exon_right, na.rm = TRUE)
    plot_data[plot_data$tx_id %in% tx_id & plot_data$exon_right == exon_right_max, "segment_start"] <- NA_real_
    plot_data[plot_data$tx_id %in% tx_id & plot_data$exon_right == exon_right_max, "segment_middle"] <- NA_real_
    plot_data[plot_data$tx_id %in% tx_id & plot_data$exon_right == exon_right_max, "segment_end"] <- NA_real_
  }

tx_id_to_name <-tx_to_gene |>
  dplyr::select(transcript_id, transcript_name, transcript_type)

  plot_data |>
    dplyr::arrange(tx_id) |>
    dplyr::left_join(tx_id_to_name, by = c("tx_id" = "transcript_id")) |>
    ggplot2::ggplot() +
    ggplot2::geom_rect(
      mapping = ggplot2::aes(
        xmin = exon_left,
        xmax = exon_right,
        ymin = 0,
        ymax = tpm,
        fill = transcript_type
      ),
      color = NA
    ) +
    ggplot2::scale_fill_manual(values = tx_type_color_names) +
    ggplot2::scale_color_manual(values = tx_type_color_names)+
    ggplot2::geom_segment(
      mapping = ggplot2::aes(
        x = segment_start,
        xend = segment_middle,
        y = 0,
        yend = tpm,
        col = transcript_type,
        alpha = 0.3
      )
    ) +
    ggplot2::geom_segment(
      mapping = ggplot2::aes(
        x = segment_middle,
        xend = segment_end,
        y = tpm,
        yend = 0,
        col = transcript_type,
        alpha = 0.3
      )
    ) +
    ggplot2::facet_wrap(
      ggplot2::vars(transcript_name),
      ncol = 1
    ) +
    ggplot2::labs(
      x = "Genomic Coordinate",
      y = "TPM"
    ) +
    ggplot2::theme_classic()
}
