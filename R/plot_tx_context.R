#' Plot Genomic Context From GTF Annotation File
#' @param exon_table Exon table generated from `prepare_exon_annotation`
#' @export
plot_tx_context <- function(exon_table) {
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

  for (tx_id in tx_id_vector) {
    tx_id_data <- plot_data[plot_data$tx_id %in% tx_id, ]
    exon_right_max <- max(tx_id_data$exon_right, na.rm = TRUE)
    plot_data[plot_data$tx_id %in% tx_id & plot_data$exon_right == exon_right_max, "segment_start"] <- NA_real_
    plot_data[plot_data$tx_id %in% tx_id & plot_data$exon_right == exon_right_max, "segment_middle"] <- NA_real_
    plot_data[plot_data$tx_id %in% tx_id & plot_data$exon_right == exon_right_max, "segment_end"] <- NA_real_
  }

  plot_data |>
    dplyr::arrange(tx_id) |>
    ggplot2::ggplot() +
    ggplot2::geom_rect(
      mapping = ggplot2::aes(
        xmin = exon_left,
        xmax = exon_right,
        ymin = 0,
        ymax = tpm
      )
    ) +
    ggplot2::geom_segment(
      mapping = ggplot2::aes(
        x = segment_start,
        xend = segment_middle,
        y = 0,
        yend = tpm
      )
    ) +
    ggplot2::geom_segment(
      mapping = ggplot2::aes(
        x = segment_middle,
        xend = segment_end,
        y = tpm,
        yend = 0
      )
    ) +
    ggplot2::facet_wrap(
      ggplot2::vars(tx_id),
      ncol = 1
    ) +
    ggplot2::labs(
      x = "Genomic Coordinate",
      y = "TPM"
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      strip.background = ggplot2::element_blank(),
      strip.text.x = ggplot2::element_blank()
    )
}
