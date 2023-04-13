plot_tx_context <- function(exon_table, custom_colors = NULL) {
  
  exon_table <- dplyr::mutate(exon_table, exon_left = as.numeric(exon_left), 
                              exon_right = as.numeric(exon_right))
  
  exon_table <- dplyr::ungroup(dplyr::arrange(dplyr::group_by(exon_table, tx_id), tx_id, exon_right))
  
  segment_end_vector <- exon_table$exon_left
  segment_end_vector <- segment_end_vector[2:length(segment_end_vector)]
  segment_end_vector <- c(segment_end_vector, NA)
  
  plot_data <- dplyr::mutate(dplyr::mutate(dplyr::mutate(exon_table, 
                                                         segment_start = exon_right, 
                                                         segment_end = segment_end_vector), 
                                           segment_middle = floor((segment_end + segment_start)/2)), 
                             tpm = 5)
  
  plot_data[nrow(plot_data), "segment_start"] <- NA_real_
  
  tx_id_vector <- plot_data$tx_id
  
  DEG_DET_table_nogene <- DEG_DET_table %>% filter(transcript_type != "gene")
  DEG_DET_table_nogene$transcript_type <- as.factor(DEG_DET_table_nogene$transcript_type)
  
  DEG_DET_table$transcript_type <- factor(DEG_DET_table$transcript_type, levels = unique(DEG_DET_table$transcript_type))
  n_levels <- length(levels(DEG_DET_table$transcript_type))
  n_colors <- if (!is.null(custom_colors)) length(custom_colors) else n_levels
  if (n_colors > n_levels) {
    n_colors <- n_levels
  }
  my_colors <- if (!is.null(custom_colors)) {
    custom_colors[1:n_colors]
  } else {
    c("#F8766D","#C77CFF", "#00BFC4", "#CD9600", "#7CAE00", "#8494FF", "#00A9FF", 
               "#FF61CC", "#0CB702", "#E68613", "#00C19A", "#ABA300", "#FF68A1")[1:n_colors]
  }
  
  for (tx_id in tx_id_vector) {
    tx_id_data <- plot_data[plot_data$tx_id %in% tx_id,]
    exon_right_max <- max(tx_id_data$exon_right, na.rm = TRUE)
    plot_data[plot_data$tx_id %in% tx_id & plot_data$exon_right == exon_right_max, "segment_start"] <- NA_real_
    plot_data[plot_data$tx_id %in% tx_id & plot_data$exon_right == exon_right_max, "segment_middle"] <- NA_real_
    plot_data[plot_data$tx_id %in% tx_id & plot_data$exon_right == exon_right_max, "segment_end"] <- NA_real_
  }
  
  tx_id_to_name <- select(tx_to_gene, transcript_id, transcript_name, transcript_type)
  
  left_join(dplyr::arrange(plot_data, tx_id), tx_id_to_name, by = c(tx_id = "transcript_id")) %>% 
    ggplot2::ggplot() + 
    ggplot2::geom_rect(mapping = ggplot2::aes(xmin = exon_left, xmax = exon_right, ymin = 0, ymax = tpm, fill = transcript_type), 
                       color = NA) + 
    ggplot2::scale_fill_manual(values = my_colors, 
                               name = "Transcript type",
                               labels = levels(DEG_DET_table$transcript_type),
                               breaks = levels(DEG_DET_table$transcript_type)) +
    ggplot2::scale_color_manual(values = my_colors, 
                                name = "Transcript type",
                                labels = levels(DEG_DET_table$transcript_type),
                                breaks = levels(DEG_DET_table$transcript_type)) +
    ggplot2::geom_segment(mapping = ggplot2::aes(x = segment_start, xend = segment_middle, y = 0, yend = tpm, col = transcript_type, alpha = 0.3)) + 
    ggplot2::geom_segment(mapping = ggplot2::aes(x = segment_middle, xend = segment_end, y = tpm, yend = 0, col = transcript_type, alpha = 0.3)) + 
    ggplot2::facet_wrap(ggplot2::vars(transcript_name), ncol = 1) + 
    ggplot2::labs(x = "Genomic Coordinate", y = "TPM") + 
    ggplot2::theme_classic()
}
