#' Plot Transcript per gene expression
#'
#' @param genes_to_plot a character vector with gene names
#'
#' @param profile_data tibble output from `prepare_profile_data`
#'
#' @return a `ggplot` object
#'
#' @export
plot_tx_expr <- function(genes_to_plot, profile_data) {
  # dependencies
  .data <- rlang::.data

  # genes_to_plot <- "ATF3"
  # profile_data <- profile_data_df

  expr_df <- profile_data

  facet_num <- length(genes_to_plot)
  facet_row_num <- round(sqrt(facet_num))

  # FIXME: this should be removed after integration with S4 object
  var <- colnames(expr_df)[2]
  var_levels <- levels(dplyr::pull(expr_df, {{ var }}))

  dodge_value <- ggplot2::position_dodge(width = 0.008)

  plot_df <- expr_df |>
    dplyr::filter(.data$parent_gene %in% genes_to_plot) |>
    dplyr::mutate(log2_mean_TPM = log2(.data$mean_TPM + 1)) |>
    dplyr::mutate(
      sd_min = dplyr::if_else(
        .data$mean_TPM - .data$SD <= 0,
        true = 0,
        false = log2((.data$mean_TPM - .data$SD) + 1)
      )
    ) |>
    dplyr::mutate(
      sd_max = log2((.data$mean_TPM + .data$SD) + 1)
    )


  transparency_vector <- c(1, 0.2)
  names(transparency_vector) <- c("Yes", "No")
  plot_object <- plot_df |>
    ggplot2::ggplot() +
    ggrepel::geom_text_repel(
      mapping = ggplot2::aes(
        x = .data[[var]], y = .data$log2_mean_TPM,
        label = dplyr::if_else(
          DE == "Yes" & .data[[var]] == var_levels[1],
          true = as.character(genename), false = ""
        )
      ),
      nudge_x = -0.2, segment.colour = "lightgray"
    ) +
    ggplot2::geom_point(
      ggplot2::aes(
        x = .data[[var]],
        y = .data$sd_max,
        fill = .data$transcript_type,
        color = .data$transcript_type
      ),
      size = 1.5, alpha = 0.4, shape = 24
    ) +
    ggplot2::geom_point(
      ggplot2::aes(
        x = .data[[var]],
        y = .data$sd_min,
        fill = .data$transcript_type,
        color = .data$transcript_type
      ),
      size = 1.5, alpha = 0.4, shape = 25
    ) +
    ggplot2::geom_line(
      mapping = ggplot2::aes(
        x = .data[[var]],
        y = log2_mean_TPM,
        group = .data$genename,
        color = .data$transcript_type,
        alpha = .data$DE
      ), size = 1.5
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        x = .data[[var]],
        ymin = .data$sd_min,
        ymax = .data$sd_max,
        group = .data$genename,
        color = .data$transcript_type
      ),
      alpha = 0.2, width = 0.01, position = dodge_value
    ) +
    ggplot2::geom_point(
      ggplot2::aes(
        x = .data[[var]],
        y = log2_mean_TPM,
        color = .data$transcript_type,
        fill = .data$transcript_type,
        alpha = .data$DE
      ),
      size = 3, shape = 21
    ) +
    ggplot2::scale_alpha_manual(values = transparency_vector) +
    ggplot2::scale_color_manual(
      values = tx_type_palette(),
      aesthetics = c("color", "fill")
    ) +
    # ggplot2::theme_bw()
    ggpubr::theme_pubr() +
    ggplot2::theme(legend.position = "right") +
    # ggplot2::facet_grid(rows = vars(transcript_name), cols = vars(genename))
    ggplot2::facet_wrap(
      facets = ~parent_gene,
      nrow = facet_row_num,
      scales = "free_y"
    )

  plot_object <- plot_object +
    ggplot2::labs(
      x = "Condition",
      y = "log2(mean(TPM) + 1)",
      color = "Transcript Biotype",
      fill = "Transcript Biotype",
      alpha = "Differentially expressed"
    )
  return(plot_object)
}
