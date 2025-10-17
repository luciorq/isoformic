#' Plot Log2 Fold-Change Results for Transcripts of Selected Genes
#'
#' Creates a bar plot of log2 fold-change values for transcripts of a selected gene,
#' differentiating transcript types and significance levels.
#'
#' @param de_data A `data.frame` or `tibble` containing combined gene and transcript differential expression results.
#'  Should contain columns for log2 fold-change, transcript type, significance, and feature symbols.
#' @param feature A character string specifying the gene name to plot.
#' @param feature_column A character string specifying the column name in `de_data` that contains gene names.
#'  Default is `"gene_name"`.
#' @param color_palette A named character vector specifying colors for different transcript types.
#'  If `NULL`, a default palette will be used.
#'
#' @return A `ggplot2` object representing the bar plot.
#'
#' @details
#' The function filters the input table for the selected gene and creates a bar plot of log2 fold-change values.
#' If all transcripts are significant, it plots without adjusting alpha transparency; otherwise, it adjusts alpha
#' based on significance. The function uses predefined colors for transcript types, which can be overridden
#' by providing `custom_colors`.
#'
#' @examples
#' # Sample data
#' de_table_long <- data.frame(
#'   feature_name = c("Transcript1", "Transcript2", "Transcript3", "GeneA"),
#'   feature_id = c("TX1", "TX2", "TX3", "GENEA"),
#'   gene_name = c("GeneA", "GeneA", "GeneA", "GeneA"),
#'   log2FC = c(1.5, -0.5, -2.0, 0.8),
#'   feature_type = c("protein_coding", "lncRNA", "retained_intron", "gene"),
#'   is_de = c("yes", "no", "yes", "yes")
#' )
#'
#' # Plot log2 fold-change for the selected gene
#' plot_obj <- plot_log2FC(
#'   de_data = de_table_long,
#'   feature = "GeneA",
#'   feature_column = "gene_name"
#' )
#'
#' # Display the plot
#' print(plot_obj)
#'
#' @export
plot_log2FC <- function(
  de_data,
  feature,
  feature_column = "gene_name",
  color_palette = NULL
) {
  .data <- rlang::.data
  if (isFALSE(feature_column %in% colnames(de_data))) {
    cli::cli_abort(
      message = c(
        x = "The specified {.var feature_column} '{feature_column}' is not found in the provided data.",
        i = "Please ensure the column name is correct and exists in the data."
      ),
      class = "isoformic_plot_log2fc_invalid_feature_column"
    )
  }
  if (rlang::is_null(color_palette)) {
    color_palette <- tx_type_palette()
  }
  x_axis_label_column <- stringr::str_replace(
    string = feature_column,
    pattern = stringr::regex("^gene_|^transcript_"),
    replacement = "feature_"
  )

  all_feature_names <- de_data[[x_axis_label_column]][
    de_data[[feature_column]] %in% feature
  ] |>
    unique()

  selected_feature_types <- de_data |>
    dplyr::filter(.data[[x_axis_label_column]] %in% all_feature_names) |>
    dplyr::pull("feature_type") |>
    unique()

  missing_types <- selected_feature_types[
    !selected_feature_types %in% unique(names(color_palette))
  ]

  if (isTRUE(length(missing_types) > 0)) {
    cli::cli_abort(
      message = c(
        x = "The provided {.var color_palette} is missing colors for the following transcript types: {paste(missing_types, collapse = ', ')}.",
        i = "Please provide colors for all transcript types present in the data."
      ),
      class = "isoformic_plot_log2fc_incomplete_color_palette"
    )
  }

  plot_obj <- ggplot2::ggplot(
    data = de_data[de_data[[feature_column]] %in% feature, ],
    mapping = ggplot2::aes(
      x = .data[[x_axis_label_column]],
      y = .data$log2FC,
      fill = .data$feature_type
    )
  )

  all_features_de <- de_data$is_de[de_data[[x_axis_label_column]] %in% feature]
  if (all(all_features_de == "yes")) {
    plot_obj <- plot_obj +
      ggplot2::geom_bar(stat = "identity")
  } else {
    plot_obj <- plot_obj +
      ggplot2::geom_bar(
        mapping = ggplot2::aes(
          x = .data[[x_axis_label_column]],
          y = .data$log2FC,
          fill = .data$feature_type,
          alpha = factor(
            .data$is_de,
            levels = c("no", "yes")
          )
        ),
        stat = "identity"
      ) +
      ggplot2::scale_alpha_manual(
        values = c(0.2, 1),
        labels = c("No", "Yes")
      ) +
      ggplot2::labs(alpha = "Differentially Expressed")
  }

  plot_obj <- plot_obj +
    ggplot2::geom_hline(yintercept = 0, linetype = "solid", color = "grey") +
    ggplot2::scale_fill_manual(values = color_palette) +
    ggplot2::labs(fill = "Feature Type") +
    ggplot2::xlab(stringr::str_to_title(stringr::str_replace(
      string = x_axis_label_column,
      pattern = "_",
      replacement = " "
    ))) +
    ggplot2::ylab(expression(~ log[2](Fold ~ Change))) +
    ggplot2::theme_bw()

  return(plot_obj)
}

# S7 Methods

plot_log2fc <- S7::new_generic("plot_log2fc", "self")

S7::method(plot_log2fc, S7::class_data.frame) <- function(
  self,
  feature,
  feature_column = "gene_name",
  color_palette = NULL
) {
  plot_obj <- plot_log2FC(
    de_data = self,
    feature = feature,
    feature_column = feature_column,
    color_palette = color_palette
  )
  return(plot_obj)
}

S7::method(plot_log2fc, IsoformicExperiment) <- function(
  self,
  feature,
  feature_column = "gene_name"
) {
  if (
    isTRUE(rlang::is_null(self@dea[["deg"]])) ||
      isTRUE(rlang::is_null(self@dea[["det"]]))
  ) {
    cli::cli_abort(
      message = c(
        x = "Differential expression results not found in {.cls IsoformicExperiment} object.",
        i = "Please run differential expression analysis and combine results before plotting."
      ),
      class = "isoformic_missing_dea_results"
    )
  }
  .data <- rlang::.data
  deg_det_table <- combine_deg_det_longer(self)
  plot_obj <- plot_log2FC(
    de_data = deg_det_table,
    feature = feature,
    feature_column = feature_column,
    color_palette = self@tx_type_palette
  )
  return(plot_obj)
}

combine_deg_det_wider <- function(isoformic_obj) {
  dplyr::left_join(
    de_tx(isoformic_obj),
    tx_annot(isoformic_obj),
    by = "transcript_id"
  ) |>
    dplyr::distinct() |>
    dplyr::left_join(
      de_gene(isoformic_obj),
      by = "gene_id",
      suffix = c("_tx", "_gene")
    )
}

combine_deg_det_longer <- function(isoformic_obj) {
  .data <- rlang::.data
  deg_df <- de_gene(isoformic_obj) |>
    dplyr::left_join(
      isoformic_obj@annot_data_genes |>
        dplyr::select(
          dplyr::any_of(c("gene_id", "gene_name"))
        ) |>
        dplyr::distinct() |>
        dplyr::collect(),
      by = "gene_id"
    ) |>
    dplyr::distinct() |>
    dplyr::mutate(
      feature_id = .data$gene_id,
      feature_name = .data$gene_name
    ) |>
    dplyr::rename(is_de = "deg_sig") |>
    dplyr::mutate(feature_type = "gene") |>
    dplyr::select(
      dplyr::any_of(c(
        "feature_id",
        "feature_name",
        "gene_id",
        "gene_name",
        "feature_type",
        "log2FC",
        "pvalue",
        "qvalue",
        "is_de"
      ))
    )
  det_df <- de_tx(isoformic_obj) |>
    dplyr::left_join(
      isoformic_obj@annot_data_transcripts |>
        dplyr::select(
          dplyr::any_of(
            c("transcript_id", "transcript_name", "transcript_type", "gene_id")
          )
        ) |>
        dplyr::distinct() |>
        dplyr::collect(),
      by = "transcript_id"
    ) |>
    dplyr::distinct() |>
    dplyr::left_join(
      isoformic_obj@annot_data_genes |>
        dplyr::select(
          dplyr::any_of(c("gene_id", "gene_name"))
        ) |>
        dplyr::distinct() |>
        dplyr::collect(),
      by = "gene_id"
    ) |>
    dplyr::distinct() |>
    # dplyr::mutate(
    #  feature_id = .data$transcript_id,
    #  feature_name = .data$transcript_name
    # ) |>
    dplyr::rename(
      feature_id = "transcript_id",
      feature_name = "transcript_name",
      feature_type = "transcript_type",
      is_de = "det_sig"
    ) |>
    dplyr::select(
      dplyr::any_of(c(
        "feature_id",
        "feature_name",
        "gene_id",
        "gene_name",
        "feature_type",
        "log2FC",
        "pvalue",
        "qvalue",
        "is_de"
      ))
    )

  dplyr::bind_rows(deg_df, det_df)
}
