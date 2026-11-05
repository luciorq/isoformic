#' Plot Log2 Fold-Changes of a Gene and Its Isoforms
#'
#' @description
#' Bar plot of the log2 fold-change of one or more genes next to the
#' fold-changes of each of their transcripts, coloured by transcript type
#' and faded when not significant. With several genes, one panel per gene
#' is drawn. The [IsoformicExperiment] method uses the stored differential
#' expression results (see [de_table()]); the data frame method accepts the
#' long table produced by [de_table()] or [join_DEG_DET()].
#'
#' @param self An [IsoformicExperiment] or a long DE table.
#' @param feature Character vector of gene names (or IDs, or transcript
#'   names, matched against `feature_column`).
#' @param feature_column Column of the DE table used to select features
#'   (default `"gene_name"`).
#' @param color_palette Named colour vector per transcript type; defaults
#'   to the object's palette or [tx_type_palette()].
#' @param ncol Number of facet columns, `NULL` for automatic.
#' @param ... Passed on to methods.
#'
#' @returns A `ggplot` object.
#'
#' @examples
#' iso <- example_isoformic()
#' plot_log2fc(iso, feature = "RBPJ")
#' plot_log2fc(iso, feature = c("RBPJ", "EGFR", "PNCK"), ncol = 1)
#'
#' # Table-based
#' plot_log2fc(de_table(iso), feature = "RBPJ")
#'
#' @name plot_log2fc
#' @rdname plot_log2fc
#' @export
plot_log2fc <- S7::new_generic(
  "plot_log2fc",
  "self",
  function(
    self,
    feature,
    feature_column = "gene_name",
    color_palette = NULL,
    ncol = NULL,
    ...
  ) {
    S7::S7_dispatch()
  }
)

S7::method(plot_log2fc, IsoformicExperiment) <- function(
  self,
  feature,
  feature_column = "gene_name",
  color_palette = NULL,
  ncol = NULL,
  ...
) {
  plot_log2fc_table(
    de_data = de_table(self),
    feature = feature,
    feature_column = feature_column,
    color_palette = color_palette %||% self@tx_type_palette,
    ncol = ncol
  )
}

# nolint start: object_name_linter.
S7::method(plot_log2fc, S7::class_data.frame) <- function(
  self,
  feature,
  feature_column = "gene_name",
  color_palette = NULL,
  ncol = NULL,
  ...
) {
  plot_log2fc_table(
    de_data = self,
    feature = feature,
    feature_column = feature_column,
    color_palette = color_palette %||% tx_type_palette(),
    ncol = ncol
  )
}
# nolint end

#' @rdname plot_log2fc
#' @param de_data A long DE table (see [de_table()]).
#' @export
# nolint start: object_name_linter.
plot_log2FC <- function(
  de_data,
  feature,
  feature_column = "gene_name",
  color_palette = NULL,
  ncol = NULL
) {
  plot_log2fc(
    de_data,
    feature = feature,
    feature_column = feature_column,
    color_palette = color_palette,
    ncol = ncol
  )
}
# nolint end

#' Draw the log2FC bar plot from a long DE table
#' @keywords internal
#' @noRd
plot_log2fc_table <- function(
  de_data,
  feature,
  feature_column,
  color_palette,
  ncol
) {
  rlang::check_required(feature)
  required <- c("feature_name", "feature_type", "log2FC", "is_de")
  missing_cols <- setdiff(required, colnames(de_data))
  if (length(missing_cols)) {
    cli::cli_abort(
      message = c(
        x = "{.arg de_data} is missing column{?s}: {.field {missing_cols}}."
      ),
      class = "isoformic_invalid_de_table"
    )
  }
  if (!feature_column %in% colnames(de_data)) {
    cli::cli_abort(
      message = c(
        x = "Column {.field {feature_column}} not found in the DE table.",
        i = "Available columns: {.val {colnames(de_data)}}."
      ),
      class = "isoformic_plot_log2fc_invalid_feature_column"
    )
  }
  data <- de_data[de_data[[feature_column]] %in% feature, ]
  # Features without a fold-change (e.g. filtered before testing) cannot be drawn
  data <- data[!is.na(data$log2FC), ]
  if (!nrow(data)) {
    cli::cli_abort(
      message = c(
        x = "No features matching {.val {feature}} in column {.field {feature_column}}."
      ),
      class = "isoformic_empty_plot_data"
    )
  }
  panel_col <- if (
    feature_column %in%
      c("gene_name", "gene_id") &&
      feature_column %in% colnames(data)
  ) {
    feature_column
  } else if ("gene_name" %in% colnames(data)) {
    "gene_name"
  } else {
    feature_column
  }
  data$panel <- factor(
    data[[panel_col]],
    levels = unique(c(intersect(feature, data[[panel_col]]), data[[panel_col]]))
  )
  data$is_de <- factor(
    ifelse(is.na(data$is_de), "no", as.character(data$is_de)),
    levels = c("yes", "no")
  )
  data$feature_type <- unknown_type(data$feature_type)
  data$feature_type <- factor(
    data$feature_type,
    levels = c("gene", setdiff(unique(data$feature_type), "gene"))
  )
  data <- data[
    order(data$panel, data$feature_type != "gene", data$feature_name),
  ]
  data$feature_name <- factor(
    data$feature_name,
    levels = unique(data$feature_name)
  )
  palette <- complete_palette(color_palette, levels(data$feature_type))

  p <- ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = .data$feature_name,
      y = .data$log2FC,
      fill = .data$feature_type,
      alpha = .data$is_de
    )
  ) +
    ggplot2::geom_hline(yintercept = 0, color = "grey60") +
    ggplot2::geom_col(width = 0.75) +
    ggplot2::scale_fill_manual(values = palette, drop = TRUE) +
    ggplot2::scale_alpha_manual(values = c(yes = 1, no = 0.3), drop = FALSE) +
    ggplot2::labs(
      x = NULL,
      y = expression(log[2](fold ~ change)),
      fill = "Feature type",
      alpha = "Differentially expressed"
    ) +
    isoformic_theme() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1)
    )
  if (nlevels(data$panel) > 1L) {
    p <- p +
      ggplot2::facet_wrap(
        ggplot2::vars(.data$panel),
        scales = "free_x",
        ncol = ncol
      )
  } else {
    p <- p + ggplot2::labs(title = levels(data$panel))
  }
  p
}
