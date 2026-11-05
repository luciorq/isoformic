#' Plot the uncertainty of isoform-level results
#'
#' @description
#' * `plot_de_confidence(x, gene_name)` shows, for the isoforms of the given
#'   genes, the median log2 fold-change over the Salmon replicates with its 95%
#'   interval, coloured by the confidence class of [de_confidence()] and with
#'   the shape marking the calls of the primary differential expression
#'   test. An interval that crosses zero, a faded colour or an ambiguous
#'   isoform tell the reader not to trust the call.
#' * `plot_isoform_correlation(x, gene_name)` draws the correlation between
#'   the replicates of every pair of isoforms of a gene. Strongly negative
#'   cells (outlined when ambiguous) mark isoforms whose reads cannot be told
#'   apart.
#' * `plot_inf_rv(x)` shows the distribution of the inferential relative
#'   variance per transcript type with the cutoff used to call isoforms
#'   uncertain.
#' * `confidence_palette()` gives the colours of the confidence classes.
#'
#' @param self An [IsoformicExperiment] with uncertainty tables (see
#'   [add_uncertainty()]).
#' @param gene_name Gene names or IDs to plot.
#' @param ncol Number of panel columns.
#' @param top_types Number of the most frequent transcript types to show.
#' @param ... Unused.
#'
#' @returns A `ggplot` object, or a named character vector for
#'   `confidence_palette()`.
#'
#' @seealso [add_uncertainty()], [de_confidence()]
#'
#' @name uncertainty-plots
#' @aliases plot_de_confidence plot_isoform_correlation plot_inf_rv
#' @aliases confidence_palette
NULL

#' @rdname uncertainty-plots
#' @export
confidence_palette <- function() {
  c(high = "#009E73", moderate = "#E69F00", low = "#D55E00")
}

#' @rdname uncertainty-plots
#' @export
plot_de_confidence <- S7::new_generic(
  "plot_de_confidence",
  "self",
  function(self, gene_name, ncol = 1L, ...) {
    S7::S7_dispatch()
  }
)

S7::method(plot_de_confidence, IsoformicExperiment) <- function(
  self,
  gene_name,
  ncol = 1L,
  ...
) {
  conf <- de_confidence(self)
  conf <- conf[conf$gene_name %in% gene_name | conf$gene_id %in% gene_name, ]
  if (!nrow(conf)) {
    cli::cli_abort(
      message = c(x = "No assessed isoforms for {.val {gene_name}}."),
      class = "isoformic_no_uncertainty"
    )
  }
  conf$label <- ifelse(
    conf$feature_type == "group",
    paste0(conf$feature_id, " (", conf$transcript_ids, ")"),
    conf$transcript_name
  )
  conf$label[is.na(conf$label)] <- conf$feature_id[is.na(conf$label)]
  conf$called <- ifelse(
    is.na(conf$n_det_sig),
    "not tested",
    ifelse(conf$n_det_sig > 0L, "significant", "not significant")
  )
  conf$label <- stats::reorder(conf$label, conf$lfc_median)
  ggplot2::ggplot(
    conf,
    ggplot2::aes(
      x = .data$lfc_median,
      y = .data$label,
      colour = .data$confidence,
      shape = .data$called
    )
  ) +
    ggplot2::geom_vline(
      xintercept = 0,
      linetype = "dashed",
      colour = "grey50"
    ) +
    ggplot2::geom_linerange(
      ggplot2::aes(xmin = .data$lfc_lo, xmax = .data$lfc_hi),
      orientation = "y",
      linewidth = 0.8
    ) +
    ggplot2::geom_point(size = 2.6) +
    ggplot2::scale_colour_manual(
      values = confidence_palette(),
      na.value = "grey60",
      drop = FALSE
    ) +
    ggplot2::scale_shape_manual(
      values = c(significant = 16, `not significant` = 1, `not tested` = 4)
    ) +
    ggplot2::facet_wrap(
      ggplot2::vars(.data$gene_name),
      ncol = ncol,
      scales = "free_y"
    ) +
    ggplot2::labs(
      x = "log2 fold-change (median and 95% interval over replicates)",
      y = NULL,
      colour = "Confidence",
      shape = "Primary test"
    ) +
    ggplot2::theme_minimal()
}

#' @rdname uncertainty-plots
#' @export
plot_isoform_correlation <- S7::new_generic(
  "plot_isoform_correlation",
  "self",
  function(self, gene_name, ncol = 1L, ...) {
    S7::S7_dispatch()
  }
)

S7::method(plot_isoform_correlation, IsoformicExperiment) <- function(
  self,
  gene_name,
  ncol = 1L,
  ...
) {
  pairs <- isoform_ambiguity(self, gene_name = gene_name)
  if (!nrow(pairs)) {
    cli::cli_abort(
      message = c(
        x = "No isoform pairs with computed correlations for {.val {gene_name}}."
      ),
      class = "isoformic_no_uncertainty"
    )
  }
  name_of <- function(name, id) ifelse(is.na(name), id, name)
  one <- data.frame(
    gene_name = pairs$gene_name,
    x = name_of(pairs$transcript_name_1, pairs$transcript_id_1),
    y = name_of(pairs$transcript_name_2, pairs$transcript_id_2),
    mean_cor = pairs$mean_cor,
    ambiguous = pairs$ambiguous
  )
  swapped <- one
  swapped$x <- one$y
  swapped$y <- one$x
  labels <- sort(unique(c(one$x, one$y)))
  genes_of <- unique(rbind(
    data.frame(gene_name = one$gene_name, label = one$x),
    data.frame(gene_name = one$gene_name, label = one$y)
  ))
  diag_rows <- data.frame(
    gene_name = genes_of$gene_name,
    x = genes_of$label,
    y = genes_of$label,
    mean_cor = 1,
    ambiguous = FALSE
  )
  d <- rbind(one, swapped, diag_rows)
  d$x <- factor(d$x, levels = labels)
  d$y <- factor(d$y, levels = rev(labels))
  ggplot2::ggplot(
    d,
    ggplot2::aes(x = .data$x, y = .data$y, fill = .data$mean_cor)
  ) +
    ggplot2::geom_tile(
      ggplot2::aes(colour = .data$ambiguous),
      linewidth = 0.8,
      width = 0.95,
      height = 0.95
    ) +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.2f", .data$mean_cor)),
      size = 3
    ) +
    ggplot2::scale_fill_gradient2(
      low = "#D55E00",
      mid = "white",
      high = "#0072B2",
      midpoint = 0,
      limits = c(-1, 1)
    ) +
    ggplot2::scale_colour_manual(
      values = c(`TRUE` = "black", `FALSE` = NA),
      guide = "none"
    ) +
    ggplot2::facet_wrap(
      ggplot2::vars(.data$gene_name),
      ncol = ncol,
      scales = "free"
    ) +
    ggplot2::labs(
      x = NULL,
      y = NULL,
      fill = "Replicate\ncorrelation",
      caption = "Outlined: ambiguous pair"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}

#' @rdname uncertainty-plots
#' @export
plot_inf_rv <- S7::new_generic(
  "plot_inf_rv",
  "self",
  function(self, top_types = 6L, ...) {
    S7::S7_dispatch()
  }
)

S7::method(plot_inf_rv, IsoformicExperiment) <- function(
  self,
  top_types = 6L,
  ...
) {
  unc <- tx_uncertainty(self)
  unc <- unc[unc$assessed %in% TRUE & !is.na(unc$mean_inf_rv), ]
  if (!nrow(unc)) {
    cli::cli_abort(
      message = c(x = "No assessed transcripts."),
      class = "isoformic_no_uncertainty"
    )
  }
  counts <- sort(table(unc$transcript_type), decreasing = TRUE)
  keep <- names(counts)[seq_len(min(top_types, length(counts)))]
  unc <- unc[unc$transcript_type %in% keep, ]
  unc$transcript_type <- factor(unc$transcript_type, levels = keep)
  cutoff <- uncertainty_settings(self)$inf_rv_cutoff
  ggplot2::ggplot(unc, ggplot2::aes(x = .data$mean_inf_rv)) +
    ggplot2::geom_histogram(bins = 50, fill = "#8da0cb") +
    ggplot2::geom_vline(xintercept = cutoff, linetype = "dashed") +
    ggplot2::scale_x_log10() +
    ggplot2::facet_wrap(
      ggplot2::vars(.data$transcript_type),
      scales = "free_y"
    ) +
    ggplot2::labs(
      x = "Mean inferential relative variance (log scale)",
      y = "Transcripts",
      caption = paste0("Dashed: uncertain from ", cutoff)
    ) +
    ggplot2::theme_minimal()
}
