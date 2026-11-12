# Expression profiles ---------------------------------------------------------

#' Expression Profiles of Genes and Their Isoforms Across Groups
#'
#' @description
#' `profile_data()` summarises the expression of selected genes and of each
#' of their isoforms per sample group (mean and standard deviation), and
#' attaches the differential expression calls. `plot_tx_expr()` draws the
#' result: one panel per gene, one line per feature (the gene plus its
#' isoforms), coloured by transcript type and faded when not differentially
#' expressed. The table-based [prepare_profile_data()] produces the same
#' table from loose tables.
#'
#' @param self An [IsoformicExperiment], or (for `plot_tx_expr()`) a
#'   profile table produced by `profile_data()`/[prepare_profile_data()].
#' @param gene_name Character vector of gene names or IDs to profile.
#' @param group Name of the `col_data` column defining the groups
#'   (x axis).
#' @param levels Levels of `group` to show, in this order. Samples of other
#'   groups are ignored, as if the object had been subset to these samples.
#'   `NULL` takes all the levels of `group` (factor levels, else order of
#'   appearance). For a profile table, the groups of the table to draw.
#' @param assay Name of the transcript assay to summarise (default
#'   `"abundance"`).
#' @param log_scale Logical, show `log2(mean + 1)` on the y axis.
#' @param label Logical, label features at the last group level.
#' @param ncol Number of facet columns, `NULL` for automatic.
#' @param color_palette Named colour vector per transcript type; defaults
#'   to the object's palette or [tx_type_palette()].
#' @param ... Passed on to methods.
#'
#' @returns `profile_data()` returns a tibble with columns `feature_id`,
#'   `feature_name`, `feature_type`, `gene_id`, `gene_name`, `group`,
#'   `mean_expr`, `sd_expr`, `n_samples` and `is_de`; the grouping variable
#'   name is stored in the `"group_var"` attribute. `plot_tx_expr()` returns
#'   a `ggplot` object.
#'
#' @examples
#' iso <- example_isoformic()
#' head(profile_data(iso, gene_name = "RBPJ", group = "condition"))
#' plot_tx_expr(iso, gene_name = c("RBPJ", "EGFR"), group = "condition")
#'
#' @name profiles
#' @rdname profiles
NULL

#' @rdname profiles
#' @export
profile_data <- S7::new_generic(
  "profile_data",
  "self",
  function(self, gene_name, group, levels = NULL, assay = "abundance", ...) {
    S7::S7_dispatch()
  }
)

S7::method(profile_data, IsoformicExperiment) <- function(
  self,
  gene_name,
  group,
  levels = NULL,
  assay = "abundance",
  ...
) {
  rlang::check_required(gene_name)
  rlang::check_required(group)
  assay <- resolve_assay_name(self, assay)
  annot <- row_data(self)
  genes <- self@annot_data_genes
  gene_ids <- genes$gene_id[
    genes$gene_name %in% gene_name | genes$gene_id %in% gene_name
  ]
  if (!length(gene_ids)) {
    cli::cli_abort(
      message = c(x = "None of {.val {gene_name}} found in the annotation."),
      class = "isoformic_gene_not_found"
    )
  }
  tx_keep <- which(annot$gene_id %in% gene_ids)
  tx_mat <- assay(self, assay)
  grp <- group_samples(self, group, colnames(tx_mat), levels = levels)
  groups <- grp$groups
  tx_mat <- tx_mat[tx_keep, grp$samples, drop = FALSE]
  gene_mat <- gene_assay(self, assay)
  gene_mat <- gene_mat[
    intersect(gene_ids, rownames(gene_mat)),
    grp$samples,
    drop = FALSE
  ]

  gene_names <- genes$gene_name[match(
    c(rownames(gene_mat), annot$gene_id[tx_keep]),
    genes$gene_id
  )]

  features <- tibble::tibble(
    feature_id = c(rownames(gene_mat), rownames(tx_mat)),
    feature_name = c(
      gene_names[seq_len(nrow(gene_mat))],
      annot$transcript_name[tx_keep]
    ),
    feature_type = c(
      rep("gene", nrow(gene_mat)),
      annot$transcript_type[tx_keep]
    ),
    gene_id = c(rownames(gene_mat), annot$gene_id[tx_keep]),
    gene_name = gene_names
  )
  mat <- rbind(gene_mat, tx_mat)
  stats <- group_stats(mat, groups)
  res <- features[rep(seq_len(nrow(features)), times = nlevels(groups)), ]
  res$group <- factor(
    rep(levels(groups), each = nrow(features)),
    levels = levels(groups)
  )
  res$mean_expr <- as.vector(stats$mean)
  res$sd_expr <- as.vector(stats$sd)
  res$n_samples <- rep(
    as.vector(table(groups)[levels(groups)]),
    each = nrow(features)
  )

  res$is_de <- "no"
  dea <- self@dea
  if (!is.null(dea[["det"]])) {
    det <- de_tx(self)
    sig_tx <- det$transcript_id[det$det_sig == "yes"]
    res$is_de[res$feature_type != "gene" & res$feature_id %in% sig_tx] <- "yes"
  }
  if (!is.null(dea[["deg"]])) {
    deg <- de_gene(self)
    sig_genes <- deg$gene_id[deg$deg_sig == "yes"]
    res$is_de[
      res$feature_type == "gene" & res$feature_id %in% sig_genes
    ] <- "yes"
  }
  res <- dplyr::arrange(
    res,
    .data$gene_name,
    .data$feature_type != "gene",
    .data$feature_name,
    .data$group
  )
  attr(res, "group_var") <- group
  res
}

#' Per-group row means and standard deviations of a matrix
#' @keywords internal
#' @noRd
group_stats <- function(mat, groups) {
  lv <- levels(groups)
  means <- matrix(
    NA_real_,
    nrow(mat),
    length(lv),
    dimnames = list(rownames(mat), lv)
  )
  sds <- means
  for (l in lv) {
    cols <- which(groups == l)
    sub <- mat[, cols, drop = FALSE]
    mu <- rowMeans(sub)
    means[, l] <- mu
    if (length(cols) > 1L) {
      sds[, l] <- sqrt(rowSums((sub - mu)^2) / (length(cols) - 1L))
    }
  }
  list(mean = means, sd = sds)
}

#' @rdname profiles
#' @export
plot_tx_expr <- S7::new_generic(
  "plot_tx_expr",
  "self",
  function(
    self,
    gene_name = NULL,
    group = NULL,
    levels = NULL,
    assay = "abundance",
    log_scale = TRUE,
    label = TRUE,
    ncol = NULL,
    color_palette = NULL,
    ...
  ) {
    S7::S7_dispatch()
  }
)

S7::method(plot_tx_expr, IsoformicExperiment) <- function(
  self,
  gene_name = NULL,
  group = NULL,
  levels = NULL,
  assay = "abundance",
  log_scale = TRUE,
  label = TRUE,
  ncol = NULL,
  color_palette = NULL,
  ...
) {
  rlang::check_required(gene_name)
  rlang::check_required(group)
  data <- profile_data(
    self,
    gene_name = gene_name,
    group = group,
    levels = levels,
    assay = assay
  )
  plot_profile_table(
    data,
    log_scale = log_scale,
    label = label,
    ncol = ncol,
    color_palette = color_palette %||% self@tx_type_palette
  )
}

# nolint start: object_name_linter.
S7::method(plot_tx_expr, S7::class_data.frame) <- function(
  self,
  gene_name = NULL,
  group = NULL,
  levels = NULL,
  assay = "abundance",
  log_scale = TRUE,
  label = TRUE,
  ncol = NULL,
  color_palette = NULL,
  ...
) {
  data <- self
  group_var <- attr(data, "group_var") %||% "group"
  if (!is.null(gene_name)) {
    data <- data[data$gene_name %in% gene_name | data$gene_id %in% gene_name, ]
  }
  if (!is.null(levels) && "group" %in% colnames(data)) {
    # The table has its groups already: `levels` chooses and orders them
    grp <- group_factor(data$group, group = group_var, levels = levels)
    data <- data[grp$keep, ]
    data$group <- grp$groups
  }
  attr(data, "group_var") <- group_var
  plot_profile_table(
    data,
    log_scale = log_scale,
    label = label,
    ncol = ncol,
    color_palette = color_palette %||% tx_type_palette()
  )
}
# nolint end

#' Draw the profile plot from a profile table
#' @keywords internal
#' @noRd
plot_profile_table <- function(data, log_scale, label, ncol, color_palette) {
  required <- c(
    "feature_name",
    "feature_type",
    "gene_name",
    "group",
    "mean_expr",
    "sd_expr",
    "is_de"
  )
  missing_cols <- setdiff(required, colnames(data))
  if (length(missing_cols)) {
    cli::cli_abort(
      message = c(
        x = "Profile data is missing column{?s}: {.field {missing_cols}}.",
        i = "Create it with {.fn profile_data} or {.fn prepare_profile_data}."
      ),
      class = "isoformic_invalid_profile_data"
    )
  }
  if (!nrow(data)) {
    cli::cli_abort(
      c(x = "No features to plot."),
      class = "isoformic_empty_plot_data"
    )
  }
  group_var <- attr(data, "group_var") %||% "group"
  if (!is.factor(data$group)) {
    data$group <- factor(data$group, levels = unique(data$group))
  }
  data$sd_expr[is.na(data$sd_expr)] <- 0
  data$lower <- pmax(data$mean_expr - data$sd_expr, 0)
  data$upper <- data$mean_expr + data$sd_expr
  if (isTRUE(log_scale)) {
    data$y <- log2(data$mean_expr + 1)
    data$ymin <- log2(data$lower + 1)
    data$ymax <- log2(data$upper + 1)
    y_lab <- expression(log[2](mean ~ expression + 1))
  } else {
    data$y <- data$mean_expr
    data$ymin <- data$lower
    data$ymax <- data$upper
    y_lab <- "Mean expression"
  }
  data$is_de <- factor(data$is_de, levels = c("yes", "no"))
  data$feature_type <- unknown_type(data$feature_type)
  data$feature_type <- factor(
    data$feature_type,
    levels = c("gene", setdiff(unique(data$feature_type), "gene"))
  )
  palette <- complete_palette(color_palette, levels(data$feature_type))
  dodge <- ggplot2::position_dodge(width = 0.15)

  p <- ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = .data$group,
      y = .data$y,
      group = .data$feature_name,
      color = .data$feature_type,
      alpha = .data$is_de
    )
  ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = .data$ymin, ymax = .data$ymax),
      width = 0.1,
      linewidth = 0.4,
      position = dodge,
      show.legend = FALSE
    ) +
    ggplot2::geom_line(linewidth = 1, position = dodge) +
    ggplot2::geom_point(
      ggplot2::aes(fill = .data$feature_type),
      shape = 21,
      size = 2.6,
      color = "white",
      position = dodge
    ) +
    ggplot2::scale_color_manual(
      values = palette,
      aesthetics = c("color", "fill"),
      drop = TRUE
    ) +
    ggplot2::scale_alpha_manual(values = c(yes = 1, no = 0.3), drop = FALSE) +
    ggplot2::facet_wrap(
      ggplot2::vars(.data$gene_name),
      scales = "free_y",
      ncol = ncol
    ) +
    ggplot2::labs(
      x = group_var,
      y = y_lab,
      color = "Feature type",
      fill = "Feature type",
      alpha = "Differentially expressed"
    ) +
    isoformic_theme()
  if (isTRUE(label)) {
    last_level <- levels(data$group)[nlevels(data$group)]
    label_data <- data[data$group == last_level, ]
    p <- p +
      ggrepel::geom_text_repel(
        data = label_data,
        ggplot2::aes(label = .data$feature_name),
        size = 3,
        direction = "y",
        hjust = 0,
        nudge_x = 0.15,
        segment.color = "grey70",
        show.legend = FALSE,
        max.overlaps = Inf
      ) +
      ggplot2::scale_x_discrete(
        expand = ggplot2::expansion(mult = c(0.15, 0.5))
      )
  }
  p
}

#' Fill a palette with grey for feature types it does not cover
#' @keywords internal
#' @noRd
complete_palette <- function(palette, types) {
  palette <- palette %||% tx_type_palette()
  types <- types[!is.na(types)]
  missing <- setdiff(types, names(palette))
  if (length(missing)) {
    palette <- c(
      palette,
      stats::setNames(rep("#9e9e9e", length(missing)), missing)
    )
  }
  palette
}

#' Replace missing transcript types with the "unknown" label
#' @keywords internal
#' @noRd
unknown_type <- function(x) {
  x <- as.character(x)
  x[is.na(x) | !nzchar(x)] <- "unknown"
  x
}

#' Shared ggplot2 theme for isoformic plots
#' @keywords internal
#' @noRd
isoformic_theme <- function(base_size = 11) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"),
      legend.position = "right",
      plot.title = ggplot2::element_text(face = "bold")
    )
}

# Table-based profile data ----------------------------------------------------

#' Prepare Expression Profile Data From Loose Tables
#'
#' Table-based counterpart of [profile_data()]: combines a transcript
#' expression table, sample metadata, a transcript-to-gene table and
#' differential expression results into the long profile table consumed by
#' [plot_tx_expr()]. Gene expression is the sum of the isoform values unless
#' `txi_gene` is supplied.
#'
#' @param txi_transcript Transcript expression: a matrix or data frame with
#'   transcript IDs as row names or in a `transcript_id` column, or a
#'   `tximport` list (its `abundance` element is used).
#' @param sample_metadata Data frame whose first column (or `sample_id`
#'   column) holds the sample IDs and that contains the `var` column.
#' @param tx_to_gene Data frame with `transcript_id`, `gene_id`,
#'   `transcript_name`, `gene_name` and `transcript_type` columns (see
#'   [make_tx_to_gene()]).
#' @param de_result_gene Gene-level results with `gene_id` or `gene_name`,
#'   `log2FC` and `qvalue`/`pvalue` columns. May be `NULL`.
#' @param de_result_transcript Transcript-level results with
#'   `transcript_id` or `transcript_name`, `log2FC` and `qvalue`/`pvalue`
#'   columns. May be `NULL`.
#' @param var Name of the grouping column in `sample_metadata`.
#' @param var_levels Levels of `var` to show, in this order; samples of
#'   other groups are ignored. `NULL` uses the factor levels or the order of
#'   appearance.
#' @param txi_gene Optional gene expression table (same formats as
#'   `txi_transcript`, with gene IDs).
#' @param gene_name Optional genes to keep (names or IDs); `NULL` keeps all.
#' @param pvalue_cutoff,lfc_cutoff,use_fdr Significance cutoffs.
#'
#' @returns A tibble in the format described in [profile_data()].
#'
#' @examples
#' iso <- example_isoformic()
#' prof <- prepare_profile_data(
#'   txi_transcript = assay(iso, "abundance"),
#'   sample_metadata = col_data(iso),
#'   tx_to_gene = tx_annot(iso),
#'   de_result_gene = de_gene(iso),
#'   de_result_transcript = de_tx(iso),
#'   var = "condition",
#'   gene_name = "RBPJ"
#' )
#' head(prof)
#'
#' @export
prepare_profile_data <- function(
  txi_transcript,
  sample_metadata,
  tx_to_gene,
  de_result_gene = NULL,
  de_result_transcript = NULL,
  var,
  var_levels = NULL,
  txi_gene = NULL,
  gene_name = NULL,
  pvalue_cutoff = 0.05,
  lfc_cutoff = 1,
  use_fdr = TRUE
) {
  rlang::check_required(var)
  tx_mat <- expression_to_matrix(txi_transcript, "transcript_id")
  tx_to_gene <- tibble::as_tibble(tx_to_gene)
  required <- c(
    "transcript_id",
    "gene_id",
    "transcript_name",
    "gene_name",
    "transcript_type"
  )
  missing_cols <- setdiff(required, colnames(tx_to_gene))
  if (length(missing_cols)) {
    cli::cli_abort(
      message = c(
        x = "{.arg tx_to_gene} is missing column{?s}: {.field {missing_cols}}."
      ),
      class = "isoformic_invalid_tx_to_gene"
    )
  }
  tx_to_gene <- dplyr::distinct(
    tx_to_gene,
    .data$transcript_id,
    .keep_all = TRUE
  )
  annot <- tx_to_gene[match(rownames(tx_mat), tx_to_gene$transcript_id), ]
  keep <- !is.na(annot$gene_id)
  tx_mat <- tx_mat[keep, , drop = FALSE]
  annot <- annot[keep, ]
  if (!is.null(gene_name)) {
    keep <- annot$gene_name %in% gene_name | annot$gene_id %in% gene_name
    tx_mat <- tx_mat[keep, , drop = FALSE]
    annot <- annot[keep, ]
  }
  if (is.null(txi_gene)) {
    gene_mat <- rowsum(tx_mat, annot$gene_id, reorder = TRUE)
  } else {
    gene_mat <- expression_to_matrix(txi_gene, "gene_id")
    gene_mat <- gene_mat[
      intersect(unique(annot$gene_id), rownames(gene_mat)),
      colnames(tx_mat),
      drop = FALSE
    ]
  }
  col_data <- standardize_col_data(sample_metadata)
  if (!var %in% colnames(col_data)) {
    cli::cli_abort(
      message = c(
        x = "Column {.field {var}} not found in {.arg sample_metadata}."
      ),
      class = "isoformic_group_not_found"
    )
  }
  samples <- intersect(colnames(tx_mat), col_data$sample_id)
  values <- col_data[[var]][match(samples, col_data$sample_id)]
  grp <- group_factor(
    values,
    group = var,
    levels = var_levels,
    arg = "var_levels"
  )
  groups <- grp$groups
  samples <- samples[grp$keep]
  tx_mat <- tx_mat[, samples, drop = FALSE]
  gene_mat <- gene_mat[, samples, drop = FALSE]

  gene_lookup <- dplyr::distinct(tx_to_gene, .data$gene_id, .data$gene_name)
  features <- tibble::tibble(
    feature_id = c(rownames(gene_mat), rownames(tx_mat)),
    feature_name = c(
      gene_lookup$gene_name[match(rownames(gene_mat), gene_lookup$gene_id)],
      annot$transcript_name
    ),
    feature_type = c(rep("gene", nrow(gene_mat)), annot$transcript_type),
    gene_id = c(rownames(gene_mat), annot$gene_id),
    gene_name = c(
      gene_lookup$gene_name[match(rownames(gene_mat), gene_lookup$gene_id)],
      annot$gene_name
    )
  )
  stats <- group_stats(rbind(gene_mat, tx_mat), groups)
  res <- features[rep(seq_len(nrow(features)), times = nlevels(groups)), ]
  res$group <- factor(
    rep(levels(groups), each = nrow(features)),
    levels = levels(groups)
  )
  res$mean_expr <- as.vector(stats$mean)
  res$sd_expr <- as.vector(stats$sd)
  res$n_samples <- rep(
    as.vector(table(groups)[levels(groups)]),
    each = nrow(features)
  )
  res$is_de <- "no"
  sig_tx <- significant_ids(
    de_result_transcript,
    "transcript",
    tx_to_gene,
    pvalue_cutoff,
    lfc_cutoff,
    use_fdr
  )
  sig_genes <- significant_ids(
    de_result_gene,
    "gene",
    tx_to_gene,
    pvalue_cutoff,
    lfc_cutoff,
    use_fdr
  )
  res$is_de[res$feature_type != "gene" & res$feature_id %in% sig_tx] <- "yes"
  res$is_de[res$feature_type == "gene" & res$feature_id %in% sig_genes] <- "yes"
  res <- dplyr::arrange(
    res,
    .data$gene_name,
    .data$feature_type != "gene",
    .data$feature_name,
    .data$group
  )
  attr(res, "group_var") <- var
  res
}

#' Coerce expression inputs (matrix, data frame, tximport list) to a matrix
#' @keywords internal
#' @noRd
expression_to_matrix <- function(x, id_col) {
  if (is.list(x) && !is.data.frame(x) && !is.null(x[["abundance"]])) {
    x <- x[["abundance"]]
  }
  if (is.data.frame(x)) {
    cand <- intersect(c(id_col, "feature_id", "id", "...1", "X1"), colnames(x))
    if (length(cand)) {
      ids <- as.character(x[[cand[1L]]])
      x <- x[, setdiff(colnames(x), cand[1L]), drop = FALSE]
      x <- as.matrix(x)
      rownames(x) <- ids
    } else {
      x <- as.matrix(x)
    }
  }
  if (!is.matrix(x) || is.null(rownames(x))) {
    cli::cli_abort(
      message = c(
        x = paste0(
          "Expression data must be a matrix with feature IDs as row names, ",
          "or a data frame with a {.field {id_col}} column."
        )
      ),
      class = "isoformic_invalid_expression"
    )
  }
  storage.mode(x) <- "double"
  x
}

#' IDs of significant features in a loose DE table
#' @keywords internal
#' @noRd
significant_ids <- function(
  de,
  level,
  tx_to_gene,
  pvalue_cutoff,
  lfc_cutoff,
  use_fdr
) {
  if (is.null(de)) {
    return(character())
  }
  de <- tibble::as_tibble(de)
  id_col <- paste0(level, "_id")
  name_col <- paste0(level, "_name")
  if (!"log2FC" %in% colnames(de) && "log2FoldChange" %in% colnames(de)) {
    de$log2FC <- de$log2FoldChange
  }
  p_col <- if (isTRUE(use_fdr) && "qvalue" %in% colnames(de)) {
    "qvalue"
  } else if ("padj" %in% colnames(de) && isTRUE(use_fdr)) {
    "padj"
  } else {
    "pvalue"
  }
  sig <- !is.na(de$log2FC) &
    abs(de$log2FC) >= lfc_cutoff &
    !is.na(de[[p_col]]) &
    de[[p_col]] <= pvalue_cutoff
  if (id_col %in% colnames(de)) {
    return(sub("\\|.*$", "", as.character(de[[id_col]][sig])))
  }
  if (name_col %in% colnames(de)) {
    lookup <- dplyr::distinct(tx_to_gene[, c(id_col, name_col)])
    return(lookup[[id_col]][lookup[[name_col]] %in% de[[name_col]][sig]])
  }
  cli::cli_abort(
    message = c(
      x = "Differential expression results need a {.field {id_col}} or {.field {name_col}} column."
    ),
    class = "isoformic_invalid_dea"
  )
}
