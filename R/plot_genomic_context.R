#' Plot the Genomic Context of a Gene's Isoforms
#'
#' @description
#' Draws the exon-intron structure of every annotated transcript of one or
#' more genes along the genome, one row per transcript, coloured by
#' transcript type, with introns as lines carrying strand arrows. This is a
#' pure `ggplot2` implementation reading the exon table of the annotation
#' directory, so it needs no genome packages. The [IsoformicExperiment]
#' method can annotate each transcript with its mean expression and
#' highlight differentially expressed isoforms.
#'
#' The data frame method takes an exon table as returned by [exon_table()]
#' or [prepare_exon_annotation()] (columns `transcript_id`, `start_pos`,
#' `end_pos`, and optionally `transcript_name`, `transcript_type`,
#' `gene_name`, `seqid`, `strand`). `plot_tx_context()` is kept as an alias
#' of that method accepting the legacy column names `tx_id`, `exon_left`
#' and `exon_right`.
#'
#' @param self An [IsoformicExperiment] or an exon table.
#' @param gene_name Character vector of gene names or IDs.
#' @param assay Name of the assay used to label transcripts with their mean
#'   expression; `NULL` disables the labels.
#' @param highlight_de Logical, outline isoforms flagged as significant in
#'   the stored transcript-level results.
#' @param highlight_poison_exons Logical, outline exons found only in
#'   unproductive isoforms of the gene (see [poison_exons()]).
#' @param color_palette Named colour vector per transcript type; defaults
#'   to the object's palette or [tx_type_palette()].
#' @param upstream_offset,downstream_offset Extra genomic span (bp) shown
#'   before and after the gene.
#' @param exon_height Height of exon boxes (0-1).
#' @param highlight Optional character vector of transcript IDs or names to
#'   outline (data frame method).
#' @param expression Optional named numeric vector (by transcript ID) used
#'   to label transcripts (data frame method).
#' @param ncol Number of facet columns when several genes are plotted.
#' @param exon_table Legacy argument of `plot_tx_context()`.
#' @param ... Passed on to methods.
#'
#' @returns A `ggplot` object.
#'
#' @examples
#' iso <- example_isoformic()
#' plot_genomic_context(iso, gene_name = "RBPJ")
#' plot_genomic_context(iso, gene_name = "EGFR", assay = NULL, highlight_de = FALSE)
#'
#' # Table-based
#' plot_genomic_context(exon_table(iso, "PNCK"))
#'
#' @name plot_genomic_context
#' @rdname plot_genomic_context
#' @export
plot_genomic_context <- S7::new_generic(
  "plot_genomic_context",
  "self",
  function(
    self,
    gene_name = NULL,
    assay = "abundance",
    highlight_de = TRUE,
    color_palette = NULL,
    expression = NULL,
    highlight = NULL,
    highlight_poison_exons = TRUE,
    upstream_offset = 0,
    downstream_offset = 0,
    exon_height = 0.6,
    ncol = 1L,
    ...
  ) {
    S7::S7_dispatch()
  }
)

S7::method(plot_genomic_context, IsoformicExperiment) <- function(
  self,
  gene_name = NULL,
  assay = "abundance",
  highlight_de = TRUE,
  color_palette = NULL,
  expression = NULL,
  highlight = NULL,
  highlight_poison_exons = TRUE,
  upstream_offset = 0,
  downstream_offset = 0,
  exon_height = 0.6,
  ncol = 1L,
  ...
) {
  if (is.null(gene_name)) {
    cli::cli_abort(
      message = c(
        x = "{.arg gene_name} is required for an {.cls IsoformicExperiment}."
      ),
      class = "isoformic_gene_not_found"
    )
  }
  exons <- exon_table(self, gene_name)
  if ((is.null(expression)) && (!is.null(assay) && length(assay_names(self)))) {
    assay <- resolve_assay_name(self, assay)
    mat <- assay(self, assay)
    ids <- intersect(unique(exons$transcript_id), rownames(mat))
    if (length(ids)) {
      expression <- rowMeans(mat[ids, , drop = FALSE])
    }
  }
  if (
    is.null(highlight) && isTRUE(highlight_de) && !is.null(self@dea[["det"]])
  ) {
    det <- de_tx(self)
    highlight <- det$transcript_id[det$det_sig == "yes"]
  }
  poison <- NULL
  if (isTRUE(highlight_poison_exons)) {
    poison <- poison_exons(self, gene_name = gene_name)
  }
  plot_exon_table(
    exons,
    color_palette = color_palette %||% self@tx_type_palette,
    expression = expression,
    highlight = highlight,
    poison = poison,
    upstream_offset = upstream_offset,
    downstream_offset = downstream_offset,
    exon_height = exon_height,
    ncol = ncol
  )
}

# nolint start: object_name_linter.
S7::method(plot_genomic_context, S7::class_data.frame) <- function(
  self,
  gene_name = NULL,
  assay = "abundance",
  highlight_de = TRUE,
  color_palette = NULL,
  expression = NULL,
  highlight = NULL,
  highlight_poison_exons = TRUE,
  upstream_offset = 0,
  downstream_offset = 0,
  exon_height = 0.6,
  ncol = 1L,
  ...
) {
  exons <- standardize_exon_table(self)
  if (!is.null(gene_name) && "gene_name" %in% colnames(exons)) {
    exons <- exons[
      exons$gene_name %in% gene_name | exons$gene_id %in% gene_name,
    ]
  }
  plot_exon_table(
    exons,
    color_palette = color_palette %||% tx_type_palette(),
    expression = expression,
    highlight = highlight,
    upstream_offset = upstream_offset,
    downstream_offset = downstream_offset,
    exon_height = exon_height,
    ncol = ncol
  )
}
# nolint end

#' @rdname plot_genomic_context
#' @export
plot_tx_context <- function(exon_table, ...) {
  plot_genomic_context(standardize_exon_table(exon_table), ...)
}

#' Accept legacy exon table column names and fill optional columns
#' @keywords internal
#' @noRd
standardize_exon_table <- function(exons) {
  exons <- tibble::as_tibble(exons)
  renames <- c(
    transcript_id = "tx_id",
    start_pos = "exon_left",
    end_pos = "exon_right",
    seqid = "chr"
  )
  for (new in names(renames)) {
    old <- renames[[new]]
    if (!new %in% colnames(exons) && old %in% colnames(exons)) {
      names(exons)[names(exons) == old] <- new
    }
  }
  required <- c("transcript_id", "start_pos", "end_pos")
  missing_cols <- setdiff(required, colnames(exons))
  if (length(missing_cols)) {
    cli::cli_abort(
      message = c(
        x = "Exon table is missing column{?s}: {.field {missing_cols}}."
      ),
      class = "isoformic_invalid_exon_table"
    )
  }
  exons$start_pos <- as.integer(exons$start_pos)
  exons$end_pos <- as.integer(exons$end_pos)
  if (!"transcript_name" %in% colnames(exons)) {
    exons$transcript_name <- exons$transcript_id
  }
  if (!"transcript_type" %in% colnames(exons)) {
    exons$transcript_type <- "unknown"
  }
  exons$transcript_type <- unknown_type(exons$transcript_type)
  if ("haplotype" %in% colnames(exons) && !all(is.na(exons$haplotype))) {
    exons$gene_name <- ifelse(
      is.na(exons$haplotype),
      exons$gene_name,
      paste0(exons$gene_name, " (", exons$haplotype, ")")
    )
  }
  if (!"gene_name" %in% colnames(exons)) {
    exons$gene_name <- if ("gene_id" %in% colnames(exons)) {
      exons$gene_id
    } else {
      "gene"
    }
  }
  if (!"gene_id" %in% colnames(exons)) {
    exons$gene_id <- exons$gene_name
  }
  if (!"seqid" %in% colnames(exons)) {
    exons$seqid <- NA_character_
  }
  if (!"strand" %in% colnames(exons)) {
    exons$strand <- "+"
  }
  exons
}

#' Draw transcript structures from a standardised exon table
#' @keywords internal
#' @noRd
plot_exon_table <- function(
  exons,
  color_palette,
  expression = NULL,
  highlight = NULL,
  poison = NULL,
  upstream_offset = 0,
  downstream_offset = 0,
  exon_height = 0.6,
  ncol = 1L
) {
  if (!nrow(exons)) {
    cli::cli_abort(
      c(x = "No exons to plot."),
      class = "isoformic_empty_plot_data"
    )
  }
  exons <- dplyr::arrange(
    exons,
    .data$gene_name,
    .data$transcript_type,
    .data$start_pos
  )
  tx <- exons |>
    dplyr::group_by(.data$gene_name, .data$transcript_id) |>
    dplyr::summarise(
      transcript_name = dplyr::first(.data$transcript_name),
      transcript_type = dplyr::first(.data$transcript_type),
      seqid = dplyr::first(.data$seqid),
      strand = dplyr::first(.data$strand),
      tx_start = min(.data$start_pos),
      tx_end = max(.data$end_pos),
      .groups = "drop"
    ) |>
    dplyr::arrange(.data$gene_name, .data$transcript_type, .data$tx_start) |>
    dplyr::mutate(y = rev(seq_len(dplyr::n())), .by = "gene_name")
  # Each gene gets its own band of rows, so that with `scales = "free"` every
  # facet shows only the y-axis labels of its own transcripts.
  gene_offset <- cumsum(c(0L, utils::head(table(tx$gene_name), -1L) + 2L))
  names(gene_offset) <- sort(unique(tx$gene_name))
  tx$y <- tx$y + unname(gene_offset[tx$gene_name])
  tx$label <- tx$transcript_name
  if (!is.null(expression)) {
    expr <- expression[tx$transcript_id]
    has <- !is.na(expr)
    tx$label[has] <- sprintf(
      "%s  (%s)",
      tx$transcript_name[has],
      format(signif(expr[has], 3), trim = TRUE)
    )
  }
  tx$highlight <- tx$transcript_id %in%
    highlight |
    tx$transcript_name %in% highlight
  exons <- dplyr::left_join(
    exons,
    tx[, c("transcript_id", "y", "highlight")],
    by = "transcript_id"
  )
  introns <- exons |>
    dplyr::arrange(.data$transcript_id, .data$start_pos) |>
    dplyr::mutate(
      next_start = dplyr::lead(.data$start_pos),
      .by = "transcript_id"
    ) |>
    dplyr::filter(
      !is.na(.data$next_start) & .data$next_start > .data$end_pos
    ) |>
    dplyr::mutate(
      x = dplyr::if_else(.data$strand == "-", .data$next_start, .data$end_pos),
      xend = dplyr::if_else(
        .data$strand == "-",
        .data$end_pos,
        .data$next_start
      )
    )
  half <- exon_height / 2
  palette <- complete_palette(color_palette, unique(exons$transcript_type))
  xlim <- tx |>
    dplyr::group_by(.data$gene_name) |>
    dplyr::summarise(
      xmin = min(.data$tx_start) - upstream_offset,
      xmax = max(.data$tx_end) + downstream_offset,
      .groups = "drop"
    )
  pad <- (xlim$xmax - xlim$xmin) * 0.02
  xlim$xmin <- xlim$xmin - pad
  xlim$xmax <- xlim$xmax + pad
  chrom <- unique(stats::na.omit(tx$seqid))
  x_label <- if (length(chrom) == 1L) {
    paste0(chrom, " position (Mb)")
  } else {
    "Genomic position (Mb)"
  }

  p <- ggplot2::ggplot() +
    ggplot2::geom_blank(data = xlim, ggplot2::aes(x = .data$xmin)) +
    ggplot2::geom_blank(data = xlim, ggplot2::aes(x = .data$xmax)) +
    ggplot2::geom_segment(
      data = introns,
      ggplot2::aes(x = .data$x, xend = .data$xend, y = .data$y, yend = .data$y),
      color = "grey45",
      linewidth = 0.4,
      arrow = grid::arrow(length = grid::unit(1.6, "mm"), type = "open")
    ) +
    # A same-colour outline keeps short exons visible when a gene spans
    # hundreds of kilobases and an exon is narrower than one pixel.
    ggplot2::geom_rect(
      data = exons,
      ggplot2::aes(
        xmin = .data$start_pos,
        xmax = .data$end_pos,
        ymin = .data$y - half,
        ymax = .data$y + half,
        fill = .data$transcript_type,
        color = .data$transcript_type
      ),
      linewidth = 0.5
    )
  if (!is.null(poison) && nrow(poison)) {
    poison_key <- paste(poison$seqid, poison$start_pos, poison$end_pos)
    exon_key <- paste(exons$seqid, exons$start_pos, exons$end_pos)
    poison_rows <- exons[exon_key %in% poison_key, ]
    if (nrow(poison_rows)) {
      p <- p +
        ggplot2::geom_rect(
          data = poison_rows,
          ggplot2::aes(
            xmin = .data$start_pos,
            xmax = .data$end_pos,
            ymin = .data$y - half - 0.05,
            ymax = .data$y + half + 0.05
          ),
          fill = NA,
          color = "black",
          linewidth = 0.8
        )
    }
  }
  if (any(tx$highlight)) {
    tx_hl <- tx[tx$highlight, ]
    p <- p +
      ggplot2::geom_rect(
        data = tx_hl,
        ggplot2::aes(
          xmin = .data$tx_start,
          xmax = .data$tx_end,
          ymin = .data$y - half - 0.08,
          ymax = .data$y + half + 0.08
        ),
        fill = NA,
        color = "black",
        linewidth = 0.4,
        linetype = "dashed"
      )
  }
  # The transcript names are the y-axis labels: they line up with their
  # track and never overlap a transcript that starts at the left edge.
  p <- p +
    ggplot2::scale_fill_manual(values = palette, name = "Transcript type") +
    ggplot2::scale_color_manual(values = palette, guide = "none") +
    ggplot2::scale_x_continuous(
      labels = scales::label_number(scale = 1e-6, accuracy = 0.001),
      expand = ggplot2::expansion(mult = c(0.01, 0.02))
    ) +
    ggplot2::scale_y_continuous(
      breaks = tx$y,
      labels = tx$label,
      minor_breaks = NULL,
      expand = ggplot2::expansion(add = 0.6)
    ) +
    ggplot2::labs(x = x_label, y = NULL) +
    isoformic_theme() +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(size = 8, color = "grey20"),
      legend.position = "bottom"
    )
  if (length(unique(tx$gene_name)) > 1L) {
    p <- p +
      ggplot2::facet_wrap(
        ggplot2::vars(.data$gene_name),
        scales = "free",
        ncol = ncol
      )
  } else {
    p <- p + ggplot2::labs(title = unique(tx$gene_name))
  }
  p
}
