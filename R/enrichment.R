# Enrichment ------------------------------------------------------------------

#' Transcript-Type Aware Gene Set Enrichment
#'
#' @description
#' Runs gene set enrichment analysis (`fgsea::fgseaMultilevel()`) separately
#' for transcripts of each type, ranking transcripts by log2 fold-change.
#' Gene sets defined on gene symbols are expanded to all transcripts of
#' their genes, so each transcript category (protein coding, retained
#' intron, nonsense-mediated decay, ...) is tested on its own and in the
#' combined `"unproductive"` category.
#'
#' The [IsoformicExperiment] method uses the stored transcript-level
#' results and the annotation; the data frame method takes a transcript
#' table with `transcript_name`, `transcript_type` and `log2FC` columns plus
#' a `tx_to_gene` table.
#'
#' `read_gmt()` reads gene sets in GMT format into a named list, and
#' `plot_enrichment()` draws a dot plot of enrichment results.
#'
#' @param self An [IsoformicExperiment] or a transcript-level DE table.
#' @param genesets_list Named list of character vectors (gene symbols), as
#'   returned by [read_gmt()].
#' @param tx_to_gene Data frame with `transcript_name` and `gene_name`
#'   columns; required by the data frame method, ignored otherwise.
#' @param pval_cutoff Keep pathways with `pval` at or below this value;
#'   `NULL` keeps all.
#' @param tx_types Named list defining the transcript categories to test.
#'   Defaults to [default_enrichment_types()], which picks
#'   `processed_transcript` or `protein_coding_CDS_not_defined` depending
#'   on the annotation version.
#' @param min_size,max_size Minimum and maximum gene set size passed to
#'   `fgsea`.
#' @param ... Passed on to `fgsea::fgseaMultilevel()`.
#' @param path Path to a GMT file.
#' @param enrichment_df Result of `run_enrichment()`.
#' @param n Number of top pathways (by adjusted p-value) to show.
#' @param experiments Categories to show; `NULL` shows all.
#' @param nes_cutoff Minimum absolute NES to show.
#' @param strip_prefix Regular expression removed from pathway names before
#'   plotting (default drops MSigDB collection prefixes such as `REACTOME_`).
#' @param label_width Pathway labels are wrapped to this many characters.
#' @param annotation_types Transcript types present in the annotation,
#'   used to choose between `processed_transcript` and
#'   `protein_coding_CDS_not_defined`.
#'
#' @returns `run_enrichment()` returns a tibble with the `fgsea` columns
#'   plus `experiment` (the transcript category). `read_gmt()` returns a
#'   named list. `plot_enrichment()` returns a `ggplot` object.
#'
#' @examples
#' iso <- example_isoformic()
#' gmt <- system.file("extdata", "reactome_subset.gmt", package = "isoformic")
#' genesets <- read_gmt(gmt)
#' length(genesets)
#'
#' if (requireNamespace("fgsea", quietly = TRUE)) {
#'   res <- run_enrichment(iso, genesets, pval_cutoff = NULL, min_size = 2)
#'   head(res)
#'   plot_enrichment(res, n = 10)
#' }
#'
#' @name enrichment
#' @rdname enrichment
NULL

#' @rdname enrichment
#' @export
run_enrichment <- S7::new_generic(
  "run_enrichment",
  "self",
  function(
    self,
    genesets_list,
    tx_to_gene = NULL,
    pval_cutoff = 0.05,
    tx_types = NULL,
    min_size = 15,
    max_size = 500,
    ...
  ) {
    S7::S7_dispatch()
  }
)

S7::method(run_enrichment, IsoformicExperiment) <- function(
  self,
  genesets_list,
  tx_to_gene = NULL,
  pval_cutoff = 0.05,
  tx_types = NULL,
  min_size = 15,
  max_size = 500,
  ...
) {
  rlang::check_required(genesets_list)
  det <- de_tx(self)
  annot <- tx_annot(self)
  det <- dplyr::left_join(
    det,
    annot[, c(
      "transcript_id",
      "transcript_name",
      "transcript_type",
      "gene_name"
    )],
    by = "transcript_id"
  )
  tx_types <- tx_types %||%
    default_enrichment_types(unique(annot$transcript_type))
  run_fgsea_by_type(
    det = det,
    rank_id = "transcript_id",
    genesets_list = genesets_list,
    tx_to_gene = annot[, c("transcript_id", "gene_name")],
    map_id = "transcript_id",
    pval_cutoff = pval_cutoff,
    tx_types = tx_types,
    min_size = min_size,
    max_size = max_size,
    ...
  )
}

# nolint start: object_name_linter.
S7::method(run_enrichment, S7::class_data.frame) <- function(
  self,
  genesets_list,
  tx_to_gene = NULL,
  pval_cutoff = 0.05,
  tx_types = NULL,
  min_size = 15,
  max_size = 500,
  ...
) {
  rlang::check_required(genesets_list)
  if (is.null(tx_to_gene)) {
    cli::cli_abort(
      message = c(
        x = "{.arg tx_to_gene} is required when {.arg self} is a data frame."
      ),
      class = "isoformic_invalid_tx_to_gene"
    )
  }
  det <- tibble::as_tibble(self)
  required <- c("transcript_name", "transcript_type", "log2FC")
  missing_cols <- setdiff(required, colnames(det))
  if (length(missing_cols)) {
    cli::cli_abort(
      message = c(
        x = "{.arg self} is missing column{?s}: {.field {missing_cols}}."
      ),
      class = "isoformic_invalid_de_table"
    )
  }
  types_in_annot <- if ("transcript_type" %in% colnames(tx_to_gene)) {
    unique(tx_to_gene$transcript_type)
  } else {
    unique(det$transcript_type)
  }
  tx_types <- tx_types %||% default_enrichment_types(types_in_annot)
  run_fgsea_by_type(
    det = det,
    rank_id = "transcript_name",
    genesets_list = genesets_list,
    tx_to_gene = tibble::as_tibble(tx_to_gene)[, c(
      "transcript_name",
      "gene_name"
    )],
    map_id = "transcript_name",
    pval_cutoff = pval_cutoff,
    tx_types = tx_types,
    min_size = min_size,
    max_size = max_size,
    ...
  )
}
# nolint end

#' @rdname enrichment
#' @export
default_enrichment_types <- function(annotation_types = NULL) {
  processed_or_cds <- "protein_coding_CDS_not_defined"
  if (
    !is.null(annotation_types) &&
      "processed_transcript" %in% annotation_types &&
      !"protein_coding_CDS_not_defined" %in% annotation_types
  ) {
    processed_or_cds <- "processed_transcript"
  }
  list(
    protein_coding = "protein_coding",
    unproductive = c(
      "retained_intron",
      processed_or_cds,
      "nonsense_mediated_decay"
    ),
    retained_intron = "retained_intron",
    processed_or_cds = processed_or_cds,
    nonsense_mediated_decay = "nonsense_mediated_decay"
  ) |>
    stats::setNames(c(
      "protein_coding",
      "unproductive",
      "retained_intron",
      processed_or_cds,
      "nonsense_mediated_decay"
    ))
}

#' Run fgsea once per transcript category
#' @keywords internal
#' @noRd
run_fgsea_by_type <- function(
  det,
  rank_id,
  genesets_list,
  tx_to_gene,
  map_id,
  pval_cutoff,
  tx_types,
  min_size,
  max_size,
  ...
) {
  rlang::check_installed("fgsea")
  if (!is.list(genesets_list) || is.null(names(genesets_list))) {
    cli::cli_abort(
      message = c(
        x = "{.arg genesets_list} must be a named list of gene symbol vectors."
      ),
      class = "isoformic_invalid_genesets"
    )
  }
  # Expand gene-symbol sets to transcript IDs once (vectorised join)
  sets_long <- tibble::tibble(
    term = rep(names(genesets_list), lengths(genesets_list)),
    gene_name = unlist(genesets_list, use.names = FALSE)
  ) |>
    dplyr::inner_join(
      dplyr::distinct(tx_to_gene[, c(map_id, "gene_name")]),
      by = "gene_name",
      relationship = "many-to-many"
    )
  tx_sets <- split(sets_long[[map_id]], sets_long$term)

  det <- det[!is.na(det$log2FC) & !is.na(det[[rank_id]]), ]
  results <- lapply(names(tx_types), function(type_name) {
    sub <- det[det$transcript_type %in% tx_types[[type_name]], ]
    if (nrow(sub) < 2L) {
      return(NULL)
    }
    ranks <- stats::setNames(sub$log2FC, sub[[rank_id]])
    ranks <- ranks[!duplicated(names(ranks))]
    ranks <- sort(ranks)
    res <- fgsea::fgseaMultilevel(
      pathways = tx_sets,
      stats = ranks,
      minSize = min_size,
      maxSize = max_size,
      eps = 0,
      ...
    )
    res <- tibble::as_tibble(as.data.frame(res))
    if (!is.null(pval_cutoff)) {
      res <- res[!is.na(res$pval) & res$pval <= pval_cutoff, ]
    }
    res$experiment <- type_name
    res
  })
  dplyr::bind_rows(results)
}

#' @rdname enrichment
#' @export
read_gmt <- function(path) {
  assert_file_exists(path)
  lines <- readLines(path, warn = FALSE)
  lines <- lines[nzchar(trimws(lines))]
  parts <- strsplit(lines, "\t", fixed = TRUE)
  sets <- lapply(parts, function(p) unique(p[-(1:2)][nzchar(p[-(1:2)])]))
  names(sets) <- vapply(parts, function(p) p[[1L]], character(1L))
  sets
}

#' @rdname enrichment
#' @export
plot_enrichment <- function(
  enrichment_df,
  n = 30,
  experiments = NULL,
  nes_cutoff = 0,
  strip_prefix = "^[A-Z]+_",
  label_width = 50
) {
  data <- tibble::as_tibble(enrichment_df)
  required <- c("pathway", "NES", "padj", "size", "experiment")
  missing_cols <- setdiff(required, colnames(data))
  if (length(missing_cols)) {
    cli::cli_abort(
      message = c(
        x = "{.arg enrichment_df} is missing column{?s}: {.field {missing_cols}}."
      ),
      class = "isoformic_invalid_enrichment"
    )
  }
  if (!is.null(experiments)) {
    data <- data[data$experiment %in% experiments, ]
  }
  data <- data[!is.na(data$NES) & abs(data$NES) >= nes_cutoff, ]
  if (!nrow(data)) {
    cli::cli_abort(
      c(x = "No pathways to plot."),
      class = "isoformic_empty_plot_data"
    )
  }
  top <- data |>
    dplyr::group_by(.data$pathway) |>
    dplyr::summarise(
      min_padj = min(.data$padj, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::arrange(.data$min_padj) |>
    dplyr::slice_head(n = n)
  data <- data[data$pathway %in% top$pathway, ]
  order_df <- data |>
    dplyr::group_by(.data$pathway) |>
    dplyr::summarise(mean_nes = mean(.data$NES), .groups = "drop") |>
    dplyr::arrange(.data$mean_nes)
  data$pathway <- factor(data$pathway, levels = order_df$pathway)
  data$experiment <- factor(
    data$experiment,
    levels = unique(enrichment_df$experiment)
  )
  limit <- max(abs(data$NES))
  pretty_label <- function(x) {
    x <- as.character(x)
    if (!is.null(strip_prefix)) {
      x <- sub(strip_prefix, "", x)
    }
    x <- gsub("_", " ", x, fixed = TRUE)
    vapply(
      strwrap(x, width = label_width, simplify = FALSE),
      paste,
      character(1L),
      collapse = "\n"
    )
  }
  ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = .data$experiment,
      y = .data$pathway,
      color = .data$NES,
      size = .data$size
    )
  ) +
    ggplot2::geom_point() +
    ggplot2::scale_color_gradient2(
      low = "#2166ac",
      mid = "grey85",
      high = "#b2182b",
      midpoint = 0,
      limits = c(-limit, limit)
    ) +
    ggplot2::scale_size_area(max_size = 7) +
    ggplot2::scale_y_discrete(labels = pretty_label) +
    ggplot2::scale_x_discrete(labels = function(x) {
      gsub("_", " ", x, fixed = TRUE)
    }) +
    ggplot2::labs(x = NULL, y = NULL, color = "NES", size = "Set size") +
    isoformic_theme() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 30, hjust = 1),
      axis.text.y = ggplot2::element_text(size = 8)
    )
}
