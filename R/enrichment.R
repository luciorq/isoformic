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
#' @section Ranking:
#' `rank_by` chooses the statistic that ranks the transcripts:
#'
#' * `"log2FC"` (default): the log2 fold-change of the stored results, a
#'   point estimate that ignores how well it is determined.
#' * `"signed_p"`: `-log10(pvalue)` with the sign of the fold-change, so that
#'   the test of the differential expression weighs the ranking. A p-value
#'   of zero gets the largest finite value of the table.
#' * `"lfc_z"`: the median fold-change over the inferential replicates
#'   divided by its standard deviation (`lfc_median / lfc_sd` of
#'   [de_confidence()]), so that a fold-change the quantification cannot
#'   determine ranks low. It needs the confidence table that
#'   `add_uncertainty(x, group = )` stores, or these two columns in a data
#'   frame (the table of [de_confidence()] works). The ranking has the
#'   contrast of that table, the last level of `add_uncertainty()` against
#'   the first; a warning (`isoformic_contrast_mismatch`) says so when most
#'   of its fold-changes have the sign opposite to the stored results.
#'   Transcripts that were not assessed, that were collapsed into a group
#'   or, for an object, that are not in its differential expression results
#'   are not ranked. A fold-change that is the same in every replicate
#'   gets the largest finite value of the table.
#'
#' @section Collections of gene sets:
#' `genesets_list` is either one collection (a named list of gene symbol
#' vectors) or several (a named list of such lists, for example
#' `list(hallmark = read_gmt(h), reactome = read_gmt(r))`). Each collection is
#' tested on its own, so that the adjusted p-values are those of the
#' collection, and the result has a `collection` column, by which
#' `plot_enrichment()` splits its panels. Every gene set and every collection
#' needs a name.
#'
#' @param self An [IsoformicExperiment] or a transcript-level DE table.
#' @param genesets_list Named list of character vectors (gene symbols), as
#'   returned by [read_gmt()], or a named list of such lists (see
#'   *Collections of gene sets*).
#' @param rank_by Statistic that ranks the transcripts: `"log2FC"`
#'   (default), `"signed_p"` or `"lfc_z"` (see *Ranking*).
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
#' @param ... Passed on to `fgsea::fgseaMultilevel()`. fgsea runs in the
#'   session unless `nproc` or `BPPARAM` is given. On Unix-like systems these
#'   fork the session: the workers share its memory, but tools that add up
#'   the resident memory of the processes count it once per worker. Give the
#'   number of workers (`nproc = 4`); the default of `BiocParallel` is all
#'   the cores of the machine but two, whatever a job scheduler allocated.
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
#'   plus `experiment` (the transcript category), `rank_by` and, for several
#'   collections, `collection`. `read_gmt()` returns a
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
    rank_by = c("log2FC", "signed_p", "lfc_z"),
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
  rank_by = c("log2FC", "signed_p", "lfc_z"),
  ...
) {
  rlang::check_required(genesets_list)
  rank_by <- rlang::arg_match(rank_by)
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
  if (identical(rank_by, "lfc_z")) {
    conf <- self@dea[["confidence"]]
    if (is.null(conf)) {
      cli::cli_abort(
        message = c(
          x = "{.code rank_by = \"lfc_z\"} needs the confidence of the fold-changes.",
          i = "Store it with {.code add_uncertainty(x, group = )}."
        ),
        class = "isoformic_no_inferential_replicates"
      )
    }
    conf <- conf[conf$feature_type %in% "transcript", ]
    idx <- match(det$transcript_id, conf$feature_id)
    det$lfc_median <- conf$lfc_median[idx]
    det$lfc_sd <- conf$lfc_sd[idx]
    # The confidence has the contrast of `add_uncertainty(group = , levels = )`,
    # which may not be the one of the stored results
    same_sign <- sign(det$lfc_median) == sign(det$log2FC)
    same_sign <- same_sign[
      !is.na(same_sign) & det$lfc_median != 0 & det$log2FC != 0
    ]
    if (length(same_sign) >= 2L && mean(same_sign) < 0.5) {
      cli::cli_warn(
        message = c(
          "!" = paste0(
            "Most fold-changes of the confidence table have the sign ",
            "opposite to the stored differential expression results."
          ),
          i = paste0(
            "{.code rank_by = \"lfc_z\"} has the contrast of ",
            "{.code add_uncertainty(x, group = , levels = )}: check its levels."
          )
        ),
        class = "isoformic_contrast_mismatch"
      )
    }
  }
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
    rank_by = rank_by,
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
  rank_by = c("log2FC", "signed_p", "lfc_z"),
  ...
) {
  rlang::check_required(genesets_list)
  rank_by <- rlang::arg_match(rank_by)
  if (is.null(tx_to_gene)) {
    cli::cli_abort(
      message = c(
        x = "{.arg tx_to_gene} is required when {.arg self} is a data frame."
      ),
      class = "isoformic_invalid_tx_to_gene"
    )
  }
  det <- tibble::as_tibble(self)
  required <- c(
    "transcript_name",
    "transcript_type",
    if (!identical(rank_by, "lfc_z")) "log2FC"
  )
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
    rank_by = rank_by,
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

#' Abort unless `genesets_list` is a named list of gene symbol vectors
#'
#' Or a named list of such lists, see `geneset_collections()`.
#' @keywords internal
#' @noRd
check_genesets <- function(genesets_list) {
  geneset_collections(genesets_list)
  invisible(genesets_list)
}

#' The collections of `genesets_list`
#'
#' One collection is a named list of gene symbol vectors; several are a named
#' list of such lists. Returns a list of collections, without names for a
#' single collection given as such.
#' @keywords internal
#' @noRd
geneset_collections <- function(genesets_list) {
  named <- function(x) {
    !is.null(names(x)) && !anyNA(names(x)) && all(nzchar(names(x)))
  }
  # A gene set: a vector, an empty set, or a list of single values (as JSON
  # and YAML readers give)
  is_set <- function(set) {
    is.null(set) ||
      (is.atomic(set) && !is.list(set)) ||
      (is.list(set) &&
        is.null(names(set)) &&
        all(vapply(
          set,
          function(v) is.atomic(v) && length(v) == 1L,
          logical(1L)
        )))
  }
  is_collection <- function(x) {
    is.list(x) && named(x) && all(vapply(x, is_set, logical(1L)))
  }
  as_collection <- function(x) {
    lapply(x, function(set) as.character(unlist(set, use.names = FALSE)))
  }
  if (is_collection(genesets_list)) {
    return(list(as_collection(genesets_list)))
  }
  if (
    is.list(genesets_list) &&
      length(genesets_list) > 0L &&
      named(genesets_list) &&
      !anyDuplicated(names(genesets_list)) &&
      all(vapply(genesets_list, is_collection, logical(1L)))
  ) {
    return(lapply(genesets_list, as_collection))
  }
  cli::cli_abort(
    message = c(
      x = "{.arg genesets_list} must be a named list of gene symbol vectors.",
      i = "Every gene set needs a name. Several collections are a named list of such lists."
    ),
    class = "isoformic_invalid_genesets"
  )
}

#' Run `run(collection)` for every collection and stack the results
#'
#' The `collection` column is added only when several collections were given.
#' @keywords internal
#' @noRd
enrich_collections <- function(genesets_list, run) {
  collections <- geneset_collections(genesets_list)
  if (is.null(names(collections))) {
    return(run(collections[[1L]]))
  }
  res <- lapply(seq_along(collections), function(i) {
    out <- run(collections[[i]])
    out$collection <- rep(names(collections)[[i]], nrow(out))
    out
  })
  dplyr::bind_rows(res)
}

#' A `fgsea` result without rows
#' @keywords internal
#' @noRd
empty_fgsea <- function() {
  tibble::tibble(
    pathway = character(),
    pval = numeric(),
    padj = numeric(),
    log2err = numeric(),
    ES = numeric(),
    NES = numeric(),
    size = integer(),
    leadingEdge = list(),
    experiment = character()
  )
}

#' Estimate over its standard deviation, for a ranking
#'
#' An estimate that is the same in every replicate (a standard deviation of
#' zero) is the most certain one: it gets the largest finite value of the
#' vector. No change in any replicate is zero.
#' @keywords internal
#' @noRd
certainty_z <- function(estimate, sd) {
  z <- estimate / sd
  z[!is.finite(z)] <- NA_real_
  certain <- sd %in% 0 & !estimate %in% c(0, NA)
  finite_z <- abs(z[!is.na(z)])
  z[certain] <- sign(estimate[certain]) * max(finite_z, 1)
  z[sd %in% 0 & estimate %in% 0] <- 0
  z
}

#' The statistic that ranks the rows of a differential expression table
#' @keywords internal
#' @noRd
enrichment_rank <- function(det, rank_by) {
  needed <- switch(
    rank_by,
    log2FC = "log2FC",
    signed_p = c("log2FC", "pvalue"),
    lfc_z = c("lfc_median", "lfc_sd")
  )
  missing_cols <- setdiff(needed, colnames(det))
  if (length(missing_cols)) {
    cli::cli_abort(
      message = c(
        x = "{.code rank_by = \"{rank_by}\"} needs the column{?s} {.field {missing_cols}}.",
        i = if (identical(rank_by, "lfc_z")) {
          "They come from {.fn de_confidence}."
        }
      ),
      class = if (identical(rank_by, "lfc_z")) {
        "isoformic_no_inferential_replicates"
      } else {
        "isoformic_invalid_de_table"
      }
    )
  }
  if (identical(rank_by, "log2FC")) {
    return(as.numeric(det$log2FC))
  }
  if (identical(rank_by, "lfc_z")) {
    return(certainty_z(as.numeric(det$lfc_median), as.numeric(det$lfc_sd)))
  }
  score <- -log10(as.numeric(det$pvalue))
  # A p-value of zero is the smallest of the table, not an infinite score
  finite <- score[is.finite(score)]
  score[is.infinite(score)] <- if (length(finite)) max(finite) else 1
  sign(as.numeric(det$log2FC)) * score
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
  rank_by = "log2FC",
  ...
) {
  rlang::check_installed("fgsea")
  check_genesets(genesets_list)
  det$rank_stat <- enrichment_rank(det, rank_by)
  if (identical(rank_by, "lfc_z") && !any(is.finite(det$rank_stat))) {
    cli::cli_abort(
      message = c(
        x = "{.code rank_by = \"lfc_z\"} needs at least two inferential replicates.",
        i = "No transcript has a standard deviation of its fold-change."
      ),
      class = "isoformic_no_inferential_replicates"
    )
  }
  det <- det[!is.na(det$rank_stat) & !is.na(det[[rank_id]]), ]
  gene_map <- dplyr::distinct(tx_to_gene[, c(map_id, "gene_name")])
  res <- enrich_collections(genesets_list, function(sets) {
    # Expand gene-symbol sets to transcript IDs once (vectorised join)
    sets_long <- tibble::tibble(
      term = rep(names(sets), lengths(sets)),
      gene_name = as.character(unlist(sets, use.names = FALSE))
    ) |>
      dplyr::inner_join(
        gene_map,
        by = "gene_name",
        relationship = "many-to-many"
      )
    tx_sets <- split(sets_long[[map_id]], sets_long$term)
    results <- lapply(names(tx_types), function(type_name) {
      sub <- det[det$transcript_type %in% tx_types[[type_name]], ]
      if (nrow(sub) < 2L) {
        return(NULL)
      }
      ranks <- stats::setNames(sub$rank_stat, sub[[rank_id]])
      ranks <- ranks[!duplicated(names(ranks))]
      ranks <- sort(ranks)
      res <- fgsea_multilevel(
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
      res$experiment <- rep(type_name, nrow(res))
      res
    })
    results <- dplyr::bind_rows(results)
    if (!ncol(results)) empty_fgsea() else results
  })
  res$rank_by <- rep(rank_by, nrow(res))
  dplyr::relocate(res, dplyr::any_of("collection"), .after = "rank_by")
}

#' Gene set enrichment on a gene table ranked by one of its columns
#'
#' The largest absolute value is kept for a gene symbol that is repeated.
#' @keywords internal
#' @noRd
gene_rank_enrichment <- function(
  table,
  rank_col,
  genesets_list,
  experiment,
  rank_by,
  min_size,
  max_size,
  ...
) {
  table <- table[!is.na(table[[rank_col]]) & !is.na(table$gene_name), ]
  table <- table[order(-abs(table[[rank_col]])), ]
  table <- table[!duplicated(table$gene_name), ]
  ranks <- sort(stats::setNames(table[[rank_col]], table$gene_name))
  res <- enrich_collections(genesets_list, function(sets) {
    res <- fgsea_multilevel(
      pathways = sets,
      stats = ranks,
      minSize = min_size,
      maxSize = max_size,
      eps = 0,
      ...
    )
    tibble::as_tibble(as.data.frame(res))
  })
  res$experiment <- rep(experiment, nrow(res))
  res$rank_by <- rep(rank_by, nrow(res))
  dplyr::relocate(res, dplyr::any_of("collection"), .after = "rank_by")
}

#' Run `fgsea::fgseaMultilevel()` serially unless asked otherwise
#'
#' Without `BPPARAM` or `nproc`, fgsea uses `BiocParallel::bpparam()`, which
#' starts socket workers (`SnowParam`) on Windows; they can hang R CMD check.
#' @keywords internal
#' @noRd
fgsea_multilevel <- function(...) {
  args <- list(...)
  if (is.null(args[["BPPARAM"]]) && is.null(args[["nproc"]])) {
    args[["BPPARAM"]] <- BiocParallel::SerialParam()
  }
  do.call(fgsea::fgseaMultilevel, args)
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
    ) +
    # The same pathway name in two collections is two results
    if ("collection" %in% colnames(data)) {
      ggplot2::facet_wrap(ggplot2::vars(.data$collection), scales = "free_y")
    }
}
