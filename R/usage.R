#' Test differential transcript usage
#'
#' @description
#' `test_usage()` tests whether the share of a gene's expression carried by
#' each isoform differs between groups of samples, with one of three
#' established frameworks:
#'
#' * `"diffsplice"`: `edgeR` quasi-likelihood fits and `diffSplice`, the
#'   pipeline of the 2025 *Nucleic Acids Research* paper below. Fast, and the
#'   one designed for counts divided by the read-ambiguity overdispersion.
#' * `"satuRn"`: quasi-binomial models of transcript usage (Gilis et al.,
#'   2021), suited to many samples and to single-cell data.
#' * `"DEXSeq"`: negative binomial models of exon (here transcript) usage
#'   (Anders et al. 2012), the slowest and often the most conservative.
#'
#' Recent benchmarks of usage tests (2024 and 2025) rank these among the best;
#' none dominates, so the backend is an argument. With inferential
#' replicates in the object (see [add_uncertainty()]), the counts are divided
#' by the overdispersion that the read ambiguity adds (see
#' [inferential_overdispersion()]) before the test (`scale = TRUE`).
#'
#' Transcripts are kept when at least `min_count` reads and a share of
#' `min_prop` of the gene fall on them in as many samples as the smallest
#' group, and genes need two such transcripts.
#'
#' The gene-level p-value combines the transcript p-values of a gene with
#' Simes' method, and both levels are adjusted with the Benjamini-Hochberg
#' procedure within each contrast, so that the three backends report the same
#' columns.
#'
#' With more than two groups, every contrast of `contrasts` is tested and the
#' results are stacked, one block per contrast. `"diffsplice"` and `"satuRn"`
#' fit all the groups once and test each contrast on that fit, so the
#' dispersions use every sample; `"DEXSeq"` compares a model with and without
#' a group effect, so it is fitted on the samples of the two groups of each
#' contrast.
#'
#' `contrasts = "omnibus"` gives one p-value per transcript for a difference
#' among any of the groups, with `contrast = "omnibus"` and no fold-change
#' (`usage_log2FC` is `NA`; look at the contrasts for the direction).
#' `"DEXSeq"` compares its models with and without the group effect on all the
#' samples. `edgeR` and `satuRn` test one contrast at a time, so for them the
#' p-value is the Simes combination of the contrasts of every level against the
#' first: the groups can only differ if some level differs from the first, and
#' those contrasts share the reference, so they are positively correlated and
#' Simes' test keeps its error rate.
#'
#' @param self An [IsoformicExperiment].
#' @param group Column of [col_data()] with the groups to compare.
#' @param levels Levels of `group` to use, reference first. Samples of other
#'   groups are ignored. `NULL` takes all the levels of `group`.
#' @param contrasts The comparisons to test: `NULL` (default) compares every
#'   level with the first one; `"pairwise"` tests every pair of levels (the
#'   later level against the earlier one); `"omnibus"` tests whether usage
#'   differs among any of the groups (see Details); otherwise a list of
#'   `c(numerator, denominator)` pairs of levels, or one such pair. With two
#'   levels the default is the single contrast second against first.
#' @param method `"diffsplice"`, `"satuRn"` or `"DEXSeq"`. The backend must be
#'   installed.
#' @param assay Counts assay to test.
#' @param scale Logical, divide the counts by the overdispersion of the
#'   replicates. `NULL` (default) does it when the object has replicates. In
#'   a read-split experiment (technical replicates of the same libraries, see
#'   the uncertainty article) this removed most false positives of `"diffsplice"`
#'   and `"DEXSeq"` and made their calls reproducible, so keep it on.
#' @param min_count,min_prop Filters described above.
#' @param workers Number of workers for `satuRn` and `DEXSeq` (parallel on
#'   Unix-like systems).
#' @param gene_name Gene names or IDs to keep in the result. The filter is
#'   applied after the test.
#' @param ... Unused.
#'
#' @returns A tibble with one row per tested transcript and contrast:
#'   `gene_id`, `gene_name`, `transcript_id`, `transcript_name`,
#'   `transcript_type`, `contrast` (`"numerator vs denominator"`),
#'   `usage_log2FC` (the change in usage on the log2 scale: for `"diffsplice"`
#'   the transcript against the other transcripts of the gene, for `"satuRn"`
#'   the log2 odds ratio, for `"DEXSeq"` the usage fold-change), `pvalue`,
#'   `qvalue`, `gene_pvalue`, `gene_qvalue`, `method` and `scaled`.
#'
#' @references
#' Dividing out quantification uncertainty enables assessment of differential
#' transcript usage with limma and edgeR (2025). *Nucleic Acids Research*,
#' doi:10.1093/nar/gkaf1305.
#'
#' Gilis J et al. (2021). satuRn: scalable analysis of differential transcript
#' usage for bulk and single-cell RNA-sequencing applications.
#' *F1000Research* 10:374.
#'
#' Anders S, Reyes A, Huber W (2012). Detecting differential usage of exons
#' from RNA-seq data. *Genome Research* 22:2008-2017.
#'
#' @seealso [test_productive_fraction()], [test_apa()], [de_confidence()]
#'
#' @examples
#' iso <- example_isoformic()
#' if (requireNamespace("edgeR", quietly = TRUE)) {
#'   res <- test_usage(iso, group = "condition", method = "diffsplice")
#'   head(res)
#'   table(res$contrast)
#' }
#'
#' @rdname test_usage
#' @export
test_usage <- S7::new_generic(
  "test_usage",
  "self",
  function(
    self,
    group,
    levels = NULL,
    contrasts = NULL,
    method = "diffsplice",
    assay = "counts",
    scale = NULL,
    min_count = 10,
    min_prop = 0.1,
    workers = 1L,
    gene_name = NULL,
    ...
  ) {
    S7::S7_dispatch()
  }
)

S7::method(test_usage, IsoformicExperiment) <- function(
  self,
  group,
  levels = NULL,
  contrasts = NULL,
  method = "diffsplice",
  assay = "counts",
  scale = NULL,
  min_count = 10,
  min_prop = 0.1,
  workers = 1L,
  gene_name = NULL,
  ...
) {
  rlang::check_required(group)
  method <- rlang::arg_match(method, c("diffsplice", "satuRn", "DEXSeq"))
  check_usage_backend(method)
  assays <- attr(self, "assays", exact = TRUE)
  if (!assay %in% names(assays)) {
    cli::cli_abort(
      message = c(
        x = "Assay {.val {assay}} not found.",
        i = "Available assays: {.val {names(assays)}}."
      ),
      class = "isoformic_assay_not_found"
    )
  }
  counts <- as.matrix(assays[[assay]])
  setup <- group_samples(
    self,
    group,
    colnames(counts),
    levels = levels,
    min_levels = 2L
  )
  counts <- counts[, setup$keep, drop = FALSE]
  groups <- setup$groups
  omnibus <- identical(contrasts, "omnibus")
  contrasts <- resolve_contrasts(
    if (omnibus) NULL else contrasts,
    levels(groups)
  )

  has_draws <- !is.null(self@metadata$inf_reps)
  scale <- scale %||% has_draws
  if (isTRUE(scale)) {
    # Estimated on the samples of the test, like the test itself
    counts <- counts /
      overdispersion_for(self, rownames(counts), samples = colnames(counts))
  }

  tx <- self@annot_data_transcripts
  gene_of <- tx$gene_id[match(rownames(counts), tx$transcript_id)]
  known <- !is.na(gene_of)
  counts <- counts[known, , drop = FALSE]
  gene_of <- gene_of[known]
  kept <- filter_usage(counts, gene_of, groups, min_count, min_prop)
  if (!any(kept)) {
    cli::cli_abort(
      message = c(
        x = "No gene has two transcripts that pass the expression filters.",
        i = "Lower {.arg min_count} or {.arg min_prop}."
      ),
      class = "isoformic_no_usage_features"
    )
  }
  counts <- counts[kept, , drop = FALSE]
  gene_of <- gene_of[kept]

  res <- if (omnibus && identical(method, "DEXSeq")) {
    fit <- usage_dexseq_fit(
      counts,
      gene_of,
      groups,
      workers,
      fold_change = FALSE
    )
    fit$contrast <- "omnibus"
    fit
  } else {
    switch(
      method,
      diffsplice = usage_diffsplice(counts, gene_of, groups, contrasts),
      satuRn = usage_saturn(counts, gene_of, groups, contrasts, workers),
      DEXSeq = usage_dexseq(counts, gene_of, groups, contrasts, workers)
    )
  }
  if (omnibus) {
    res <- omnibus_usage(res)
    contrasts <- list(omnibus = "omnibus")
  }
  res$gene_id <- gene_of[match(res$transcript_id, rownames(counts))]
  res <- dplyr::bind_rows(lapply(split(res, res$contrast), adjust_usage))
  res <- res |>
    dplyr::left_join(
      tx[, c("transcript_id", "transcript_name", "transcript_type")],
      by = "transcript_id"
    ) |>
    dplyr::left_join(
      self@annot_data_genes[, c("gene_id", "gene_name")],
      by = "gene_id"
    ) |>
    dplyr::mutate(
      contrast = factor(.data$contrast, levels = names(contrasts)),
      method = method,
      scaled = isTRUE(scale)
    ) |>
    dplyr::select(
      "gene_id",
      "gene_name",
      "transcript_id",
      "transcript_name",
      "transcript_type",
      "contrast",
      "usage_log2FC",
      "pvalue",
      "qvalue",
      "gene_pvalue",
      "gene_qvalue",
      "method",
      "scaled"
    ) |>
    dplyr::arrange(.data$contrast, .data$pvalue) |>
    dplyr::mutate(contrast = as.character(.data$contrast))
  filter_genes(res, gene_name)
}

#' Packages needed by each usage-test backend
#' @keywords internal
#' @noRd
check_usage_backend <- function(method) {
  pkgs <- switch(
    method,
    diffsplice = c("edgeR", "limma"),
    satuRn = c("satuRn", "SummarizedExperiment", "BiocParallel"),
    DEXSeq = c("DEXSeq", "BiocParallel")
  )
  rlang::check_installed(
    pkgs,
    reason = paste0("to test transcript usage with ", method, ".")
  )
}

#' Transcripts with enough reads and enough share of their gene
#'
#' A transcript is kept when it has at least `min_count` reads and at least a
#' share of `min_prop` of the gene in as many samples as the smallest group,
#' and its gene keeps at least two transcripts.
#' @keywords internal
#' @noRd
filter_usage <- function(counts, gene_of, groups, min_count, min_prop) {
  n_min <- min(table(groups))
  totals <- rowsum(counts, gene_of)
  share <- counts / totals[match(gene_of, rownames(totals)), , drop = FALSE]
  share[is.nan(share)] <- 0
  keep <- rowSums(counts >= min_count & share >= min_prop) >= n_min
  n_kept <- table(gene_of[keep])
  keep & gene_of %in% names(n_kept)[n_kept >= 2L]
}

#' Combine the p-values of a gene with Simes' method
#' @keywords internal
#' @noRd
simes_pvalue <- function(p) {
  p <- sort(p[!is.na(p)])
  if (!length(p)) {
    return(NA_real_)
  }
  min(length(p) * p / seq_along(p))
}

#' The comparisons to test, as a named list of `c(numerator, denominator)`
#' @keywords internal
#' @noRd
resolve_contrasts <- function(contrasts, levels) {
  if (is.null(contrasts)) {
    pairs <- lapply(levels[-1L], function(l) c(l, levels[[1L]]))
  } else if (identical(contrasts, "pairwise")) {
    combos <- utils::combn(levels, 2L)
    pairs <- lapply(seq_len(ncol(combos)), function(k) rev(combos[, k]))
  } else {
    pairs <- if (is.character(contrasts)) list(contrasts) else contrasts
    valid <- is.list(pairs) &&
      length(pairs) > 0L &&
      all(vapply(
        pairs,
        function(p) {
          is.character(p) &&
            length(p) == 2L &&
            all(p %in% levels) &&
            p[[1L]] != p[[2L]]
        },
        logical(1L)
      ))
    if (!valid) {
      cli::cli_abort(
        message = c(
          x = "{.arg contrasts} must be {.val pairwise} or pairs of distinct levels.",
          i = "Available levels: {.val {levels}}."
        ),
        class = "isoformic_invalid_contrasts"
      )
    }
  }
  names(pairs) <- vapply(pairs, paste, character(1L), collapse = " vs ")
  if (anyDuplicated(names(pairs))) {
    pairs <- pairs[!duplicated(names(pairs))]
  }
  pairs
}

#' One omnibus p-value per transcript from the contrasts against the reference
#'
#' The Simes combination of the contrast p-values of each transcript; results
#' that already are omnibus (DEXSeq) pass through.
#' @keywords internal
#' @noRd
omnibus_usage <- function(res) {
  if (all(res$contrast == "omnibus")) {
    res$usage_log2FC <- NA_real_
    return(res)
  }
  by_tx <- split(res$pvalue, res$transcript_id)
  tibble::tibble(
    transcript_id = names(by_tx),
    contrast = "omnibus",
    usage_log2FC = NA_real_,
    pvalue = vapply(by_tx, simes_pvalue, numeric(1L), USE.NAMES = FALSE)
  )
}

#' Benjamini-Hochberg and Simes gene p-values within one contrast
#' @keywords internal
#' @noRd
adjust_usage <- function(res) {
  res$qvalue <- stats::p.adjust(res$pvalue, method = "BH")
  genes <- tibble::tibble(gene_id = unique(res$gene_id))
  genes$gene_pvalue <- vapply(
    genes$gene_id,
    function(g) simes_pvalue(res$pvalue[res$gene_id == g]),
    numeric(1L)
  )
  genes$gene_qvalue <- stats::p.adjust(genes$gene_pvalue, method = "BH")
  dplyr::left_join(res, genes, by = "gene_id")
}

#' Workers for the Bioconductor backends
#' @keywords internal
#' @noRd
usage_bpparam <- function(workers) {
  if (workers > 1L && .Platform$OS.type == "unix") {
    BiocParallel::MulticoreParam(workers)
  } else {
    BiocParallel::SerialParam()
  }
}

#' Contrast vector over the levels of a `~ 0 + group` design
#' @keywords internal
#' @noRd
contrast_vector <- function(pair, levels) {
  v <- numeric(length(levels))
  v[match(pair[[1L]], levels)] <- 1
  v[match(pair[[2L]], levels)] <- -1
  v
}

#' edgeR quasi-likelihood fits and diffSplice
#'
#' One fit with a coefficient per group; each contrast is tested on it.
#' @keywords internal
#' @noRd
usage_diffsplice <- function(counts, gene_of, groups, contrasts) {
  d <- edgeR::DGEList(
    counts,
    group = groups,
    genes = data.frame(GeneID = gene_of, row.names = rownames(counts))
  )
  d <- edgeR::normLibSizes(d)
  design <- stats::model.matrix(~ 0 + groups)
  colnames(design) <- paste0("group", seq_len(nlevels(groups)))
  d <- edgeR::estimateDisp(d, design)
  fit <- edgeR::glmQLFit(d, design)
  dplyr::bind_rows(lapply(names(contrasts), function(nm) {
    splice <- edgeR::diffSpliceDGE(
      fit,
      contrast = contrast_vector(contrasts[[nm]], levels(groups)),
      geneid = "GeneID",
      verbose = FALSE
    )
    top <- edgeR::topSpliceDGE(splice, test = "exon", number = Inf)
    tibble::tibble(
      transcript_id = rownames(top),
      contrast = nm,
      usage_log2FC = top$logFC,
      pvalue = top$P.Value
    )
  }))
}

#' satuRn quasi-binomial models of transcript usage
#'
#' One fit with a coefficient per group; the contrasts are the columns of one
#' contrast matrix.
#' @keywords internal
#' @noRd
usage_saturn <- function(counts, gene_of, groups, contrasts, workers) {
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts),
    rowData = data.frame(
      isoform_id = rownames(counts),
      gene_id = gene_of,
      row.names = rownames(counts)
    ),
    colData = data.frame(group = groups, row.names = colnames(counts))
  )
  se <- satuRn::fitDTU(
    object = se,
    formula = ~ 0 + group,
    BPPARAM = usage_bpparam(workers),
    verbose = FALSE
  )
  keys <- paste0("contrast", seq_along(contrasts))
  contrast <- vapply(
    contrasts,
    contrast_vector,
    numeric(nlevels(groups)),
    levels = levels(groups)
  )
  contrast <- matrix(
    contrast,
    ncol = length(contrasts),
    dimnames = list(paste0("group", levels(groups)), keys)
  )
  se <- satuRn::testDTU(
    object = se,
    contrasts = contrast,
    diagplot1 = FALSE,
    diagplot2 = FALSE,
    sort = FALSE
  )
  dplyr::bind_rows(lapply(seq_along(contrasts), function(k) {
    res <- as.data.frame(
      SummarizedExperiment::rowData(se)[[paste0("fitDTUResult_", keys[[k]])]]
    )
    tibble::tibble(
      transcript_id = rownames(res),
      contrast = names(contrasts)[[k]],
      usage_log2FC = res$estimates / log(2),
      pvalue = res$pval
    )
  }))
}

#' DEXSeq negative binomial models of transcript usage
#'
#' DEXSeq tests a model with a group effect against one without, so each
#' contrast is fitted on the samples of its two groups.
#' @keywords internal
#' @noRd
usage_dexseq <- function(counts, gene_of, groups, contrasts, workers) {
  dplyr::bind_rows(lapply(names(contrasts), function(nm) {
    pair <- contrasts[[nm]]
    in_pair <- groups %in% pair
    res <- usage_dexseq_fit(
      counts[, in_pair, drop = FALSE],
      gene_of,
      factor(as.character(groups[in_pair]), levels = rev(pair)),
      workers
    )
    res$contrast <- nm
    res
  }))
}

#' DEXSeq on the given groups: the model with a group effect against the one
#' without; with two groups the fold-change is the second level against the
#' first
#' @keywords internal
#' @noRd
usage_dexseq_fit <- function(
  counts,
  gene_of,
  groups,
  workers,
  fold_change = TRUE
) {
  bp <- usage_bpparam(workers)
  dx <- suppressMessages(suppressWarnings(DEXSeq::DEXSeqDataSet(
    countData = round(counts),
    sampleData = data.frame(condition = groups, row.names = colnames(counts)),
    design = ~ sample + exon + condition:exon,
    featureID = rownames(counts),
    groupID = gene_of
  )))
  dx <- DEXSeq::estimateSizeFactors(dx)
  dx <- DEXSeq::estimateDispersions(dx, quiet = TRUE, BPPARAM = bp)
  dx <- DEXSeq::testForDEU(dx, reducedModel = ~ sample + exon, BPPARAM = bp)
  if (isTRUE(fold_change)) {
    dx <- DEXSeq::estimateExonFoldChanges(
      dx,
      fitExpToVar = "condition",
      BPPARAM = bp
    )
  }
  res <- DEXSeq::DEXSeqResults(dx)
  lfc <- grep("^log2fold", colnames(res), value = TRUE)
  tibble::tibble(
    transcript_id = as.character(res[["featureID"]]),
    usage_log2FC = if (length(lfc)) as.numeric(res[[lfc[[1L]]]]) else NA_real_,
    pvalue = as.numeric(res[["pvalue"]])
  )
}
