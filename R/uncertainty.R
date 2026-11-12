#' Quantify and report the uncertainty of isoform-level results
#'
#' @description
#' An isoform-level result is only as good as the quantification behind it.
#' Reads that map to several isoforms of a gene can be shared out between
#' them in more than one way, so the estimate of each isoform carries an
#' *inferential* uncertainty that is unrelated to biological variability.
#' Salmon measures it by drawing replicates from the posterior (Gibbs
#' sampling, `--numGibbsSamples`) and by writing the equivalence classes
#' (`--dumpEq`). These functions turn that information into columns that can
#' be filtered on and reported.
#'
#' Every isoform gets four kinds of evidence:
#'
#' * **Quantification uncertainty**: the inferential relative variance
#'   (`mean_inf_rv`, Zhu et al. 2019, as in `fishpond`), the variance of the
#'   replicates in excess of what counting noise alone would give, relative
#'   to the expression. Isoforms with a high value are not well determined by
#'   the reads. `uncertain` flags `mean_inf_rv >= inf_rv_cutoff`.
#' * **Ambiguity between isoforms**: the correlation between the replicates
#'   of isoforms of the same gene. A strong negative correlation means the
#'   reads cannot tell the two isoforms apart and every read taken from one
#'   is given to the other. Pairs at or below `cor_cutoff` are *ambiguous*
#'   and are joined into `ambiguity_group`s that can be tested as one
#'   feature. This is the idea of `mmcollapse` (Turro et al. 2014), with a
#'   fixed cutoff.
#' * **Read support**: the share of the fragments of an isoform that map to
#'   it alone (`unique_fraction`, from the equivalence classes).
#' * **Stability of the fold-change**: the contrast between two groups is
#'   repeated on every replicate. `lfc_lo`/`lfc_hi` is the 95% interval of
#'   the log2 fold-change over replicates, `sign_agreement` the share of
#'   replicates that agree on its direction and `frac_draws_sig` the share in
#'   which the isoform is called differentially expressed.
#'
#' `de_confidence()` combines the last three into one class per isoform:
#' `high` when the fold-change is stable (`sign_agreement >= min_agreement`),
#' the isoform is not `uncertain` and has no ambiguous sibling; `moderate`
#' when exactly one of these fails; `low` when two or more fail. After
#' `add_uncertainty(identifiability = TRUE)` a fourth condition counts: the
#' equivalence classes must identify the isoform (`identifiable`). The class
#' says how far a call on the isoform can be trusted, so read it for the
#' isoforms flagged by the primary test (`n_det_sig`, `confident_de`): an
#' isoform without a real change is `moderate` or `low` simply because its
#' fold-change is not stable.
#'
#' * `add_uncertainty(x, paths)` reads the replicates, equivalence classes and
#'   run metadata from the Salmon output directories `paths` (named by sample
#'   ID, in any order), computes the tables below and stores them in the
#'   object, so that they are written by [write_isoformic()]. When `group` is
#'   given it also stores the fold-change stability of that contrast.
#' * `tx_uncertainty(x)`: per-transcript table, for the transcripts whose
#'   draws are in the object (the expressed ones).
#' * `isoform_ambiguity(x)`: per-pair table of replicate correlations.
#' * `de_confidence(x, group)`: per-feature table of fold-change stability
#'   and confidence class; with `collapse_groups = TRUE`, ambiguous isoforms
#'   are summed into one feature before testing.
#' * `uncertainty_report(x)`: a list of summary tables for a methods or
#'   results section.
#' * `uncertainty_methods(x)`: a ready-to-edit paragraph that states the
#'   settings used.
#'
#' The replicates themselves are kept in memory only (in
#' `metadata(x)$inf_reps`) and are not written by [write_isoformic()]; the
#' tables derived from them are. Draws are kept for the transcripts with an
#' estimated mean count of at least `min_count`; the others are treated as
#' unexpressed.
#'
#' @details
#' The thresholds are heuristics, not universal constants: `inf_rv_cutoff`
#' (0.5), `cor_cutoff` (-0.5) and `min_agreement` (0.95) are conservative
#' starting points. With 30 draws, sampling noise alone puts the inferential
#' relative variance of a well-determined, highly expressed transcript near
#' 0.1, so cutoffs much below 0.5 flag noise; poorly determined isoforms
#' typically reach values of several units. Report the numbers next to the
#' classes, and rerun with other cutoffs if the data call for it. The replicates measure the
#' uncertainty of the quantification only; they do not replace biological
#' replicates and cannot make a noisy contrast significant.
#'
#' What is not measured. The replicates describe how much an estimate moves
#' under the quantification model. They do not find the transcripts whose
#' abundance is not identifiable, that is, where several abundance vectors
#' explain the reads equally well (Zheng et al. 2022 report that only about
#' half of the transcripts flagged by replicate variance are also
#' non-identifiable in that sense), and they do not find an estimate that is
#' wrong because the annotation lacks the transcript the reads came from
#' (Ma and Kingsford 2019). An ambiguous pair is a statement about the
#' replicates, not a proof of non-identifiability.
#'
#' How many replicates. The inferential relative variance and the
#' correlations are computed per sample, from the replicates of that sample.
#' The overdispersion used to scale counts (see [inferential_overdispersion()])
#' pools all samples, and its precision depends on the total
#' `n_samples * (n_replicates - 1)`: Baldoni, Chen and Smyth (2024) recommend
#' at least 60 and prefer 200 to 300, and find that Gibbs samples resolve
#' low counts better than bootstrap replicates. `add_uncertainty()` says so
#' when the total is below 60.
#'
#' Fold-changes are computed on counts scaled to a common library size, as
#' `log2(mean_b + pseudocount) - log2(mean_a + pseudocount)` with a pseudocount
#' of 5 (the `swish` convention), and the per-replicate test is a Welch
#' t-test on the log scale with Benjamini-Hochberg adjustment across features.
#' This is a screening statistic for stability, not a replacement for the
#' primary test (`swish`, `DESeq2`, ...), whose results are shown next to it.
#' When the sample table (`col_data`) has a `library_size` column, with the
#' total estimated count of the full quantification, it is used for the scaling.
#' Set it when the object holds only a subset of the transcripts, so that the
#' subset does not change the scaling.
#'
#' @param self An [IsoformicExperiment].
#' @param paths Salmon output directories (see [read_inferential_replicates()]),
#'   named by sample ID. `NULL` uses the replicates already stored in the
#'   object.
#' @param group Name of a `col_data()` column with the sample groups.
#' @param levels The two group levels to contrast, reference first. Samples
#'   of other groups are left out of the contrast. `NULL` uses the levels of
#'   `group`, which must then be exactly two.
#' @param gene_name Optional gene names or IDs to keep.
#' @param n_draws Number of replicates to read; `NULL` reads all of them.
#' @param min_count Minimum mean estimated count for a transcript to be
#'   assessed.
#' @param inf_rv_cutoff Inferential relative variance from which an isoform is
#'   called uncertain. `NULL` uses the value the tables were computed with.
#' @param cor_cutoff Mean replicate correlation at or below which two isoforms
#'   of a gene are called ambiguous: a number, or `"adaptive"` to derive it
#'   from the data as `mmcollapse` does (minus the 97.5th percentile of the
#'   largest mean correlation of each isoform, which treats the positive
#'   correlations as noise). The value used is stored and stated by
#'   `uncertainty_methods()`.
#' @param min_agreement Share of replicates that must agree on the direction
#'   of a fold-change for it to be called stable. `NULL` uses the stored
#'   value.
#' @param pseudocount Pseudocount of the log2 fold-change (5, as in `swish`).
#' @param equivalence Logical, read the equivalence classes when present, or
#'   the table returned by [equivalence_support()], to attach read support when
#'   the draws are supplied without `paths`.
#' @param identifiability Logical, also compute which transcripts the
#'   equivalence classes can tell apart ([equivalence_identifiability()];
#'   minutes per sample, needs `paths` and, for human samples, the `Matrix`
#'   package), or the table that function returned. `de_confidence()` then
#'   counts a non-identifiable transcript as one more reason against the
#'   call.
#' @param run_qc Logical, add the Salmon run metadata (see
#'   [salmon_run_qc()]) to `col_data()`.
#' @param keep_draws Logical, keep the replicates in the object (needed to
#'   recompute the tables and for [productive_fraction()] uncertainty).
#' @param groups External groups of transcripts that the reads cannot tell
#'   apart, replacing the groups built from the correlation of the replicates:
#'   a table with `transcript_id` and `group_id`, or the path to Terminus output
#'   (see [read_terminus_groups()]), or the cut of TreeTerminus trees (see
#'   [cut_transcript_trees()]). The confidence classes and
#'   `collapse_groups` then follow these groups; the pairs of
#'   [isoform_ambiguity()] keep describing the correlation.
#' @param collapse_groups Logical, sum ambiguous isoforms into one feature
#'   before testing the contrast.
#' @param ambiguous_only Logical, `isoform_ambiguity()` returns only the
#'   ambiguous pairs.
#' @param ... Unused.
#'
#' @returns
#' * `add_uncertainty()`: the object with `dea$uncertainty`, `dea$ambiguity`
#'   and, when `group` is given, `dea$confidence` filled in.
#' * `tx_uncertainty()`: a tibble with `transcript_id`, `gene_id`, names and
#'   type, `mean_count`, `mean_inf_rv`, `max_inf_rv`, `mean_inf_cv`,
#'   `uncertain`, `min_sibling_cor`, `n_ambiguous_partners`,
#'   `ambiguity_group` and, when equivalence classes were read, `n_classes`,
#'   `n_unique_fragments`, `unique_fraction` and `low_information` (`TRUE`
#'   for an uncertain isoform without a single fragment of its own, the
#'   isoforms `mmcollapse` selects for collapsing). With `identifiability`,
#'   also `identifiable_fraction` and `identifiable`.
#' * `isoform_ambiguity()`: a tibble of isoform pairs of the same gene with
#'   `mean_cor`, `min_cor`, `n_samples`, `collapse_gain` and `ambiguous`.
#'   `collapse_gain` is the score Terminus groups on (Sarkar et al. 2020):
#'   the inferential relative variance of the sum of the two isoforms minus
#'   the mean of their own, averaged over samples. The more negative, the
#'   more is gained by testing the pair as one feature.
#' * `de_confidence()`: a tibble with one row per feature: `feature_id`,
#'   `feature_type` (`"transcript"` or `"group"`), gene and transcript names,
#'   the interval and stability of the fold-change, `mean_inf_rv`,
#'   `ambiguous`, `identifiable` (only when the identifiability was
#'   computed; `NA` for groups), `confidence` (an ordered factor), `reason`,
#'   the existing
#'   differential expression result (`det_log2FC`, `det_qvalue`,
#'   `n_det_sig`) and `confident_de`, which is `TRUE` when a significant
#'   isoform has `high` confidence.
#' * `uncertainty_report()`: a list with `settings`, `replicates` (replicate
#'   type, number of draws, degrees of freedom over the study and, when the
#'   draws are in the object, the median overdispersion, which is close to 1
#'   for long reads), `run`, `transcripts` (by transcript type), `genes`,
#'   `de` and `switches` tibbles.
#' * `uncertainty_methods()`: a character string.
#'
#' @references
#' Zhu, A., Srivastava, A., Ibrahim, J.G., Patro, R., Love, M.I. (2019).
#' Nonparametric expression analysis using inferential replicate counts.
#' *Nucleic Acids Research* 47(18), e105.
#'
#' Sarkar, H., Srivastava, A., Bravo, H.C., Love, M.I., Patro, R. (2020).
#' Terminus enables the discovery of data-driven, robust transcript groups
#' from RNA-seq data. *Bioinformatics* 36(Suppl 1), i102-i110.
#'
#' Turro, E., Astle, W.J., Tavare, S. (2014). Flexible analysis of RNA-seq
#' data using mixed effects models. *Bioinformatics* 30(2), 180-188.
#'
#' Zheng, H., Ma, C., Kingsford, C. (2022). Deriving ranges of optimal
#' estimated transcript expression due to nonidentifiability. *Journal of
#' Computational Biology* 29(2), 121-139.
#'
#' Ma, C., Kingsford, C. (2019). Detecting, categorizing, and correcting
#' coverage anomalies of RNA-seq quantification. *Cell Systems* 9(6),
#' 589-599.
#'
#' Baldoni, P.L., Chen, L., Smyth, G.K. (2024). Faster and more accurate
#' assessment of differential transcript expression with Gibbs sampling and
#' edgeR v4. *NAR Genomics and Bioinformatics* 6(4), lqae151.
#'
#' @seealso [read_inferential_replicates()], [plot_de_confidence()],
#'   [isoform_switches()]
#'
#' @examples
#' quant_dirs <- list.files(
#'   system.file("extdata", "mini_quants", package = "isoformic"),
#'   full.names = TRUE
#' )
#' names(quant_dirs) <- sub("_quant$", "", basename(quant_dirs))
#' annot <- system.file("extdata", "mini_annotation", package = "isoformic")
#' col_data <- data.frame(
#'   sample_id = names(quant_dirs),
#'   condition = rep(c("control", "treatment"), each = 3)
#' )
#' iso <- read_salmon(quant_dirs, annot_path = annot, col_data = col_data)
#'
#' # With your own data, `read_salmon(dirs, annot, inferential = TRUE)` or
#' # `add_uncertainty(iso, paths = dirs)` read the replicates that
#' # `salmon quant --numGibbsSamples 30` wrote; here bundled draws are used
#' draws <- as_inferential_replicates(read_parquet(
#'   system.file("extdata", "mini_draws.parquet", package = "isoformic")
#' ))
#' iso@metadata$inf_reps <- draws
#' iso <- add_uncertainty(iso, group = "condition")
#' head(tx_uncertainty(iso))
#' isoform_ambiguity(iso, ambiguous_only = TRUE)
#' head(de_confidence(iso))
#' uncertainty_report(iso)$transcripts
#'
#' @name uncertainty
#' @aliases add_uncertainty tx_uncertainty isoform_ambiguity de_confidence
#' @aliases uncertainty_report uncertainty_methods
NULL

uncertainty_defaults <- function() {
  list(
    min_count = 1,
    inf_rv_cutoff = 0.5,
    cor_cutoff = -0.5,
    min_agreement = 0.95,
    pseudocount = 5,
    shift = 0.01
  )
}

# Generics ---------------------------------------------------------------------

#' @rdname uncertainty
#' @export
add_uncertainty <- S7::new_generic(
  "add_uncertainty",
  "self",
  function(
    self,
    paths = NULL,
    group = NULL,
    levels = NULL,
    n_draws = NULL,
    min_count = 1,
    inf_rv_cutoff = 0.5,
    cor_cutoff = -0.5,
    min_agreement = 0.95,
    equivalence = TRUE,
    identifiability = FALSE,
    run_qc = TRUE,
    keep_draws = TRUE,
    collapse_groups = FALSE,
    groups = NULL,
    ...
  ) {
    S7::S7_dispatch()
  }
)

S7::method(add_uncertainty, IsoformicExperiment) <- function(
  self,
  paths = NULL,
  group = NULL,
  levels = NULL,
  n_draws = NULL,
  min_count = 1,
  inf_rv_cutoff = 0.5,
  cor_cutoff = -0.5,
  min_agreement = 0.95,
  equivalence = TRUE,
  identifiability = FALSE,
  run_qc = TRUE,
  keep_draws = TRUE,
  collapse_groups = FALSE,
  groups = NULL,
  ...
) {
  if (is.null(group)) {
    # `levels` without `group` is an error
    group_samples(self, NULL, character(), levels = levels)
  }
  groups <- resolve_groups(groups)
  samples <- object_samples(self)
  equiv <- NULL
  tool <- NULL
  if (!is.null(paths)) {
    dirs <- resolve_quant_dirs(paths)
    tool <- dirs$kind[[1L]]
    if (!setequal(dirs$sample_id, samples)) {
      cli::cli_abort(
        message = c(
          x = "The sample IDs of {.arg paths} do not match the samples of the object.",
          i = "Name {.arg paths} with {.val {samples}}."
        ),
        class = "isoformic_sample_mismatch"
      )
    }
    dirs <- dirs[match(samples, dirs$sample_id), ]
    named <- stats::setNames(dirs$dir, dirs$sample_id)
    keep_tx <- NULL
    if ("counts" %in% assay_names(self)) {
      counts <- assay(self, "counts")
      keep_tx <- rownames(counts)[rowMeans(counts) >= min_count]
    }
    reps <- read_inferential_replicates(
      named,
      transcripts = keep_tx,
      n_draws = n_draws
    )
    if (is.data.frame(equivalence)) {
      equiv <- equivalence
    } else if (isTRUE(equivalence)) {
      equiv <- rlang::try_fetch(
        equivalence_support(named, transcripts = rownames(reps[[1L]])),
        isoformic_no_equivalence_classes = function(cnd) {
          cli::cli_inform(c(
            i = "No equivalence classes found; run {.code salmon quant --dumpEq} to add read support."
          ))
          NULL
        }
      )
    }
    if (isTRUE(identifiability)) {
      identifiability <- equivalence_identifiability(
        named,
        transcripts = rownames(reps[[1L]])
      )
    }
    if (isTRUE(run_qc)) {
      qc <- salmon_run_qc(named)
      qc_cols <- setdiff(names(qc), "sample_id")
      qc_prefix <- paste0(attr(qc, "tool") %||% "salmon", "_")
      names(qc)[match(qc_cols, names(qc))] <- ifelse(
        startsWith(qc_cols, qc_prefix),
        qc_cols,
        paste0(qc_prefix, qc_cols)
      )
      col_data <- as.data.frame(self@col_data)
      col_data <- col_data[,
        setdiff(colnames(col_data), colnames(qc)[-1L]),
        drop = FALSE
      ]
      self@col_data <- dplyr::left_join(col_data, qc, by = "sample_id")
    }
  } else {
    reps <- self@metadata$inf_reps
    equiv <- if (is.data.frame(equivalence)) {
      equivalence
    } else {
      stored_equivalence(self)
    }
    if (isTRUE(identifiability)) {
      cli::cli_abort(
        message = c(
          x = "{.code identifiability = TRUE} needs {.arg paths}.",
          i = "Or pass the table of {.fn equivalence_identifiability}."
        ),
        class = "isoformic_identifiability_needs_paths"
      )
    }
  }
  equiv <- join_identifiability(equiv, identifiability)
  arr <- draws_view(self, reps)
  settings <- list(
    min_count = min_count,
    inf_rv_cutoff = inf_rv_cutoff,
    cor_cutoff = cor_cutoff,
    min_agreement = min_agreement,
    pseudocount = uncertainty_defaults()$pseudocount,
    shift = uncertainty_defaults()$shift
  )
  tables <- build_uncertainty(self, arr, equiv, settings, groups)
  settings$cor_cutoff <- tables$cor_cutoff
  dea <- self@dea
  dea[["uncertainty"]] <- tables$transcripts
  dea[["ambiguity"]] <- tables$pairs
  dea[["confidence"]] <- NULL
  metadata <- self@metadata
  metadata$uncertainty_n_draws <- draws_dim(arr)[[3L]]
  metadata$uncertainty_resampling_df <- resampling_df(
    draws_dim(arr)[[2L]],
    draws_dim(arr)[[3L]]
  )
  metadata$uncertainty_tool <- tool %||% metadata$uncertainty_tool
  metadata$uncertainty_ambiguity_source <- if (is.null(groups)) {
    "correlation"
  } else {
    attr(groups, "source") %||% "external"
  }
  metadata$uncertainty_tree_cutoff <- attr(groups, "cutoff")
  metadata$uncertainty_replicate_type <- attr(reps, "samp_type") %||%
    metadata$uncertainty_replicate_type %||%
    NA_character_
  metadata$uncertainty_min_count <- min_count
  metadata$uncertainty_inf_rv_cutoff <- inf_rv_cutoff
  metadata$uncertainty_cor_cutoff <- tables$cor_cutoff
  metadata$uncertainty_cor_cutoff_rule <- if (is.numeric(cor_cutoff)) {
    "fixed"
  } else {
    "adaptive"
  }
  metadata$uncertainty_min_agreement <- min_agreement
  metadata$uncertainty_n_assessed <- sum(tables$transcripts$assessed)
  inform_resampling_df(
    metadata$uncertainty_resampling_df,
    metadata$uncertainty_replicate_type
  )
  if (isTRUE(keep_draws)) {
    metadata$inf_reps <- reps
  } else {
    metadata$inf_reps <- NULL
  }
  self@dea <- dea
  self@metadata <- metadata
  clear_cache(self)
  if (!is.null(group)) {
    conf <- compute_de_confidence(
      self,
      arr,
      tables,
      group = group,
      levels = levels,
      collapse_groups = collapse_groups,
      settings = settings
    )
    dea <- self@dea
    dea[["confidence"]] <- conf$table
    self@dea <- dea
    metadata <- self@metadata
    metadata$uncertainty_group <- group
    metadata$uncertainty_levels <- paste(conf$levels, collapse = "|")
    metadata$uncertainty_collapsed <- collapse_groups
    self@metadata <- metadata
  }
  self
}

#' @rdname uncertainty
#' @export
tx_uncertainty <- S7::new_generic(
  "tx_uncertainty",
  "self",
  function(self, gene_name = NULL, ...) {
    S7::S7_dispatch()
  }
)

S7::method(tx_uncertainty, IsoformicExperiment) <- function(
  self,
  gene_name = NULL,
  ...
) {
  res <- resolve_uncertainty(self)$transcripts
  tx <- self@annot_data_transcripts[, c(
    "transcript_id",
    "transcript_name",
    "transcript_type"
  )]
  genes <- self@annot_data_genes[, c("gene_id", "gene_name")]
  res <- res |>
    dplyr::left_join(tx, by = "transcript_id") |>
    dplyr::left_join(genes, by = "gene_id") |>
    dplyr::relocate(
      "gene_name",
      "transcript_name",
      "transcript_type",
      .after = "gene_id"
    )
  filter_genes(res, gene_name)
}

#' @rdname uncertainty
#' @export
isoform_ambiguity <- S7::new_generic(
  "isoform_ambiguity",
  "self",
  function(self, gene_name = NULL, ambiguous_only = FALSE, ...) {
    S7::S7_dispatch()
  }
)

S7::method(isoform_ambiguity, IsoformicExperiment) <- function(
  self,
  gene_name = NULL,
  ambiguous_only = FALSE,
  ...
) {
  res <- resolve_uncertainty(self)$pairs
  tx <- self@annot_data_transcripts[, c("transcript_id", "transcript_name")]
  res$transcript_name_1 <- tx$transcript_name[match(
    res$transcript_id_1,
    tx$transcript_id
  )]
  res$transcript_name_2 <- tx$transcript_name[match(
    res$transcript_id_2,
    tx$transcript_id
  )]
  res <- dplyr::left_join(
    res,
    self@annot_data_genes[, c("gene_id", "gene_name")],
    by = "gene_id"
  ) |>
    dplyr::relocate("gene_name", .after = "gene_id") |>
    dplyr::arrange(.data$mean_cor)
  if (isTRUE(ambiguous_only)) {
    res <- res[res$ambiguous %in% TRUE, ]
  }
  filter_genes(res, gene_name)
}

#' @rdname uncertainty
#' @export
de_confidence <- S7::new_generic(
  "de_confidence",
  "self",
  function(
    self,
    group = NULL,
    levels = NULL,
    collapse_groups = FALSE,
    gene_name = NULL,
    pseudocount = 5,
    inf_rv_cutoff = NULL,
    min_agreement = NULL,
    ...
  ) {
    S7::S7_dispatch()
  }
)

S7::method(de_confidence, IsoformicExperiment) <- function(
  self,
  group = NULL,
  levels = NULL,
  collapse_groups = FALSE,
  gene_name = NULL,
  pseudocount = 5,
  inf_rv_cutoff = NULL,
  min_agreement = NULL,
  ...
) {
  if (is.null(group)) {
    group_samples(self, NULL, character(), levels = levels)
    res <- self@dea[["confidence"]]
    if (is.null(res)) {
      cli::cli_abort(
        message = c(
          x = "The object has no stored differential expression confidence.",
          i = "Compute it with {.code de_confidence(x, group = )} or {.code add_uncertainty(x, group = )}."
        ),
        class = "isoformic_no_uncertainty"
      )
    }
    return(filter_genes(res, gene_name))
  }
  settings <- uncertainty_settings(self)
  settings$inf_rv_cutoff <- inf_rv_cutoff %||% settings$inf_rv_cutoff
  settings$min_agreement <- min_agreement %||% settings$min_agreement
  settings$pseudocount <- pseudocount
  arr <- draws_view(self)
  tables <- resolve_uncertainty(self, arr = arr, settings = settings)
  conf <- compute_de_confidence(
    self,
    arr,
    tables,
    group = group,
    levels = levels,
    collapse_groups = collapse_groups,
    settings = settings
  )
  filter_genes(conf$table, gene_name)
}

# Computation ------------------------------------------------------------------

#' Sample IDs of the assays of an object
#' @keywords internal
#' @noRd
object_samples <- function(self) {
  assays <- attr(self, "assays", exact = TRUE)
  if (!length(assays)) {
    cli::cli_abort(
      message = c(x = "The object has no assays."),
      class = "isoformic_no_assays"
    )
  }
  colnames(assays[[1L]])
}

#' Stored settings of the uncertainty tables, with package defaults
#' @keywords internal
#' @noRd
uncertainty_settings <- function(self) {
  defaults <- uncertainty_defaults()
  meta <- self@metadata
  pick <- function(name) {
    value <- suppressWarnings(as.numeric(meta[[paste0("uncertainty_", name)]]))
    if (length(value) != 1L || is.na(value)) defaults[[name]] else value
  }
  list(
    min_count = pick("min_count"),
    inf_rv_cutoff = pick("inf_rv_cutoff"),
    cor_cutoff = pick("cor_cutoff"),
    min_agreement = pick("min_agreement"),
    pseudocount = defaults$pseudocount,
    shift = defaults$shift
  )
}

#' Equivalence support already stored in the uncertainty table
#' @keywords internal
#' @noRd
stored_equivalence <- function(self) {
  unc <- self@dea[["uncertainty"]]
  value_cols <- c(
    "n_classes",
    "n_fragments",
    "n_unique_fragments",
    "unique_fraction",
    "n_samples",
    "identifiable_fraction",
    "identifiable"
  )
  if (
    is.null(unc) ||
      !any(c("unique_fraction", "identifiable") %in% colnames(unc))
  ) {
    return(NULL)
  }
  cols <- intersect(value_cols, colnames(unc))
  has_value <- rowSums(!is.na(unc[, cols, drop = FALSE])) > 0L
  unc[has_value, c("transcript_id", cols)]
}

#' Add the identifiability of the transcripts to the read support table
#'
#' `identifiability` is the table of `equivalence_identifiability()`; anything
#' else leaves `equiv` as it is (its stored columns are kept).
#' @keywords internal
#' @noRd
join_identifiability <- function(equiv, identifiability) {
  if (!is.data.frame(identifiability)) {
    return(equiv)
  }
  cols <- c("transcript_id", "identifiable_fraction", "identifiable")
  if (!all(cols %in% colnames(identifiability))) {
    cli::cli_abort(
      message = c(
        x = "{.arg identifiability} must be the table of {.fn equivalence_identifiability}.",
        i = "It needs the columns {.val {cols}}."
      ),
      class = "isoformic_invalid_identifiability"
    )
  }
  ident <- tibble::as_tibble(identifiability[, cols])
  if (is.null(equiv)) {
    return(ident)
  }
  equiv <- equiv[, setdiff(colnames(equiv), cols[-1L]), drop = FALSE]
  dplyr::full_join(equiv, ident, by = "transcript_id")
}

#' Inferential replicates as a view of transcripts x samples x draws
#'
#' The replicates stay the list of matrices they are stored as: the view
#' holds the transcripts (those of the first replicate that the annotation
#' knows) and the samples of the object, with the rows and columns of each
#' replicate that hold them. Its consumers read one replicate
#' (`draw_slice()`) or one sample (`draws_sample()`) at a time, so the dense
#' array of all the replicates is never built.
#' @keywords internal
#' @noRd
draws_view <- function(self, inf_reps = NULL) {
  reps <- inf_reps %||% self@metadata$inf_reps
  if (is.null(reps) || !length(reps)) {
    cli::cli_abort(
      message = c(
        x = "The object has no inferential replicates.",
        i = "Attach them with {.code add_uncertainty(x, paths = )} or {.code read_salmon(inferential = TRUE)}."
      ),
      class = "isoformic_no_inferential_replicates"
    )
  }
  samp_type <- attr(reps, "samp_type")
  reps <- lapply(as.list(reps), function(m) {
    if (is.matrix(m)) m else as.matrix(m)
  })
  first <- reps[[1L]]
  samples <- object_samples(self)
  if (is.null(rownames(first)) || !all(samples %in% colnames(first))) {
    cli::cli_abort(
      message = c(
        x = paste0(
          "Each inferential replicate must be a matrix with transcript IDs ",
          "as row names and the sample IDs of the object as column names."
        )
      ),
      class = "isoformic_invalid_inf_reps"
    )
  }
  tx_known <- self@annot_data_transcripts$transcript_id
  ids <- rownames(first)[rownames(first) %in% tx_known]
  if (!length(ids)) {
    cli::cli_abort(
      message = c(
        x = "None of the replicate transcripts are present in the annotation."
      ),
      class = "isoformic_no_annotated_transcripts"
    )
  }
  # Replicates with the dimnames of the first one share its indices
  first_rows <- match(ids, rownames(first))
  first_cols <- match(samples, colnames(first))
  rows <- vector("list", length(reps))
  cols <- vector("list", length(reps))
  for (k in seq_along(reps)) {
    m <- reps[[k]]
    if (identical(rownames(m), rownames(first))) {
      rows[[k]] <- first_rows
    } else {
      rows[[k]] <- match(ids, rownames(m))
      if (anyNA(rows[[k]])) {
        cli::cli_abort(
          message = c(
            x = "Inferential replicates differ in their transcripts."
          ),
          class = "isoformic_invalid_inf_reps"
        )
      }
    }
    cols[[k]] <- if (identical(colnames(m), colnames(first))) {
      first_cols
    } else {
      match(samples, colnames(m))
    }
    if (anyNA(cols[[k]])) {
      cli::cli_abort(
        message = c(x = "Inferential replicates differ in their samples."),
        class = "isoformic_invalid_inf_reps"
      )
    }
  }
  structure(
    list(reps = reps, ids = ids, samples = samples, rows = rows, cols = cols),
    samp_type = samp_type,
    class = "isoformic_draws"
  )
}

#' Inferential replicates as an array of transcripts x samples x draws
#'
#' A dense copy of all the replicates: use `draws_view()` unless the array
#' itself is needed.
#' @keywords internal
#' @noRd
draws_array <- function(self, inf_reps = NULL) {
  view <- draws_view(self, inf_reps)
  size <- draws_dim(view)
  arr <- array(
    NA_real_,
    dim = size,
    dimnames = list(view$ids, view$samples, NULL)
  )
  for (k in seq_len(size[[3L]])) {
    arr[seq_len(size[[1L]]), seq_len(size[[2L]]), k] <- draw_slice(view, k)
  }
  attr(arr, "samp_type") <- attr(view, "samp_type")
  arr
}

#' Numbers of transcripts, samples and draws of a view or an array of draws
#' @keywords internal
#' @noRd
draws_dim <- function(x) {
  if (inherits(x, "isoformic_draws")) {
    return(c(length(x$ids), length(x$samples), length(x$reps)))
  }
  dim(x)
}

#' @keywords internal
#' @noRd
draws_ids <- function(x) {
  if (inherits(x, "isoformic_draws")) x$ids else rownames(x)
}

#' @keywords internal
#' @noRd
draws_samples <- function(x) {
  if (inherits(x, "isoformic_draws")) x$samples else dimnames(x)[[2L]]
}

#' The draws of the samples `keep` (logical, positions or IDs)
#' @keywords internal
#' @noRd
draws_select_samples <- function(x, keep) {
  if (!inherits(x, "isoformic_draws")) {
    return(x[, keep, , drop = FALSE])
  }
  pos <- if (is.character(keep)) {
    match(keep, x$samples)
  } else {
    seq_along(x$samples)[keep]
  }
  x$samples <- x$samples[pos]
  x$cols <- lapply(x$cols, function(idx) idx[pos])
  x
}

#' One draw as a transcript x sample matrix of doubles
#' @keywords internal
#' @noRd
draw_slice <- function(arr, k) {
  if (inherits(arr, "isoformic_draws")) {
    m <- arr$reps[[k]]
    rows <- arr$rows[[k]]
    cols <- arr$cols[[k]]
    if (
      length(rows) != nrow(m) ||
        length(cols) != ncol(m) ||
        !identical(rows, seq_len(nrow(m))) ||
        !identical(cols, seq_len(ncol(m)))
    ) {
      m <- m[rows, cols, drop = FALSE]
    }
    if (!is.double(m)) {
      storage.mode(m) <- "double"
    }
    return(m)
  }
  n_tx <- dim(arr)[1L]
  matrix(arr[seq_len(n_tx), seq_len(dim(arr)[2L]), k], nrow = n_tx)
}

#' The draws of one sample as a transcript x draw matrix
#'
#' `rows` are positions in the transcripts of the draws, `s` the position of
#' the sample.
#' @keywords internal
#' @noRd
draws_sample <- function(x, rows, s) {
  if (!inherits(x, "isoformic_draws")) {
    return(matrix(x[rows, s, ], nrow = length(rows)))
  }
  out <- matrix(0, nrow = length(rows), ncol = length(x$reps))
  for (k in seq_along(x$reps)) {
    out[, k] <- x$reps[[k]][x$rows[[k]][rows], x$cols[[k]][[s]]]
  }
  out
}

#' Per transcript and sample mean and variance over the draws
#' @keywords internal
#' @noRd
draw_moments <- function(arr) {
  size <- draws_dim(arr)
  slice_moments(
    function(k) draw_slice(arr, k),
    n_draws = size[[3L]],
    dimnames = list(draws_ids(arr), draws_samples(arr))
  )
}

#' Mean and variance over the draws of the matrices that `slice(k)` returns
#' @keywords internal
#' @noRd
slice_moments <- function(slice, n_draws, dimnames) {
  s1 <- matrix(0, length(dimnames[[1L]]), length(dimnames[[2L]]))
  s2 <- s1
  for (k in seq_len(n_draws)) {
    x <- unname(slice(k))
    s1 <- s1 + x
    s2 <- s2 + x * x
  }
  mean_draw <- s1 / n_draws
  var_draw <- if (n_draws > 1L) {
    pmax(s2 - n_draws * mean_draw^2, 0) / (n_draws - 1L)
  } else {
    mean_draw * 0
  }
  dimnames(mean_draw) <- dimnames(var_draw) <- dimnames
  list(mean = mean_draw, var = var_draw, n = n_draws)
}

#' Degrees of freedom of the replicates over the whole study
#'
#' `n_samples * (n_draws - 1)`, the quantity on which the precision of the
#' pooled overdispersion depends (Baldoni, Chen and Smyth 2024).
#' @keywords internal
#' @noRd
resampling_df <- function(n_samples, n_draws) {
  as.integer(n_samples) * (as.integer(n_draws) - 1L)
}

#' Say so when the replicates are too few for a precise overdispersion
#' @keywords internal
#' @noRd
inform_resampling_df <- function(
  df,
  replicate_type = NA_character_,
  min_df = 60L
) {
  if (is.null(df) || is.na(df) || df >= min_df) {
    return(invisible(FALSE))
  }
  cli::cli_inform(
    message = c(
      "!" = "The replicates give {df} degrees of freedom over the study (samples x (replicates - 1)).",
      i = "At least {min_df} are recommended for the overdispersion used to scale counts (200 to 300 preferred).",
      i = if (identical(replicate_type, "bootstrap")) {
        "Gibbs samples ({.code --numGibbsSamples}) resolve low counts better than bootstrap replicates."
      }
    ),
    class = "isoformic_few_inferential_replicates"
  )
  invisible(TRUE)
}

#' Inferential relative variance (Zhu et al. 2019)
#'
#' Same formula as `fishpond::computeInfRV()` and Terminus (Sarkar et al.
#' 2020): `max(var - mean, 0) / (mean + pseudocount) + shift`, so a transcript
#' whose draws vary no more than counts would has the value `shift`.
#' @keywords internal
#' @noRd
inferential_rv <- function(moments, pseudocount = 5, shift = 0.01) {
  pmax(moments$var - moments$mean, 0) / (moments$mean + pseudocount) + shift
}

#' Correlation of the draws of isoforms of the same gene
#'
#' Pairs are generated within genes only. The correlation across draws is
#' computed per sample from standardised draws and averaged over the samples
#' where both isoforms are expressed and vary. `collapse_gain` is the Terminus
#' score (Sarkar et al. 2020) averaged over the same samples: the InfRV of the
#' sum of the two isoforms minus the mean InfRV of the two. Negative values
#' mean that the sum is better determined than its members.
#' @keywords internal
#' @noRd
inferential_correlation <- function(
  arr,
  sample_means,
  gene_id,
  assessed,
  min_count,
  pseudocount = uncertainty_defaults()$pseudocount,
  shift = uncertainty_defaults()$shift
) {
  empty <- tibble::tibble(
    gene_id = character(),
    transcript_id_1 = character(),
    transcript_id_2 = character(),
    mean_cor = numeric(),
    min_cor = numeric(),
    n_samples = integer(),
    collapse_gain = numeric()
  )
  d <- draws_dim(arr)[[3L]]
  eligible <- which(assessed & !is.na(gene_id))
  if (d < 3L || length(eligible) < 2L) {
    return(empty)
  }
  ord <- eligible[order(gene_id[eligible])]
  g <- gene_id[ord]
  runs <- rle(g)
  end_of <- rep(cumsum(runs$lengths), runs$lengths)
  remaining <- end_of - seq_along(g)
  if (sum(remaining) == 0L) {
    return(empty)
  }
  first <- rep(seq_along(g), remaining)
  second <- first + sequence(remaining)
  n_pairs <- length(first)
  r_sum <- numeric(n_pairs)
  r_n <- integer(n_pairs)
  r_min <- rep(Inf, n_pairs)
  gain_sum <- numeric(n_pairs)
  chunk <- 50000L
  for (s in seq_len(draws_dim(arr)[[2L]])) {
    x <- draws_sample(arr, ord, s)
    mu <- rowMeans(x)
    centred <- x - mu
    sds <- sqrt(rowSums(centred^2) / (d - 1L))
    rv <- inferential_rv(list(mean = mu, var = sds^2), pseudocount, shift)
    z <- centred / sds
    z[sds == 0 | sample_means[ord, s] < min_count, ] <- NA_real_
    for (start in seq.int(1L, n_pairs, by = chunk)) {
      idx <- start:min(start + chunk - 1L, n_pairs)
      r <- rowSums(
        z[first[idx], , drop = FALSE] * z[second[idx], , drop = FALSE]
      ) /
        (d - 1L)
      ok <- !is.na(r)
      hit <- idx[ok]
      r_sum[hit] <- r_sum[hit] + r[ok]
      r_n[hit] <- r_n[hit] + 1L
      r_min[hit] <- pmin(r_min[hit], r[ok])
      # InfRV of the sum of the pair minus the mean InfRV of its members
      a <- first[hit]
      b <- second[hit]
      sum_moments <- list(
        mean = mu[a] + mu[b],
        var = pmax(sds[a]^2 + sds[b]^2 + 2 * r[ok] * sds[a] * sds[b], 0)
      )
      gain_sum[hit] <- gain_sum[hit] +
        inferential_rv(sum_moments, pseudocount, shift) -
        (rv[a] + rv[b]) / 2
    }
  }
  keep <- r_n > 0L
  ids <- draws_ids(arr)
  tibble::tibble(
    gene_id = g[first][keep],
    transcript_id_1 = ids[ord[first]][keep],
    transcript_id_2 = ids[ord[second]][keep],
    mean_cor = r_sum[keep] / r_n[keep],
    min_cor = r_min[keep],
    n_samples = r_n[keep],
    collapse_gain = gain_sum[keep] / r_n[keep]
  )
}

#' Connected components of the ambiguous pairs, named per gene
#' @keywords internal
#' @noRd
assign_ambiguity_groups <- function(pairs) {
  edges <- pairs[pairs$ambiguous %in% TRUE, ]
  none <- tibble::tibble(
    transcript_id = character(),
    ambiguity_group = character()
  )
  if (!nrow(edges)) {
    return(none)
  }
  ends <- c(edges$transcript_id_1, edges$transcript_id_2)
  nodes <- unique(ends)
  ia <- match(edges$transcript_id_1, nodes)
  ib <- match(edges$transcript_id_2, nodes)
  label <- seq_along(nodes)
  repeat {
    low <- pmin(label[ia], label[ib])
    best <- tapply(c(low, low), c(ia, ib), min)
    node <- as.integer(names(best))
    updated <- label
    updated[node] <- pmin(updated[node], as.vector(best))
    updated <- updated[updated]
    if (identical(updated, label)) {
      break
    }
    label <- updated
  }
  gene <- c(edges$gene_id, edges$gene_id)[match(nodes, ends)]
  rank <- stats::ave(label, gene, FUN = function(x) match(x, sort(unique(x))))
  tibble::tibble(
    transcript_id = nodes,
    ambiguity_group = paste0(gene, "_amb", rank)
  )
}

#' Build the per-transcript uncertainty table and the pair table
#' @keywords internal
#' @noRd
build_uncertainty <- function(
  self,
  arr,
  equiv,
  settings,
  external_groups = NULL
) {
  defaults <- uncertainty_defaults()
  tx <- self@annot_data_transcripts
  ids <- draws_ids(arr)
  gene_id <- tx$gene_id[match(ids, tx$transcript_id)]
  mom <- draw_moments(arr)
  mean_count <- unname(rowMeans(mom$mean))
  assessed <- mean_count >= settings$min_count
  rv <- inferential_rv(mom, defaults$pseudocount, defaults$shift)
  cv <- sqrt(mom$var) / mom$mean
  cv[!is.finite(cv)] <- NA_real_
  mean_cv <- unname(rowMeans(cv, na.rm = TRUE))
  mean_cv[is.nan(mean_cv)] <- NA_real_
  mean_rv <- unname(rowMeans(rv))
  max_rv <- unname(apply(rv, 1L, max))
  mean_rv[!assessed] <- NA_real_
  max_rv[!assessed] <- NA_real_
  mean_cv[!assessed] <- NA_real_
  pairs <- inferential_correlation(
    arr,
    mom$mean,
    gene_id,
    assessed,
    settings$min_count
  )
  cor_cutoff <- resolve_cor_cutoff(settings$cor_cutoff, pairs)
  pairs$ambiguous <- pairs$mean_cor <= cor_cutoff
  groups <- assign_ambiguity_groups(pairs)
  long <- rbind(
    data.frame(
      id = pairs$transcript_id_1,
      cor = pairs$mean_cor,
      amb = pairs$ambiguous
    ),
    data.frame(
      id = pairs$transcript_id_2,
      cor = pairs$mean_cor,
      amb = pairs$ambiguous
    )
  )
  min_cor <- tapply(long$cor, long$id, min)
  n_amb <- tapply(as.integer(long$amb), long$id, sum)
  res <- tibble::tibble(
    transcript_id = ids,
    gene_id = gene_id,
    mean_count = mean_count,
    assessed = assessed,
    mean_inf_rv = mean_rv,
    max_inf_rv = max_rv,
    mean_inf_cv = mean_cv,
    uncertain = ifelse(assessed, mean_rv >= settings$inf_rv_cutoff, NA),
    min_sibling_cor = as.vector(min_cor[ids]),
    n_ambiguous_partners = as.integer(ifelse(
      is.na(n_amb[ids]),
      0L,
      n_amb[ids]
    )),
    ambiguity_group = groups$ambiguity_group[match(ids, groups$transcript_id)]
  )
  res$n_ambiguous_partners[!assessed] <- NA_integer_
  if (!is.null(external_groups)) {
    # Groups from another tool replace the groups of correlated replicates
    res$ambiguity_group <- external_groups$group_id[match(
      ids,
      external_groups$transcript_id
    )]
  }
  if (!is.null(equiv)) {
    equiv <- equiv[, intersect(
      c(
        "transcript_id",
        "n_classes",
        "n_fragments",
        "n_unique_fragments",
        "unique_fraction",
        "identifiable_fraction",
        "identifiable"
      ),
      colnames(equiv)
    )]
    res <- dplyr::left_join(res, equiv, by = "transcript_id")
    if ("unique_fraction" %in% colnames(res)) {
      # No fragment of its own and a high spread of the replicates: the
      # selection of mmcollapse (Turro et al. 2014)
      res$low_information <- res$unique_fraction == 0 & res$uncertain
    }
  }
  list(transcripts = res, pairs = pairs, cor_cutoff = cor_cutoff)
}

#' Correlation cutoff, fixed or derived from the data
#'
#' `"adaptive"` is the rule of mmcollapse (Turro et al. 2014): the positive
#' correlations between isoforms are taken as noise, and a pair is ambiguous
#' when its correlation is more negative than minus the 97.5th percentile of
#' the largest mean correlation of each isoform.
#' @keywords internal
#' @noRd
resolve_cor_cutoff <- function(cor_cutoff, pairs, prob = 0.975) {
  fallback <- uncertainty_defaults()$cor_cutoff
  if (
    is.numeric(cor_cutoff) && length(cor_cutoff) == 1L && !is.na(cor_cutoff)
  ) {
    return(cor_cutoff)
  }
  if (!identical(cor_cutoff, "adaptive")) {
    cli::cli_abort(
      message = c(
        x = "{.arg cor_cutoff} must be one number or {.val adaptive}."
      ),
      class = "isoformic_invalid_cor_cutoff"
    )
  }
  top <- tapply(
    c(pairs$mean_cor, pairs$mean_cor),
    c(pairs$transcript_id_1, pairs$transcript_id_2),
    max
  )
  value <- if (length(top)) {
    unname(stats::quantile(top, prob, names = FALSE, na.rm = TRUE))
  } else {
    NA_real_
  }
  if (is.na(value) || value <= 0) {
    cli::cli_inform(
      message = c(
        "!" = "No adaptive correlation cutoff could be derived from the replicates.",
        i = "Using {fallback}."
      ),
      class = "isoformic_no_adaptive_cutoff"
    )
    return(fallback)
  }
  -min(value, 1)
}

#' Stored uncertainty tables, or computed from the draws
#' @keywords internal
#' @noRd
resolve_uncertainty <- function(self, arr = NULL, settings = NULL) {
  stored <- self@dea[["uncertainty"]]
  if (!is.null(stored) && !is.null(self@dea[["ambiguity"]])) {
    return(list(transcripts = stored, pairs = self@dea[["ambiguity"]]))
  }
  if (is.null(self@metadata$inf_reps) && is.null(arr)) {
    cli::cli_abort(
      message = c(
        x = "The object has no uncertainty tables.",
        i = "Compute them with {.code add_uncertainty(x, paths = )}."
      ),
      class = "isoformic_no_uncertainty"
    )
  }
  arr <- arr %||% draws_view(self)
  build_uncertainty(self, arr, NULL, settings %||% uncertainty_settings(self))
}

#' Row-wise quantiles of a matrix
#' @keywords internal
#' @noRd
row_quantiles <- function(m, probs) {
  q <- apply(m, 1L, stats::quantile, probs = probs, names = FALSE, na.rm = TRUE)
  matrix(q, nrow = length(probs), dimnames = list(NULL, rownames(m)))
}

#' Welch t statistics of a matrix of log values between two groups
#' @keywords internal
#' @noRd
welch_matrix <- function(x, in_b) {
  xa <- x[, !in_b, drop = FALSE]
  xb <- x[, in_b, drop = FALSE]
  na <- ncol(xa)
  nb <- ncol(xb)
  ma <- rowMeans(xa)
  mb <- rowMeans(xb)
  diff <- mb - ma
  if (na < 2L || nb < 2L) {
    return(list(lfc = diff, p = rep(NA_real_, nrow(x))))
  }
  va <- rowSums((xa - ma)^2) / (na - 1L)
  vb <- rowSums((xb - mb)^2) / (nb - 1L)
  se2 <- va / na + vb / nb
  stat <- diff / sqrt(se2)
  df <- se2^2 / ((va / na)^2 / (na - 1L) + (vb / nb)^2 / (nb - 1L))
  p <- 2 * stats::pt(-abs(stat), df)
  p[!is.finite(p) | se2 <= 0] <- NA_real_
  list(lfc = diff, p = p)
}

#' Repeat a two-group contrast on every draw and summarise its stability
#' @keywords internal
#' @noRd
compute_de_confidence <- function(
  self,
  arr,
  tables,
  group,
  levels,
  collapse_groups,
  settings
) {
  # The contrast uses the samples of the two levels; the uncertainty tables
  # in `tables` describe all the samples of the object
  setup <- two_group_setup(self, group, draws_samples(arr), levels)
  lv <- base::levels(setup$groups)
  arr <- draws_select_samples(arr, setup$keep)
  in_b <- as.character(setup$groups) == lv[2L]
  samples <- setup$samples
  d <- draws_dim(arr)[[3L]]
  tx_table <- tables$transcripts
  ids <- draws_ids(arr)
  amb_group <- tx_table$ambiguity_group[match(ids, tx_table$transcript_id)]
  # A group is collapsed when more than one of its members is in the draws; a
  # member left alone (after `x[i, ]`) stays a transcript, flagged as ambiguous
  collapsible <- isTRUE(collapse_groups) &
    !is.na(amb_group) &
    amb_group %in% amb_group[duplicated(amb_group)]
  feature_id <- ifelse(collapsible, amb_group, ids)
  # Features: transcripts, or ambiguous groups summed over their members, in
  # the order of their first transcript. One draw at a time.
  features <- unique(feature_id)
  # Transcripts on several rows of the replicates are summed as well
  feature_slice <- if (anyDuplicated(feature_id) > 0L) {
    function(k) {
      unname(rowsum(unname(draw_slice(arr, k)), feature_id, reorder = FALSE))
    }
  } else {
    function(k) unname(draw_slice(arr, k))
  }
  mom <- slice_moments(feature_slice, d, list(features, samples))
  mean_count <- unname(rowMeans(mom$mean))
  rv <- inferential_rv(
    mom,
    uncertainty_defaults()$pseudocount,
    uncertainty_defaults()$shift
  )
  mean_rv <- unname(rowMeans(rv))
  assessed <- mean_count >= settings$min_count
  # Library-size scaling: the `library_size` column of the sample table (for
  # objects holding a subset of the transcripts), else the counts assay, else
  # the draws
  libsize <- NULL
  sample_table <- col_data(self)
  if ("library_size" %in% colnames(sample_table)) {
    libsize <- sample_table$library_size[match(samples, sample_table$sample_id)]
  } else if ("counts" %in% assay_names(self)) {
    libsize <- colSums(assay(self, "counts"))[samples]
  }
  if (is.null(libsize) || anyNA(libsize) || any(libsize <= 0)) {
    libsize <- colSums(mom$mean)
  }
  scale <- mean(libsize) / libsize
  cutoffs <- de_cutoffs(self)
  lfc <- matrix(
    NA_real_,
    nrow = length(features),
    ncol = d,
    dimnames = list(features, NULL)
  )
  sig <- matrix(FALSE, nrow = length(features), ncol = d)
  for (k in seq_len(d)) {
    x <- feature_slice(k)
    x <- log2(sweep(x, 2L, scale, "*") + settings$pseudocount)
    w <- welch_matrix(x, in_b)
    lfc[, k] <- w$lfc
    padj <- stats::p.adjust(w$p, method = "BH")
    sig[, k] <- !is.na(padj) &
      padj < cutoffs$pvalue &
      abs(w$lfc) >= cutoffs$log2fc
  }
  med <- unname(apply(lfc, 1L, stats::median))
  q <- row_quantiles(lfc, c(0.025, 0.975))
  agreement <- unname(rowMeans(sign(lfc) == sign(med)))
  effect <- unname(rowMeans(
    sign(lfc) == sign(med) & abs(lfc) >= cutoffs$log2fc
  ))
  stable <- agreement >= settings$min_agreement
  member_tx <- split(ids, factor(feature_id, levels = features))
  is_group <- unname(
    lengths(member_tx) > 1L & isTRUE(collapse_groups)
  )
  gene_of <- self@annot_data_transcripts$gene_id[
    match(
      vapply(member_tx, function(x) x[[1L]], character(1L)),
      self@annot_data_transcripts$transcript_id
    )
  ]
  # Collapsed groups are not transcript IDs, so they are never matched here
  ambiguous <- features %in% ids[!is.na(amb_group)] & !is_group
  uncertain <- mean_rv >= settings$inf_rv_cutoff
  # Structural identifiability, when add_uncertainty() computed it. A group
  # is not assessed: summing its members is the remedy
  has_ident <- "identifiable" %in% colnames(tx_table)
  identifiable <- if (has_ident) {
    tx_table$identifiable[match(features, tx_table$transcript_id)]
  } else {
    rep(NA, length(features))
  }
  identifiable[is_group] <- NA
  not_identifiable <- identifiable %in% FALSE
  n_fail <- as.integer(!stable) +
    as.integer(uncertain) +
    as.integer(ambiguous) +
    as.integer(not_identifiable)
  confidence <- c("high", "moderate", "low")[pmin(n_fail, 2L) + 1L]
  confidence[!assessed] <- NA_character_
  reason <- vapply(
    seq_along(features),
    function(i) {
      why <- c(
        if (!stable[i]) "unstable fold-change",
        if (uncertain[i]) "uncertain quantification",
        if (ambiguous[i]) "ambiguous with a sibling isoform",
        if (not_identifiable[i]) "not identifiable from the equivalence classes"
      )
      if (length(why)) paste(why, collapse = "; ") else NA_character_
    },
    character(1L)
  )
  tx <- self@annot_data_transcripts
  first_tx <- vapply(member_tx, function(x) x[[1L]], character(1L))
  res <- tibble::tibble(
    feature_id = features,
    feature_type = ifelse(is_group, "group", "transcript"),
    gene_id = gene_of,
    transcript_ids = vapply(member_tx, paste, character(1L), collapse = ";"),
    transcript_name = ifelse(
      is_group,
      NA_character_,
      tx$transcript_name[match(first_tx, tx$transcript_id)]
    ),
    transcript_type = ifelse(
      is_group,
      NA_character_,
      tx$transcript_type[match(first_tx, tx$transcript_id)]
    ),
    mean_count = mean_count,
    lfc_median = med,
    lfc_lo = q[1L, ],
    lfc_hi = q[2L, ],
    lfc_sd = unname(apply(lfc, 1L, stats::sd)),
    sign_agreement = agreement,
    frac_draws_effect = effect,
    frac_draws_sig = unname(rowMeans(sig)),
    mean_inf_rv = mean_rv,
    ambiguous = ambiguous,
    confidence = factor(
      confidence,
      levels = c("high", "moderate", "low"),
      ordered = TRUE
    ),
    reason = reason
  )
  if (has_ident) {
    res <- tibble::add_column(
      res,
      identifiable = identifiable,
      .after = "ambiguous"
    )
  }
  res$gene_name <- self@annot_data_genes$gene_name[match(
    res$gene_id,
    self@annot_data_genes$gene_id
  )]
  res <- dplyr::relocate(res, "gene_name", .after = "gene_id")
  res[
    !assessed,
    c(
      "lfc_median",
      "lfc_lo",
      "lfc_hi",
      "lfc_sd",
      "sign_agreement",
      "frac_draws_effect",
      "frac_draws_sig"
    )
  ] <-
    NA_real_
  # Existing differential expression results next to the stability
  det <- self@dea[["det"]]
  res$det_log2FC <- NA_real_
  res$det_qvalue <- NA_real_
  res$n_det_sig <- NA_integer_
  if (!is.null(det)) {
    flagged <- de_tx(self)
    is_sig <- flagged$det_sig == "yes"
    pos <- match(flagged$transcript_id, ids)
    ok <- !is.na(pos) & is_sig
    res$n_det_sig <- tabulate(
      match(feature_id[pos[ok]], features),
      length(features)
    )
    single <- res$feature_type == "transcript"
    idx <- match(res$feature_id[single], flagged$transcript_id)
    res$det_log2FC[single] <- flagged$log2FC[idx]
    if ("qvalue" %in% colnames(flagged)) {
      res$det_qvalue[single] <- flagged$qvalue[idx]
    }
  }
  res$confident_de <- ifelse(
    is.na(res$n_det_sig),
    NA,
    res$n_det_sig > 0L & res$confidence %in% "high"
  )
  res <- res[order(-abs(res$lfc_median), na.last = TRUE), ]
  list(table = res, levels = lv)
}

# Reporting --------------------------------------------------------------------

#' @rdname uncertainty
#' @export
uncertainty_report <- S7::new_generic(
  "uncertainty_report",
  "self",
  function(self, ...) {
    S7::S7_dispatch()
  }
)

S7::method(uncertainty_report, IsoformicExperiment) <- function(self, ...) {
  tables <- resolve_uncertainty(self)
  unc <- tables$transcripts
  tx <- self@annot_data_transcripts[, c("transcript_id", "transcript_type")]
  unc$transcript_type <- tx$transcript_type[match(
    unc$transcript_id,
    tx$transcript_id
  )]
  settings <- uncertainty_settings(self)
  run_cols <- grep("^(salmon|oarfish)_", colnames(self@col_data), value = TRUE)
  run <- if (length(run_cols)) {
    self@col_data[, c("sample_id", run_cols)]
  } else {
    NULL
  }
  assessed <- unc[unc$assessed %in% TRUE, ]
  by_type <- assessed |>
    dplyr::group_by(.data$transcript_type) |>
    dplyr::summarise(
      n_assessed = dplyr::n(),
      n_uncertain = sum(.data$uncertain, na.rm = TRUE),
      pct_uncertain = 100 * mean(.data$uncertain, na.rm = TRUE),
      median_inf_rv = stats::median(.data$mean_inf_rv, na.rm = TRUE),
      n_ambiguous = sum(!is.na(.data$ambiguity_group)),
      median_unique_fraction = if ("unique_fraction" %in% colnames(assessed)) {
        stats::median(.data$unique_fraction, na.rm = TRUE)
      } else {
        NA_real_
      },
      .groups = "drop"
    ) |>
    dplyr::arrange(dplyr::desc(.data$n_assessed))
  by_gene <- assessed |>
    dplyr::group_by(.data$gene_id) |>
    dplyr::summarise(
      n_assessed = dplyr::n(),
      n_uncertain = sum(.data$uncertain, na.rm = TRUE),
      n_ambiguous_groups = dplyr::n_distinct(
        .data$ambiguity_group,
        na.rm = TRUE
      ),
      .groups = "drop"
    ) |>
    dplyr::left_join(
      self@annot_data_genes[, c("gene_id", "gene_name")],
      by = "gene_id"
    ) |>
    dplyr::relocate("gene_name", .after = "gene_id") |>
    dplyr::arrange(
      dplyr::desc(.data$n_ambiguous_groups),
      dplyr::desc(.data$n_uncertain)
    )
  conf <- self@dea[["confidence"]]
  de <- NULL
  switches <- NULL
  if (!is.null(conf) && !all(is.na(conf$n_det_sig))) {
    sig <- conf[conf$feature_type == "transcript" & conf$n_det_sig %in% 1L, ]
    de <- tibble::tibble(
      confidence = factor(
        c("high", "moderate", "low"),
        levels = c("high", "moderate", "low"),
        ordered = TRUE
      ),
      n_significant = vapply(
        c("high", "moderate", "low"),
        function(x) sum(sig$confidence %in% x),
        integer(1L),
        USE.NAMES = FALSE
      ),
      n_unstable = vapply(
        c("high", "moderate", "low"),
        function(x) {
          sum(
            sig$confidence %in% x & sig$sign_agreement < settings$min_agreement
          )
        },
        integer(1L),
        USE.NAMES = FALSE
      ),
      n_uncertain = vapply(
        c("high", "moderate", "low"),
        function(x) {
          sum(sig$confidence %in% x & sig$mean_inf_rv >= settings$inf_rv_cutoff)
        },
        integer(1L),
        USE.NAMES = FALSE
      ),
      n_ambiguous = vapply(
        c("high", "moderate", "low"),
        function(x) sum(sig$confidence %in% x & sig$ambiguous),
        integer(1L),
        USE.NAMES = FALSE
      )
    )
    if (!is.null(self@dea[["det"]]) && !is.null(self@dea[["deg"]])) {
      all_sw <- isoform_switches(self)
      conf_sw <- isoform_switches(self, confident_only = TRUE)
      switches <- tibble::tibble(
        set = c("all significant isoforms", "confident isoforms only"),
        n_genes = c(nrow(all_sw), nrow(conf_sw)),
        n_opposite_direction = c(
          sum(all_sw$switch_type == "opposite_direction"),
          sum(conf_sw$switch_type == "opposite_direction")
        ),
        n_det_not_deg = c(
          sum(all_sw$switch_type == "det_not_deg"),
          sum(conf_sw$switch_type == "det_not_deg")
        )
      )
    }
  }
  meta <- self@metadata
  overdispersion <- NA_real_
  if (!is.null(meta$inf_reps)) {
    od <- compute_overdispersion(draws_view(self))
    overdispersion <- stats::median(od$overdispersion[od$df > 0L])
  }
  replicates <- tibble::tibble(
    replicate_type = as.character(
      meta$uncertainty_replicate_type %||% NA_character_
    ),
    n_draws = as.integer(meta$uncertainty_n_draws %||% NA_integer_),
    resampling_df = as.integer(meta$uncertainty_resampling_df %||% NA_integer_),
    median_overdispersion = overdispersion
  )
  list(
    settings = settings,
    replicates = replicates,
    run = run,
    transcripts = by_type,
    genes = by_gene,
    de = de,
    switches = switches
  )
}

#' @rdname uncertainty
#' @export
uncertainty_methods <- S7::new_generic(
  "uncertainty_methods",
  "self",
  function(self, ...) {
    S7::S7_dispatch()
  }
)

S7::method(uncertainty_methods, IsoformicExperiment) <- function(self, ...) {
  meta <- self@metadata
  settings <- uncertainty_settings(self)
  n_draws <- meta$uncertainty_n_draws %||% NA
  type <- meta$uncertainty_replicate_type %||% NA_character_
  type <- if (identical(type, "gibbs")) {
    "Gibbs samples"
  } else if (identical(type, "bootstrap")) {
    "bootstrap replicates"
  } else {
    "inferential replicates"
  }
  tool_label <- switch(
    self@metadata$uncertainty_tool %||% "salmon",
    oarfish = "Oarfish",
    "Salmon"
  )
  versions <- NULL
  if ("salmon_version" %in% colnames(self@col_data)) {
    versions <- unique(stats::na.omit(self@col_data$salmon_version))
  }
  sentences <- c(
    paste0(
      "Transcript abundances were quantified with ",
      tool_label,
      if (length(versions)) {
        paste0(" ", paste(versions, collapse = "/"))
      } else {
        ""
      },
      " and ",
      n_draws,
      " ",
      type,
      " per sample were used to measure the inferential uncertainty of each isoform",
      if (!is.null(meta$uncertainty_resampling_df)) {
        paste0(
          " (",
          meta$uncertainty_resampling_df,
          " degrees of freedom over the study)"
        )
      } else {
        ""
      },
      "."
    ),
    paste0(
      "For every transcript, the inferential relative variance (InfRV; ",
      "Zhu et al. 2019) was averaged over samples, and transcripts with an ",
      "InfRV of at least ",
      settings$inf_rv_cutoff,
      " were called uncertain."
    ),
    if (identical(meta$uncertainty_ambiguity_source, "treeterminus")) {
      paste0(
        "Isoforms were placed in the transcript trees of TreeTerminus, and every ",
        "isoform with an InfRV of at least ",
        meta$uncertainty_tree_cutoff %||% settings$inf_rv_cutoff,
        " was replaced by its lowest ancestor below that value (nested choices ",
        "were merged into the larger node); an isoform in such a group was ",
        "called ambiguous."
      )
    } else if (identical(meta$uncertainty_ambiguity_source, "terminus")) {
      paste0(
        "Isoforms that the reads cannot tell apart were grouped with Terminus, ",
        "and an isoform in a group was called ambiguous."
      )
    } else if (identical(meta$uncertainty_ambiguity_source, "external")) {
      paste0(
        "Isoforms that the reads cannot tell apart were grouped with an external ",
        "grouping, and an isoform in a group was called ambiguous."
      )
    } else {
      paste0(
        "The correlation between the replicates of isoforms of the same gene ",
        "was averaged over samples, and isoform pairs with a correlation of ",
        signif(settings$cor_cutoff, 3),
        " or lower were called ambiguous and grouped (as in mmcollapse, Turro ",
        "et al. 2014, ",
        if (identical(meta$uncertainty_cor_cutoff_rule, "adaptive")) {
          paste0(
            "with the cutoff set to minus the 97.5th percentile of the largest ",
            "correlation of each isoform)."
          )
        } else {
          "with a fixed cutoff)."
        }
      )
    },
    paste0(
      "The log2 fold-change of each isoform was recomputed on every replicate ",
      "(counts scaled to a common library size, pseudocount 5), and a ",
      "fold-change was called stable when at least ",
      format(100 * settings$min_agreement),
      "% of the replicates agreed on its direction."
    ),
    paste0(
      "Differentially expressed isoforms were assigned high confidence when ",
      "the fold-change was stable, the quantification was not uncertain",
      if ("identifiable" %in% colnames(self@dea[["uncertainty"]])) {
        paste0(
          ", no sibling isoform was ambiguous and the equivalence classes of ",
          "every sample identified the isoform (its abundance is fixed by ",
          "the class counts; Zheng et al. 2022)"
        )
      } else {
        " and no sibling isoform was ambiguous"
      },
      ", moderate confidence when one of ",
      "these conditions failed and low confidence when two or more failed."
    ),
    paste0(
      "Only transcripts with a mean estimated count of at least ",
      settings$min_count,
      " were assessed."
    )
  )
  text <- paste(sentences, collapse = " ")
  if (!is.null(meta$uncertainty_group)) {
    lv <- strsplit(meta$uncertainty_levels, "|", fixed = TRUE)[[1L]]
    text <- paste0(
      text,
      " The contrast was ",
      lv[2L],
      " versus ",
      lv[1L],
      " in the ",
      meta$uncertainty_group,
      " variable."
    )
  }
  text
}

#' Put the quantification uncertainty into the counts
#'
#' @description
#' Reads that map to several isoforms make the counts of an isoform more
#' variable than a Poisson count, by an amount that differs between
#' transcripts. Baldoni et al. (2024) estimate this overdispersion for every
#' transcript from the Salmon replicates and divide it out of the counts, so
#' that tests built for counts (`edgeR`, `limma`, `DESeq2`) keep their error
#' rate without a model of their own for the read ambiguity.
#'
#' * `inferential_overdispersion()` computes the overdispersion of every
#'   transcript from the replicates stored in the object, as
#'   `edgeR::catchSalmon()` does: the Pearson statistic of a quasi-Poisson fit
#'   to the replicates of each sample, pooled over the samples, moderated
#'   towards the median with a prior of `prior_df` degrees of freedom and
#'   floored at one. The estimates use the transcripts whose replicates are
#'   in the object (those with a mean count of at least `min_count` in
#'   [add_uncertainty()]), so the prior differs from the one `catchSalmon()`
#'   computes over every transcript of the index; transcripts without
#'   replicates get the prior.
#' * `scale_counts()` divides an assay by the overdispersion and stores the
#'   result as a new assay, ready for a count-based test.
#'
#' The overdispersion complements the confidence classes of
#' [de_confidence()]: it changes the test, while the classes describe how far
#' each call can be trusted.
#'
#' @param self An [IsoformicExperiment] with inferential replicates (see
#'   [add_uncertainty()]).
#' @param prior_df Degrees of freedom of the prior (3 in `edgeR`).
#' @param assay Name of the counts assay to scale.
#' @param name Name of the new assay.
#' @param ... Unused.
#'
#' @returns `inferential_overdispersion()`: a tibble with `transcript_id`,
#'   `overdispersion` and `df` (degrees of freedom of the estimate), with the
#'   prior in the `prior` attribute. `scale_counts()`: the object with the
#'   scaled assay added.
#'
#' @references
#' Baldoni PL, Chen Y, Hediyeh-Zadeh S, Liao Y, Dong X, Ritchie ME, Shi W,
#' Smyth GK (2024). Dividing out quantification uncertainty allows efficient
#' assessment of differential transcript expression with edgeR. *Nucleic Acids
#' Research* 52(3):e13.
#'
#' @seealso [add_uncertainty()], [de_confidence()]
#'
#' @examples
#' extdata <- system.file("extdata", package = "isoformic")
#' quant_dirs <- list.files(file.path(extdata, "mini_quants"), full.names = TRUE)
#' names(quant_dirs) <- sub("_quant$", "", basename(quant_dirs))
#' iso <- read_salmon(
#'   quant_dirs,
#'   annot_path = file.path(extdata, "mini_annotation"),
#'   col_data = data.frame(
#'     sample_id = names(quant_dirs),
#'     condition = rep(c("control", "treatment"), each = 3)
#'   )
#' )
#' iso@metadata$inf_reps <- as_inferential_replicates(read_parquet(
#'   file.path(extdata, "mini_draws.parquet")
#' ))
#' iso <- scale_counts(iso)
#' assay_names(iso)
#' head(inferential_overdispersion(iso))
#'
#' @name inferential-overdispersion
#' @aliases inferential_overdispersion scale_counts
NULL

#' @rdname inferential-overdispersion
#' @export
inferential_overdispersion <- S7::new_generic(
  "inferential_overdispersion",
  "self",
  function(self, prior_df = 3, ...) {
    S7::S7_dispatch()
  }
)

S7::method(inferential_overdispersion, IsoformicExperiment) <- function(
  self,
  prior_df = 3,
  ...
) {
  compute_overdispersion(draws_view(self), prior_df = prior_df)
}

#' Overdispersion of every transcript of an array of draws
#' @keywords internal
#' @noRd
compute_overdispersion <- function(arr, prior_df = 3) {
  moments <- draw_moments(arr)
  n_draw <- moments$n
  if (n_draw < 2L) {
    cli::cli_abort(
      message = c(x = "At least two inferential replicates are needed."),
      class = "isoformic_no_inferential_replicates"
    )
  }
  expressed <- moments$mean > 0
  pearson <- ifelse(expressed, moments$var * (n_draw - 1L) / moments$mean, 0)
  df <- rowSums(expressed) * (n_draw - 1L)
  phi <- rowSums(pearson) / df
  has_df <- df > 0L
  prior <- NA_real_
  if (any(has_df)) {
    prior <- max(
      1,
      stats::median(phi[has_df]) /
        stats::qf(0.5, df1 = stats::median(df[has_df]), df2 = prior_df)
    )
    phi[has_df] <- (prior_df * prior + df[has_df] * phi[has_df]) /
      (prior_df + df[has_df])
    phi <- pmax(phi, 1)
    phi[!has_df] <- prior
  }
  res <- tibble::tibble(
    transcript_id = rownames(moments$mean),
    overdispersion = unname(phi),
    df = unname(df)
  )
  attr(res, "prior") <- prior
  res
}

#' @rdname inferential-overdispersion
#' @export
scale_counts <- S7::new_generic(
  "scale_counts",
  "self",
  function(self, assay = "counts", name = paste0(assay, "_scaled"), ...) {
    S7::S7_dispatch()
  }
)

S7::method(scale_counts, IsoformicExperiment) <- function(
  self,
  assay = "counts",
  name = paste0(assay, "_scaled"),
  ...
) {
  check_string(assay)
  check_string(name)
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
  assays[[name]] <- counts / overdispersion_for(self, rownames(counts))
  self@assays <- assays
  self
}

#' Overdispersion of the given transcripts, the prior for those without draws
#' @keywords internal
#' @noRd
overdispersion_for <- function(self, transcript_ids, samples = NULL) {
  arr <- draws_view(self)
  if (!is.null(samples)) {
    arr <- draws_select_samples(arr, samples)
  }
  od <- compute_overdispersion(arr)
  factor <- od$overdispersion[match(transcript_ids, od$transcript_id)]
  factor[is.na(factor)] <- attr(od, "prior") %||% 1
  factor
}
