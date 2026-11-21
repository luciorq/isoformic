# Productive versus unproductive transcription -------------------------------

#' Productivity of Transcripts and Genes
#'
#' @description
#' The biotype split used by the workflow (protein coding versus retained
#' intron, nonsense-mediated decay, ...) is a proxy. These functions make
#' productivity a computed, graded property of every transcript and
#' summarise it per gene:
#'
#' * `transcript_productivity()` classifies each assayed transcript as
#'   `"productive"` (complete ORF, no predicted NMD, same ORF as the
#'   canonical isoform), `"altered_productive"` (complete ORF, no NMD, but
#'   a different ORF than the canonical isoform: alternative start, in-frame
#'   change, frameshift or stop in novel sequence), `"unproductive"` (no
#'   ORF, premature stop codon by the 50-nt rule, no stop codon, an
#'   unproductive biotype such as `nonsense_mediated_decay`,
#'   `retained_intron`, `protein_coding_CDS_not_defined`,
#'   `protein_coding_LoF`, or a longest ORF without coding potential),
#'   `"noncoding"` (transcripts of genes without any coding isoform) or
#'   `"uncertain"` (incomplete CDS annotation without
#'   a premature stop, a credible ORF unrelated to the canonical one, FANTOM
#'   `uncertain_coding`).
#'   The `reason` column explains each call and `source` tells whether it
#'   came from the sequence (`"sequence"`, when the annotation has
#'   transcript sequences, see [extract_transcript_sequences()] and
#'   [transcript_orfs()]), the annotation (`"annotation"`) or external
#'   evidence (`"evidence"`). Evidence columns take precedence, in this
#'   order: a `productivity_evidence` column added with
#'   [add_annotation_columns()], a validated NMD exon ([add_nmd_exons()]),
#'   a Ribo-seq ORF matching the transcript's ORF ([add_riboseq_orfs()]),
#'   an APPRIS principal isoform ([add_appris()]).
#' * `productive_fraction()` is the share of a gene's expression carried by
#'   productive (and altered productive) isoforms, per sample or per group,
#'   with the difference between the last and the first group level. A
#'   gene whose expression is flat while its productive fraction drops is a
#'   post-transcriptional regulation call that gene-level differential
#'   expression misses. With inferential replicates (`inf_reps`, Salmon
#'   bootstraps or Gibbs samples) it reports the uncertainty of each
#'   fraction, and with a `reference` ([gtex_reference()]) the difference to
#'   the healthy-tissue baseline.
#' * `productivity_fractions()` splits the gene's expression over all
#'   productivity classes and, for unproductive isoforms, over the
#'   mechanism (`nmd`, `retained_intron`, `no_orf`, `non_stop`, `other`).
#' * `test_productive_fraction()` tests the change in productive fraction
#'   between two groups on counts, with a quasi-binomial model per gene
#'   (pooled proportions per group, Pearson over-dispersion, t test on the
#'   log odds ratio, BH adjustment). A gene with fewer than 5 reads in its
#'   smaller class (productive or not), pooled over both groups, carries
#'   too little information on usage: it gets a `log_odds_ratio` of 0 and no
#'   statistic or p-value, and does not count in the BH adjustment. A gene
#'   without reads in one group gets `NA` there too.
#' * `run_productivity_enrichment()` ranks genes by the change in
#'   productive fraction between the last and the first group level
#'   (`rank_by = "delta"`), by that change over its standard deviation
#'   (`rank_by = "delta_z"`) or by the statistic of
#'   `test_productive_fraction()` (`rank_by = "statistic"`), and runs gene
#'   set enrichment (`fgsea`) on gene symbols, which captures regulated
#'   unproductive splicing (AS-NMD) directly. It takes an
#'   [IsoformicExperiment], a `productive_fraction()` table or a
#'   `test_productive_fraction()` table, and one or several gene set
#'   collections.
#' * `poison_exons()` lists exons that occur only in unproductive isoforms
#'   of genes that also have productive isoforms and that do not overlap
#'   any exon of a productive isoform (cassette poison / NMD exons, as
#'   opposed to retained introns or alternative splice sites);
#'   [plot_genomic_context()] outlines them with `highlight_poison_exons`.
#' * `benchmark_productivity()` measures the rules against the annotation
#'   biotypes: recall of `nonsense_mediated_decay` and `retained_intron`
#'   transcripts by the sequence/structure rules alone, and the share of
#'   `protein_coding` and canonical transcripts called productive.
#' * `productivity_palette()` gives colours for the classes.
#'
#' @details
#' **Low-count genes.** A productive fraction computed from a few reads is
#' mostly sampling noise: a gene with 4 reads moves by 0.25 when one read
#' changes isoform. `delta_z` does not catch this, because the inferential
#' replicates measure how the reads are split between isoforms, not how many
#' reads were sampled. `min_gene_count` removes genes whose mean count (the
#' sum of the counts of their isoforms in the `"counts"` assay) is below it in
#' any group, before the statistics and the BH adjustment; with
#' `rank_by = "statistic"`, before the ranking. The filter is off by default
#' (0); 10 reads is a reasonable start. The rule asks for enough reads in
#' every group, because a fraction is compared between groups: it is
#' stricter than the transcript filter of [test_usage()] (`min_count` reads
#' in as many samples as the smallest group), and `min_expr` sets single
#' isoform values to zero rather than removing genes.
#' `rank_by = "statistic"` accounts for the reads of every gene and for the
#' variation between its samples. Its dispersion is estimated for each gene
#' alone, on `n - 2` degrees of freedom and not moderated, so with three
#' samples per group the ranking is still noisy. The t statistic is put on
#' the normal scale through its p-value
#' (`sign(statistic) * qnorm(pvalue / 2, lower.tail = FALSE)`), because the
#' degrees of freedom differ between genes when samples without reads drop
#' out. Its sign is that of the pooled change (`delta_productive_fraction`
#' of `test_productive_fraction()`), which can differ from the mean of the
#' per-sample fractions in `productive_fraction()` when the samples differ
#' much in depth. Genes with fewer than 5 reads in one class get no
#' statistic and are left out of the ranking, so that a change in the
#' expression of the whole gene does not rank it.
#'
#' **Uncertainty of the fractions.** With inferential replicates,
#' `productive_fraction()` computes the productive fraction of every gene and
#' sample in each draw. The draws are counts; for another `assay` each draw
#' is first put on the scale of that assay, multiplied by `assay / counts` of
#' the same transcript and sample (for TPM, the inverse of the effective
#' length times the scaling of the sample), so that the fractions of the
#' draws and of the assay, and `min_expr`, are on one scale. Where the point
#' count or value is zero (Gibbs draws often give reads to such cells) the
#' factor is estimated: from the `length` assay when the assay follows it
#' (TPM of [read_salmon()] and [read_oarfish()]), otherwise from a transcript
#' effect plus a sample effect on the log scale fitted on the other cells,
#' and for a transcript without any count from its annotated length.
#' `productive_fraction_sd` is the sample standard deviation of a gene's
#' fraction over the draws, around the mean of the draws, averaged over the
#' samples kept; `delta_productive_fraction_sd` is the standard deviation
#' over the draws of the change between the group means, and `delta_z` is
#' `delta_productive_fraction` divided by it. In every draw, only the samples
#' in which the assay gives the gene a fraction are used. Both measure the
#' spread of the draws, not their distance to the point estimate: the mean of
#' Gibbs draws can differ from the EM estimate of the assay.
#'
#' @param self An [IsoformicExperiment]. `run_productivity_enrichment()` also
#'   takes the table of `productive_fraction()` for two or more groups, or of
#'   `test_productive_fraction()` for `rank_by = "statistic"`, so that one
#'   table serves several gene set collections; the arguments that build the
#'   table (`group`, `levels`, `assay`, `include_altered`, `min_expr`,
#'   `min_gene_count`, `inf_reps`) are then ignored.
#' @param gene_name Optional genes (names or IDs) to restrict the result to.
#' @param nmd_rule Distance (nt) from the stop codon to the last junction
#'   above which NMD is predicted; default 50.
#' @param use_biotype Logical, also use unproductive biotypes and FANTOM
#'   coding status (default `TRUE`).
#' @param use_sequence Logical, use the sequence-derived ORFs of
#'   [transcript_orfs()] when the annotation has sequences (default `TRUE`).
#' @param use_evidence Logical, apply the evidence columns (default `TRUE`).
#' @param coding_cutoff Coding probability below which a longest ORF (one
#'   that could not be projected from the canonical start) is not believed
#'   (default 0.364).
#' @param include_altered Logical, count `altered_productive` isoforms as
#'   productive in `productive_fraction()` (default `TRUE`).
#' @param assay Transcript assay to use (default `"abundance"`;
#'   `"counts"` for the test).
#' @param group Name of a `col_data` column defining groups; `NULL` for
#'   per-sample values. `run_productivity_enrichment()` requires it for an
#'   [IsoformicExperiment] and ignores it for a `productive_fraction()` table.
#' @param levels Levels of `group` to use, in this order; the change is
#'   `last - first`. Samples of other groups are ignored, as if the object
#'   had been subset to these samples, so every returned value comes from
#'   the samples kept. `NULL` takes all the levels of `group`. The tests
#'   (`test_productive_fraction()`, `test_apa()`) need exactly two levels.
#' @param min_expr Minimum assay value for an isoform to contribute
#'   (default 0): isoforms below it count as zero.
#' @param min_gene_count Minimum mean count of a gene, the sum of its
#'   isoforms in the `"counts"` assay, in every group (over all the samples
#'   kept without `group`); genes below it are dropped. Needs a `"counts"`
#'   assay whatever `assay` is. The default 0 keeps every gene. See Details.
#' @param inf_reps Optional list of transcript matrices (inferential
#'   replicates with the transcripts in rows and the samples of the assay in
#'   columns, as counts, as [read_inferential_replicates()] returns them);
#'   `NULL` uses `x@metadata$inf_reps` when present and `FALSE` skips the
#'   replicates, which leaves out the standard deviations and `delta_z` and
#'   is much faster. For an `assay` other than `"counts"` the draws are put on
#'   its scale with the `"counts"` assay (see Details); an object without a
#'   `"counts"` assay takes the draws to be on the scale of `assay`.
#'   `run_productivity_enrichment()` reads the replicates only for
#'   `rank_by = "delta_z"`.
#' @param rank_by Statistic that ranks the genes in
#'   `run_productivity_enrichment()`: `"delta"` (default), the change in
#'   productive fraction, `"delta_z"`, that change divided by its
#'   standard deviation over the inferential replicates, which moves the
#'   genes with an uncertain change towards the middle of the ranking, or
#'   `"statistic"`, the test of `test_productive_fraction()` on the
#'   `"counts"` assay between the first and the last level of `levels`, as
#'   a signed normal quantile of its p-value (see Details); `assay` and
#'   `inf_reps` are then not used.
#'   `"delta_z"` needs inferential replicates; a change that is the same in
#'   every replicate (`delta_z` is `NA` in the `productive_fraction()`
#'   table) is ranked with the largest absolute `delta_z` of the table (at
#'   least 1), with the sign of the change. It reflects the
#'   quantification uncertainty only, not the variation between samples or
#'   the number of reads; `"statistic"` reflects both (see Details).
#' @param reference Optional tibble from [gtex_reference()], or a GTEx
#'   tissue name, giving the healthy-tissue productive fraction per gene.
#' @param as_matrix Logical, return the genes x samples matrix of
#'   productive fractions instead of the summary table.
#' @param genesets_list Named list of gene symbol vectors (see [read_gmt()]),
#'   or a named list of such lists: each collection is then tested on its own
#'   and the result has a `collection` column (see [run_enrichment()]).
#' @param min_size,max_size Gene set size limits passed to `fgsea`.
#' @param ... Used only by `run_productivity_enrichment()`, which passes them
#'   to `fgsea::fgseaMultilevel()`; the other functions do not use them
#'   (reserved for methods). Each `fgsea` call runs with
#'   `BiocParallel::SerialParam()` unless `nproc` or `BPPARAM` is given;
#'   prefer `workers` (see *Workers* in [run_enrichment()]).
#' @param workers Number of worker processes that run the `fgsea` calls of
#'   `run_productivity_enrichment()` (one per gene set collection) in
#'   parallel: a whole number, or `"auto"` for the CPUs available to the
#'   session. Defaults to the `isoformic.workers` option, or 1 (serial).
#'   Needs the `mirai` package above 1; see *Workers* in [run_enrichment()],
#'   also for the random seeds.
#'
#' @returns
#' * `transcript_productivity()`: one row per transcript with
#'   `transcript_id`, `transcript_name`, `transcript_type`, `gene_id`,
#'   `gene_name`, `productivity` (a factor with the levels of
#'   `productivity_palette()`), `reason`, `source`, `has_cds`,
#'   `cds_complete`, `orf_length`, `nmd_predicted`, `stop_to_last_junction`,
#'   `orf_fraction_of_canonical`, `is_canonical` and, with sequences, the
#'   ORF columns behind the call (`orf_source`, `nmd_class`,
#'   `frame_vs_canonical`, `coding_probability`, `n_uorf`).
#' * `productive_fraction()`: one row per coding gene with
#'   `gene_id`, `gene_name`, `mean_gene_expr`, `productive_fraction` (mean
#'   over the samples kept), with `group` one `productive_fraction_<level>`
#'   column per level and `delta_productive_fraction`, with inferential
#'   replicates `productive_fraction_sd` (the standard deviation over the
#'   draws, see Details) and, with `group`,
#'   `delta_productive_fraction_sd` and `delta_z`, and with `reference` the
#'   `reference_*` columns and `delta_vs_reference`. A genes x samples
#'   matrix with `as_matrix = TRUE`.
#' * `productivity_fractions()`: one row per gene and class present
#'   (`productivity`, a factor), with `fraction` and, with `group`, the
#'   `fraction_<level>` and `delta_fraction` columns.
#' * `test_productive_fraction()`: one row per coding gene with
#'   `mean_gene_expr`, the two `productive_fraction_<level>` columns,
#'   `delta_productive_fraction`, `log_odds_ratio`, `se`, `statistic`,
#'   `dispersion`, `df`, `pvalue` and `padj`, ordered by `pvalue`.
#' * `run_productivity_enrichment()`: the `fgsea` result with the columns
#'   `experiment` (`"productive_fraction"`), `rank_by` and, for a list of
#'   collections, `collection`.
#' * `poison_exons()`: one row per exon with `gene_id`, `gene_name`,
#'   `seqid`, `start_pos`, `end_pos`, `strand`, `n_unproductive`,
#'   `transcript_ids` and `exon_length`.
#' * `benchmark_productivity()`: a list with the `summary` (`metric`, `n`,
#'   `value`, `source`) and `confusion` (`transcript_type`, `productivity`,
#'   `n`) tibbles.
#' * `productivity_palette()`: a named character vector of colours.
#'
#' @examples
#' iso <- example_isoformic()
#' prod <- transcript_productivity(iso, gene_name = "RBPJ")
#' prod[, c("transcript_name", "transcript_type", "productivity", "reason", "source")]
#' head(productive_fraction(iso, group = "condition"))
#' head(productivity_fractions(iso, group = "condition"))
#' head(test_productive_fraction(iso, group = "condition"))
#' poison_exons(iso, gene_name = "RBPJ")
#' benchmark_productivity(iso)$summary
#'
#' @name productivity
#' @rdname productivity
NULL

#' @rdname productivity
#' @export
productivity_palette <- function() {
  c(
    productive = "#4daf4a",
    altered_productive = "#ff7f00",
    unproductive = "#e41a1c",
    noncoding = "#999999",
    uncertain = "#cccccc"
  )
}

unproductive_biotypes <- function() {
  c(
    "nonsense_mediated_decay",
    "non_stop_decay",
    "retained_intron",
    "protein_coding_CDS_not_defined",
    "processed_transcript",
    # Not translated in the reference genome because of a variant
    "protein_coding_LoF"
  )
}

#' @rdname productivity
#' @export
transcript_productivity <- S7::new_generic(
  "transcript_productivity",
  "self",
  function(
    self,
    gene_name = NULL,
    nmd_rule = 50,
    use_biotype = TRUE,
    use_sequence = TRUE,
    use_evidence = TRUE,
    coding_cutoff = 0.364,
    ...
  ) {
    S7::S7_dispatch()
  }
)

S7::method(transcript_productivity, IsoformicExperiment) <- function(
  self,
  gene_name = NULL,
  nmd_rule = 50,
  use_biotype = TRUE,
  use_sequence = TRUE,
  use_evidence = TRUE,
  coding_cutoff = 0.364,
  ...
) {
  with_sequence <- isTRUE(use_sequence) && has_sequences(self)
  signature <- paste(
    annot_full_signature(self),
    if (with_sequence) {
      annot_signature(self@annot_path, "sequences")
    } else {
      "noseq"
    },
    rlang::hash(row_names(self)),
    nmd_rule,
    use_biotype,
    use_evidence,
    coding_cutoff
  )
  res <- object_cache_fetch(
    self,
    "transcript_productivity",
    signature,
    function() {
      classify_productivity(
        feats = transcript_features(self, nmd_rule = nmd_rule),
        annot = row_data(self),
        genes = self@annot_data_genes,
        orfs = if (with_sequence) {
          transcript_orfs(self, nmd_rule = nmd_rule)
        } else {
          NULL
        },
        use_biotype = use_biotype,
        use_evidence = use_evidence,
        coding_cutoff = coding_cutoff
      )
    }
  )
  filter_genes(res, gene_name)
}

#' Classify transcripts from their features, ORFs, biotypes and evidence
#' @keywords internal
#' @noRd
classify_productivity <- function(
  feats,
  annot,
  genes,
  orfs = NULL,
  use_biotype = TRUE,
  use_evidence = TRUE,
  coding_cutoff = 0.364
) {
  n <- nrow(feats)
  idx <- match(feats$transcript_id, annot$transcript_id)
  biotype <- feats$transcript_type
  annot_col <- function(col) {
    if (col %in% colnames(annot)) annot[[col]][idx] else rep(NA, n)
  }
  coding_status <- as.character(annot_col("coding_status"))
  gene_type <- genes$gene_type[match(feats$gene_id, genes$gene_id)]
  calls <- list(
    class = rep("productive", n),
    reason = rep("complete ORF identical to the canonical isoform", n),
    source = rep(if (is.null(orfs)) "annotation" else "sequence", n)
  )
  set <- function(calls, mask, cls, why, src = NULL) {
    mask <- mask & !is.na(mask)
    calls$class[mask] <- cls
    calls$reason[mask] <- if (length(why) == 1L) why else why[mask]
    if (!is.null(src)) {
      calls$source[mask] <- src
    }
    calls
  }

  if (is.null(orfs)) {
    gene_coding <- stats::ave(feats$has_cds, feats$gene_id, FUN = any) |
      gene_type %in% "protein_coding" |
      coding_status %in% "coding"
    calls <- set(
      calls,
      !gene_coding,
      "noncoding",
      "gene without coding isoforms"
    )
    calls <- set(
      calls,
      gene_coding & !feats$has_cds,
      "unproductive",
      "no ORF annotated in a coding gene"
    )
    calls <- set(
      calls,
      gene_coding & feats$has_cds & !feats$cds_complete,
      "uncertain",
      "incomplete CDS annotation"
    )
    calls <- set(
      calls,
      gene_coding & feats$has_cds & feats$cds_complete & feats$nmd_predicted,
      "unproductive",
      "stop codon more than the NMD rule distance upstream of the last junction"
    )
    altered <- gene_coding &
      feats$has_cds &
      feats$cds_complete &
      !feats$nmd_predicted &
      feats$cds_differs_from_canonical %in% TRUE
    calls <- set(
      calls,
      altered,
      "altered_productive",
      "complete ORF differing from the canonical isoform"
    )
    nmd_predicted <- feats$nmd_predicted
    extra <- NULL
  } else {
    o <- orfs[match(feats$transcript_id, orfs$transcript_id), ]
    credible <- o$orf_source %in%
      c("annotation", "projected", "annotation_partial") |
      (o$orf_source %in% "longest" & o$coding_probability > coding_cutoff)
    gene_coding <- stats::ave(
      o$orf_source %in% c("annotation", "projected"),
      feats$gene_id,
      FUN = any
    ) |
      gene_type %in% "protein_coding" |
      coding_status %in% "coding"
    no_seq <- !o$has_sequence
    calls <- set(
      calls,
      !gene_coding,
      "noncoding",
      "gene without coding isoforms"
    )
    calls <- set(
      calls,
      gene_coding & no_seq & !feats$has_cds,
      "unproductive",
      "no ORF annotated in a coding gene (no sequence)",
      "annotation"
    )
    calls <- set(
      calls,
      gene_coding & no_seq & feats$has_cds & !feats$cds_complete,
      "uncertain",
      "incomplete CDS annotation (no sequence)",
      "annotation"
    )
    calls <- set(
      calls,
      gene_coding &
        no_seq &
        feats$has_cds &
        feats$cds_complete &
        feats$nmd_predicted,
      "unproductive",
      "stop codon more than the NMD rule distance upstream of the last junction (no sequence)",
      "annotation"
    )
    calls <- set(
      calls,
      gene_coding & !no_seq & o$orf_source %in% "none",
      "unproductive",
      "no ORF found in a coding gene"
    )
    calls <- set(
      calls,
      gene_coding & !no_seq & o$orf_source %in% "longest" & !credible,
      "unproductive",
      paste0(
        "canonical start codon not projectable; longest ORF without coding potential (probability ",
        formatC(o$coding_probability, digits = 2, format = "f"),
        ")"
      )
    )
    can_len <- o$orf_length[match(o$canonical_transcript_id, o$transcript_id)]
    start_same <- o$start_same_as_canonical %in% TRUE
    stop_same <- o$stop_same_as_canonical %in% TRUE
    same_length <- !is.na(can_len) & o$orf_length == can_len
    altered_reason <- rep("", n)
    altered_reason[
      o$orf_source %in% "longest"
    ] <- "credible ORF not projectable from the canonical start"
    altered_reason[
      !start_same & !(o$orf_source %in% "longest")
    ] <- "alternative start codon"
    altered_reason[
      start_same & stop_same & !same_length
    ] <- "in-frame internal change of the ORF"
    altered_reason[
      start_same & !stop_same & o$frame_vs_canonical %in% "in_frame"
    ] <-
      "in-frame change of the ORF end"
    altered_reason[start_same & o$frame_vs_canonical %in% "frameshift"] <-
      "frameshift relative to the canonical ORF"
    altered_reason[start_same & o$frame_vs_canonical %in% "novel_sequence"] <-
      "stop codon in sequence absent from the canonical isoform"
    altered <- gene_coding &
      credible &
      o$nmd_class %in% c("none", "long_3utr") &
      !(start_same & stop_same & same_length)
    calls <- set(calls, altered, "altered_productive", altered_reason)
    long_utr <- gene_coding & credible & o$nmd_class %in% "long_3utr"
    calls$reason[long_utr] <- paste0(calls$reason[long_utr], "; long 3' UTR")
    # A credible but unprojectable ORF is only an altered protein when it ends
    # at the canonical stop codon (N-terminal variant); otherwise its relation
    # to the gene's protein is unknown
    calls <- set(
      calls,
      altered & o$orf_source %in% "longest" & !stop_same,
      "uncertain",
      "credible ORF not projectable from the canonical start and not sharing its stop codon"
    )
    calls <- set(
      calls,
      gene_coding &
        o$orf_source %in% "annotation_partial" &
        !(o$nmd_class %in% c("ptc", "ptc_start_proximal", "non_stop")),
      "uncertain",
      "incomplete CDS annotation"
    )
    calls <- set(
      calls,
      gene_coding & credible & is.na(o$nmd_class) & o$end_incomplete %in% TRUE,
      "uncertain",
      "incomplete 3' end, no stop codon reachable"
    )
    calls <- set(
      calls,
      gene_coding & credible & o$nmd_class %in% "ptc_start_proximal",
      "unproductive",
      "premature stop codon close to the start codon (NMD escape by re-initiation possible)"
    )
    calls <- set(
      calls,
      gene_coding & credible & o$nmd_class %in% "non_stop",
      "unproductive",
      "no in-frame stop codon (non-stop decay)"
    )
    calls <- set(
      calls,
      gene_coding & credible & o$nmd_class %in% "ptc",
      "unproductive",
      "premature stop codon more than the NMD rule distance upstream of the last junction"
    )
    nmd_predicted <- o$nmd_class %in% c("ptc", "ptc_start_proximal")
    extra <- tibble::tibble(
      orf_source = o$orf_source,
      nmd_class = o$nmd_class,
      frame_vs_canonical = o$frame_vs_canonical,
      coding_probability = o$coding_probability,
      n_uorf = o$n_uorf
    )
  }

  if (isTRUE(use_biotype)) {
    calls <- set(
      calls,
      biotype %in% unproductive_biotypes() & calls$class != "noncoding",
      "unproductive",
      "unproductive biotype",
      "annotation"
    )
    calls <- set(
      calls,
      coding_status %in% "uncertain" | biotype %in% "uncertain_coding",
      "uncertain",
      "FANTOM uncertain coding status",
      "annotation"
    )
    calls <- set(
      calls,
      coding_status %in%
        "nonCoding" &
        gene_coding &
        calls$class == "productive",
      "unproductive",
      "FANTOM non-coding status in a coding gene",
      "annotation"
    )
  }
  if (isTRUE(use_evidence)) {
    calls <- set(
      calls,
      annot_col("appris_principal") %in% TRUE & calls$class != "noncoding",
      "productive",
      "APPRIS principal isoform",
      "evidence"
    )
    calls <- set(
      calls,
      annot_col("riboseq_main_orf") %in% TRUE & calls$class != "noncoding",
      "productive",
      "ORF translated in Ribo-seq (GENCODE Ribo-seq ORFs)",
      "evidence"
    )
    calls <- set(
      calls,
      annot_col("nmd_validated_exon") %in% TRUE,
      "unproductive",
      "contains an NMD-inducing exon validated by NMD inhibition",
      "evidence"
    )
    ev <- tolower(as.character(annot_col("productivity_evidence")))
    ok <- ev %in% names(productivity_palette())
    if (any(ok)) {
      calls <- set(calls, ok, ev[ok], "external evidence", "evidence")
    }
  }
  res <- tibble::tibble(
    transcript_id = feats$transcript_id,
    transcript_name = feats$transcript_name,
    transcript_type = biotype,
    gene_id = feats$gene_id,
    gene_name = feats$gene_name,
    productivity = factor(calls$class, levels = names(productivity_palette())),
    reason = calls$reason,
    source = calls$source,
    has_cds = feats$has_cds,
    cds_complete = feats$cds_complete,
    orf_length = if (is.null(orfs)) feats$orf_length else o$orf_length,
    nmd_predicted = nmd_predicted,
    stop_to_last_junction = if (is.null(orfs)) {
      feats$stop_to_last_junction
    } else {
      o$stop_to_last_junction
    },
    orf_fraction_of_canonical = feats$orf_fraction_of_canonical,
    is_canonical = feats$is_canonical
  )
  if (!is.null(extra)) {
    res <- dplyr::bind_cols(res, extra)
  }
  res
}

#' Mechanism of an unproductive call, from its reason
#' @keywords internal
#' @noRd
productivity_detail <- function(prod) {
  class <- as.character(prod$productivity)
  reason <- prod$reason
  biotype <- prod$transcript_type
  detail <- class
  unp <- class == "unproductive"
  detail[unp] <- "other"
  detail[unp & grepl("premature stop|NMD rule|NMD-inducing", reason)] <- "nmd"
  detail[
    unp &
      grepl("unproductive biotype", reason, fixed = TRUE) &
      biotype %in% "nonsense_mediated_decay"
  ] <- "nmd"
  detail[
    unp &
      grepl("unproductive biotype", reason, fixed = TRUE) &
      biotype %in% "retained_intron"
  ] <- "retained_intron"
  detail[
    unp & grepl("no ORF|without coding potential|CDS_not_defined", reason)
  ] <- "no_orf"
  detail[
    unp &
      grepl("unproductive biotype", reason, fixed = TRUE) &
      biotype %in% "protein_coding_CDS_not_defined"
  ] <- "no_orf"
  detail[unp & grepl("non-stop", reason, fixed = TRUE)] <- "non_stop"
  detail
}

#' Expression matrix with the productivity flags of its rows
#' @keywords internal
#' @noRd
productivity_matrix <- function(self, assay, min_expr, include_altered) {
  assay <- resolve_assay_name(self, assay)
  prod <- transcript_productivity(self)
  tx_mat <- tx_assay(self, assay)
  if (isTRUE(min_expr > 0)) {
    tx_mat[tx_mat < min_expr] <- 0
  }
  idx <- match(rownames(tx_mat), prod$transcript_id)
  productive_classes <- c(
    "productive",
    if (isTRUE(include_altered)) "altered_productive"
  )
  list(
    tx_mat = tx_mat,
    prod = prod,
    is_prod = as.character(prod$productivity[idx]) %in% productive_classes,
    gene_ids = prod$gene_id[idx],
    coding_genes = unique(prod$gene_id[!prod$productivity %in% "noncoding"])
  )
}

#' Genes that pass `min_gene_count`, `NULL` when no filter applies
#'
#' The gene count is the sum of the counts of its transcripts in the
#' `"counts"` assay over `sample_ids`; its mean must reach `min_gene_count`
#' in every level of `groups` (over all the samples when `groups` is `NULL`).
#' @keywords internal
#' @noRd
count_filtered_genes <- function(
  self,
  prod,
  sample_ids,
  groups,
  min_gene_count
) {
  if (
    !is.numeric(min_gene_count) ||
      length(min_gene_count) != 1L ||
      is.na(min_gene_count) ||
      min_gene_count < 0
  ) {
    cli::cli_abort(
      message = c(
        x = "{.arg min_gene_count} must be a single number of at least 0."
      ),
      class = "isoformic_invalid_argument"
    )
  }
  if (min_gene_count == 0) {
    return(NULL)
  }
  if (!"counts" %in% assay_names(self)) {
    cli::cli_abort(
      message = c(
        x = "{.arg min_gene_count} needs a {.val counts} assay.",
        i = "Available assays: {.val {assay_names(self)}}. Use {.code min_gene_count = 0} to skip the filter."
      ),
      class = "isoformic_assay_not_found"
    )
  }
  counts <- tx_assay(self, "counts")[, sample_ids, drop = FALSE]
  if (!is.double(counts)) {
    # Integer sums overflow to NA
    storage.mode(counts) <- "double"
  }
  gene_ids <- prod$gene_id[match(rownames(counts), prod$transcript_id)]
  known <- !is.na(gene_ids)
  total <- rowsum(counts[known, , drop = FALSE], gene_ids[known])
  if (is.null(groups)) {
    groups <- factor(rep("all", ncol(total)))
  }
  pass <- rep(TRUE, nrow(total))
  for (lv in levels(groups)) {
    lv_mean <- rowMeans(total[, groups == lv, drop = FALSE])
    pass <- pass & !is.na(lv_mean) & lv_mean >= min_gene_count
  }
  rownames(total)[pass]
}

#' @rdname productivity
#' @export
productive_fraction <- S7::new_generic(
  "productive_fraction",
  "self",
  function(
    self,
    group = NULL,
    levels = NULL,
    assay = "abundance",
    include_altered = TRUE,
    as_matrix = FALSE,
    gene_name = NULL,
    min_expr = 0,
    min_gene_count = 0,
    inf_reps = NULL,
    reference = NULL,
    ...
  ) {
    S7::S7_dispatch()
  }
)

S7::method(productive_fraction, IsoformicExperiment) <- function(
  self,
  group = NULL,
  levels = NULL,
  assay = "abundance",
  include_altered = TRUE,
  as_matrix = FALSE,
  gene_name = NULL,
  min_expr = 0,
  min_gene_count = 0,
  inf_reps = NULL,
  reference = NULL,
  ...
) {
  pm <- productivity_matrix(self, assay, min_expr, include_altered)
  grp <- group_samples(self, group, colnames(pm$tx_mat), levels = levels)
  groups <- grp$groups
  # Genes with too few reads in some group are dropped from the result
  expressed <- count_filtered_genes(
    self,
    pm$prod,
    grp$samples,
    groups,
    min_gene_count
  )
  # Only the transcripts of the genes reported and the samples kept are
  # summed, in the assay and in every draw
  coding <- pm$gene_ids %in% pm$coding_genes
  if (!is.null(expressed)) {
    coding <- coding & pm$gene_ids %in% expressed
  }
  if (!is.null(gene_name) && !isTRUE(as_matrix)) {
    coding <- coding & pm$gene_ids %in% gene_ids_named(self, gene_name)
  }
  tx_mat <- pm$tx_mat[coding, grp$samples, drop = FALSE]
  # The gene order of rowsum()
  gene_levels <- sort(
    unique(pm$gene_ids[coding]),
    na.last = TRUE,
    method = "quick"
  )
  gene_idx <- match(pm$gene_ids[coding], gene_levels)
  is_prod <- pm$is_prod[coding]
  total <- gene_sums(tx_mat, gene_idx, gene_levels)
  frac <- gene_sums(
    tx_mat[is_prod, , drop = FALSE],
    gene_idx[is_prod],
    gene_levels
  ) /
    total
  # A total that is not finite gives no fraction
  frac[!is.finite(frac) | !is.finite(total)] <- NA_real_
  if (isTRUE(as_matrix)) {
    return(frac)
  }
  genes <- self@annot_data_genes[, c("gene_id", "gene_name")]
  res <- tibble::tibble(
    gene_id = rownames(frac) %||% character(),
    gene_name = genes$gene_name[match(rownames(frac), genes$gene_id)],
    mean_gene_expr = as.vector(rowMeans(total)),
    productive_fraction = as.vector(rowMeans(frac, na.rm = TRUE))
  )
  res$productive_fraction[is.nan(res$productive_fraction)] <- NA_real_
  cols <- if (is.null(group)) {
    list()
  } else {
    group_fraction_columns(frac, groups, "productive_fraction")
  }
  has_delta <- "delta_productive_fraction" %in% names(cols)
  reps <- resolve_inf_reps(self, inf_reps)
  delta_sd <- NULL
  if (!is.null(reps) && "productive_fraction_sd" %in% names(cols)) {
    cli::cli_abort(
      message = c(
        x = paste0(
          "A level of {.arg group} named {.val sd} gives the column ",
          "{.field productive_fraction_sd}, which holds the standard ",
          "deviation over the inferential replicates."
        ),
        i = "Rename the level, or use {.code inf_reps = FALSE}."
      ),
      class = "isoformic_group_levels"
    )
  }
  if (!is.null(reps)) {
    # Fitted on every transcript of the object, so that the factors of a gene
    # do not depend on the other genes asked for
    scale <- draw_assay_scale(
      self,
      resolve_assay_name(self, assay),
      rownames(tx_mat),
      colnames(tx_mat),
      fit_ids = rownames(pm$tx_mat)
    )
    rep_sd <- replicate_fraction_sd(
      reps,
      frac = frac,
      tx_mat = tx_mat,
      gene_idx = gene_idx,
      gene_levels = gene_levels,
      is_prod = is_prod,
      min_expr = min_expr,
      scale = scale,
      groups = if (has_delta) groups
    )
    res$productive_fraction_sd <- rep_sd$sd
    delta_sd <- rep_sd$delta_sd
  }
  if (!is.null(group)) {
    for (nm in names(cols)) {
      res[[nm]] <- cols[[nm]]
    }
    if (!is.null(delta_sd)) {
      res$delta_productive_fraction_sd <- delta_sd
      res$delta_z <- res$delta_productive_fraction /
        res$delta_productive_fraction_sd
      res$delta_z[!is.finite(res$delta_z)] <- NA_real_
      # No change in any draw is a z-score of zero, not a missing one
      no_change <- res$delta_productive_fraction %in%
        0 &
        res$delta_productive_fraction_sd %in% 0
      res$delta_z[no_change] <- 0
    }
  }
  if (has_delta) {
    res <- dplyr::arrange(
      res,
      dplyr::desc(abs(.data$delta_productive_fraction))
    )
  }
  if (!is.null(reference)) {
    if (is.character(reference)) {
      reference <- gtex_reference(
        self,
        tissue = reference,
        include_altered = include_altered
      )
    }
    ref <- reference[, intersect(
      c(
        "gene_id",
        "reference_tissue",
        "reference_productive_fraction",
        "reference_dominant_productivity"
      ),
      colnames(reference)
    )]
    res <- dplyr::left_join(res, ref, by = "gene_id")
    last_col <- if (!is.null(group)) {
      paste0("productive_fraction_", levels(groups)[nlevels(groups)])
    } else {
      "productive_fraction"
    }
    res$delta_vs_reference <- res[[last_col]] -
      res$reference_productive_fraction
  }
  filter_genes(res, gene_name)
}

#' Standard deviations of the productive fractions over the replicates
#'
#' One draw at a time: only the running moments are kept, never the
#' fractions of all the draws. The draws are put on the scale of the assay
#' (`scale`, see `draw_assay_scale()`), and a draw fraction counts only where
#' the assay gives the gene a fraction in that sample (`frac`), so that the
#' spread describes the samples behind the reported values. `sd` is the mean
#' over the samples of the sample standard deviation of each gene's fraction
#' over the draws (around the mean of the draws, not around the point
#' estimate); `delta_sd` is the standard deviation over the draws of the
#' change between the last and the first level of `groups` (`NULL` without
#' `groups`).
#' @keywords internal
#' @noRd
replicate_fraction_sd <- function(
  reps,
  frac,
  tx_mat,
  gene_idx,
  gene_levels,
  is_prod,
  min_expr,
  scale = NULL,
  groups = NULL
) {
  no_point <- is.na(frac)
  cell_stats <- running_moments(length(frac))
  delta_stats <- running_moments(nrow(frac))
  read_draw <- inf_reps_reader(reps)
  for (k in seq_along(reps)) {
    draw_frac <- draw_productive_fraction(
      read_draw(k),
      tx_ids = rownames(tx_mat),
      sample_ids = colnames(tx_mat),
      gene_idx = gene_idx,
      gene_levels = gene_levels,
      is_prod = is_prod,
      min_expr = min_expr,
      scale = scale
    )
    draw_frac[no_point] <- NA_real_
    cell_stats <- running_moments_add(cell_stats, as.vector(draw_frac))
    if (!is.null(groups)) {
      delta_stats <- running_moments_add(
        delta_stats,
        group_fraction_columns(
          draw_frac,
          groups,
          "productive_fraction"
        )$delta_productive_fraction
      )
    }
  }
  cell_sd <- matrix(running_moments_sd(cell_stats), nrow(frac), ncol(frac))
  sd <- as.vector(rowMeans(cell_sd, na.rm = TRUE))
  sd[is.nan(sd)] <- NA_real_
  list(
    sd = sd,
    delta_sd = if (!is.null(groups)) running_moments_sd(delta_stats)
  )
}

#' Factors that put count draws on the scale of an assay
#'
#' The inferential replicates are counts (Salmon `NumReads`, oarfish read
#' counts). For another assay, each draw is multiplied by `assay / counts`
#' of the same transcript and sample, the factor that turns the point count
#' into the point value of the assay (for TPM, the scaling of the sample over
#' the effective length). Where the count or the assay value is zero, the
#' factor is estimated:
#'
#' 1. With a `length` assay that the factors follow (TPM: `assay * length /
#'    counts` is one number per sample), from the length of the transcript
#'    in that sample.
#' 2. Otherwise as a transcript effect plus a sample effect on the log
#'    scale: the sample effects are the medians over the transcripts seen in
#'    every sample, the transcript effect the mean over its cells. A
#'    transcript without any cell takes its effect from its annotated length
#'    (`tx_length`, minus a fragment length fitted on the other transcripts)
#'    when the effects follow the lengths, and the median transcript effect
#'    when they do not.
#'
#' The terms shared by the transcripts (the sample terms and the fragment
#' length) are fitted on at most 5000 evenly spaced transcripts of `fit_ids`
#' (by default `tx_ids`), so that the factors of `tx_ids` do not depend on
#' which of the `fit_ids` were asked for. Returns a `tx_ids` x `sample_ids`
#' matrix, or `NULL` (the draws are used as they are) for the `"counts"`
#' assay and for an object without a `"counts"` assay.
#' @keywords internal
#' @noRd
draw_assay_scale <- function(
  self,
  assay,
  tx_ids,
  sample_ids,
  fit_ids = tx_ids
) {
  nms <- assay_names(self)
  if (identical(assay, "counts") || !"counts" %in% nms) {
    return(NULL)
  }
  if (length(fit_ids) > 5000L) {
    fit_ids <- fit_ids[round(seq(1, length(fit_ids), length.out = 5000L))]
  }
  as_double <- function(name, ids) {
    m <- tx_assay(self, name)
    m <- m[match(ids, rownames(m)), sample_ids, drop = FALSE]
    storage.mode(m) <- "double"
    m
  }
  log_ratio <- function(ids) {
    v <- as_double(assay, ids)
    n <- as_double("counts", ids)
    out <- log(v / n)
    out[!(n > 0 & v > 0 & is.finite(out))] <- NA_real_
    dimnames(out) <- list(ids, sample_ids)
    out
  }
  lr_tx <- log_ratio(tx_ids)
  scale <- exp(lr_tx)
  if (!anyNA(scale)) {
    return(scale)
  }
  lr_fit <- log_ratio(fit_ids)
  with_len <- "length" %in% nms && !identical(assay, "length")
  if (with_len) {
    per_len <- lr_fit + log(as_double("length", fit_ids))
    per_len[!is.finite(per_len)] <- NA_real_
    sample_k <- apply(per_len, 2L, stats::median, na.rm = TRUE)
    resid <- abs(per_len - rep(sample_k, each = nrow(per_len)))
    if (
      all(is.finite(sample_k)) &&
        isTRUE(stats::median(resid, na.rm = TRUE) < 0.01)
    ) {
      len <- as_double("length", tx_ids)
      fill <- is.na(scale) & is.finite(len) & len > 0
      scale[fill] <- (rep(exp(sample_k), each = nrow(scale)) / len)[fill]
      if (!anyNA(scale)) {
        return(scale)
      }
    }
  }
  # Additive fit on the log scale: sample effects from the transcripts seen
  # in every sample, then each transcript effect from its own cells
  sample_eff <- numeric(length(sample_ids))
  complete <- rowSums(is.na(lr_fit)) == 0L
  if (any(complete)) {
    centred <- lr_fit[complete, , drop = FALSE] -
      rowMeans(lr_fit[complete, , drop = FALSE])
    sample_eff <- apply(centred, 2L, stats::median)
  }
  tx_effect <- function(lr) {
    rowMeans(lr - rep(sample_eff, each = nrow(lr)), na.rm = TRUE)
  }
  tx_eff <- tx_effect(lr_tx)
  unseen <- !is.finite(tx_eff)
  if (any(unseen)) {
    tx_eff[unseen] <- unseen_tx_effect(
      self,
      fit_ids,
      tx_effect(lr_fit),
      tx_ids[unseen]
    )
  }
  fitted <- exp(outer(tx_eff, sample_eff, "+"))
  fill <- is.na(scale)
  scale[fill] <- fitted[fill]
  scale
}

#' Transcript effects of `draw_assay_scale()` for transcripts without a cell
#'
#' `fit_eff` holds the effects of the transcripts `fit_ids` (`NaN` for the
#' ones without a cell). The effect of each of `target_ids` comes from its
#' annotated length when the known effects follow minus the log of the
#' effective length, approximated as `tx_length - fragment` with at least a
#' tenth of `tx_length` (`fragment` from 0 to 500 nt, fitted on the known
#' effects); it is the median of the known effects otherwise.
#' @keywords internal
#' @noRd
unseen_tx_effect <- function(self, fit_ids, fit_eff, target_ids) {
  known <- is.finite(fit_eff)
  typical <- if (any(known)) stats::median(fit_eff[known]) else 0
  out <- rep(typical, length(target_ids))
  rd <- row_data(self)
  if (!"tx_length" %in% colnames(rd)) {
    return(out)
  }
  length_of <- function(ids) {
    as.double(rd$tx_length[match(ids, rd$transcript_id)])
  }
  fit_len <- length_of(fit_ids)
  known <- known & is.finite(fit_len) & fit_len > 0
  if (sum(known) < 10L) {
    return(out)
  }
  fit_len <- fit_len[known]
  fit_eff <- fit_eff[known]
  eff_len <- function(len, fragment) pmax(len - fragment, len / 10)
  spread <- function(x) stats::median(abs(x - stats::median(x)))
  fragments <- seq(0, 500, by = 10)
  scores <- vapply(
    fragments,
    function(f) spread(fit_eff + log(eff_len(fit_len, f))),
    numeric(1L)
  )
  # Use the lengths only when they explain the transcript effects
  if (min(scores) >= spread(fit_eff)) {
    return(out)
  }
  best <- fragments[[which.min(scores)]]
  offset <- stats::median(fit_eff + log(eff_len(fit_len, best)))
  from_len <- offset - log(eff_len(length_of(target_ids), best))
  out[is.finite(from_len)] <- from_len[is.finite(from_len)]
  out
}

#' Gene IDs of genes given by name or ID
#' @keywords internal
#' @noRd
gene_ids_named <- function(self, gene_name) {
  genes <- self@annot_data_genes
  unique(c(
    genes$gene_id[genes$gene_name %in% gene_name],
    gene_name[gene_name %in% genes$gene_id]
  ))
}

#' Inferential replicates to use, `NULL` when there are none
#'
#' `inf_reps = FALSE` switches the stored replicates off. The stored draws
#' are returned through `object_inf_reps()` (not copied):
#' `replicate_fraction_sd()` reads them one at a time through
#' `inf_reps_reader()`, which applies the selection of a subset, and
#' `draw_productive_fraction()` reads the rows and columns it needs from
#' each of them. Zero-row draws (a subset that keeps only the moments) give
#' `NULL`, as an object without replicates.
#' @keywords internal
#' @noRd
resolve_inf_reps <- function(self, inf_reps) {
  if (isFALSE(inf_reps)) {
    return(NULL)
  }
  reps <- inf_reps %||% object_inf_reps(self)
  if (is.null(reps)) {
    return(NULL)
  }
  if (!is.list(reps) || is.data.frame(reps) || !length(reps)) {
    cli::cli_abort(
      message = c(
        x = "{.arg inf_reps} must be a non-empty list of transcript matrices, or {.code FALSE}."
      ),
      class = "isoformic_invalid_inf_reps"
    )
  }
  if (!inf_reps_has_rows(reps)) {
    # Zero-row draws (only the moments of a subset): no draws to resample
    return(NULL)
  }
  reps
}

#' Sum the rows of `mat` per gene into a genes x samples matrix
#'
#' `gene_idx` gives, for each row, the position of its gene in `gene_levels`.
#' Genes without rows get zero.
#' @keywords internal
#' @noRd
gene_sums <- function(mat, gene_idx, gene_levels) {
  out <- matrix(
    0,
    nrow = length(gene_levels),
    ncol = ncol(mat),
    dimnames = list(gene_levels, colnames(mat))
  )
  if (nrow(mat)) {
    sums <- rowsum(mat, gene_idx, reorder = FALSE)
    out[as.integer(rownames(sums)), ] <- sums
  }
  out
}

#' Productive fraction per gene and sample of one inferential replicate
#'
#' Transcripts that the draw does not have count as zero ([add_uncertainty()]
#' keeps the draws of the expressed transcripts only). `scale` (rows
#' `tx_ids`, columns `sample_ids`, from `draw_assay_scale()`) puts the draw
#' on the scale of the assay before `min_expr` applies.
#' @keywords internal
#' @noRd
draw_productive_fraction <- function(
  draw,
  tx_ids,
  sample_ids,
  gene_idx,
  gene_levels,
  is_prod,
  min_expr = 0,
  scale = NULL
) {
  if (!is.matrix(draw)) {
    draw <- as.matrix(draw)
  }
  # A draw without rows has no transcript: every fraction is missing
  if (
    (nrow(draw) > 0L && is.null(rownames(draw))) ||
      !all(sample_ids %in% colnames(draw))
  ) {
    cli::cli_abort(
      message = c(
        x = "Each inferential replicate must be a matrix with the transcripts and samples of the assay in its dimnames."
      ),
      class = "isoformic_invalid_inf_reps"
    )
  }
  rows <- match(tx_ids, rownames(draw))
  found <- !is.na(rows)
  sub <- draw[rows[found], sample_ids, drop = FALSE]
  if (!is.double(sub)) {
    # Integer sums overflow to NA
    storage.mode(sub) <- "double"
  }
  if (!is.null(scale)) {
    sub <- sub * scale[found, , drop = FALSE]
  }
  if (isTRUE(min_expr > 0)) {
    sub[sub < min_expr] <- 0
  }
  prod_rows <- is_prod[found]
  idx <- gene_idx[found]
  total <- gene_sums(sub, idx, gene_levels)
  frac <- gene_sums(
    sub[prod_rows, , drop = FALSE],
    idx[prod_rows],
    gene_levels
  ) /
    total
  frac[!is.finite(frac) | !is.finite(total)] <- NA_real_
  frac
}

#' Running mean and variance of vectors that may have missing values
#'
#' Welford's update, per element, over the values that are not missing.
#' `running_moments_sd()` gives the sample standard deviation (`NA` with
#' fewer than two values).
#' @keywords internal
#' @noRd
running_moments <- function(n) {
  list(n = numeric(n), mean = numeric(n), m2 = numeric(n))
}

#' @keywords internal
#' @noRd
running_moments_add <- function(stats, x) {
  ok <- !is.na(x)
  n <- stats$n[ok] + 1
  d <- x[ok] - stats$mean[ok]
  mean_new <- stats$mean[ok] + d / n
  stats$m2[ok] <- stats$m2[ok] + d * (x[ok] - mean_new)
  stats$mean[ok] <- mean_new
  stats$n[ok] <- n
  stats
}

#' @keywords internal
#' @noRd
running_moments_sd <- function(stats) {
  out <- sqrt(pmax(stats$m2, 0) / (stats$n - 1))
  out[stats$n < 2] <- NA_real_
  out
}

#' @rdname productivity
#' @export
productivity_fractions <- S7::new_generic(
  "productivity_fractions",
  "self",
  function(
    self,
    group = NULL,
    levels = NULL,
    assay = "abundance",
    min_expr = 0,
    gene_name = NULL,
    ...
  ) {
    S7::S7_dispatch()
  }
)

S7::method(productivity_fractions, IsoformicExperiment) <- function(
  self,
  group = NULL,
  levels = NULL,
  assay = "abundance",
  min_expr = 0,
  gene_name = NULL,
  ...
) {
  pm <- productivity_matrix(self, assay, min_expr, include_altered = TRUE)
  grp <- group_samples(self, group, colnames(pm$tx_mat), levels = levels)
  groups <- grp$groups
  pm$tx_mat <- pm$tx_mat[, grp$samples, drop = FALSE]
  idx <- match(rownames(pm$tx_mat), pm$prod$transcript_id)
  detail <- productivity_detail(pm$prod)[idx]
  classes <- c(
    "productive",
    "altered_productive",
    "nmd",
    "retained_intron",
    "no_orf",
    "non_stop",
    "other",
    "noncoding",
    "uncertain"
  )
  total <- rowsum(pm$tx_mat, pm$gene_ids)
  genes <- self@annot_data_genes[, c("gene_id", "gene_name")]
  out <- lapply(classes, function(cls) {
    frac <- rowsum(pm$tx_mat * (detail == cls), pm$gene_ids) / total
    frac[!is.finite(frac)] <- NA_real_
    res <- tibble::tibble(
      gene_id = rownames(frac) %||% character(),
      gene_name = genes$gene_name[match(rownames(frac), genes$gene_id)],
      productivity = cls,
      fraction = as.vector(rowMeans(frac, na.rm = TRUE))
    )
    res$fraction[is.nan(res$fraction)] <- NA_real_
    if (!is.null(groups)) {
      cols <- group_fraction_columns(frac, groups, "fraction")
      for (nm in names(cols)) {
        res[[nm]] <- cols[[nm]]
      }
    }
    res
  })
  res <- dplyr::bind_rows(out)
  res <- res[!is.na(res$fraction) & res$fraction > 0, ]
  res$productivity <- factor(res$productivity, levels = classes)
  res <- dplyr::arrange(res, .data$gene_id, .data$productivity)
  filter_genes(res, gene_name)
}

#' @rdname productivity
#' @export
test_productive_fraction <- S7::new_generic(
  "test_productive_fraction",
  "self",
  function(
    self,
    group,
    levels = NULL,
    assay = "counts",
    include_altered = TRUE,
    min_expr = 0,
    min_gene_count = 0,
    gene_name = NULL,
    ...
  ) {
    S7::S7_dispatch()
  }
)

S7::method(test_productive_fraction, IsoformicExperiment) <- function(
  self,
  group,
  levels = NULL,
  assay = "counts",
  include_altered = TRUE,
  min_expr = 0,
  min_gene_count = 0,
  gene_name = NULL,
  ...
) {
  rlang::check_required(group)
  pm <- productivity_matrix(self, assay, min_expr, include_altered)
  total <- rowsum(pm$tx_mat, pm$gene_ids)
  productive <- rowsum(pm$tx_mat * pm$is_prod, pm$gene_ids)
  keep <- rownames(total) %in% pm$coding_genes
  total <- total[keep, , drop = FALSE]
  productive <- productive[keep, , drop = FALSE]
  setup <- two_group_setup(self, group, colnames(total), levels)
  total <- total[, setup$keep, drop = FALSE]
  productive <- productive[, setup$keep, drop = FALSE]
  # Genes with too few reads are left out before the test, so that they do
  # not enter the BH adjustment
  expressed <- count_filtered_genes(
    self,
    pm$prod,
    setup$samples,
    setup$groups,
    min_gene_count
  )
  if (!is.null(expressed)) {
    keep_genes <- rownames(total) %in% expressed
    total <- total[keep_genes, , drop = FALSE]
    productive <- productive[keep_genes, , drop = FALSE]
  }
  test <- two_group_proportion_test(productive, total, setup$groups)
  genes <- self@annot_data_genes[, c("gene_id", "gene_name")]
  res <- dplyr::bind_cols(
    tibble::tibble(
      gene_id = rownames(total) %||% character(),
      gene_name = genes$gene_name[match(rownames(total), genes$gene_id)],
      mean_gene_expr = as.vector(rowMeans(total))
    ),
    test
  )
  names(res) <- sub("^proportion_", "productive_fraction_", names(res))
  names(res)[names(res) == "delta"] <- "delta_productive_fraction"
  res <- dplyr::arrange(res, .data$pvalue)
  filter_genes(res, gene_name)
}

#' @rdname productivity
#' @export
run_productivity_enrichment <- S7::new_generic(
  "run_productivity_enrichment",
  "self",
  function(
    self,
    genesets_list,
    group = NULL,
    levels = NULL,
    assay = "abundance",
    min_size = 15,
    max_size = 500,
    rank_by = c("delta", "delta_z", "statistic"),
    include_altered = TRUE,
    min_expr = 0,
    min_gene_count = 0,
    inf_reps = NULL,
    workers = getOption("isoformic.workers", 1L),
    ...
  ) {
    S7::S7_dispatch()
  }
)

S7::method(run_productivity_enrichment, IsoformicExperiment) <- function(
  self,
  genesets_list,
  group = NULL,
  levels = NULL,
  assay = "abundance",
  min_size = 15,
  max_size = 500,
  rank_by = c("delta", "delta_z", "statistic"),
  include_altered = TRUE,
  min_expr = 0,
  min_gene_count = 0,
  inf_reps = NULL,
  workers = getOption("isoformic.workers", 1L),
  ...
) {
  rlang::check_installed("fgsea")
  rlang::check_required(genesets_list)
  if (is.null(group)) {
    cli::cli_abort(
      message = c(
        x = "{.arg group} is required to rank genes by change in productive fraction."
      ),
      class = "isoformic_group_levels"
    )
  }
  rank_by <- rlang::arg_match(rank_by)
  check_genesets(genesets_list)
  if (identical(rank_by, "statistic")) {
    if (!identical(assay, "abundance") || !is.null(inf_reps)) {
      cli::cli_warn(
        message = c(
          "!" = paste(
            "{.code rank_by = \"statistic\"} tests the {.val counts} assay",
            "and does not use {.arg assay} or {.arg inf_reps}."
          )
        ),
        class = "isoformic_ignored_argument"
      )
    }
    # The test needs counts and two levels: the first and the last, the two
    # that the change in productive fraction compares
    lv <- levels(
      group_samples(
        self,
        group,
        col_names(self),
        levels = levels,
        min_levels = 2L
      )$groups
    )
    tested <- test_productive_fraction(
      self,
      group = group,
      levels = lv[c(1L, length(lv))],
      include_altered = include_altered,
      min_expr = min_expr,
      min_gene_count = min_gene_count
    )
    res <- productivity_enrichment(
      tested,
      genesets_list,
      rank_by = rank_by,
      min_size = min_size,
      max_size = max_size,
      workers = workers,
      ...
    )
    return(res)
  }
  pf <- productive_fraction(
    self,
    group = group,
    levels = levels,
    assay = assay,
    include_altered = include_altered,
    min_expr = min_expr,
    min_gene_count = min_gene_count,
    # The change itself does not need the replicates: they are not read
    inf_reps = if (identical(rank_by, "delta")) FALSE else inf_reps
  )
  if (!"delta_productive_fraction" %in% colnames(pf)) {
    cli::cli_abort(
      message = c(
        x = "{.arg group} must have at least two levels to rank genes by change in productive fraction."
      ),
      class = "isoformic_group_levels"
    )
  }
  productivity_enrichment(
    pf,
    genesets_list,
    rank_by = rank_by,
    min_size = min_size,
    max_size = max_size,
    workers = workers,
    ...
  )
}

# nolint start: object_name_linter.
S7::method(run_productivity_enrichment, S7::class_data.frame) <- function(
  self,
  genesets_list,
  group = NULL,
  levels = NULL,
  assay = "abundance",
  min_size = 15,
  max_size = 500,
  rank_by = c("delta", "delta_z", "statistic"),
  include_altered = TRUE,
  min_expr = 0,
  min_gene_count = 0,
  inf_reps = NULL,
  workers = getOption("isoformic.workers", 1L),
  ...
) {
  rlang::check_installed("fgsea")
  rlang::check_required(genesets_list)
  rank_by <- rlang::arg_match(rank_by)
  check_genesets(genesets_list)
  pf <- tibble::as_tibble(self)
  missing_cols <- setdiff(
    c("gene_name", "delta_productive_fraction"),
    colnames(pf)
  )
  if (length(missing_cols)) {
    cli::cli_abort(
      message = c(
        x = "{.arg self} is missing column{?s}: {.field {missing_cols}}.",
        i = "Pass the result of {.fn productive_fraction} for two or more groups, or of {.fn test_productive_fraction}."
      ),
      class = "isoformic_invalid_productive_fraction"
    )
  }
  needed <- c("delta_z", "delta_productive_fraction_sd")
  if (identical(rank_by, "delta_z") && !all(needed %in% colnames(pf))) {
    missing_cols <- setdiff(needed, colnames(pf))
    cli::cli_abort(
      message = c(
        x = "{.code rank_by = \"delta_z\"} needs the column{?s} {.field {missing_cols}} in {.arg self}.",
        i = "Compute the table with {.fn productive_fraction} on an object with inferential replicates."
      ),
      class = "isoformic_no_inferential_replicates"
    )
  }
  productivity_enrichment(
    pf,
    genesets_list,
    rank_by = rank_by,
    min_size = min_size,
    max_size = max_size,
    workers = workers,
    ...
  )
}
# nolint end

#' Gene set enrichment on a `productive_fraction()` table
#' @keywords internal
#' @noRd
productivity_enrichment <- function(
  pf,
  genesets_list,
  rank_by,
  min_size,
  max_size,
  workers = 1L,
  ...
) {
  rank_col <- switch(
    rank_by,
    delta = "delta_productive_fraction",
    delta_z = "delta_z",
    statistic = "signed_z"
  )
  if (identical(rank_by, "statistic")) {
    missing_cols <- setdiff(c("statistic", "pvalue"), colnames(pf))
    if (length(missing_cols)) {
      cli::cli_abort(
        message = c(
          x = "{.code rank_by = \"statistic\"} needs the column{?s} {.field {missing_cols}}.",
          i = "Pass the result of {.fn test_productive_fraction}, or the object itself."
        ),
        class = "isoformic_invalid_productive_fraction"
      )
    }
    # The degrees of freedom differ between genes: compare them through the
    # p-value, on the normal scale
    pf$signed_z <- signed_normal_quantile(pf$statistic, pf$pvalue)
  }
  if (
    identical(rank_by, "delta_z") &&
      !all(c("delta_z", "delta_productive_fraction_sd") %in% colnames(pf))
  ) {
    cli::cli_abort(
      message = c(
        x = "{.code rank_by = \"delta_z\"} needs inferential replicates.",
        i = "Pass them with {.arg inf_reps} or store them with {.fn add_uncertainty}."
      ),
      class = "isoformic_no_inferential_replicates"
    )
  }
  if (identical(rank_by, "delta_z")) {
    # A change that is the same in every draw is the most certain one, not a
    # missing one: it gets the largest finite |z| of the table
    pf$delta_z <- certainty_z(
      pf$delta_productive_fraction,
      pf$delta_productive_fraction_sd
    )
    if (!any(is.finite(pf$delta_z))) {
      cli::cli_abort(
        message = c(
          x = "{.code rank_by = \"delta_z\"} needs at least two inferential replicates.",
          i = "No gene has a standard deviation of its change in productive fraction."
        ),
        class = "isoformic_no_inferential_replicates"
      )
    }
  }
  gene_rank_enrichment(
    pf,
    rank_col = rank_col,
    genesets_list = genesets_list,
    experiment = "productive_fraction",
    rank_by = rank_by,
    min_size = min_size,
    max_size = max_size,
    workers = workers,
    ...
  )
}

#' Signed normal quantile of a two-sided p-value
#'
#' A p-value of zero is the smallest of the vector, not an infinite score.
#' @keywords internal
#' @noRd
signed_normal_quantile <- function(statistic, pvalue) {
  z <- stats::qnorm(pvalue / 2, lower.tail = FALSE)
  finite <- z[is.finite(z)]
  z[z %in% Inf] <- if (length(finite)) max(finite) else 1
  sign(statistic) * z
}

#' @rdname productivity
#' @export
poison_exons <- S7::new_generic(
  "poison_exons",
  "self",
  function(self, gene_name = NULL, ...) {
    S7::S7_dispatch()
  }
)

S7::method(poison_exons, IsoformicExperiment) <- function(
  self,
  gene_name = NULL,
  ...
) {
  prod <- transcript_productivity(self, gene_name = gene_name)
  exons <- exons_of_transcripts(
    self,
    prod$transcript_id,
    columns = c("transcript_id", "seqid", "start_pos", "end_pos", "strand")
  )
  exons$productivity <- as.character(prod$productivity[match(
    exons$transcript_id,
    prod$transcript_id
  )])
  exons$gene_id <- prod$gene_id[match(exons$transcript_id, prod$transcript_id)]
  summary <- poison_exon_candidates(exons)
  gene_has_productive <- stats::setNames(
    stats::ave(
      prod$productivity %in% c("productive", "altered_productive"),
      prod$gene_id,
      FUN = any
    ),
    prod$gene_id
  )
  productive_exons <- exons[
    exons$productivity %in% c("productive", "altered_productive"),
    c("gene_id", "start_pos", "end_pos")
  ]
  names(productive_exons) <- c("gene_id", "p_start", "p_end")
  productive_exons <- dplyr::distinct(productive_exons)
  candidates <- summary |>
    dplyr::filter(gene_has_productive[.data$gene_id])
  overlapping <- candidates |>
    dplyr::inner_join(
      productive_exons,
      by = "gene_id",
      relationship = "many-to-many"
    ) |>
    dplyr::filter(
      .data$start_pos <= .data$p_end,
      .data$end_pos >= .data$p_start
    ) |>
    dplyr::distinct(.data$gene_id, .data$start_pos, .data$end_pos)
  res <- candidates |>
    dplyr::anti_join(overlapping, by = c("gene_id", "start_pos", "end_pos")) |>
    dplyr::left_join(
      self@annot_data_genes[, c("gene_id", "gene_name")],
      by = "gene_id"
    ) |>
    dplyr::relocate("gene_name", .after = "gene_id") |>
    dplyr::mutate(exon_length = .data$end_pos - .data$start_pos + 1L) |>
    dplyr::select(-"n_other") |>
    dplyr::arrange(.data$gene_id, .data$start_pos)
  res
}

#' Exons used only by unproductive transcripts, one row per exon
#'
#' Exons are grouped by gene and coordinates, in the order in which they
#' first appear. Only the groups whose transcripts are all unproductive are
#' kept, and only those get the list of their transcripts. The counts come
#' from one grouping of all rows, without a closure per group.
#' @keywords internal
#' @noRd
poison_exon_candidates <- function(exons) {
  cols <- c("gene_id", "seqid", "start_pos", "end_pos", "strand")
  gid <- dplyr::group_indices(dplyr::group_by(
    exons[, cols],
    dplyr::across(dplyr::all_of(cols))
  ))
  n <- max(c(0L, gid))
  unproductive <- exons$productivity == "unproductive"
  n_unproductive <- tabulate(gid[unproductive %in% TRUE], n)
  n_other <- tabulate(gid[unproductive %in% FALSE], n)
  # A missing class makes the counts of the old summary missing
  has_na <- tabulate(gid[is.na(unproductive)], n) > 0L
  keep_group <- n_unproductive > 0L & n_other == 0L & !has_na
  # Groups in order of first appearance, as summarise(.by = ) gives them
  first <- which(!duplicated(gid))
  first <- first[keep_group[gid[first]]]
  res <- exons[first, cols]
  res$n_unproductive <- n_unproductive[gid[first]]
  res$n_other <- n_other[gid[first]]
  rows <- keep_group[gid]
  ids <- split(exons$transcript_id[rows], factor(gid[rows], gid[first]))
  res$transcript_ids <- vapply(
    ids,
    function(x) paste(sort(unique(x)), collapse = ";"),
    character(1L),
    USE.NAMES = FALSE
  )
  res
}

#' Exons of some transcripts of an object
#'
#' From the annotation cache when the exon table is there, otherwise only the
#' rows of the transcripts and the requested columns are read.
#' @keywords internal
#' @noRd
exons_of_transcripts <- function(self, transcript_id, columns) {
  annot_path <- prop_raw(self, "annot_path")
  if (is.null(annot_path) || annot_table_cached(self, "exons")) {
    exons <- self@annot_data_exons
    return(exons[exons$transcript_id %in% transcript_id, columns])
  }
  read_annotation_table(
    annot_path,
    "exons",
    columns = columns,
    filter = list(transcript_id = transcript_id)
  )
}

#' @rdname productivity
#' @export
benchmark_productivity <- S7::new_generic(
  "benchmark_productivity",
  "self",
  function(self, nmd_rule = 50, use_sequence = TRUE, ...) {
    S7::S7_dispatch()
  }
)

S7::method(benchmark_productivity, IsoformicExperiment) <- function(
  self,
  nmd_rule = 50,
  use_sequence = TRUE,
  ...
) {
  rules <- transcript_productivity(
    self,
    nmd_rule = nmd_rule,
    use_biotype = FALSE,
    use_sequence = use_sequence,
    use_evidence = FALSE
  )
  benchmark_table(rules)
}

#' @keywords internal
#' @noRd
benchmark_table <- function(rules) {
  cls <- as.character(rules$productivity)
  type <- rules$transcript_type
  rate <- function(mask, target) {
    if (!any(mask, na.rm = TRUE)) {
      return(NA_real_)
    }
    mean(cls[mask %in% TRUE] %in% target)
  }
  summary <- tibble::tibble(
    metric = c(
      "nmd_biotype_called_unproductive",
      "retained_intron_called_unproductive",
      "protein_coding_called_productive",
      "canonical_called_productive",
      "lncrna_called_noncoding"
    ),
    n = c(
      sum(type %in% "nonsense_mediated_decay"),
      sum(type %in% "retained_intron"),
      sum(type %in% "protein_coding"),
      sum(rules$is_canonical %in% TRUE & !cls %in% "noncoding"),
      sum(type %in% "lncRNA")
    ),
    value = c(
      rate(type %in% "nonsense_mediated_decay", "unproductive"),
      rate(type %in% "retained_intron", "unproductive"),
      rate(type %in% "protein_coding", c("productive", "altered_productive")),
      rate(rules$is_canonical %in% TRUE & !cls %in% "noncoding", "productive"),
      rate(type %in% "lncRNA", "noncoding")
    ),
    source = unique(rules$source)[1L]
  )
  confusion <- rules |>
    dplyr::count(.data$transcript_type, .data$productivity, name = "n") |>
    dplyr::arrange(.data$transcript_type, .data$productivity)
  list(summary = summary, confusion = confusion)
}
