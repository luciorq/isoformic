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
#'   `retained_intron`, `protein_coding_CDS_not_defined`, or a longest ORF
#'   without coding potential), `"noncoding"` (transcripts of genes without
#'   any coding isoform) or `"uncertain"` (incomplete CDS annotation without
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
#'   log odds ratio, BH adjustment).
#' * `run_productivity_enrichment()` ranks genes by the change in
#'   productive fraction between two groups and runs gene set enrichment
#'   (`fgsea`) on gene symbols, which captures regulated unproductive
#'   splicing (AS-NMD) directly.
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
#' @param self An [IsoformicExperiment].
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
#'   per-sample values.
#' @param levels Order of the group levels; the change is `last - first`.
#'   For the tests (`test_productive_fraction()`, `test_apa()`), the two
#'   groups to compare; samples of other groups are ignored.
#' @param min_expr Minimum assay value for an isoform to contribute
#'   (default 0): isoforms below it count as zero.
#' @param inf_reps Optional list of transcript matrices (inferential
#'   replicates with the same dimnames as the assay); `NULL` uses
#'   `metadata(x)$inf_reps` when present.
#' @param reference Optional tibble from [gtex_reference()], or a GTEx
#'   tissue name, giving the healthy-tissue productive fraction per gene.
#' @param as_matrix Logical, return the genes x samples matrix of
#'   productive fractions instead of the summary table.
#' @param genesets_list Named list of gene symbol vectors (see [read_gmt()]).
#' @param min_size,max_size Gene set size limits passed to `fgsea`.
#' @param ... Passed on to methods or to `fgsea::fgseaMultilevel()`.
#'
#' @returns Tibbles, a matrix for `productive_fraction(as_matrix = TRUE)`,
#'   a named character vector for `productivity_palette()`, a list with
#'   `summary` and `confusion` tibbles for `benchmark_productivity()`.
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
    "processed_transcript"
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
  tx_mat <- assay(self, assay)
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
  inf_reps = NULL,
  reference = NULL,
  ...
) {
  pm <- productivity_matrix(self, assay, min_expr, include_altered)
  fraction_of <- function(tx_mat) {
    total <- rowsum(tx_mat, pm$gene_ids)
    frac <- rowsum(tx_mat * pm$is_prod, pm$gene_ids) / total
    frac[!is.finite(frac)] <- NA_real_
    frac[rownames(frac) %in% pm$coding_genes, , drop = FALSE]
  }
  frac <- fraction_of(pm$tx_mat)
  if (isTRUE(as_matrix)) {
    return(frac)
  }
  total <- rowsum(pm$tx_mat, pm$gene_ids)[rownames(frac), , drop = FALSE]
  genes <- self@annot_data_genes[, c("gene_id", "gene_name")]
  res <- tibble::tibble(
    gene_id = rownames(frac) %||% character(),
    gene_name = genes$gene_name[match(rownames(frac), genes$gene_id)],
    mean_gene_expr = as.vector(rowMeans(total)),
    productive_fraction = as.vector(rowMeans(frac, na.rm = TRUE))
  )
  res$productive_fraction[is.nan(res$productive_fraction)] <- NA_real_
  reps <- resolve_inf_reps(self, inf_reps, pm$tx_mat)
  rep_fracs <- NULL
  if (!is.null(reps)) {
    rep_fracs <- lapply(reps, function(m) {
      if (isTRUE(min_expr > 0)) {
        m[m < min_expr] <- 0
      }
      fraction_of(m)[rownames(frac), , drop = FALSE]
    })
    per_sample_sd <- sqrt(
      Reduce(
        `+`,
        lapply(rep_fracs, function(f) {
          (f - frac)^2
        })
      ) /
        length(rep_fracs)
    )
    res$productive_fraction_sd <- as.vector(rowMeans(
      per_sample_sd,
      na.rm = TRUE
    ))
    res$productive_fraction_sd[is.nan(res$productive_fraction_sd)] <- NA_real_
  }
  if (!is.null(group)) {
    groups <- sample_groups(self, group, colnames(frac), levels = levels)
    cols <- group_fraction_columns(frac, groups, "productive_fraction")
    for (nm in names(cols)) {
      res[[nm]] <- cols[[nm]]
    }
    if (!is.null(rep_fracs) && "delta_productive_fraction" %in% names(cols)) {
      rep_delta <- vapply(
        rep_fracs,
        function(f) {
          group_fraction_columns(
            f,
            groups,
            "productive_fraction"
          )$delta_productive_fraction
        },
        numeric(nrow(frac))
      )
      rep_delta <- matrix(rep_delta, nrow = nrow(frac))
      res$delta_productive_fraction_sd <- as.vector(apply(
        rep_delta,
        1L,
        stats::sd,
        na.rm = TRUE
      ))
      res$delta_z <- res$delta_productive_fraction /
        res$delta_productive_fraction_sd
      res$delta_z[!is.finite(res$delta_z)] <- NA_real_
    }
    if ("delta_productive_fraction" %in% names(cols)) {
      res <- dplyr::arrange(
        res,
        dplyr::desc(abs(.data$delta_productive_fraction))
      )
    }
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

#' Inferential replicates as a list of matrices aligned with `tx_mat`
#' @keywords internal
#' @noRd
resolve_inf_reps <- function(self, inf_reps, tx_mat) {
  reps <- inf_reps %||% self@metadata$inf_reps
  if (is.null(reps)) {
    return(NULL)
  }
  if (!is.list(reps) || !length(reps)) {
    cli::cli_abort(
      message = c(
        x = "{.arg inf_reps} must be a non-empty list of transcript matrices."
      ),
      class = "isoformic_invalid_inf_reps"
    )
  }
  lapply(reps, function(m) {
    m <- as.matrix(m)
    missing_tx <- setdiff(rownames(tx_mat), rownames(m))
    if (length(missing_tx) && all(colnames(tx_mat) %in% colnames(m))) {
      # add_uncertainty() keeps the draws of expressed transcripts only
      m <- rbind(
        m[, colnames(tx_mat), drop = FALSE],
        matrix(
          0,
          nrow = length(missing_tx),
          ncol = ncol(tx_mat),
          dimnames = list(missing_tx, colnames(tx_mat))
        )
      )
    }
    if (
      !all(rownames(tx_mat) %in% rownames(m)) ||
        !all(colnames(tx_mat) %in% colnames(m))
    ) {
      cli::cli_abort(
        message = c(
          x = "Each inferential replicate must have the transcripts and samples of the assay in its dimnames."
        ),
        class = "isoformic_invalid_inf_reps"
      )
    }
    m[rownames(tx_mat), colnames(tx_mat), drop = FALSE]
  })
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
  groups <- if (!is.null(group)) {
    sample_groups(self, group, colnames(total), levels = levels)
  } else {
    NULL
  }
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
    group,
    levels = NULL,
    assay = "abundance",
    min_size = 15,
    max_size = 500,
    ...
  ) {
    S7::S7_dispatch()
  }
)

S7::method(run_productivity_enrichment, IsoformicExperiment) <- function(
  self,
  genesets_list,
  group,
  levels = NULL,
  assay = "abundance",
  min_size = 15,
  max_size = 500,
  ...
) {
  rlang::check_installed("fgsea")
  rlang::check_required(genesets_list)
  rlang::check_required(group)
  pf <- productive_fraction(self, group = group, levels = levels, assay = assay)
  if (!"delta_productive_fraction" %in% colnames(pf)) {
    cli::cli_abort(
      message = c(
        x = "{.arg group} must have at least two levels to rank genes by change in productive fraction."
      ),
      class = "isoformic_group_levels"
    )
  }
  pf <- pf[!is.na(pf$delta_productive_fraction) & !is.na(pf$gene_name), ]
  pf <- pf[!duplicated(pf$gene_name), ]
  ranks <- sort(stats::setNames(pf$delta_productive_fraction, pf$gene_name))
  res <- fgsea::fgseaMultilevel(
    pathways = genesets_list,
    stats = ranks,
    minSize = min_size,
    maxSize = max_size,
    eps = 0,
    ...
  )
  res <- tibble::as_tibble(as.data.frame(res))
  res$experiment <- "productive_fraction"
  res
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
  exons <- self@annot_data_exons
  exons <- exons[
    exons$transcript_id %in% prod$transcript_id,
    c("transcript_id", "seqid", "start_pos", "end_pos", "strand")
  ]
  exons$productivity <- as.character(prod$productivity[match(
    exons$transcript_id,
    prod$transcript_id
  )])
  exons$gene_id <- prod$gene_id[match(exons$transcript_id, prod$transcript_id)]
  summary <- exons |>
    dplyr::summarise(
      n_unproductive = sum(.data$productivity == "unproductive"),
      n_other = sum(.data$productivity != "unproductive"),
      transcript_ids = paste(sort(unique(.data$transcript_id)), collapse = ";"),
      .by = c("gene_id", "seqid", "start_pos", "end_pos", "strand")
    )
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
    dplyr::filter(
      .data$n_unproductive > 0L,
      .data$n_other == 0L,
      gene_has_productive[.data$gene_id]
    )
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
