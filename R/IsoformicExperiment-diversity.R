# Transcriptomic diversity summaries -----------------------------------------

#' Transcriptomic Diversity Summaries
#'
#' @description
#' These generics quantify how expression is distributed across the isoforms
#' of each gene, which is the question an [IsoformicExperiment] is built to
#' answer. All computations use matrix operations on the cached assays and
#' scale to full transcriptome experiments.
#'
#' * `isoform_fraction(x)`: matrix (transcripts x samples) of the fraction
#'   of its gene's expression carried by each transcript. `NA` where the gene
#'   is not expressed in a sample.
#' * `tx_diversity(x)`: one row per gene summarising isoform usage:
#'   number of annotated, assayed and expressed isoforms, number of
#'   transcript types, the dominant isoform and its mean fraction, and the
#'   *effective number of isoforms* (`exp(Shannon entropy)` of the mean
#'   isoform fractions, 1 = a single isoform carries all expression).
#' * `dominant_isoform(x, group)`: the dominant isoform of each gene,
#'   overall or per level of a `col_data` column, flagging genes whose
#'   dominant isoform changes between groups (`dominant_switch`).
#' * `isoform_switches(x)`: relates transcript- and gene-level differential
#'   expression per gene: counts of significant isoforms up and down, the
#'   gene-level call and a `switch_type` classification
#'   (`"opposite_direction"`, `"det_not_deg"`, `"consistent"`).
#' * `exon_diversity(x, gene_name)`: per gene, the number of distinct exons,
#'   how many are constitutive (shared by all isoforms) versus alternative,
#'   and the mean number of exons per isoform.
#' * `haplotype_summary(x, group)`: for haplotype-resolved annotations
#'   (see [merge_annotations()]), expression of each gene per haplotype, the
#'   haplotype fraction (allelic imbalance), the dominant isoform of each
#'   haplotype (compared through `source_transcript_id`) and whether the
#'   dominant isoform differs between haplotypes
#'   (`haplotype_isoform_switch`). Genes are matched across haplotypes by
#'   `source_gene_id`, falling back to `gene_name`.
#'
#' `dominant_isoform()` also reports `dominant_is_canonical`, whether the
#' dominant isoform is the MANE Select / Ensembl canonical / RefSeq Select
#' transcript of its gene. When the object carries uncertainty tables (see
#' [add_uncertainty()]) it adds `dominant_inf_rv`, the inferential relative
#' variance of the dominant isoform, and `dominant_ambiguous`, whether it
#' cannot be told apart from a sibling isoform. `isoform_switches()` adds
#' `n_det_confident`, the number of significant isoforms with high
#' confidence, once [de_confidence()] results are stored.
#'
#' @param self An [IsoformicExperiment].
#' @param assay Name of the transcript assay to use (default `"abundance"`,
#'   falling back to the first assay).
#' @param min_expr Minimum assay value for a transcript to count as
#'   expressed in a sample.
#' @param min_samples Minimum number of samples in which a transcript must
#'   be expressed to count as expressed.
#' @param group Name of a `col_data` column defining sample groups, or
#'   `NULL` for all samples together.
#' @param gene_name Optional character vector of gene names or IDs to
#'   restrict the result to.
#' @param confident_only Logical. Count only the significant isoforms with
#'   high confidence (see [de_confidence()]) in `isoform_switches()`.
#' @param ... Passed on to methods.
#'
#' @returns A matrix for `isoform_fraction()`, tibbles otherwise.
#'
#' @examples
#' iso <- example_isoformic()
#' isoform_fraction(iso)[1:3, 1:3]
#' head(tx_diversity(iso))
#' dominant_isoform(iso, group = "condition")
#' isoform_switches(iso)
#' exon_diversity(iso, gene_name = c("RBPJ", "EGFR"))
#'
#' @name diversity
#' @rdname diversity
NULL

#' @rdname diversity
#' @export
isoform_fraction <- S7::new_generic(
  "isoform_fraction",
  "self",
  function(self, assay = "abundance", ...) {
    S7::S7_dispatch()
  }
)

S7::method(isoform_fraction, IsoformicExperiment) <- function(
  self,
  assay = "abundance",
  ...
) {
  assay <- resolve_assay_name(self, assay)
  signature <- paste(assay, assays_signature(self), annot_full_signature(self))
  object_cache_fetch(self, "isoform_fraction", signature, function() {
    tx_mat <- assay(self, assay)
    gene_mat <- gene_assay(self, assay)
    gene_ids <- gene_ids_of(self, rownames(tx_mat))
    gene_tot <- gene_mat[gene_ids, , drop = FALSE]
    frac <- tx_mat / gene_tot
    frac[!is.finite(frac)] <- NA_real_
    dimnames(frac) <- dimnames(tx_mat)
    frac
  })
}

resolve_assay_name <- function(self, assay) {
  nms <- assay_names(self)
  if (!length(nms)) {
    cli::cli_abort(
      message = c(x = "The object has no assays."),
      class = "isoformic_no_assays"
    )
  }
  if (is.null(assay) || !assay %in% nms) {
    if (!is.null(assay) && !identical(assay, "abundance")) {
      cli::cli_abort(
        message = c(
          x = "Assay {.val {assay}} not found.",
          i = "Available assays: {.val {nms}}."
        ),
        class = "isoformic_assay_not_found"
      )
    }
    assay <- nms[1L]
  }
  assay
}

#' @rdname diversity
#' @export
tx_diversity <- S7::new_generic(
  "tx_diversity",
  "self",
  function(
    self,
    assay = "abundance",
    min_expr = 1,
    min_samples = 1L,
    gene_name = NULL,
    ...
  ) {
    S7::S7_dispatch()
  }
)

S7::method(tx_diversity, IsoformicExperiment) <- function(
  self,
  assay = "abundance",
  min_expr = 1,
  min_samples = 1L,
  gene_name = NULL,
  ...
) {
  assay <- resolve_assay_name(self, assay)
  tx_mat <- assay(self, assay)
  frac <- isoform_fraction(self, assay = assay)
  annot <- row_data(self)
  gene_id <- annot$gene_id

  expressed <- rowSums(tx_mat >= min_expr) >= min_samples
  mean_frac <- rowMeans(frac, na.rm = TRUE)
  mean_frac[is.nan(mean_frac)] <- NA_real_
  mean_expr <- rowMeans(tx_mat)

  # Per-gene aggregates with rowsum() (C-level, sorted by gene_id)
  n_assayed <- rowsum(rep(1L, length(gene_id)), gene_id)
  genes <- rownames(n_assayed)
  n_expressed <- rowsum(as.integer(expressed), gene_id)
  gene_expr <- rowsum(mean_expr, gene_id)
  eff <- effective_number_by_group(mean_frac, gene_id, genes)
  dominant <- first_by_group(gene_id, mean_frac, genes)
  type_tbl <- type_summary_by_group(
    annot$transcript_type[expressed],
    gene_id[expressed],
    genes
  )

  n_annot <- dplyr::count(
    self@annot_data_transcripts,
    .data$gene_id,
    name = "n_tx_annotated"
  )
  gene_annot <- self@annot_data_genes[, c("gene_id", "gene_name", "gene_type")]
  res <- tibble::tibble(
    gene_id = genes,
    n_tx_assayed = as.vector(n_assayed),
    n_tx_expressed = as.vector(n_expressed),
    n_tx_types = type_tbl$n_types,
    tx_types = type_tbl$types,
    mean_gene_expr = as.vector(gene_expr),
    dominant_transcript_id = annot$transcript_id[dominant],
    dominant_transcript_name = annot$transcript_name[dominant],
    dominant_transcript_type = annot$transcript_type[dominant],
    dominant_fraction = unname(mean_frac[dominant]),
    effective_isoforms = eff
  ) |>
    dplyr::left_join(n_annot, by = "gene_id") |>
    dplyr::left_join(gene_annot, by = "gene_id") |>
    dplyr::relocate(
      "gene_name",
      "gene_type",
      "n_tx_annotated",
      .after = "gene_id"
    ) |>
    dplyr::arrange(
      dplyr::desc(.data$effective_isoforms),
      dplyr::desc(.data$n_tx_expressed)
    )
  filter_genes(res, gene_name)
}

#' Index of the row with the largest value within each group
#'
#' Returns, for each group in `groups` (sorted order), the index into the
#' original vectors of the maximal `values` (NA last). Fully vectorised.
#'
#' @keywords internal
#' @noRd
first_by_group <- function(group, values, groups = sort(unique(group))) {
  ord <- order(group, -values, na.last = TRUE)
  first <- ord[!duplicated(group[ord])]
  first[match(groups, group[first])]
}

#' Effective number of isoforms per group: exp(Shannon entropy) of shares
#' @keywords internal
#' @noRd
effective_number_by_group <- function(p, group, groups) {
  ok <- !is.na(p) & p > 0
  p <- p[ok]
  g <- group[ok]
  if (!length(p)) {
    return(rep(NA_real_, length(groups)))
  }
  totals <- rowsum(p, g)
  share <- p / as.vector(totals)[match(g, rownames(totals))]
  h <- rowsum(-share * log(share), g)
  res <- exp(as.vector(h))[match(groups, rownames(h))]
  res
}

#' Number of distinct transcript types and their collapsed labels per group
#' @keywords internal
#' @noRd
type_summary_by_group <- function(types, group, groups) {
  if (!length(types)) {
    return(list(
      n_types = rep(0L, length(groups)),
      types = rep("", length(groups))
    ))
  }
  pairs <- unique(data.frame(
    group = group,
    type = types,
    stringsAsFactors = FALSE
  ))
  pairs <- pairs[order(pairs$group, pairs$type), ]
  n_types <- rowsum(rep(1L, nrow(pairs)), pairs$group)
  collapsed <- vapply(
    split(pairs$type, pairs$group),
    paste,
    character(1L),
    collapse = ";"
  )
  n_types_out <- as.vector(n_types)[match(groups, rownames(n_types))]
  n_types_out[is.na(n_types_out)] <- 0L
  types_out <- unname(collapsed[match(groups, names(collapsed))])
  types_out[is.na(types_out)] <- ""
  list(n_types = n_types_out, types = types_out)
}

filter_genes <- function(res, gene_name) {
  if (is.null(gene_name)) {
    return(res)
  }
  res[res$gene_name %in% gene_name | res$gene_id %in% gene_name, ]
}

#' @rdname diversity
#' @export
dominant_isoform <- S7::new_generic(
  "dominant_isoform",
  "self",
  function(
    self,
    group = NULL,
    assay = "abundance",
    gene_name = NULL,
    min_expr = 0,
    ...
  ) {
    S7::S7_dispatch()
  }
)

S7::method(dominant_isoform, IsoformicExperiment) <- function(
  self,
  group = NULL,
  assay = "abundance",
  gene_name = NULL,
  min_expr = 0,
  ...
) {
  assay <- resolve_assay_name(self, assay)
  tx_mat <- assay(self, assay)
  annot <- row_data(self)
  gene_ids <- annot$gene_id
  if (isTRUE(min_expr > 0)) {
    # Isoforms below the floor neither dominate nor count towards the total
    tx_mat[tx_mat < min_expr] <- 0
    gene_mat <- rowsum(tx_mat, gene_ids, reorder = TRUE)
  } else {
    gene_mat <- gene_assay(self, assay)
  }
  group_means <- group_row_means(self, tx_mat, group)
  gene_means <- group_row_means(self, gene_mat, group)
  frac <- group_means / gene_means[gene_ids, , drop = FALSE]
  frac[!is.finite(frac)] <- NA_real_

  levels_g <- colnames(frac)
  n_tx <- nrow(frac)
  gene_long <- rep(gene_ids, times = length(levels_g))
  group_long <- rep(levels_g, each = n_tx)
  frac_long <- as.vector(frac)
  key <- paste(group_long, gene_long, sep = "\r")
  genes_sorted <- sort(unique(gene_ids))
  keys_sorted <- paste(
    rep(levels_g, each = length(genes_sorted)),
    rep(genes_sorted, times = length(levels_g)),
    sep = "\r"
  )
  dominant <- first_by_group(key, frac_long, keys_sorted)
  tx_idx <- ((dominant - 1L) %% n_tx) + 1L
  gene_expr <- rowsum(as.vector(group_means), key)

  canonical <- if ("is_canonical" %in% colnames(annot)) {
    annot$is_canonical
  } else {
    rep(NA, nrow(annot))
  }
  res <- tibble::tibble(
    gene_id = rep(genes_sorted, times = length(levels_g)),
    group = rep(levels_g, each = length(genes_sorted)),
    dominant_transcript_id = annot$transcript_id[tx_idx],
    dominant_transcript_name = annot$transcript_name[tx_idx],
    dominant_transcript_type = annot$transcript_type[tx_idx],
    dominant_is_canonical = canonical[tx_idx],
    dominant_fraction = unname(frac_long[dominant]),
    gene_expr = as.vector(gene_expr)[match(keys_sorted, rownames(gene_expr))]
  )
  if (!is.null(group)) {
    valid <- !is.na(res$dominant_fraction)
    n_distinct <- rowsum(
      as.integer(
        !duplicated(paste(res$gene_id, res$dominant_transcript_id)[valid])
      ),
      res$gene_id[valid]
    )
    switch_flag <- as.vector(n_distinct)[match(
      res$gene_id,
      rownames(n_distinct)
    )] >
      1L
    res$dominant_switch <- ifelse(is.na(switch_flag), FALSE, switch_flag)
    names(res)[names(res) == "group"] <- group
  } else {
    res$group <- NULL
  }
  unc <- self@dea[["uncertainty"]]
  if (!is.null(unc)) {
    pos <- match(res$dominant_transcript_id, unc$transcript_id)
    res$dominant_inf_rv <- unc$mean_inf_rv[pos]
    res$dominant_ambiguous <- ifelse(
      is.na(pos),
      NA,
      !is.na(unc$ambiguity_group[pos])
    )
  }
  res <- res |>
    dplyr::left_join(
      self@annot_data_genes[, c("gene_id", "gene_name", "gene_type")],
      by = "gene_id"
    ) |>
    dplyr::relocate("gene_name", "gene_type", .after = "gene_id")
  if (!is.null(group)) {
    res <- dplyr::arrange(
      res,
      dplyr::desc(.data$dominant_switch),
      .data$gene_id
    )
  }
  filter_genes(res, gene_name)
}

#' Row means per sample group via an indicator matrix (one BLAS call)
#' @keywords internal
#' @noRd
group_row_means <- function(self, mat, group = NULL) {
  if (is.null(group)) {
    res <- matrix(
      rowMeans(mat),
      ncol = 1L,
      dimnames = list(rownames(mat), "all")
    )
    return(res)
  }
  groups <- sample_groups(self, group, colnames(mat))
  indicator <- stats::model.matrix(~ 0 + groups)
  colnames(indicator) <- levels(groups)
  sums <- mat %*% indicator
  sweep(sums, 2L, colSums(indicator), "/")
}

#' Factor of sample groups aligned with assay columns
#'
#' With `allow_na = TRUE`, samples whose value is outside `levels` (or missing)
#' get `NA` instead of an error, so that a subset of the groups can be used.
#' @keywords internal
#' @noRd
sample_groups <- function(
  self,
  group,
  sample_ids,
  levels = NULL,
  allow_na = FALSE
) {
  check_string(group)
  col_data <- self@col_data
  if (!group %in% colnames(col_data)) {
    cli::cli_abort(
      message = c(
        x = "Column {.field {group}} not found in {.fn col_data}.",
        i = "Available columns: {.val {colnames(col_data)}}."
      ),
      class = "isoformic_group_not_found"
    )
  }
  values <- col_data[[group]][match(sample_ids, col_data$sample_id)]
  if (is.null(levels)) {
    levels <- if (is.factor(values)) {
      levels(values)
    } else {
      unique(as.character(values))
    }
  }
  groups <- factor(as.character(values), levels = levels)
  if (anyNA(groups) && !isTRUE(allow_na)) {
    cli::cli_abort(
      message = c(
        x = "Some samples have a {.field {group}} value outside {.arg levels} or missing."
      ),
      class = "isoformic_group_levels"
    )
  }
  groups
}

#' @rdname diversity
#' @export
isoform_switches <- S7::new_generic(
  "isoform_switches",
  "self",
  function(self, gene_name = NULL, confident_only = FALSE, ...) {
    S7::S7_dispatch()
  }
)

S7::method(isoform_switches, IsoformicExperiment) <- function(
  self,
  gene_name = NULL,
  confident_only = FALSE,
  ...
) {
  det <- de_tx(self)
  tx <- self@annot_data_transcripts
  idx <- match(det$transcript_id, tx$transcript_id)
  gene_id <- tx$gene_id[idx]
  tx_type <- tx$transcript_type[idx]
  sig <- det$det_sig == "yes"
  conf <- self@dea[["confidence"]]
  conf_high <- NULL
  if (!is.null(conf)) {
    conf <- conf[conf$feature_type == "transcript", ]
    conf_high <- conf$confidence[match(det$transcript_id, conf$feature_id)] %in%
      "high"
  }
  if (isTRUE(confident_only)) {
    if (is.null(conf_high)) {
      cli::cli_abort(
        message = c(
          x = "{.arg confident_only} needs the confidence of the differential expression calls.",
          i = "Compute it with {.code add_uncertainty(x, paths = , group = )} or {.fn de_confidence}."
        ),
        class = "isoformic_no_uncertainty"
      )
    }
    sig <- sig & conf_high
  }

  n_tested <- rowsum(rep(1L, length(gene_id)), gene_id)
  genes <- rownames(n_tested)
  n_sig <- rowsum(as.integer(sig), gene_id)
  n_up <- rowsum(as.integer(sig & det$log2FC > 0), gene_id)
  n_down <- rowsum(as.integer(sig & det$log2FC < 0), gene_id)
  abs_lfc <- abs(det$log2FC)
  abs_lfc[is.na(abs_lfc)] <- -Inf
  max_idx <- first_by_group(gene_id, abs_lfc, genes)
  sig_types <- type_summary_by_group(tx_type[sig], gene_id[sig], genes)

  res <- tibble::tibble(
    gene_id = genes,
    n_tx_tested = as.vector(n_tested),
    n_det_sig = as.vector(n_sig),
    n_det_up = as.vector(n_up),
    n_det_down = as.vector(n_down),
    det_sig_types = sig_types$types,
    max_abs_log2FC = abs_lfc[max_idx]
  )
  if (!is.null(conf_high)) {
    n_conf <- rowsum(as.integer(det$det_sig == "yes" & conf_high), gene_id)
    res$n_det_confident <- as.vector(n_conf)
    res <- dplyr::relocate(res, "n_det_confident", .after = "n_det_sig")
  }
  res <- dplyr::filter(res, .data$n_det_sig > 0L)

  if (!is.null(self@dea[["deg"]])) {
    deg <- de_gene(self)[, c("gene_id", "log2FC", "deg_sig")]
    names(deg) <- c("gene_id", "gene_log2FC", "deg_sig")
    res <- dplyr::left_join(res, deg, by = "gene_id")
  } else {
    res$gene_log2FC <- NA_real_
    res$deg_sig <- NA_character_
  }
  res <- res |>
    dplyr::mutate(
      switch_type = dplyr::case_when(
        .data$n_det_up > 0L & .data$n_det_down > 0L ~ "opposite_direction",
        is.na(.data$deg_sig) | .data$deg_sig != "yes" ~ "det_not_deg",
        TRUE ~ "consistent"
      )
    ) |>
    dplyr::left_join(
      self@annot_data_genes[, c("gene_id", "gene_name", "gene_type")],
      by = "gene_id"
    ) |>
    dplyr::relocate("gene_name", "gene_type", .after = "gene_id") |>
    dplyr::arrange(
      factor(
        .data$switch_type,
        levels = c("opposite_direction", "det_not_deg", "consistent")
      ),
      dplyr::desc(.data$n_det_sig),
      dplyr::desc(.data$max_abs_log2FC)
    )
  filter_genes(res, gene_name)
}

#' @rdname diversity
#' @export
exon_diversity <- S7::new_generic(
  "exon_diversity",
  "self",
  function(self, gene_name = NULL, ...) {
    S7::S7_dispatch()
  }
)

S7::method(exon_diversity, IsoformicExperiment) <- function(
  self,
  gene_name = NULL,
  ...
) {
  genes <- self@annot_data_genes[, c("gene_id", "gene_name", "gene_type")]
  tx <- self@annot_data_transcripts[, c("transcript_id", "gene_id", "n_exons")]
  if (!is.null(gene_name)) {
    gene_ids <- genes$gene_id[
      genes$gene_name %in% gene_name | genes$gene_id %in% gene_name
    ]
    tx <- tx[tx$gene_id %in% gene_ids, ]
  }
  exons <- self@annot_data_exons[, c("transcript_id", "start_pos", "end_pos")]
  if (!is.null(gene_name)) {
    exons <- exons[exons$transcript_id %in% tx$transcript_id, ]
  }
  # Aggregate over (potentially millions of) exon rows in DuckDB
  con <- duckdb_connect()
  duckdb::duckdb_register(con, "tx", as.data.frame(tx))
  duckdb::duckdb_register(con, "exons", as.data.frame(exons))
  res <- DBI::dbGetQuery(
    con,
    "WITH ex AS (
       SELECT DISTINCT t.gene_id, e.transcript_id, e.start_pos, e.end_pos
       FROM exons e JOIN tx t USING (transcript_id)
     ),
     ntx AS (
       SELECT gene_id, count(*)::INTEGER AS n_tx,
              avg(n_exons) AS mean_exons_per_tx,
              max(n_exons)::INTEGER AS max_exons_per_tx
       FROM tx GROUP BY gene_id
     ),
     ek AS (
       SELECT gene_id, start_pos, end_pos, count(*)::INTEGER AS n_tx_with
       FROM ex GROUP BY gene_id, start_pos, end_pos
     ),
     agg AS (
       SELECT ek.gene_id,
              count(*)::INTEGER AS n_exons,
              sum(CASE WHEN ek.n_tx_with = ntx.n_tx THEN 1 ELSE 0 END)::INTEGER AS n_constitutive,
              sum(CASE WHEN ek.n_tx_with < ntx.n_tx THEN 1 ELSE 0 END)::INTEGER AS n_alternative
       FROM ek JOIN ntx USING (gene_id)
       GROUP BY ek.gene_id
     )
     SELECT ntx.gene_id, ntx.n_tx, ntx.mean_exons_per_tx, ntx.max_exons_per_tx,
            coalesce(agg.n_exons, 0) AS n_exons,
            coalesce(agg.n_constitutive, 0) AS n_constitutive,
            coalesce(agg.n_alternative, 0) AS n_alternative
     FROM ntx LEFT JOIN agg USING (gene_id)"
  )
  tibble::as_tibble(res) |>
    dplyr::left_join(genes, by = "gene_id") |>
    dplyr::relocate("gene_name", "gene_type", .after = "gene_id") |>
    dplyr::arrange(dplyr::desc(.data$n_alternative), dplyr::desc(.data$n_exons))
}

#' @rdname diversity
#' @export
haplotype_summary <- S7::new_generic(
  "haplotype_summary",
  "self",
  function(
    self,
    group = NULL,
    assay = "abundance",
    gene_name = NULL,
    ...
  ) {
    S7::S7_dispatch()
  }
)

S7::method(haplotype_summary, IsoformicExperiment) <- function(
  self,
  group = NULL,
  assay = "abundance",
  gene_name = NULL,
  ...
) {
  assay <- resolve_assay_name(self, assay)
  annot <- row_data(self)
  haps <- annotation_haplotypes(annot)
  if (length(haps) < 2L) {
    cli::cli_abort(
      message = c(
        x = "The annotation is not haplotype-resolved (fewer than two haplotype labels).",
        i = "Build one with {.fn merge_annotations}."
      ),
      class = "isoformic_not_haplotype_resolved"
    )
  }
  tx_mat <- assay(self, assay)
  group_means <- group_row_means(self, tx_mat, group)
  levels_g <- colnames(group_means)
  gene_key <- dplyr::coalesce(
    annot$source_gene_id,
    annot$gene_name,
    annot$gene_id
  )
  hap <- annot$haplotype
  keep <- !is.na(hap)
  tx_key <- dplyr::coalesce(
    annot$source_transcript_id,
    annot$transcript_name,
    annot$transcript_id
  )

  rows <- lapply(levels_g, function(lv) {
    expr <- group_means[keep, lv]
    key <- paste(gene_key[keep], hap[keep], sep = "\r")
    hap_expr <- rowsum(expr, key)
    gene_total <- rowsum(expr, gene_key[keep])
    parts <- strsplit(rownames(hap_expr), "\r", fixed = TRUE)
    gk <- vapply(parts, `[[`, character(1L), 1L)
    hp <- vapply(parts, `[[`, character(1L), 2L)
    frac_within <- expr / as.vector(hap_expr)[match(key, rownames(hap_expr))]
    frac_within[!is.finite(frac_within)] <- NA_real_
    dom <- first_by_group(key, frac_within, rownames(hap_expr))
    tibble::tibble(
      gene_key = gk,
      haplotype = hp,
      group = lv,
      mean_expr = as.vector(hap_expr),
      haplotype_fraction = as.vector(hap_expr) /
        as.vector(gene_total)[match(gk, rownames(gene_total))],
      dominant_transcript_id = annot$transcript_id[keep][dom],
      dominant_source_transcript_id = tx_key[keep][dom],
      dominant_transcript_type = annot$transcript_type[keep][dom],
      dominant_fraction = unname(frac_within[dom])
    )
  })
  res <- dplyr::bind_rows(rows)
  res$haplotype_fraction[!is.finite(res$haplotype_fraction)] <- NA_real_
  switches <- res |>
    dplyr::filter(!is.na(.data$dominant_fraction)) |>
    dplyr::summarise(
      haplotype_isoform_switch = dplyr::n_distinct(
        .data$dominant_source_transcript_id
      ) >
        1L,
      .by = c("gene_key", "group")
    )
  res <- dplyr::left_join(res, switches, by = c("gene_key", "group"))
  res$haplotype_isoform_switch[is.na(res$haplotype_isoform_switch)] <- FALSE
  names_lookup <- dplyr::distinct(
    tibble::tibble(gene_key = gene_key, gene_name = annot$gene_name),
    .data$gene_key,
    .keep_all = TRUE
  )
  res <- dplyr::left_join(res, names_lookup, by = "gene_key") |>
    dplyr::relocate("gene_name", .after = "gene_key")
  if (is.null(group)) {
    res$group <- NULL
  } else {
    names(res)[names(res) == "group"] <- group
  }
  res <- dplyr::arrange(
    res,
    dplyr::desc(.data$haplotype_isoform_switch),
    .data$gene_key,
    .data$haplotype
  )
  if (!is.null(gene_name)) {
    res <- res[res$gene_name %in% gene_name | res$gene_key %in% gene_name, ]
  }
  res
}

#' Samples of the requested groups, for a test between groups
#'
#' Returns the logical vector of the columns to keep and their group factor.
#' Samples outside `levels` are dropped, so a contrast between some of the
#' groups is possible. `n_levels` is the number of levels required: `2` for a
#' two-group test, `NULL` for two or more.
#' @keywords internal
#' @noRd
group_setup <- function(
  self,
  group,
  sample_ids,
  levels = NULL,
  n_levels = NULL
) {
  groups <- sample_groups(
    self,
    group,
    sample_ids,
    levels = levels,
    allow_na = TRUE
  )
  observed <- unique(as.character(groups[!is.na(groups)]))
  groups <- factor(
    as.character(groups),
    levels = intersect(levels(groups), observed)
  )
  ok <- if (is.null(n_levels)) {
    nlevels(groups) >= 2L
  } else {
    nlevels(groups) == n_levels
  }
  if (!ok) {
    wanted <- if (is.null(n_levels)) {
      "at least two groups"
    } else {
      "exactly two groups"
    }
    cli::cli_abort(
      message = c(
        x = "The test needs {wanted}; {.field {group}} has {.val {levels(groups)}}.",
        i = "Choose them with {.arg levels}, reference first."
      ),
      class = "isoformic_group_levels"
    )
  }
  keep <- !is.na(groups)
  list(keep = keep, groups = droplevels(groups[keep]))
}

#' Samples of exactly two groups, for a two-group test
#' @keywords internal
#' @noRd
two_group_setup <- function(self, group, sample_ids, levels = NULL) {
  group_setup(self, group, sample_ids, levels = levels, n_levels = 2L)
}
