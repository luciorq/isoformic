# lib_null_methods.R: the null-split methods M1-M11 on plain matrices
# (sourced, not run).
#
# Purpose: one implementation of the pre-registered methods (prereg.yaml,
# `methods`) for the per-split loops of the null analyses (placenta WP-A,
# mouse WP-C1), on transcript x sample matrices, with no IsoformicExperiment
# subsetting. Every method calls the frozen package code where the package
# has it (filter_usage(), usage_diffsplice(), usage_saturn(),
# usage_dexseq_fit(), welch_matrix(), inferential_rv(),
# overdispersion_from_moments()) or edgeR / fishpond directly:
#   M1  edgeR QL, transcript counts           nm_edger_dte(counts, groups)
#   M2  edgeR QL, divided counts              nm_edger_dte(counts, groups, od)
#   M6  fishpond swish on the draws           nm_swish(arr, counts, lengths, groups, seed)
#   M7  M1 calls kept when class == high      nm_confidence() + M1
#   M8  diffSplice usage, raw counts          nm_usage(counts, gene_of, groups)
#   M9  diffSplice usage, divided counts      nm_usage(..., od = od)
#   M10 satuRn / DEXSeq usage                 nm_usage(..., backend = )
#   M11 edgeR QL on gene sums                 nm_gene_dte(counts, gene_of, groups)
# M3-M5 (limma-voom, DESeq2) are not used by WP-C and are not defined here.
#
# Filters: DTE methods test the universe rows kept by a design-free
# edgeR::filterByExpr() (WP-A rule, see below); usage methods use the package
# filter_usage() on all transcripts (min_count 10, min_prop 0.1, as
# test_usage() and 22_null_dtu.R); swish runs on the rows M1 tests, then
# fishpond::labelKeep() defaults (as 21_null_dte.R).
#
# Simes: adjust_usage() of the package loops over genes (quadratic in the
# number of genes); nm_simes() is the same formula vectorised and
# 41_mouse_null.R checks it against adjust_usage() on one split.
#
# Inputs: none (functions only). Sources null_common.R (WP-A); requires the
# snapshot loaded (ue_setup()).
# Env vars: none.

# DTE ------------------------------------------------------------------------------
# The WP-A implementations (null_common.R: edger_ql(), swish_test(),
# two_group_design()) are sourced, so the mouse null runs the same code as
# the placenta null. As in 21_null_dte.R: the tested rows are the universe
# (`rows`, U) kept by edgeR::filterByExpr() without design or group (the
# filter does not see the labels), with the library sizes of all
# transcripts; divided counts (M2) are filtered the same way on the divided
# counts; BH over the tested rows.
source("data-raw/uncertainty_evidence/null_common.R")

nm_design_free_keep <- function(counts, lib) {
  suppressMessages(edgeR::filterByExpr(
    edgeR::DGEList(counts, lib.size = lib),
    group = NULL
  ))
}

nm_in_b <- function(groups) {
  groups <- factor(groups)
  groups == levels(groups)[2L]
}

nm_edger_dte <- function(counts, groups, od = NULL, rows = NULL) {
  x <- if (is.null(od)) counts else counts / od[rownames(counts)]
  lib <- colSums(x)
  if (!is.null(rows)) {
    x <- x[rownames(x) %in% rows, , drop = FALSE]
  }
  keep <- nm_design_free_keep(x, lib)
  r <- edger_ql(x[keep, , drop = FALSE], nm_in_b(groups))
  data.frame(
    feature_id = r$id,
    logFC = r$log2FC,
    pvalue = r$pvalue,
    qvalue = stats::p.adjust(r$pvalue, "BH")
  )
}

# Gene level (prereg M11 = WP-A M11 raw): gene sums of all transcripts
nm_gene_dte <- function(counts, gene_of, groups) {
  nm_edger_dte(rowsum(counts, gene_of, reorder = TRUE), groups)
}

# Usage ----------------------------------------------------------------------------

# Gene-level Simes combination (simes_pvalue() of the package, vectorised)
nm_simes <- function(p, gene) {
  ok <- !is.na(p)
  p <- p[ok]
  gene <- gene[ok]
  o <- order(gene, p)
  p <- p[o]
  gene <- gene[o]
  n <- stats::ave(p, gene, FUN = length)
  r <- stats::ave(p, gene, FUN = seq_along)
  v <- n * p / r
  m <- tapply(v, gene, min)
  data.frame(gene_id = names(m), gene_pvalue = as.vector(m))
}

nm_usage <- function(
  counts,
  gene_of,
  groups,
  od = NULL,
  backend = c("diffsplice", "saturn", "dexseq"),
  min_count = 10,
  min_prop = 0.1,
  workers = 1L
) {
  backend <- match.arg(backend)
  if (!is.null(od)) {
    counts <- counts / od[rownames(counts)]
  }
  groups <- factor(groups)
  kept <- filter_usage(counts, gene_of, groups, min_count, min_prop)
  counts <- counts[kept, , drop = FALSE]
  g <- gene_of[kept]
  contrasts <- resolve_contrasts(NULL, levels(groups))
  res <- switch(
    backend,
    diffsplice = usage_diffsplice(counts, g, groups, contrasts),
    saturn = usage_saturn(counts, g, groups, contrasts, workers),
    dexseq = usage_dexseq_fit(
      counts,
      g,
      factor(as.character(groups), levels = rev(levels(groups))),
      workers,
      fold_change = FALSE
    )
  )
  res <- as.data.frame(res)
  res$gene_id <- g[match(res$transcript_id, rownames(counts))]
  res$qvalue <- stats::p.adjust(res$pvalue, "BH")
  genes <- nm_simes(res$pvalue, res$gene_id)
  genes$gene_qvalue <- stats::p.adjust(genes$gene_pvalue, "BH")
  list(transcripts = res, genes = genes)
}

# Swish ----------------------------------------------------------------------------

# WP-A swish_test(): arr transcript x sample x draw, counts and effective
# lengths of the same rows and columns (scaleInfReps with length correction,
# labelKeep, swish with 100 permutations, seed per split)
nm_swish <- function(arr, counts, lengths, groups, seed, nperms = 100L) {
  r <- swish_test(counts, arr, lengths, nm_in_b(groups), seed, nperms)
  data.frame(
    feature_id = r$id,
    logFC = r$log2FC,
    pvalue = r$pvalue,
    qvalue = r$swish_qvalue
  )
}

# Confidence class (compute_de_confidence(), transcripts, no collapse) -------------

# arr: transcript x sample x draw of the samples `s` (columns of arr);
# in_b: logical, second group among `s`; lib: library sizes of `s` (counts
# column sums); pair_cor: data frame with transcript_id_1, transcript_id_2
# and one column of per-sample draw correlation per sample name (NA where
# masked), as written by 40_mouse_prepare.R. A transcript is ambiguous when a
# sibling pair, with both members assessed (mean count >= 1 over `s`), has a
# mean correlation over `s` (non-NA samples) <= -0.5.
nm_confidence <- function(
  arr,
  s,
  in_b,
  lib,
  ids,
  pair_cor,
  min_count = 1,
  inf_rv_cutoff = 0.5,
  cor_cutoff = -0.5,
  min_agreement = 0.95,
  pc = 5
) {
  rn <- dimnames(arr)[[1L]]
  sub <- arr[, s, , drop = FALSE]
  n <- dim(sub)[3L]
  mu <- rowMeans(sub, dims = 2L)
  v <- pmax((rowSums(sub^2, dims = 2L) - n * mu^2) / (n - 1L), 0)
  mean_count <- rowMeans(mu)
  rv <- rowMeans(inferential_rv(list(mean = mu, var = v)))
  assessed <- mean_count >= min_count
  names(assessed) <- rn
  sn <- dimnames(arr)[[2L]][s]
  pr <- as.matrix(pair_cor[, sn, drop = FALSE])
  okp <- assessed[pair_cor$transcript_id_1] %in%
    TRUE &
    assessed[pair_cor$transcript_id_2] %in% TRUE &
    rowSums(!is.na(pr)) > 0
  mc <- rowMeans(pr[okp, , drop = FALSE], na.rm = TRUE)
  amb <- okp
  amb[okp] <- mc <= cor_cutoff
  amb_ids <- unique(c(
    pair_cor$transcript_id_1[amb],
    pair_cor$transcript_id_2[amb]
  ))
  idx <- match(ids, rn)
  scale <- mean(lib) / lib
  lfc <- matrix(NA_real_, length(idx), n)
  for (k in seq_len(n)) {
    x <- log2(
      sweep(matrix(sub[idx, , k], nrow = length(idx)), 2L, scale, "*") + pc
    )
    lfc[, k] <- welch_matrix(x, in_b)$lfc
  }
  med <- apply(lfc, 1L, stats::median)
  agree <- rowMeans(sign(lfc) == sign(med))
  stable <- agree >= min_agreement
  uncertain <- rv[idx] >= inf_rv_cutoff
  ambiguous <- ids %in% amb_ids
  n_fail <- as.integer(!stable) + as.integer(uncertain) + as.integer(ambiguous)
  cls <- c("high", "moderate", "low")[pmin(n_fail, 2L) + 1L]
  cls[!assessed[idx]] <- NA_character_
  data.frame(
    transcript_id = ids,
    mean_count = mean_count[idx],
    mean_inf_rv = rv[idx],
    sign_agreement = agree,
    stable = stable,
    uncertain = uncertain,
    ambiguous = ambiguous,
    confidence = cls
  )
}

# Jackknife with Monte Carlo split sampling -----------------------------------------

# Delete-one-unit jackknife of a statistic of a set of splits. `full` holds
# the splits of the complete design; `units` a list (one element per deleted
# unit) of the splits sampled for that unit (m of a pool of `pool` distinct
# splits). `stat` maps a set of splits to a number. The Monte Carlo variance
# of each unit estimate (inner bootstrap over its m splits, times the finite
# population correction 1 - m / pool) is removed from the jackknife sum.
# `log` puts the interval on the log scale (ratios).
# B = 0 skips the correction (the interval is the uncorrected one). Units
# that hold all `pool` splits have no Monte Carlo variance and are not
# bootstrapped (exact jackknife). The verification of WP-C showed that the
# inner bootstrap removes only a part of the Monte Carlo variance of
# heavy-tailed, nonlinear statistics (adjusted rate ratios: 31% of the sum
# of squares removed while 28-split unit estimates were 4 times as dispersed
# as the exact 126-split ones), so 41_mouse_null.R now uses exact units
# where the fits exist and B = 0 (uncorrected, conservative) elsewhere.
nm_jackknife <- function(
  full,
  units,
  stat,
  pool,
  B = 200L,
  log = FALSE,
  seed = 1L
) {
  tr <- if (log) base::log else identity
  inv <- if (log) base::exp else identity
  est <- stat(full)
  th <- vapply(units, function(u) tr(stat(u)), numeric(1))
  set.seed(seed)
  mc_var <- vapply(
    units,
    function(u) {
      # units are data frames of split keys (rows) or vectors
      m <- if (is.data.frame(u)) nrow(u) else length(u)
      if (B <= 0L || m >= pool) {
        return(0)
      }
      pick <- function(i) if (is.data.frame(u)) u[i, , drop = FALSE] else u[i]
      b <- replicate(B, tr(stat(pick(sample.int(m, m, replace = TRUE)))))
      stats::var(b[is.finite(b)]) * max(0, 1 - m / pool)
    },
    numeric(1)
  )
  n <- length(units)
  if (!all(is.finite(th)) || !is.finite(tr(est))) {
    # A unit estimate is 0 or undefined on the log scale (for example no
    # p < 0.01 call at all): no interval
    return(data.frame(
      estimate = est,
      se = NA_real_,
      se_uncorrected = NA_real_,
      ci_lo = NA_real_,
      ci_hi = NA_real_,
      ci_lo_uncorrected = NA_real_,
      ci_hi_uncorrected = NA_real_,
      n_units = n,
      mc_share = NA_real_,
      mc_corrected = NA
    ))
  }
  mc_var[!is.finite(mc_var)] <- 0
  ss <- sum((th - mean(th))^2)
  ss_corr <- ss - (1 - 1 / n) * sum(mc_var)
  se_raw <- sqrt((n - 1) / n * ss)
  # When the estimated Monte Carlo share reaches the whole jackknife sum of
  # squares the corrected variance is not estimable (it would be 0 or
  # negative); the uncorrected (conservative) interval is reported instead.
  corrected <- ss_corr > 0 && sum(mc_var) > 0
  se <- if (corrected) sqrt((n - 1) / n * ss_corr) else se_raw
  data.frame(
    estimate = est,
    se = se,
    se_uncorrected = se_raw,
    ci_lo = inv(tr(est) - 1.96 * se),
    ci_hi = inv(tr(est) + 1.96 * se),
    ci_lo_uncorrected = inv(tr(est) - 1.96 * se_raw),
    ci_hi_uncorrected = inv(tr(est) + 1.96 * se_raw),
    n_units = n,
    mc_share = if (ss > 0) (1 - 1 / n) * sum(mc_var) / ss else NA_real_,
    mc_corrected = corrected
  )
}

# Confidence class in two steps, for loops over many splits of the same
# samples: nm_conf_design() computes what does not depend on the labels
# (InfRV, assessed, ambiguity, the log2 scaled draws of `ids`), and
# nm_conf_split() the draw stability of one labelling. The result equals
# nm_confidence() (checked in 41_mouse_null.R); welch_matrix()$lfc is the
# difference of the group means of the log2 values.
# `mean`, `var` (optional): per-sample moments of the draws of every
# transcript (rows named, columns = the samples `s`); when given, InfRV and
# the `assessed` rule come from them and `arr` needs only the rows `ids`
# (the moments are those of the same draws, so the result is identical).
nm_conf_design <- function(
  arr,
  s,
  lib,
  ids,
  pair_cor,
  min_count = 1,
  inf_rv_cutoff = 0.5,
  cor_cutoff = -0.5,
  pc = 5,
  mean = NULL,
  var = NULL
) {
  n <- dim(arr)[3L]
  if (!is.null(mean)) {
    rn <- rownames(mean)
    mu <- mean
    v <- var
  } else {
    rn <- dimnames(arr)[[1L]]
    mu <- matrix(0, length(rn), length(s))
    s2 <- mu
    for (k in seq_len(n)) {
      x <- arr[, s, k]
      mu <- mu + x
      s2 <- s2 + x * x
    }
    mu <- mu / n
    v <- pmax((s2 - n * mu^2) / (n - 1L), 0)
    rm(s2)
  }
  mean_count <- rowMeans(mu)
  rv <- rowMeans(inferential_rv(list(mean = mu, var = v)))
  assessed <- stats::setNames(mean_count >= min_count, rn)
  sn <- dimnames(arr)[[2L]][s]
  pr <- as.matrix(pair_cor[, sn, drop = FALSE])
  okp <- assessed[pair_cor$transcript_id_1] %in%
    TRUE &
    assessed[pair_cor$transcript_id_2] %in% TRUE &
    rowSums(!is.na(pr)) > 0
  amb <- okp
  amb[okp] <- rowMeans(pr[okp, , drop = FALSE], na.rm = TRUE) <= cor_cutoff
  amb_ids <- unique(c(
    pair_cor$transcript_id_1[amb],
    pair_cor$transcript_id_2[amb]
  ))
  idx <- match(ids, rn)
  ia <- match(ids, dimnames(arr)[[1L]])
  stopifnot(!anyNA(idx), !anyNA(ia))
  scale <- base::mean(lib) / lib
  logd <- array(0, c(length(ia), length(s), n))
  for (k in seq_len(n)) {
    logd[,, k] <- log2(
      sweep(matrix(arr[ia, s, k], nrow = length(ia)), 2L, scale, "*") + pc
    )
  }
  list(
    ids = ids,
    mean_count = mean_count[idx],
    mean_inf_rv = rv[idx],
    uncertain = rv[idx] >= inf_rv_cutoff,
    ambiguous = ids %in% amb_ids,
    assessed = assessed[idx],
    logd = logd,
    samples = sn
  )
}

nm_conf_split <- function(cd, in_b, min_agreement = 0.95) {
  n <- dim(cd$logd)[3L]
  lfc <- matrix(NA_real_, length(cd$ids), n)
  for (k in seq_len(n)) {
    x <- cd$logd[,, k]
    lfc[, k] <- rowMeans(x[, in_b, drop = FALSE]) -
      rowMeans(x[, !in_b, drop = FALSE])
  }
  med <- apply(lfc, 1L, stats::median)
  agree <- rowMeans(sign(lfc) == sign(med))
  stable <- agree >= min_agreement
  n_fail <- as.integer(!stable) +
    as.integer(cd$uncertain) +
    as.integer(cd$ambiguous)
  cls <- c("high", "moderate", "low")[pmin(n_fail, 2L) + 1L]
  cls[!cd$assessed] <- NA_character_
  data.frame(
    transcript_id = cd$ids,
    sign_agreement = agree,
    stable = stable,
    uncertain = cd$uncertain,
    ambiguous = cd$ambiguous,
    confidence = cls
  )
}
