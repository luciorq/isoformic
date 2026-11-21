# 49_replication_summarise.R: the WP-C headline tables (mouse null, read
# halves) and their comparison with the placenta results.
#
# Purpose: collect the numbers WP-C must report, each with its 95% CI and
# its source table, into one long table, and check whether the placenta
# pattern replicates: same sign and overlapping 95% CIs for the matching
# statistic. Nothing is re-estimated here; every value is copied from the
# Parquet tables of 41_mouse_null.R and 44_rhalf_replicate.R (and of WP-A /
# WP-B for the placenta, when they exist).
#
# Inputs:
#   $ISOFORMIC_UE_OUT/null_mouse/results_mouse_null.parquet, m10_splits.parquet,
#     splits.parquet, prepare/validation*.parquet       (40, 41)
#   $ISOFORMIC_UE_OUT/rhalf/replication/summary_*.parquet, split_check.parquet,
#     quant_times.parquet, summary_meta.parquet          (42, 43, 44)
#   Placenta (optional; comparisons are skipped when absent):
#     PLACENTA_NULL_SUMMARY  WP-A summary directory (29_null_summarise.R),
#       default $ISOFORMIC_UE_OUT/null_placenta/summary
#     PLACENTA_REPLICATION_FLAGS  default
#       $ISOFORMIC_UE_OUT/placenta_audit/replication/summary_flags.parquet (WP-B)
# Env vars: ISOFORMIC_UE_OUT, ISOFORMIC_UE_PKG (ue_common.R) and the two above.
# Outputs: $ISOFORMIC_UE_OUT/wpc_summary/wpc_headlines.parquet (+ .tsv),
#   wpc_replication_check.parquet (+ .tsv), wpc_rule_RA.parquet (+ .tsv; the
#   pre-registered rule R_A on the mouse null, all and batch-balanced
#   splits), session_info.parquet
#
# Resources (estimate): seconds, < 1 GB.
# Run from the repository root, inside a memory cap of about 4 GB (see
# README.md):
#   pixi run --frozen Rscript data-raw/uncertainty_evidence/49_replication_summarise.R

source("data-raw/uncertainty_evidence/ue_common.R")
ue <- ue_setup("49_replication_summarise")
root <- ue_out_root()
out <- file.path(root, "wpc_summary")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
rp <- function(...) {
  f <- file.path(root, ...)
  if (!file.exists(f)) {
    return(NULL)
  }
  x <- as.data.frame(read_parquet(f))
  x[,
    setdiff(colnames(x), c("git_sha", "snapshot_sha256", "seed", "script")),
    drop = FALSE
  ]
}
parse_ci <- function(x) {
  v <- suppressWarnings(as.numeric(strsplit(gsub("\\[|\\]", "", x), ", ")[[1]]))
  if (length(v) != 2L) c(NA_real_, NA_real_) else v
}
rows <- list()
add <- function(
  arm,
  what,
  estimate,
  ci_lo = NA_real_,
  ci_hi = NA_real_,
  source,
  note = ""
) {
  # jarl-ignore undesirable_operator: appends to a result list of the enclosing scope
  rows[[length(rows) + 1L]] <<- data.frame(
    arm = arm,
    what = what,
    estimate = estimate,
    ci_lo = ci_lo,
    ci_hi = ci_hi,
    source = source,
    note = note
  )
}

# Mouse null ------------------------------------------------------------------------------
mres <- rp("null_mouse", "results_mouse_null.parquet")
msrc <- "null_mouse/results_mouse_null.parquet"
if (!is.null(mres)) {
  for (i in seq_len(nrow(mres))) {
    r <- mres[i, ]
    add(
      "mouse_null",
      paste(r$method, r$statistic),
      r$estimate,
      r$ci_lo,
      r$ci_hi,
      sprintf("%s (method %s, statistic %s)", msrc, r$method, r$statistic),
      paste(r$note, r$jackknife, sep = "; ")
    )
  }
}
m10 <- rp("null_mouse", "m10_splits.parquet")
if (!is.null(m10)) {
  for (col in c("saturn_q05_U", "dexseq_q05_U", "M9_q05_U", "M8_q05_U")) {
    v <- m10[[col]]
    # t interval over the 10 splits (splits share animals: indicative only)
    se <- stats::sd(v) / sqrt(length(v))
    add(
      "mouse_null",
      paste("M10 splits: false DTU genes per split,", col),
      mean(v),
      mean(v) - stats::qt(0.975, length(v) - 1) * se,
      mean(v) + stats::qt(0.975, length(v) - 1) * se,
      sprintf("null_mouse/m10_splits.parquet column %s", col),
      "10 random 5 vs 5 splits; t interval treats splits as independent (they share animals), indicative only"
    )
  }
  add(
    "mouse_null",
    "M10 DEXSeq minutes per split (median)",
    stats::median(m10$dexseq_min),
    NA,
    NA,
    "null_mouse/m10_splits.parquet column dexseq_min"
  )
}
sp <- rp("null_mouse", "splits.parquet")
if (!is.null(sp)) {
  f <- sp[sp$design == "full", ]
  add(
    "mouse_null",
    "full splits: max |log2 depth difference| between groups",
    max(abs(f$log2_depth_diff)),
    NA,
    NA,
    "null_mouse/splits.parquet column log2_depth_diff (design full)"
  )
  add(
    "mouse_null",
    "full splits: max |mapping-rate difference| (percentage points)",
    max(abs(f$maprate_diff)),
    NA,
    NA,
    "null_mouse/splits.parquet column maprate_diff (design full)"
  )
}

# Read halves -------------------------------------------------------------------------------
rdir <- c("rhalf", "replication")
fl <- rp(rdir[1], rdir[2], "summary_flags.parquet")
if (!is.null(fl)) {
  for (i in seq_len(nrow(fl))) {
    r <- fl[i, ]
    key <- sprintf("%s %s %s %s", r$halving, r$series, r$method, r$flag)
    src <- sprintf(
      "rhalf/replication/summary_flags.parquet (halving %s, series %s, method %s, flag %s)",
      r$halving,
      r$series,
      r$method,
      r$flag
    )
    ci <- parse_ci(r$rate_flagged_ci)
    add(
      "read_halves",
      paste(key, "replication of flagged calls"),
      r$rate_flagged,
      ci[1],
      ci[2],
      paste(src, "rate_flagged")
    )
    ci <- parse_ci(r$rate_unflagged_ci)
    add(
      "read_halves",
      paste(key, "replication of unflagged calls"),
      r$rate_unflagged,
      ci[1],
      ci[2],
      paste(src, "rate_unflagged")
    )
    ci <- parse_ci(r$deficit_ci)
    add(
      "read_halves",
      paste(key, "deficit (unflagged - flagged)"),
      r$deficit,
      ci[1],
      ci[2],
      paste(src, "deficit")
    )
    ci <- parse_ci(r$adjusted_or_ci)
    add(
      "read_halves",
      paste(key, "adjusted OR of replication"),
      r$adjusted_or,
      ci[1],
      ci[2],
      paste(src, "adjusted_or"),
      "gene-clustered bootstrap"
    )
    ci <- parse_ci(r$delta_auc_ci)
    add(
      "read_halves",
      paste(key, "delta-AUC within series"),
      r$delta_auc,
      ci[1],
      ci[2],
      paste(src, "delta_auc")
    )
    add(
      "read_halves",
      paste(key, "percent of calls flagged"),
      r$pct_flagged,
      NA,
      NA,
      paste(src, "pct_flagged"),
      sprintf("%d calls, %d genes", r$n_calls, r$n_genes)
    )
    if (!is.null(r$n_untestable_v)) {
      add(
        "read_halves",
        paste(
          key,
          "calls not tested in the validation half (its filterByExpr)"
        ),
        r$n_untestable_v,
        NA,
        NA,
        paste(src, "n_untestable_v"),
        sprintf(
          "%d of them flagged; they count as not replicated in the rates above",
          r$n_untestable_v_flagged
        )
      )
      ci <- parse_ci(r$rate_flagged_testable_ci)
      add(
        "read_halves",
        paste(
          key,
          "replication of flagged calls tested in the validation half"
        ),
        r$rate_flagged_testable,
        ci[1],
        ci[2],
        paste(src, "rate_flagged_testable")
      )
      ci <- parse_ci(r$rate_unflagged_testable_ci)
      add(
        "read_halves",
        paste(
          key,
          "replication of unflagged calls tested in the validation half"
        ),
        r$rate_unflagged_testable,
        ci[1],
        ci[2],
        paste(src, "rate_unflagged_testable")
      )
      ci <- parse_ci(r$deficit_testable_ci)
      add(
        "read_halves",
        paste(key, "deficit among calls tested in the validation half"),
        r$deficit_testable,
        ci[1],
        ci[2],
        paste(src, "deficit_testable")
      )
    }
  }
}
c4 <- rp(rdir[1], rdir[2], "summary_c4.parquet")
if (!is.null(c4)) {
  for (i in seq_len(nrow(c4))) {
    r <- c4[i, ]
    ci <- parse_ci(r$ci)
    add(
      "read_halves",
      sprintf(
        "%s %s %s deficit(primary) - deficit(count-matched C4)",
        r$halving,
        r$series,
        r$method
      ),
      r$estimate,
      ci[1],
      ci[2],
      "rhalf/replication/summary_c4.parquet"
    )
    if (!is.null(r$n_unflagged_both)) {
      add(
        "read_halves",
        sprintf(
          "%s %s %s unflagged calls shared by the primary flag and C4",
          r$halving,
          r$series,
          r$method
        ),
        r$n_unflagged_both,
        NA,
        NA,
        "rhalf/replication/summary_c4.parquet column n_unflagged_both",
        sprintf(
          "unflagged: primary %d (%d replicated), C4 %d (%d replicated); share of calls with the same flag %.3f",
          r$n_unflagged_primary,
          r$n_unflagged_primary_replicated,
          r$n_unflagged_c4,
          r$n_unflagged_c4_replicated,
          r$share_calls_same_flag
        )
      )
    }
  }
}
pl <- rp(rdir[1], rdir[2], "summary_placebo.parquet")
if (!is.null(pl)) {
  for (i in seq_len(nrow(pl))) {
    r <- pl[i, ]
    add(
      "read_halves",
      sprintf(
        "%s %s %s primary flag %s vs placebo within count x |z| cells",
        r$halving,
        r$series,
        r$method,
        r$statistic
      ),
      r$observed,
      r$placebo_2.5,
      r$placebo_97.5,
      "rhalf/replication/summary_placebo.parquet",
      sprintf(
        "CI columns: 2.5%% and 97.5%% of 200 placebo flags (mean %.3f); one-sided empirical p = %.3f",
        r$placebo_mean,
        r$empirical_p_one_sided
      )
    )
  }
}
tr <- rp(rdir[1], rdir[2], "summary_transfer.parquet")
if (!is.null(tr)) {
  for (i in seq_len(nrow(tr))) {
    r <- tr[i, ]
    ci <- parse_ci(r$ci)
    add(
      "read_halves",
      sprintf("delta-AUC transferred %s -> %s", r$train, r$test),
      r$delta_auc_transfer,
      ci[1],
      ci[2],
      "rhalf/replication/summary_transfer.parquet",
      "models fitted on the first series, evaluated on the second; gene bootstrap of the second"
    )
  }
}
pm <- rp(rdir[1], rdir[2], "summary_perm.parquet")
if (!is.null(pm)) {
  for (k in unique(paste(pm$halving, pm$series))) {
    x <- pm[paste(pm$halving, pm$series) == k, ]
    tru <- x[x$true_partition, ]
    oth <- x[!x$true_partition, ]
    add(
      "read_halves",
      paste(k, "replication, true labels (mean of 2 directions)"),
      mean(tru$replication),
      NA,
      NA,
      "rhalf/replication/summary_perm.parquet (true_partition)"
    )
    add(
      "read_halves",
      paste(
        k,
        "replication, 9 other 3 vs 3 partitions (mean; range in CI columns)"
      ),
      mean(oth$replication),
      min(oth$replication),
      max(oth$replication),
      "rhalf/replication/summary_perm.parquet (not true_partition)",
      "CI columns hold the min and max over the 18 permuted directed analyses, not a CI"
    )
    add(
      "read_halves",
      paste(k, "deficit, 9 other partitions (mean; range in CI columns)"),
      mean(oth$deficit, na.rm = TRUE),
      min(oth$deficit, na.rm = TRUE),
      max(oth$deficit, na.rm = TRUE),
      "rhalf/replication/summary_perm.parquet",
      "range, not a CI"
    )
  }
}
sc <- rp(rdir[1], rdir[2], "split_check.parquet")
if (!is.null(sc)) {
  add(
    "read_halves",
    "split check: runs with names agreeing in both halves and both mates",
    sum(sc$names_agree_A & sc$names_agree_B),
    NA,
    NA,
    "rhalf/split_check.tsv",
    sprintf("of %d runs", nrow(sc))
  )
  add(
    "read_halves",
    "split check: read pairs in both halves (all runs)",
    sum(sc$spots_in_both_halves),
    NA,
    NA,
    "rhalf/split_check.tsv"
  )
  add(
    "read_halves",
    "split check: runs whose halves cover the original exactly",
    sum(sc$halves_cover_original),
    NA,
    NA,
    "rhalf/split_check.tsv"
  )
  add(
    "read_halves",
    "split check: pairs in half A (per run; identical seed)",
    stats::median(sc$n_A_1),
    min(sc$n_A_1),
    max(sc$n_A_1),
    "rhalf/split_check.tsv column n_A_1",
    "CI columns: min and max over runs"
  )
}
qt <- rp(rdir[1], rdir[2], "quant_times.parquet")
if (!is.null(qt)) {
  ok <- qt[as.character(qt$exit_status) == "0", ]
  add(
    "read_halves",
    "Salmon wall time per half (s, median of verified runs)",
    stats::median(ok$wall_s),
    min(ok$wall_s),
    max(ok$wall_s),
    "rhalf/quant_times.tsv (exit_status 0)",
    sprintf(
      "%d verified runs; %d failed attempts (memory watchdog); CI columns: min and max",
      nrow(ok),
      sum(as.character(qt$exit_status) != "0")
    )
  )
  add(
    "read_halves",
    "Salmon peak RSS per half (GB, max)",
    max(ok$max_rss_kb) / 1024^2,
    NA,
    NA,
    "rhalf/quant_times.tsv"
  )
}

# Pre-registered rule R_A on the mouse null (prereg.yaml decision_rules):
# (1) the M1 adjusted rate ratio of the top tertile has a 95% CI above 1,
# (2) M2 or M6 reduces it (log RR difference with a CI above 0), (3) M11
# shows no comparable excess (its adjusted RR CI includes 1, or lies below
# the M1 estimate). Evaluated on all splits and on the batch-balanced ones.
if (!is.null(mres)) {
  ra <- list()
  for (sfx in c("", "_batch_balanced")) {
    g <- function(m, st) {
      mres[mres$method == m & mres$statistic == paste0(st, sfx), ]
    }
    m1 <- g("M1", "adjusted_rate_ratio_high_low")
    d2 <- g("M1-M2", "log_adjusted_rr_difference")
    d6 <- g("M1-M6", "log_adjusted_rr_difference")
    m11 <- g("M11", "adjusted_rate_ratio_high_low")
    if (!nrow(m1) || !nrow(m11)) {
      next
    }
    c1 <- isTRUE(m1$ci_lo > 1)
    c2 <- (nrow(d2) && isTRUE(d2$ci_lo > 0)) ||
      (nrow(d6) && isTRUE(d6$ci_lo > 0))
    c3 <- isTRUE(m11$ci_lo <= 1 && m11$ci_hi >= 1) ||
      isTRUE(m11$ci_hi < m1$estimate)
    ra[[length(ra) + 1L]] <- data.frame(
      splits = if (nzchar(sfx, keepNA = TRUE)) "batch_balanced" else "all",
      m1_adj_rr = m1$estimate,
      m1_ci_lo = m1$ci_lo,
      m1_ci_hi = m1$ci_hi,
      m1_jackknife = m1$jackknife,
      m1_m2_diff = if (nrow(d2)) d2$estimate else NA,
      m1_m2_ci_lo = if (nrow(d2)) d2$ci_lo else NA,
      m1_m2_ci_hi = if (nrow(d2)) d2$ci_hi else NA,
      m1_m6_diff = if (nrow(d6)) d6$estimate else NA,
      m1_m6_ci_lo = if (nrow(d6)) d6$ci_lo else NA,
      m1_m6_ci_hi = if (nrow(d6)) d6$ci_hi else NA,
      m11_adj_rr = m11$estimate,
      m11_ci_lo = m11$ci_lo,
      m11_ci_hi = m11$ci_hi,
      criterion1 = c1,
      criterion2 = c2,
      criterion3 = c3,
      R_A_passes = c1 && c2 && c3
    )
    add(
      "mouse_null",
      sprintf(
        "R_A (prereg) on %s splits: 1 = passes",
        if (nzchar(sfx, keepNA = TRUE)) "batch-balanced" else "all"
      ),
      as.numeric(c1 && c2 && c3),
      NA,
      NA,
      "wpc_summary/wpc_rule_RA.tsv",
      sprintf(
        "criterion 1 (M1 adjusted RR CI > 1) %s; 2 (M2 or M6 reduces it, CI > 0) %s; 3 (no comparable M11 excess) %s",
        c1,
        c2,
        c3
      )
    )
  }
  if (length(ra)) {
    ra <- do.call(rbind, ra)
    ue_write(ra, file.path(out, "wpc_rule_RA.parquet"), ue)
    utils::write.table(
      ra,
      file.path(out, "wpc_rule_RA.tsv"),
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )
  }
}
cb <- rp("null_mouse", "control_batches.parquet")
if (!is.null(cb)) {
  for (b in unique(cb$batch)) {
    x <- cb[cb$batch == b, ]
    add(
      "mouse_null",
      sprintf(
        "control animals in %s: fragment length mean (median; min and max in CI columns)",
        b
      ),
      stats::median(x$frag_length_mean),
      min(x$frag_length_mean),
      max(x$frag_length_mean),
      "null_mouse/control_batches.parquet column frag_length_mean",
      paste(x$run, collapse = ", ")
    )
  }
}

heads <- do.call(rbind, rows)
ue_write(heads, file.path(out, "wpc_headlines.parquet"), ue)
utils::write.table(
  heads,
  file.path(out, "wpc_headlines.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# Replication of the placenta pattern ------------------------------------------------------------
# Placenta null (WP-A, 29_null_summarise.R): q2_tertile_rate_ratio (crude
# p < 0.01 rate ratio, top vs bottom InfRV tertile, per design size; the
# 5v5 design matches the mouse splits), q3_glm_rate_ratios (adjusted rate
# ratio of the top tertile, the pre-registered model "prereg" on expression
# decile and isoform class, the same covariates as the mouse fit; 7v7
# splits; amendment A6) and q3_h2_method_interaction (model "prereg"). WP-A numbering: M9a = prereg M7 (class high).
# Placenta CI: WP-A's corrected jackknife (jk_lo, jk_hi); where WP-A's
# Monte Carlo correction collapsed it to a point (jk_lo == jk_hi), its
# uncorrected jackknife interval (jk_lo_raw, jk_hi_raw) is used and the row
# says so. Placenta replication (WP-B): placenta_audit/replication/summary_flags.
pl_dir <- Sys.getenv(
  "PLACENTA_NULL_SUMMARY",
  file.path(root, "null_placenta", "summary")
)
pl_flags <- Sys.getenv(
  "PLACENTA_REPLICATION_FLAGS",
  file.path(root, "placenta_audit", "replication", "summary_flags.parquet")
)
rpf <- function(f) if (file.exists(f)) as.data.frame(read_parquet(f)) else NULL
q2 <- rpf(file.path(pl_dir, "q2_tertile_rate_ratio.parquet"))
q3 <- rpf(file.path(pl_dir, "q3_glm_rate_ratios.parquet"))
h2 <- rpf(file.path(pl_dir, "q3_h2_method_interaction.parquet"))
pf <- rpf(pl_flags)
cmp <- list()
check <- function(
  what,
  m_est,
  m_lo,
  m_hi,
  p_est,
  p_lo,
  p_hi,
  src_m,
  src_p,
  note = ""
) {
  same_sign <- if (anyNA(c(m_est, p_est))) NA else sign(m_est) == sign(p_est)
  overlap <- if (anyNA(c(m_lo, m_hi, p_lo, p_hi))) {
    NA
  } else {
    (m_lo <= p_hi && p_lo <= m_hi)
  }
  # jarl-ignore undesirable_operator: appends to a result list of the enclosing scope
  cmp[[length(cmp) + 1L]] <<- data.frame(
    what = what,
    other = m_est,
    other_ci_lo = m_lo,
    other_ci_hi = m_hi,
    placenta = p_est,
    placenta_ci_lo = p_lo,
    placenta_ci_hi = p_hi,
    same_sign = same_sign,
    ci_overlap = overlap,
    replicates = isTRUE(same_sign) && isTRUE(overlap),
    source_other = src_m,
    source_placenta = src_p,
    note = note
  )
}
pl_ci <- function(r) {
  if (is.na(r$jk_lo) || is.na(r$jk_hi) || r$jk_lo == r$jk_hi) {
    list(
      lo = r$jk_lo_raw,
      hi = r$jk_hi_raw,
      note = "placenta CI: WP-A uncorrected jackknife (its corrected CI collapsed to a point)"
    )
  } else {
    list(
      lo = r$jk_lo,
      hi = r$jk_hi,
      note = "placenta CI: WP-A corrected jackknife"
    )
  }
}
mrow <- function(m, st) {
  if (is.null(mres)) NULL else mres[mres$method == m & mres$statistic == st, ]
}
# 1. Crude p < 0.01 rate ratio, top vs bottom tertile (log scale)
for (pair in list(
  c("M1", "M1"),
  c("M2", "M2"),
  c("M6", "M6"),
  c("M7", "M9a")
)) {
  r <- mrow(pair[1], "p01_ratio_high_low")
  p <- if (is.null(q2)) {
    NULL
  } else {
    q2[
      q2$design_set == "5v5" & q2$method == pair[2] & q2$threshold == "n_p01",
    ]
  }
  if (is.null(r) || !nrow(r)) {
    next
  }
  if (!is.null(p) && nrow(p)) {
    ci <- pl_ci(p[1, ])
    check(
      sprintf(
        "crude p<0.01 rate ratio top/bottom InfRV tertile, 5v5, %s (log scale)",
        pair[1]
      ),
      log(r$estimate),
      log(r$ci_lo),
      log(r$ci_hi),
      log(p$estimate[1]),
      log(ci$lo),
      log(ci$hi),
      paste(msrc, pair[1], "p01_ratio_high_low"),
      sprintf(
        "null_placenta/summary/q2_tertile_rate_ratio.parquet 5v5 %s n_p01",
        pair[2]
      ),
      ci$note
    )
  }
}
# 2. Adjusted rate ratio of the top tertile (log scale); also over the
# batch-balanced mouse splits
for (m in c("M1", "M2", "M6", "M1@_batch_balanced")) {
  sfx <- if (grepl("@", m, fixed = TRUE)) sub(".*@", "", m) else ""
  m <- sub("@.*", "", m)
  r <- mrow(m, paste0("adjusted_rate_ratio_high_low", sfx))
  p <- if (is.null(q3)) {
    NULL
  } else {
    q3[q3$method == m & q3$model == "prereg" & q3$term == "infrv_tertilehigh", ]
  }
  if (is.null(r) || !nrow(r) || is.null(p) || !nrow(p)) {
    next
  }
  ci <- pl_ci(p[1, ])
  check(
    sprintf(
      "adjusted rate ratio top/bottom tertile, %s%s (log scale; mouse 5v5 quasi-Poisson decile + isoform class; placenta 7v7 quasi-Poisson decile + isoform class)",
      m,
      sfx
    ),
    log(r$estimate),
    log(r$ci_lo),
    log(r$ci_hi),
    log(p$estimate[1]),
    log(ci$lo),
    log(ci$hi),
    paste(msrc, m, paste0("adjusted_rate_ratio_high_low", sfx)),
    sprintf(
      "null_placenta/summary/q3_glm_rate_ratios.parquet %s prereg infrv_tertilehigh",
      m
    ),
    ci$note
  )
}
# 3. Method interaction: log RR(M1) - log RR(alternative); placenta gives RR(alt)/RR(M1)
for (alt in c("M2", "M6")) {
  r <- mrow(paste0("M1-", alt), "log_adjusted_rr_difference")
  p <- if (is.null(h2)) {
    NULL
  } else {
    h2[h2$reference == "M1" & h2$alternative == alt & h2$model == "prereg", ]
  }
  if (is.null(r) || !nrow(r) || is.null(p) || !nrow(p)) {
    next
  }
  ci <- pl_ci(p[1, ])
  check(
    sprintf(
      "reduction of the top-tertile adjusted RR by %s (log RR(M1) - log RR(%s); > 0 = reduces)",
      alt,
      alt
    ),
    r$estimate,
    r$ci_lo,
    r$ci_hi,
    -log(p$estimate[1]),
    -log(ci$hi),
    -log(ci$lo),
    paste(msrc, paste0("M1-", alt), "log_adjusted_rr_difference"),
    sprintf(
      "null_placenta/summary/q3_h2_method_interaction.parquet reference M1 alternative %s model prereg (sign flipped, log)",
      alt
    ),
    ci$note
  )
}
# 4. Flags against non-replication: read halves vs placenta half-splits
flag_map <- c(
  flag_primary = "flag_primary",
  flag_unstable = "flag_unstable",
  flag_infrv_tertile = "flag_infrv_top_tertile",
  flag_ambiguous = "flag_ambiguous",
  flag_c4 = "flag_c4"
)
if (!is.null(fl)) {
  for (i in which(fl$flag %in% names(flag_map) & fl$method == "M1")) {
    r <- fl[i, ]
    p <- if (!is.null(pf)) pf[pf$flag == flag_map[[r$flag]], ] else NULL
    has <- !is.null(p) && nrow(p) >= 1L
    src_p <- if (has) {
      sprintf(
        "placenta_audit/replication/summary_flags.parquet flag %s",
        flag_map[[r$flag]]
      )
    } else {
      paste("not available:", pl_flags)
    }
    ci_d <- parse_ci(r$deficit_ci)
    pci_d <- if (has) parse_ci(p$deficit_ci[1]) else c(NA, NA)
    check(
      sprintf(
        "read halves %s %s vs placenta half-splits: %s deficit (unflagged - flagged replication)",
        r$halving,
        r$series,
        r$flag
      ),
      r$deficit,
      ci_d[1],
      ci_d[2],
      if (has) p$deficit[1] else NA,
      pci_d[1],
      pci_d[2],
      "rhalf/replication/summary_flags.parquet",
      src_p
    )
    ci_m <- parse_ci(r$matched_deficit_ci)
    pci_m <- if (has) parse_ci(p$matched_deficit_ci[1]) else c(NA, NA)
    check(
      sprintf(
        "read halves %s %s vs placenta half-splits: %s matched deficit (count x |z| cells)",
        r$halving,
        r$series,
        r$flag
      ),
      r$matched_deficit,
      ci_m[1],
      ci_m[2],
      if (has) p$matched_deficit[1] else NA,
      pci_m[1],
      pci_m[2],
      "rhalf/replication/summary_flags.parquet",
      src_p
    )
    ci_or <- parse_ci(r$adjusted_or_ci)
    pci_or <- if (has) parse_ci(p$adjusted_or_ci[1]) else c(NA, NA)
    check(
      sprintf(
        "read halves %s %s vs placenta half-splits: %s adjusted OR (log scale)",
        r$halving,
        r$series,
        r$flag
      ),
      log(r$adjusted_or),
      log(ci_or[1]),
      log(ci_or[2]),
      if (has) log(p$adjusted_or[1]) else NA,
      log(pci_or[1]),
      log(pci_or[2]),
      "rhalf/replication/summary_flags.parquet",
      src_p,
      "read-half ORs can be degenerate (quasi-separation: almost every unflagged call replicates)"
    )
  }
}
if (length(cmp)) {
  cmp <- do.call(rbind, cmp)
  ue_write(cmp, file.path(out, "wpc_replication_check.parquet"), ue)
  utils::write.table(
    cmp,
    file.path(out, "wpc_replication_check.tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
}
ue_write_session(out, ue)
message(
  "WP-C summary: ",
  nrow(heads),
  " headline rows; ",
  if (is.data.frame(cmp)) nrow(cmp) else 0,
  " comparison rows"
)
