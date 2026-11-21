# WP-B (placenta audit), part 4: the MUST REPORT table of the package and the
# verdict on decision rule R_B.
#
# Purpose: collect, from the Parquet outputs of 30_audit_prevalence.R,
# 31_halfsplit_replication.R and 32_story.R, every number the WP-B plan asks
# for into one long table (item, statistic, value, CI, exact source file and
# column), add the non-circular link to the WP-A null (flag rate of real calls
# against the flag rate of WP-A null calls, read from WP-A's per-split
# tables), and write the R_B verdict. No new statistics are estimated here
# except the WP-A link.
#
# Inputs: $ISOFORMIC_UE_OUT/placenta_audit/{prevalence,replication,story}/;
#   WP-A per-split tables $ISOFORMIC_UE_OUT/null_placenta/dte/main/tx/*.parquet
#   and null_placenta/splits/splits.parquet (optional: the link is skipped
#   when they are absent).
# Env vars: ISOFORMIC_UE_OUT, ISOFORMIC_UE_PKG, ISOFORMIC_UE_THREADS,
#   ISOFORMIC_UE_PILOT.
# Also reports the corrections of the WP-B verification: q-matched gene
#   references, reasons for lost conclusions, the null link within
#   expression deciles, the stricter control C4b, sign concordance as the
#   outcome, C1 per number of relabelled samples, the H3 rescue share and the
#   genuinely held-out replication of the story isoforms.
# Outputs: placenta_audit/must_report.parquet, placenta_audit/wpa_link.parquet,
#   placenta_audit/verdict_RB.parquet (copies of the replication verdict).
# Resources (estimate): seconds to a minute, < 2 GB.
# Run from the repository root, inside a memory cap of about 4 GB (see
# README.md):
#   pixi run --frozen Rscript data-raw/uncertainty_evidence/39_audit_summarise.R

source("data-raw/uncertainty_evidence/ue_common.R")
ue <- ue_setup("39_audit_summarise")
source("data-raw/uncertainty_evidence/38_audit_common.R")
seed <- wpb_seed(ue, 900L)
out_dir <- wpb_dir()
say <- wpb_logger(wpb_dir("logs", "39_audit_summarise.log"))
P <- function(...) file.path(out_dir, ...)
rd <- function(...) wpb_read(P(...))
rows <- list()
add <- function(
  item,
  statistic,
  value,
  ci = NA_character_,
  source,
  note = NA_character_
) {
  # jarl-ignore undesirable_operator: appends to a result list of the enclosing scope
  rows[[length(rows) + 1L]] <<- data.frame(
    item = item,
    statistic = statistic,
    value = as.numeric(value),
    ci = as.character(ci),
    source = source,
    note = note
  )
}

# (1) Prevalence per family -----------------------------------------------------------
sf <- rd("prevalence", "summary_family.parquet")
for (i in seq_len(nrow(sf))) {
  r <- sf[i, ]
  add(
    paste0("1_prevalence_", r$family),
    r$statistic,
    r$estimate,
    paste0("gene-cluster ", r$ci_gene_cluster, "; draw ", r$ci_draw),
    "placenta_audit/prevalence/summary_family.parquet: estimate, ci_gene_cluster, ci_draw",
    "sensitivity of conclusions to quantification uncertainty (self-consistency), not an error rate"
  )
  add(
    paste0("1_prevalence_", r$family),
    paste0(r$statistic, "__p_matched_gene_reference"),
    r$ref_gene_matched,
    paste0("gene-cluster ", r$ref_ci_gene_cluster, "; draw ", r$ref_ci_draw),
    "placenta_audit/prevalence/summary_family.parquet: ref_gene_matched, ref_ci_gene_cluster, ref_ci_draw"
  )
  add(
    paste0("1_prevalence_", r$family),
    paste0(r$statistic, "__q_matched_gene_reference"),
    r$ref_gene_q_matched,
    paste0(
      "gene-cluster ",
      r$ref_q_ci_gene_cluster,
      "; draw ",
      r$ref_q_ci_draw
    ),
    "placenta_audit/prevalence/summary_family.parquet: ref_gene_q_matched, ref_q_ci_gene_cluster, ref_q_ci_draw",
    sprintf(
      "matched on naive BH q (largest naive p of the conclusions %.3g); the families have different BH cut-offs",
      r$max_naive_p_of_conclusions
    )
  )
}
lr <- rd("prevalence", "summary_lost_reason.parquet")
for (i in seq_len(nrow(lr))) {
  it <- paste0(
    "1_lost_reason_",
    lr$family[i],
    "_",
    gsub(" ", "_", lr$alternative[i], fixed = TRUE)
  )
  add(
    it,
    "pct_lost_not_tested",
    lr$pct_lost_not_tested[i],
    NA,
    "placenta_audit/prevalence/summary_lost_reason.parquet: pct_lost_not_tested",
    sprintf(
      "%d of %d lost (of %d conclusions) not tested by the alternative",
      as.integer(lr$n_lost_not_tested[i]),
      as.integer(lr$n_lost[i]),
      as.integer(lr$n_conclusions[i])
    )
  )
  add(
    it,
    "pct_lost_tested_not_called",
    lr$pct_lost_tested_not_called[i],
    NA,
    "placenta_audit/prevalence/summary_lost_reason.parquet: pct_lost_tested_not_called",
    sprintf(
      "%d tested and not called",
      as.integer(lr$n_lost_tested_not_called[i])
    )
  )
}
for (fam in unique(sf$family)) {
  r <- sf[sf$family == fam, ][1, ]
  add(
    paste0("1_prevalence_", fam),
    "n_gained_divided_or_scaled",
    r$gained_divided,
    NA,
    "placenta_audit/prevalence/summary_family.parquet: gained_divided"
  )
  add(
    paste0("1_prevalence_", fam),
    "n_gained_swish",
    r$gained_swish,
    NA,
    "placenta_audit/prevalence/summary_family.parquet: gained_swish"
  )
}
ch <- rd("prevalence", "summary_changed_q001.parquet")
for (i in seq_len(nrow(ch))) {
  add(
    paste0("1_prevalence_", ch$family[i]),
    "pct_of_changed_with_naive_q_lt_0.001",
    ch$pct_of_changed_with_naive_q_lt_0.001[i],
    NA,
    "placenta_audit/prevalence/summary_changed_q001.parquet",
    ch$definition[i]
  )
}

# Non-circular link: real calls against WP-A null calls -------------------------------
nl <- rd("prevalence", "null_link.parquet")
for (i in seq_len(nrow(nl))) {
  add(
    "1_null_link_wpb_20_splits",
    nl$set[i],
    nl$pct_flagged[i],
    sprintf("Wilson [%.4g, %.4g]", nl$wilson_lo[i], nl$wilson_hi[i]),
    "placenta_audit/prevalence/null_link.parquet: pct_flagged",
    sprintf("n = %d; %s", nl$n_calls[i], nl$note[i])
  )
}
nds <- rd("prevalence", "null_link_decile_summary.parquet")
add(
  "1_null_link_wpb_by_expression_decile",
  "pct_flagged_null_standardised_to_real_deciles",
  nds$pct_flagged_null_standardised_to_real_deciles,
  NA,
  "placenta_audit/prevalence/null_link_decile_summary.parquet",
  sprintf(
    "real %.2f%%, null crude %.2f%%; null rate above real in %d of %d deciles (null_link_by_decile.parquet)",
    nds$pct_flagged_real,
    nds$pct_flagged_null_crude,
    as.integer(nds$n_deciles_null_flag_rate_above_real),
    as.integer(nds$n_deciles_compared)
  )
)
add(
  "1_null_link_wpb_by_expression_decile",
  "MH_OR_flagged_null_vs_real",
  nds$mh_or_flagged_null_vs_real,
  paste("MH", nds$mh_or_ci),
  "placenta_audit/prevalence/null_link_decile_summary.parquet",
  nds$note
)
wpa_dir <- file.path(ue_out_root(), "null_placenta", "dte", "main", "tx")
wpa_files <- list.files(wpa_dir, pattern = "\\.parquet$", full.names = TRUE)
if (length(wpa_files) && !wpb_pilot()) {
  sp <- wpb_read(file.path(
    ue_out_root(),
    "null_placenta",
    "splits",
    "splits.parquet"
  ))
  x <- wpb_query(sprintf(
    "SELECT split_id, transcript_id, M1_p, M9_assessed, M9_high FROM read_parquet([%s], union_by_name = true)",
    paste(sprintf("'%s'", wpa_files), collapse = ",")
  ))
  # WP-A computes M9 (the classes) only on some splits: keep those
  has_m9 <- tapply(!is.na(x$M9_assessed), x$split_id, any)
  n_m9 <- sum(has_m9)
  x <- x[!is.na(x$M1_p) & x$split_id %in% names(has_m9)[has_m9], ]
  x$M1_q <- stats::ave(x$M1_p, x$split_id, FUN = function(p) {
    stats::p.adjust(p, "BH")
  })
  x$size <- sp$size[match(x$split_id, sp$split_id)]
  link <- do.call(
    rbind,
    lapply(c("all", sort(unique(x$size))), function(s) {
      d <- if (s == "all") x else x[x$size == s, ]
      d <- d[d$M9_assessed %in% TRUE, ]
      one <- function(sel, set) {
        n <- sum(sel)
        k <- sum(sel & !d$M9_high)
        ci <- if (n) stats::prop.test(k, n)$conf.int * 100 else c(NA, NA)
        data.frame(
          size = s,
          set = set,
          n_splits = length(unique(d$split_id[sel])),
          n_calls = n,
          pct_flagged = if (n) 100 * k / n else NA_real_,
          ci_lo = ci[1],
          ci_hi = ci[2]
        )
      }
      rbind(
        one(d$M1_p < 0.001, "WP-A null M1 calls p < 0.001"),
        one(d$M1_q < 0.05, "WP-A null M1 calls BH q < 0.05 (per split)"),
        one(rep(TRUE, nrow(d)), "WP-A all tested transcripts (base rate)")
      )
    })
  )
  link$n_split_files <- length(wpa_files)
  link$n_split_files_with_m9 <- n_m9
  link$note <- "flag = !M9_high (WP-A matrix form of the package classes, draws of the split samples); CI: Wilson, calls pooled over splits (ignores clustering)"
  ue_write(link, P("wpa_link.parquet"), ue, seed)
  print(link)
  for (i in seq_len(nrow(link))) {
    add(
      "1_null_link_wpa",
      paste0(link$set[i], " [", link$size[i], "]"),
      link$pct_flagged[i],
      sprintf("Wilson [%.4g, %.4g]", link$ci_lo[i], link$ci_hi[i]),
      "placenta_audit/wpa_link.parquet: pct_flagged (from null_placenta/dte/main/tx)",
      sprintf(
        "n = %d calls in %d splits; %d WP-A split files read, %d with M9",
        link$n_calls[i],
        link$n_splits[i],
        length(wpa_files),
        n_m9
      )
    )
  }
}

# (2)-(7) Replication ------------------------------------------------------------------------
fl <- rd("replication", "summary_flags.parquet")
src <- "placenta_audit/replication/summary_flags.parquet"
for (i in seq_len(nrow(fl))) {
  r <- fl[i, ]
  it <- paste0("2_replication_", r$flag)
  add(
    it,
    "pct_flagged",
    r$pct_flagged,
    NA,
    paste0(src, ": pct_flagged"),
    sprintf("n = %d calls, %d flagged", r$n_calls, r$n_flagged)
  )
  add(
    it,
    "replication_rate_flagged",
    r$rate_flagged,
    NA,
    paste0(src, ": rate_flagged")
  )
  add(
    it,
    "replication_rate_unflagged",
    r$rate_unflagged,
    NA,
    paste0(src, ": rate_unflagged")
  )
  add(
    it,
    "deficit_unflagged_minus_flagged",
    r$deficit,
    r$deficit_ci,
    paste0(src, ": deficit, deficit_ci"),
    paste("split range (not a CI):", r$deficit_split_range_not_ci)
  )
  add(
    it,
    "matched_deficit",
    r$matched_deficit,
    r$matched_deficit_ci,
    paste0(src, ": matched_deficit")
  )
  add(
    it,
    "adjusted_OR_replication",
    r$adjusted_or,
    r$adjusted_or_ci,
    paste0(src, ": adjusted_or"),
    paste("split range (not a CI):", r$adjusted_or_split_range_not_ci)
  )
  add(
    paste0("4_delta_auc_", r$flag),
    "delta_AUC_nonreplication",
    r$delta_auc,
    r$delta_auc_ci,
    paste0(src, ": delta_auc"),
    sprintf("AUC of count + |z| alone: %.4f", r$auc_count_z)
  )
  add(
    paste0("5_ppv_cost_", r$flag),
    "PPV_nonreplication",
    r$ppv_nonreplication,
    NA,
    paste0(src, ": ppv_nonreplication")
  )
  add(
    paste0("5_ppv_cost_", r$flag),
    "sensitivity_nonreplication",
    r$sensitivity_nonreplication,
    NA,
    paste0(src, ": sensitivity_nonreplication")
  )
  add(
    paste0("5_ppv_cost_", r$flag),
    "replicated_calls_removed_n",
    r$replicated_removed_n,
    NA,
    paste0(src, ": replicated_removed_n")
  )
  add(
    paste0("5_ppv_cost_", r$flag),
    "replicated_calls_removed_pct",
    r$replicated_removed_pct,
    NA,
    paste0(src, ": replicated_removed_pct")
  )
  add(
    paste0("5_ppv_cost_", r$flag),
    "nonreplicated_calls_removed_n",
    r$nonreplicated_removed_n,
    NA,
    paste0(src, ": nonreplicated_removed_n")
  )
}
pl <- rd("replication", "summary_placebo.parquet")
for (i in seq_len(nrow(pl))) {
  add(
    "2_placebo_flag",
    paste0(pl$statistic[i], "_placebo_mean"),
    pl$placebo_mean[i],
    sprintf(
      "placebo 2.5-97.5%% [%.4g, %.4g]",
      pl$placebo_2.5[i],
      pl$placebo_97.5[i]
    ),
    "placenta_audit/replication/summary_placebo.parquet",
    sprintf(
      "observed %.4g; one-sided empirical p %.4g (200 placebo draws)",
      pl$observed[i],
      pl$empirical_p_one_sided[i]
    )
  )
}
if (file.exists(P("replication", "summary_c1.parquet"))) {
  c1 <- rd("replication", "summary_c1.parquet")
  add(
    "2_C1_permutation",
    "replication_rate_observed",
    c1$replication_rate_observed,
    NA,
    "placenta_audit/replication/summary_c1.parquet"
  )
  add(
    "2_C1_permutation",
    "replication_rate_label_permuted",
    c1$replication_rate_perm,
    paste(
      "range over permutations (not a CI):",
      c1$replication_rate_perm_range_not_ci
    ),
    "placenta_audit/replication/summary_c1.parquet",
    sprintf("%d permutations", c1$n_permutations)
  )
  add(
    "2_C1_permutation",
    "deficit_primary_label_permuted",
    c1$deficit_perm_pooled,
    paste("range (not a CI):", c1$deficit_perm_range_not_ci),
    "placenta_audit/replication/summary_c1.parquet"
  )
  add(
    "2_C1_permutation",
    "excess_deficit_over_C1",
    c1$excess_deficit_over_c1,
    c1$excess_deficit_ci,
    "placenta_audit/replication/summary_c1.parquet",
    c1$note
  )
  c1s <- rd("replication", "summary_c1_strata.parquet")
  for (i in seq_len(nrow(c1s))) {
    add(
      "2_C1_permutation_by_fixed_InfRV_tertile",
      paste0("excess_replication_rate_", c1s$fixed_infrv_tertile[i]),
      c1s$excess_rate_over_c1[i],
      NA,
      "placenta_audit/replication/summary_c1_strata.parquet",
      sprintf(
        "observed %.4f (n %d) vs permuted %.4f (n %d); flagged %.1f%%",
        c1s$rate_observed[i],
        as.integer(c1s$n_observed[i]),
        c1s$rate_perm[i],
        as.integer(c1s$n_perm[i]),
        c1s$pct_flagged_observed[i]
      )
    )
  }
}
c4 <- rd("replication", "summary_c4.parquet")
add(
  "3_C4_count_matched",
  "deficit_primary_minus_deficit_C4",
  c4$estimate,
  c4$ci,
  "placenta_audit/replication/summary_c4.parquet: estimate, ci",
  sprintf(
    "deficit primary %.4f, deficit C4 %.4f",
    c4$deficit_primary,
    c4$deficit_c4
  )
)
c4b <- rd("replication", "summary_c4b.parquet")
add(
  "3_C4b_count_z_predicted",
  "deficit_primary_minus_deficit_C4b",
  c4b$estimate,
  c4b$ci,
  "placenta_audit/replication/summary_c4b.parquet: estimate, ci",
  sprintf(
    "deficit primary %.4f, deficit C4b %.4f (not preregistered; %s)",
    c4b$deficit_primary,
    c4b$deficit_c4b,
    c4b$note
  )
)
sct <- rd("replication", "summary_sign_concordance.parquet")
for (i in seq_len(nrow(sct))) {
  add(
    paste0("2_sign_concordance_outcome_", sct$flag[i]),
    "adjusted_OR_sign_concordance",
    sct$adjusted_or_sign_concordance[i],
    sct$ci[i],
    "placenta_audit/replication/summary_sign_concordance.parquet",
    sprintf(
      "sign concordance flagged %.4f vs unflagged %.4f",
      sct$sign_concordance_flagged[i],
      sct$sign_concordance_unflagged[i]
    )
  )
}
if (file.exists(P("replication", "summary_c1_by_relabelled.parquet"))) {
  cb <- rd("replication", "summary_c1_by_relabelled.parquet")
  for (i in seq_len(nrow(cb))) {
    add(
      "2_C1_by_n_relabelled",
      paste0("deficit_n_relabelled_", cb$n_relabelled[i]),
      cb$deficit[i],
      NA,
      "placenta_audit/replication/summary_c1_by_relabelled.parquet",
      sprintf(
        "%d permutations, %d true PE of 4 in the permuted PE group, residual contrast %.2f, replication %.4f, flagged %.1f%%",
        as.integer(cb$n_permutations[i]),
        as.integer(cb$true_pe_in_permuted_pe_group[i]),
        cb$mean_residual_contrast[i],
        cb$replication_rate[i],
        cb$pct_flagged[i]
      )
    )
  }
}
c3 <- rd("replication", "summary_c3.parquet")
add(
  "6_C3_gene_level",
  "slope_log_odds_replication_per_doubling_max_InfRV",
  c3$slope_log_odds_per_doubling_max_infrv,
  c3$ci,
  "placenta_audit/replication/summary_c3.parquet",
  sprintf(
    "n = %d gene calls; replication %.4f (transcripts %.4f); low vs high max-InfRV tertile %.4f vs %.4f",
    as.integer(c3$n_gene_calls),
    c3$replication_rate_gene,
    c3$replication_rate_transcript,
    c3$rate_low_max_infrv_tertile,
    c3$rate_high_max_infrv_tertile
  )
)
if (file.exists(P("replication", "summary_h3.parquet"))) {
  h3 <- rd("replication", "summary_h3.parquet")
  add(
    "7_H3_groups",
    "group_replication_rate",
    h3$group_replication_rate,
    h3$group_rate_ci_wilson_no_clustering,
    "placenta_audit/replication/summary_h3.parquet",
    sprintf(
      "%d groups of %d flagged non-replicating calls in a group (of %d flagged non-replicating)",
      as.integer(h3$n_groups_tested),
      as.integer(h3$n_of_them_in_a_group),
      as.integer(h3$n_flagged_nonreplicating_calls)
    )
  )
  add(
    "7_H3_groups",
    "pct_of_flagged_nonreplicating_calls_rescued_by_group",
    h3$pct_of_flagged_nonreplicating_calls_rescued,
    NA,
    "placenta_audit/replication/summary_h3.parquet",
    h3$note
  )
  add(
    "7_H3_groups",
    "pct_selected_groups_discovery_p_lt_0.1",
    h3$pct_selected_groups_discovery_p_lt_0.1,
    NA,
    "placenta_audit/replication/summary_h3.parquet"
  )
  add(
    "7_H3_groups",
    "all_groups_with_a_call_replication_rate",
    h3$all_groups_with_a_call_replication_rate,
    NA,
    "placenta_audit/replication/summary_h3.parquet",
    sprintf("n = %d groups", as.integer(h3$n_groups_with_any_call))
  )
  add(
    "7_H3_groups",
    "unflagged_single_transcript_replication_rate",
    h3$unflagged_single_rate,
    NA,
    "placenta_audit/replication/summary_h3.parquet",
    sprintf(
      "n = %d; biased reference (the groups are selected on non-replication)",
      as.integer(h3$n_unflagged_single)
    )
  )
}
sec <- rd("replication", "summary_secondary.parquet")
for (i in seq_len(nrow(sec))) {
  add(
    "2_secondary_outcomes",
    sec$metric[i],
    sec$pooled_or_mean[i],
    paste("split range (not a CI):", sec$split_range_not_ci[i]),
    "placenta_audit/replication/summary_secondary.parquet",
    sprintf("flagged %.4g, unflagged %.4g", sec$flagged[i], sec$unflagged[i])
  )
}
jm <- rd("replication", "summary_joint_model.parquet")
for (i in seq_len(nrow(jm))) {
  add(
    "2_joint_model",
    paste0("adjusted_OR_", jm$flag[i]),
    jm$adjusted_or_joint[i],
    jm$ci[i],
    "placenta_audit/replication/summary_joint_model.parquet"
  )
}
ab <- rd("replication", "summary_ablation.parquet")
add(
  "2_ablation_15_draws",
  "primary_flag_agreement_30_vs_15",
  ab$primary_flag_agreement_30_vs_15,
  NA,
  "placenta_audit/replication/summary_ablation.parquet"
)
if (file.exists(P("replication", "summary_dtu.parquet"))) {
  dt <- rd("replication", "summary_dtu.parquet")
  add(
    "2_DTU_replication",
    "rate_flagged",
    dt$rate_flagged,
    NA,
    "placenta_audit/replication/summary_dtu.parquet",
    sprintf(
      "n = %d calls, %.1f%% flagged, unflagged rate %.4f",
      as.integer(dt$n_calls),
      dt$pct_flagged,
      dt$rate_unflagged
    )
  )
}
vd <- rd("replication", "summary_verdict_RB.parquet")
ue_write(vd, P("verdict_RB.parquet"), ue, seed)
add(
  "verdict_R_B",
  "passed",
  as.numeric(vd$passed),
  NA,
  "placenta_audit/replication/summary_verdict_RB.parquet",
  sprintf(
    "adjusted OR %.4f %s (below 1: %s); deficit(primary) - deficit(C4) %.4f %s (above 0: %s)",
    vd$adjusted_or,
    vd$adjusted_or_ci,
    vd$or_ci_below_1,
    vd$deficit_minus_c4,
    vd$deficit_minus_c4_ci,
    vd$beats_count_matched
  )
)

# (8) Story -----------------------------------------------------------------------------------
if (file.exists(P("story", "flt1_isoforms.parquet"))) {
  fl1 <- rd("story", "flt1_isoforms.parquet")
  for (i in seq_len(nrow(fl1))) {
    r <- fl1[i, ]
    add(
      "8_story_FLT1",
      paste0(r$transcript_name, "_naive_log2FC"),
      r$naive_log2FC,
      sprintf("[%.3f, %.3f]", r$naive_ci_lo, r$naive_ci_hi),
      "placenta_audit/story/flt1_isoforms.parquet",
      sprintf(
        "%s, %s; mean count %.0f; q %.3g; swish log2FC %.3f q %.3g; sign_agree %.2f; sig_frac %.2f; class %s; replicated when called %s of %d pairs; cohort-sign support in V halves %.2f (not held out); null p<0.001 in %s splits",
        r$form,
        r$transcript_type,
        r$mean_count,
        r$naive_q,
        r$swish_log2FC,
        r$swish_q,
        r$sign_agree,
        r$sig_frac,
        r$confidence,
        format(round(r$replicated_when_called, 3)),
        as.integer(r$n_called_discovery),
        r$cohort_sign_support_in_v_halves,
        r$null_n_p_lt_0.001
      )
    )
  }
  gs <- rd("story", "story_genes_gene_level.parquet")
  for (i in seq_len(nrow(gs))) {
    r <- gs[i, ]
    add(
      "8_story_gene_level",
      paste0(r$gene_name, "_naive_log2FC"),
      r$naive_log2FC,
      sprintf("[%.3f, %.3f]", r$naive_ci_lo, r$naive_ci_hi),
      "placenta_audit/story/story_genes_gene_level.parquet",
      sprintf(
        "q %.3g; swish log2FC %.3f q %.3g; sign_agree %.2f; sig_frac %.2f; replicated when called %s of %d pairs; cohort-sign support in V halves %.2f (not held out)",
        r$naive_q,
        r$swish_log2FC,
        r$swish_q,
        r$sign_agree,
        r$sig_frac,
        format(round(r$replicated_when_called, 3)),
        as.integer(r$n_called_discovery),
        r$cohort_sign_support_in_v_halves
      )
    )
  }
  for (f in c(
    "story_genes_isoforms.parquet",
    "rule_selected_isoforms.parquet"
  )) {
    if (!file.exists(P("story", f))) {
      next
    }
    si <- rd("story", f)
    si <- si[(si$naive_q < 0.05) %in% TRUE, ]
    for (i in seq_len(nrow(si))) {
      r <- si[i, ]
      add(
        "8_story_isoforms_called",
        paste0(r$transcript_name, "_replicated_when_called"),
        r$replicated_when_called,
        NA,
        paste0(
          "placenta_audit/story/",
          f,
          ": replicated_when_called, n_called_discovery"
        ),
        sprintf(
          "called in %d discovery halves; naive log2FC %.3f, swish %.3f; mean count %.0f; class %s; sign_agree %.2f",
          as.integer(r$n_called_discovery),
          r$naive_log2FC,
          r$swish_log2FC,
          r$mean_count,
          r$confidence,
          r$sign_agree
        )
      )
    }
  }
  sm <- rd("story", "flt1_soluble_membrane.parquet")
  for (i in seq_len(nrow(sm))) {
    add(
      "8_story_FLT1_soluble_membrane",
      sm$quantity[i],
      sm$point[i],
      sprintf("draw range [%.3f, %.3f]", sm$draw_min[i], sm$draw_max[i]),
      "placenta_audit/story/flt1_soluble_membrane.parquet",
      sprintf("same sign in all 30 draws: %s", sm$same_sign_in_all_draws[i])
    )
  }
}

mr <- do.call(rbind, rows)
ue_write(mr, P("must_report.parquet"), ue, seed)
options(width = 250)
print(mr[, c("item", "statistic", "value", "ci")], right = FALSE)
ue_write_session(out_dir, ue)
say("must_report.parquet written: %d rows", nrow(mr))
