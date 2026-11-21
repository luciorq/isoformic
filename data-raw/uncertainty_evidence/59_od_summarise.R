# WP-D, step 3: headline numbers of the overdispersion prior package
#
# Purpose: collect the numbers that WP-D must report from the tables of
# 50_od_prior.R and 51_od_prior_impact.R into one table, each with its 95%
# interval, the interval method, the source table and a tag ("verified":
# computed by these scripts from the data; "estimate": an approximation
# such as a percentile of a resampling distribution that depends on the
# number of resamples). Also writes the proposal for the maintainer as text
# (no package file is changed).
#
# Inputs: $ISOFORMIC_UE_OUT/od_prior/*.parquet (50_ to 54_ scripts).
#
# Intervals: the primary interval of a sample-level statistic is point +/-
# 1.96 x the delete-one-block jackknife SE (53_od_prior_jackknife.R). The
# percentile interval of the sample bootstrap is reported in the note only:
# the prior falls as the number of distinct samples grows, resamples with
# replacement have fewer distinct samples, so that bootstrap is biased upward
# (99.8% of placenta resamples lie above the point estimate).
# Transcript-level intervals (genome-wide priors of one dataset, draw-count
# ratios) come from the gene-clustered bootstrap of
# 54_od_prior_gene_bootstrap.R. The transcript bootstrap of 50_od_prior.R
# resamples isoforms of the same gene independently and is about 6 times too
# narrow; it is kept in the notes as a lower bound only.
# DTE call changes (51_od_prior_impact.R) are reported for the common tested
# set (both methods on the same transcripts, one BH correction): the
# separate-filter rows mix real changes with transcripts that only one run
# tested (n_lost_untested, n_gained_untested).
# Seeds are read from the tables (boot_seed, panel_seed, seed columns) and
# listed in the summary, so that no seed is typed by hand.
# Env vars: ISOFORMIC_UE_OUT (required), ISOFORMIC_OD_OUT; ISOFORMIC_V50_DIR
#   and ISOFORMIC_UE_RLIB only for the optional extra R library.
# Outputs: od_prior/summary.parquet, od_prior/summary.tsv (same rows, for
#   reading), od_prior/proposal.txt.
# Resources (estimate): seconds, under 2 GB.
# Run from the repository root, inside a memory cap of about 4 GB (see
# README.md):
#   pixi run --frozen Rscript data-raw/uncertainty_evidence/59_od_summarise.R

source("data-raw/uncertainty_evidence/ue_common.R")
ue_use_rlib()
ue <- ue_setup("59_od_summarise")
out <- Sys.getenv("ISOFORMIC_OD_OUT", file.path(ue_out_root(), "od_prior"))
rd <- function(name) read_parquet(file.path(out, paste0(name, ".parquet")))

priors <- rd("priors")
boot <- rd("bootstrap_replicates")
splits <- rd("splits_prior")
panels <- rd("panels_prior")
fdf <- rd("prior_df_fit")
jkt <- rd("jackknife")
bbias <- rd("bootstrap_bias")
gb <- rd("gene_bootstrap")
gq <- function(d, st) {
  r <- gb[gb$dataset == d & gb$statistic == st, ]
  stopifnot(nrow(r) == 1L)
  r
}
gci <- "gene-clustered bootstrap (genes resampled), percentile"
jq <- function(d, st) {
  r <- jkt[jkt$dataset == d & jkt$statistic == st, ]
  stopifnot(nrow(r) == 1L)
  r
}
jci <- "point +/- 1.96 x delete-one-block jackknife SE"
cs <- rd("catchsalmon_check")
dct <- rd("draw_count_transcripts")
imp_phi <- rd("impact_phi")
imp_calls <- rd("impact_calls")
imp_pl <- rd("impact_placenta_priors")
fold <- rd("fold_map")
fpr <- rd("filter_priors")
byd <- rd("impact_phi_by_df")
fprior <- function(d, floor) {
  fpr$prior[fpr$dataset == d & fpr$count_floor %in% floor]
}

rows <- list()
add <- function(
  item,
  value,
  lo = NA,
  hi = NA,
  ci = "",
  tag = "verified",
  source = "",
  note = ""
) {
  # jarl-ignore undesirable_operator: appends to a result list of the enclosing scope
  rows[[length(rows) + 1L]] <<- data.frame(
    item = item,
    value = value,
    ci_lo = lo,
    ci_hi = hi,
    ci_method = ci,
    tag = tag,
    source = source,
    note = note
  )
}
pr <- function(d, v) priors[priors$dataset == d & priors$variant == v, ]
bq <- function(d, v, k) {
  x <- boot$value[boot$dataset == d & boot$variant == v & boot$kind == k]
  c(stats::quantile(x, c(0.025, 0.975), names = FALSE), length(x))
}
wilson <- function(x, n) {
  z <- 1.959964
  p <- x / n
  c0 <- (p + z^2 / (2 * n)) / (1 + z^2 / n)
  h <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / (1 + z^2 / n)
  c(c0 - h, c0 + h)
}

# (1) Placenta genome-wide prior
a <- pr("placenta", "a_all_targets")
sb <- bq("placenta", "a_all_targets", "sample_block")
tb <- bq("placenta", "a_all_targets", "transcript")
ja <- jq("placenta", "prior a_all_targets")
bb <- bbias[bbias$variant == "a_all_targets", ]
add(
  "(1) placenta prior, all 654,828 targets, 25 runs",
  a$prior,
  ja$ci_lo,
  ja$ci_hi,
  paste0(jci, " (23 blocks)"),
  tag = "verified (point); estimate (CI)",
  source = "priors.parquet a_all_targets; jackknife.parquet",
  note = sprintf(
    paste(
      "jackknife SE %.4f; bootstrap SE %.4f (normal interval %.4f-%.4f);",
      "biased percentile bootstrap %.4f-%.4f (mean %.4f, %.1f%% above point);",
      "unfloored %.4f; median phi %.4f; median df %g; %d transcripts with df > 0"
    ),
    ja$jackknife_se,
    bb$boot_se,
    a$prior - 1.959964 * bb$boot_se,
    a$prior + 1.959964 * bb$boot_se,
    sb[1],
    sb[2],
    bb$boot_mean,
    100 * bb$frac_boot_above_point,
    a$prior_unfloored,
    a$median_phi,
    a$median_df,
    a$n_transcripts_df
  )
)
ga <- gq("placenta", "prior a_all_targets")
add(
  "(1) placenta prior, gene-clustered bootstrap",
  a$prior,
  ga$ci_lo,
  ga$ci_hi,
  sprintf("%s, %d resamples", gci, ga$n_boot),
  tag = "estimate",
  source = "gene_bootstrap.parquet",
  note = sprintf(
    paste(
      "bootstrap SD %.4f over %d genes; the transcript bootstrap of",
      "50_od_prior.R (%.4f-%.4f, SD %.4f) ignores the correlation of isoforms",
      "of a gene and is a lower bound only"
    ),
    ga$boot_sd,
    ga$n_genes,
    tb[1],
    tb[2],
    stats::sd(boot$value[
      boot$dataset == "placenta" &
        boot$variant == "a_all_targets" &
        boot$kind == "transcript"
    ])
  )
)

# (2) Variants
for (v in unique(priors$variant[priors$dataset == "placenta"])) {
  p <- pr("placenta", v)
  ci <- c(NA, NA, NA)
  jst <- c(
    a_all_targets = "prior a_all_targets",
    b_primary_chromosomes = "prior b_primary_chromosomes",
    c_haplotype_folded = "prior c_haplotype_folded (S0 map)"
  )
  if (v %in% names(jst)) {
    jr <- jq("placenta", jst[[v]])
    ci <- c(jr$ci_lo, jr$ci_hi)
  }
  add(
    paste0("(2) placenta prior, ", v),
    p$prior,
    ci[1],
    ci[2],
    if (!is.na(ci[1])) jci else "",
    source = "priors.parquet",
    note = sprintf(
      "n samples %d; snapshot function %.6f; %s",
      p$n_samples,
      p$prior_snapshot_function,
      p$note
    )
  )
}
for (sz in sort(unique(splits$size))) {
  s <- splits[splits$size == sz & splits$source != "leave_one_unit_out", ]
  add(
    sprintf(
      "(2e) placenta prior, %dv%d null splits of control units, split range",
      sz,
      sz
    ),
    stats::median(s$prior),
    stats::quantile(s$prior, 0.025, names = FALSE),
    stats::quantile(s$prior, 0.975, names = FALSE),
    sprintf("split range (2.5-97.5%% of %d splits)", nrow(s)),
    source = "splits_prior.parquet",
    note = unique(s$source)
  )
}
loo <- splits[splits$source == "leave_one_unit_out", ]
# WP-A 7v7_louo: 15 left-out units x 14 splits = 210 splits; the prior of a
# split depends only on its 14 units, so each of the 15 values counts 14 times
v210 <- rep(loo$prior, each = 14L)
add(
  "(2e) placenta prior, the 210 WP-A-style 7v7 splits (7v7_louo), split range",
  stats::median(v210),
  stats::quantile(v210, 0.025, names = FALSE),
  stats::quantile(v210, 0.975, names = FALSE),
  "split range (2.5-97.5% of 210 splits = 15 distinct 14-unit sets x 14)",
  source = "splits_prior.parquet leave_one_unit_out",
  note = sprintf(
    "min %.4f, max %.4f over the 15 distinct sets",
    min(loo$prior),
    max(loo$prior)
  )
)

# (3) Knockdown panel prior vs genome-wide prior
for (s in c("GSE129540", "GSE75489")) {
  pb <- pr(s, "panel_bundled_draws")
  pu <- pr(s, "panel_unrounded_moments")
  g <- pr(s, "genome_all_targets")
  add(
    sprintf("(3) %s panel prior (bundled draws, exact)", s),
    pb$prior_snapshot_function,
    source = "priors.parquet panel_bundled_draws",
    note = sprintf(
      "%s; unfloored %.6f; same panel from unrounded moments %.6f",
      pb$note,
      pb$prior_unfloored,
      pu$prior
    )
  )
  pg <- bq(s, "panel_unrounded_moments", "panel_gene")
  add(
    sprintf("(3) %s panel prior, panel-gene bootstrap", s),
    pu$prior,
    pg[1],
    pg[2],
    sprintf("bootstrap of the panel genes, %d resamples", pg[3]),
    source = "bootstrap_replicates.parquet panel_gene",
    tag = "estimate"
  )
  gs <- bq(s, "genome_all_targets", "sample")
  gt <- bq(s, "genome_all_targets", "transcript")
  gg <- gq(s, "prior genome_all_targets")
  add(
    sprintf("(3) %s genome-wide prior (quants_sub, 6 runs)", s),
    g$prior,
    gg$ci_lo,
    gg$ci_hi,
    sprintf("%s, %d resamples", gci, gg$n_boot),
    tag = "verified (point); estimate (CI)",
    source = "priors.parquet genome_all_targets; gene_bootstrap.parquet",
    note = sprintf(
      paste(
        "transcript bootstrap %.4f-%.4f (lower bound only);",
        "percentile sample bootstrap %.4f-%.4f (biased, see header); unfloored %.4f"
      ),
      gt[1],
      gt[2],
      gs[1],
      gs[2],
      g$prior_unfloored
    )
  )
  rs <- bq(s, "ratio_panel_over_genome", "sample")
  jr <- jq(s, "ratio panel (unrounded moments) / genome-wide prior")
  jg <- jq(s, "genome-wide prior")
  add(
    sprintf("(3) %s genome-wide prior, jackknife interval", s),
    g$prior,
    jg$ci_lo,
    jg$ci_hi,
    paste0(jci, " (6 runs; rough: each fit leaves out one of 6 samples)"),
    tag = "estimate (rough)",
    source = "jackknife.parquet"
  )
  add(
    sprintf("(3) HEADLINE %s ratio panel prior / genome-wide prior", s),
    pb$prior_snapshot_function / g$prior,
    jr$ci_lo,
    jr$ci_hi,
    paste0(
      jci,
      " (6 runs; rough: each fit leaves out one of 6 samples), of the",
      " unrounded-panel / genome ratio"
    ),
    tag = "verified (point); estimate (rough CI)",
    source = "priors.parquet; jackknife.parquet",
    note = sprintf(
      "panel-gene bootstrap of the panel prior %.3f-%.3f; percentile sample bootstrap of the ratio %.3f-%.3f (biased)",
      pg[1],
      pg[2],
      rs[1],
      rs[2]
    )
  )
}

# (4) Random panels
pctl <- function(x, v) mean(x <= v)
for (d in unique(panels$dataset)) {
  for (ty in setdiff(
    unique(panels$panel_type[panels$dataset == d]),
    "knockdown_panel_genes"
  )) {
    x <- panels$prior[panels$dataset == d & panels$panel_type == ty]
    gw <- if (d == "placenta") {
      pr(d, "a_all_targets")$prior
    } else {
      pr(d, "genome_all_targets")$prior
    }
    add(
      sprintf("(4) %s %s panel priors, median", d, ty),
      stats::median(x),
      stats::quantile(x, 0.025, names = FALSE),
      stats::quantile(x, 0.975, names = FALSE),
      sprintf("2.5-97.5%% of %d panels", length(x)),
      tag = "estimate",
      source = "panels_prior.parquet",
      note = paste0(
        sprintf(
          "genome-wide %.4f; median panel / genome %.3f; panel seed %s",
          gw,
          stats::median(x) / gw,
          paste(
            unique(panels$panel_seed[
              panels$dataset == d & panels$panel_type == ty
            ]),
            collapse = ","
          )
        ),
        if (ty == "matched_39") {
          paste(
            "; the 20 top-DE genes of the matched panels come from a DE fit",
            "with all 17 controls, so both members of each technical A/B pair",
            "count (affects only which genes are picked)"
          )
        } else {
          ""
        }
      )
    )
    k <- sum(abs(x / gw - 1) <= 0.10)
    w <- wilson(k, length(x))
    add(
      sprintf("(4) %s %s P(panel prior within +/-10%% of genome-wide)", d, ty),
      k / length(x),
      w[1],
      w[2],
      "Wilson 95% CI",
      tag = "estimate",
      source = "panels_prior.parquet"
    )
    gf <- fprior(d, 1)
    k1 <- sum(abs(x / gf - 1) <= 0.10)
    w1 <- wilson(k1, length(x))
    add(
      sprintf(
        "(4) %s %s P(panel prior within +/-10%% of the genome-wide prior under the panel filter, mean NumReads >= 1)",
        d,
        ty
      ),
      k1 / length(x),
      w1[1],
      w1[2],
      "Wilson 95% CI",
      tag = "estimate",
      source = "panels_prior.parquet; filter_priors.parquet",
      note = sprintf(
        "genome-wide prior with mean NumReads >= 1: %.4f; median panel / it %.3f",
        gf,
        stats::median(x) / gf
      )
    )
    actual <- if (d == "placenta") {
      panels$prior[
        panels$dataset == d & panels$panel_type == "knockdown_panel_genes"
      ]
    } else {
      pr(d, "panel_bundled_draws")$prior_snapshot_function
    }
    add(
      sprintf("(4) %s %s percentile of the actual panel prior", d, ty),
      pctl(x, actual),
      tag = "estimate",
      source = "panels_prior.parquet",
      note = sprintf(
        "actual panel prior %.4f (%s)",
        actual,
        if (d == "placenta") {
          "the knockdown panel genes on placenta"
        } else {
          "bundled panel"
        }
      )
    )
  }
}

# (4b) Decomposition: expression filter vs gene selection
for (i in seq_len(nrow(fpr))) {
  r <- fpr[i, ]
  add(
    sprintf("(4b) %s genome-wide prior, %s", r$dataset, r$rule),
    r$prior,
    source = "filter_priors.parquet",
    note = sprintf(
      "median phi %.4f; median df %g; %d transcripts with df > 0",
      r$median_phi,
      r$median_df,
      r$n_transcripts
    )
  )
}
for (s in c("GSE129540", "GSE75489")) {
  pb <- pr(s, "panel_bundled_draws")$prior_snapshot_function
  g0 <- fprior(s, 0)
  g1 <- fprior(s, 1)
  add(
    sprintf(
      "(4b) %s panel/genome ratio split: filter part (genome >=1 / genome all)",
      s
    ),
    g1 / g0,
    source = "filter_priors.parquet"
  )
  add(
    sprintf(
      "(4b) %s panel/genome ratio split: gene-selection part (panel / genome >=1)",
      s
    ),
    pb / g1,
    source = "priors.parquet; filter_priors.parquet"
  )
}

# (5) Draw count
r100 <- pr("GSE129540_g100_runs", "g100_100_draws")$prior
r30 <- pr("GSE129540_g100_runs", "g100_first_30_draws")$prior
rsub <- pr("GSE129540_g100_runs", "quants_sub_30_draws")$prior
c1 <- bq("GSE129540_g100_runs", "ratio_100_over_first30", "transcript")
c2 <- bq("GSE129540_g100_runs", "ratio_100_over_sub30", "transcript")
g1 <- gq("GSE129540_g100_runs", "ratio prior 100 / first 30 draws")
g2 <- gq("GSE129540_g100_runs", "ratio prior 100 draws / quants_sub 30 draws")
add(
  "(5) prior ratio 100 / first 30 draws (3 g100 runs)",
  r100 / r30,
  g1$ci_lo,
  g1$ci_hi,
  sprintf("%s, %d resamples", gci, g1$n_boot),
  tag = "verified (point); estimate (CI)",
  source = "priors.parquet GSE129540_g100_runs; gene_bootstrap.parquet",
  note = sprintf(
    paste(
      "priors %.4f / %.4f; median per-transcript phi ratio %.4f;",
      "transcript bootstrap %.4f-%.4f (lower bound only)"
    ),
    r100,
    r30,
    dct$median_transcript_phi_ratio[dct$comparison == "100_over_first30"],
    c1[1],
    c1[2]
  )
)
add(
  "(5) prior ratio 100 draws / independent 30-draw runs (quants_sub)",
  r100 / rsub,
  g2$ci_lo,
  g2$ci_hi,
  sprintf("%s, %d resamples", gci, g2$n_boot),
  tag = "verified (point); estimate (CI)",
  source = "priors.parquet; gene_bootstrap.parquet",
  note = sprintf(
    paste(
      "priors %.4f / %.4f; median per-transcript phi ratio %.4f;",
      "transcript bootstrap %.4f-%.4f (lower bound only)"
    ),
    r100,
    rsub,
    dct$median_transcript_phi_ratio[dct$comparison == "100_over_sub30"],
    c2[1],
    c2[2]
  )
)
p15 <- pr("placenta", "draws_first_15")$prior
c3 <- bq("placenta", "ratio_15_over_30", "sample_block")
j15 <- jq("placenta", "ratio prior first 15 / 30 draws")
add(
  "(5) placenta prior ratio first 15 / 30 draws",
  p15 / a$prior,
  j15$ci_lo,
  j15$ci_hi,
  jci,
  tag = "verified (point); estimate (CI)",
  source = "priors.parquet draws_first_15; jackknife.parquet",
  note = sprintf("percentile sample bootstrap %.4f-%.4f (biased)", c3[1], c3[2])
)

# (6) Prior df
jd <- jq("placenta", "fitFDistRobustly df2 (prior df)")
js <- jq("placenta", "fitFDistRobustly scale")
fnr <- fdf$value[fdf$statistic == "nonrobust_df2"]
add(
  "(6) prior df (df2), moment fit (limma::fitFDistRobustly df2 = fitFDist df2), placenta all targets",
  jd$point,
  jd$ci_lo,
  jd$ci_hi,
  jci,
  tag = "estimate",
  source = "jackknife.parquet; prior_df_fit.parquet",
  note = sprintf(
    paste(
      "edgeR/isoformic fix the prior df at 3; the robust fit returns the",
      "non-robust df2 (%.6f vs %.6f), only its scale differs; the F model of",
      "Pearson phi over Gibbs draws is approximate, so 3.30 vs 3 is an",
      "empirical check, not a test"
    ),
    jd$point,
    fnr
  )
)
add(
  "(6) limma::fitFDistRobustly scale, placenta all targets",
  js$point,
  js$ci_lo,
  js$ci_hi,
  jci,
  tag = "estimate",
  source = "jackknife.parquet; prior_df_fit.parquet",
  note = sprintf("the package prior (median rule) is %.4f", a$prior)
)
for (st in c(
  "robust_df2",
  "robust_scale",
  "full_df_df2",
  "full_df_scale",
  "nonrobust_df2"
)) {
  f <- fdf[fdf$statistic == st, ]
  add(
    paste0("(6) limma::fitFDist ", st, " (placenta, all targets), bootstrap"),
    f$value,
    f$boot_lo,
    f$boot_hi,
    sprintf(
      "percentile sample bootstrap (blocks), %d resamples; biased, see header",
      f$n_boot
    ),
    source = "prior_df_fit.parquet",
    tag = "estimate",
    note = "full_df_*: transcripts with a nonzero draw mean in all 25 runs; the bootstrap changes that set"
  )
}

# (7) Impact
for (i in seq_len(nrow(imp_phi))) {
  r <- imp_phi[i, ]
  add(
    paste0("(7) phi change: ", r$context),
    r$median_ratio,
    r$min_ratio,
    r$max_ratio,
    "min-max over transcripts",
    source = "impact_phi.parquet",
    note = sprintf(
      "priors %.4f -> %.4f; %d of %d transcripts change > 10%%, %d > 1%%",
      r$prior_reference,
      r$prior_alternative,
      r$n_changed_gt_10pct,
      r$n_transcripts,
      r$n_changed_gt_1pct
    )
  )
}
for (i in seq_len(nrow(imp_calls))) {
  r <- imp_calls[i, ]
  add(
    paste0(
      "(7) calls (",
      r$tested_set,
      " tested set): ",
      r$context,
      " | ",
      r$comparison
    ),
    r$n_shared,
    source = "impact_calls.parquet",
    note = sprintf(
      paste(
        "value = shared calls; reference %d calls of %d tested; alternative %d",
        "of %d; tested by both %d; lost %d (%d not tested by the alternative);",
        "gained %d (%d not tested by the reference)%s"
      ),
      r$n_calls_reference,
      r$n_tested_reference,
      r$n_calls_alternative,
      r$n_tested_alternative,
      r$n_tested_both,
      r$n_lost,
      r$n_lost_untested,
      r$n_gained,
      r$n_gained_untested,
      if (r$tested_set == "separate") {
        "; separate filterByExpr and BH per method: see the common row"
      } else {
        "; same transcripts and one BH correction for both methods"
      }
    )
  )
}

for (i in seq_len(nrow(byd))) {
  r <- byd[i, ]
  add(
    sprintf(
      "(7b) placenta phi change by subset: genome prior -> %s, %s",
      r$alternative,
      r$subset
    ),
    r$n_changed_gt_10pct,
    source = "impact_phi_by_df.parquet",
    note = sprintf(
      "of %d transcripts; > 1%%: %d; median ratio %.4f; max ratio %.4f; priors %.4f -> %.4f",
      r$n_transcripts,
      r$n_changed_gt_1pct,
      r$median_ratio,
      r$max_ratio,
      r$prior_genome,
      r$prior_alternative
    )
  )
}

# (8) catchSalmon
for (i in seq_len(nrow(cs))) {
  add(
    paste0("(8) catchSalmon max relative difference of phi, ", cs$dataset[i]),
    cs$max_rel_diff_phi[i],
    source = "catchsalmon_check.parquet",
    note = sprintf(
      "prior: package %.10f, catchSalmon %.10f, relative difference %.2e",
      cs$prior_package[i],
      cs$prior_catchsalmon[i],
      cs$rel_diff_prior[i]
    )
  )
}
m <- pr("mouse_bnst", "genome_all_targets")
if (nrow(m)) {
  mb <- bq("mouse_bnst", "genome_all_targets", "sample")
  add(
    "(9) mouse BNST genome-wide prior (19 runs)",
    m$prior,
    NA,
    NA,
    "",
    source = "priors.parquet mouse_bnst",
    note = sprintf(
      "%s; percentile sample bootstrap %.4f-%.4f (biased upward, see header)",
      m$note,
      mb[1],
      mb[2]
    )
  )
}
fr <- table(fold$rule)
add(
  "folding map (31-mer sensitivity): non-primary targets folded",
  sum(fr[c("duplicate_cluster", "kmer_match")], na.rm = TRUE),
  source = "fold_map.parquet",
  note = paste(names(fr), fr, sep = "=", collapse = "; ")
)
add(
  "folding map (S0, pre-registered): non-primary targets folded",
  sum(fold$fold_id_s0 != fold$transcript_id),
  source = "fold_map.parquet fold_id_s0",
  note = sprintf("%d non-primary targets", sum(!fold$primary))
)
if ("agree_with_s0" %in% names(fold)) {
  np <- !fold$primary
  add(
    "folding map: agreement with the S0 store map (non-primary targets)",
    mean(fold$agree_with_s0[np]),
    source = "fold_map.parquet"
  )
}

# Seeds, read from the tables
sd_b <- unique(boot[, c("dataset", "variant", "kind", "boot_seed")])
for (i in seq_len(nrow(sd_b))) {
  add(
    sprintf(
      "seed: bootstrap_replicates %s %s %s",
      sd_b$dataset[i],
      sd_b$variant[i],
      sd_b$kind[i]
    ),
    sd_b$boot_seed[i],
    source = "bootstrap_replicates.parquet boot_seed"
  )
}
sd_g <- unique(gb[, c("dataset", "statistic", "boot_seed")])
for (i in seq_len(nrow(sd_g))) {
  add(
    sprintf("seed: gene_bootstrap %s %s", sd_g$dataset[i], sd_g$statistic[i]),
    sd_g$boot_seed[i],
    source = "gene_bootstrap.parquet boot_seed"
  )
}
sd_p <- unique(panels[
  !is.na(panels$panel_seed),
  c("dataset", "panel_type", "panel_seed")
])
for (i in seq_len(nrow(sd_p))) {
  add(
    sprintf("seed: panels %s %s", sd_p$dataset[i], sd_p$panel_type[i]),
    sd_p$panel_seed[i],
    source = "panels_prior.parquet panel_seed"
  )
}
for (src in c("splits_prior", "prior_df_boot")) {
  add(
    paste0("seed: ", src),
    unique(rd(src)$seed),
    source = paste0(src, ".parquet seed")
  )
}

summ <- dplyr::bind_rows(rows)
ue_write(summ, file.path(out, "summary.parquet"), ue)
utils::write.table(
  summ,
  file.path(out, "summary.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
options(width = 200)
print(
  summ[, c("item", "value", "ci_lo", "ci_hi", "tag")],
  digits = 4,
  right = FALSE
)

# Proposal for the maintainer (text only; no package file is changed) -------
v <- function(item) summ$value[summ$item == item][1]
pl1 <- fprior("placenta", 1)
pl100 <- fprior("placenta", 100)
# Calls on the common tested set (51_od_prior_impact.R)
crow <- function(ctx, cmp) {
  r <- imp_calls[
    imp_calls$tested_set == "common" &
      startsWith(imp_calls$context, ctx) &
      imp_calls$comparison == cmp,
  ]
  stopifnot(nrow(r) == 1L)
  r
}
kd_common <- imp_calls[
  imp_calls$tested_set == "common" & grepl("^GSE", imp_calls$context),
]
usage_rows <- imp_calls[
  grepl("^GSE", imp_calls$context) &
    grepl("test_usage", imp_calls$context, fixed = TRUE),
]
pl_med <- crow(
  "placenta 8 PE vs 15 control units: DTE edgeR QL, divided counts",
  "genome prior -> panel_median"
)
pl_raw <- crow(
  "placenta 8 PE vs 15 control units: DTE edgeR QL",
  "raw counts -> divided counts, genome-wide prior"
)
ratio129 <- pr("GSE129540", "panel_bundled_draws")$prior /
  pr("GSE129540", "genome_all_targets")$prior
ratio754 <- pr("GSE75489", "panel_bundled_draws")$prior /
  pr("GSE75489", "genome_all_targets")$prior
filt129 <- fprior("GSE129540", 1) / pr("GSE129540", "genome_all_targets")$prior
filt754 <- fprior("GSE75489", 1) / pr("GSE75489", "genome_all_targets")$prior
proposal <- c(
  "WP-D proposal (not applied; for the maintainer)",
  "",
  sprintf(
    paste(
      "1. The bundled knockdown draws give a panel prior of %.3f (GSE129540) and",
      "%.3f (GSE75489): %.2f and %.2f times the genome-wide prior of the same",
      "runs under the catchSalmon convention (every target; %.3f and %.3f).",
      "On a log scale %.0f%% and %.0f%% of that gap is the panel's own expression",
      "filter (mean NumReads >= 1; genome-wide prior under that filter %.3f and",
      "%.3f, factors %.2f and %.2f); the rest (factors %.2f and %.2f) is gene",
      "selection (top-DE and ESRP-target, multi-isoform genes). The prior also",
      "falls as samples are added, so it is a convention, not a property of the",
      "data. The difference has a negligible effect on results: moderated phi",
      "changes by at most %.1f%%; on the common tested set %s, and test_usage",
      "calls change by at most %d transcripts, because every",
      "panel transcript has df >= 145 (6 runs x 29) >> prior df 3."
    ),
    pr("GSE129540", "panel_bundled_draws")$prior,
    pr("GSE75489", "panel_bundled_draws")$prior,
    ratio129,
    ratio754,
    pr("GSE129540", "genome_all_targets")$prior,
    pr("GSE75489", "genome_all_targets")$prior,
    100 * log(filt129) / log(ratio129),
    100 * log(filt754) / log(ratio754),
    fprior("GSE129540", 1),
    fprior("GSE75489", 1),
    filt129,
    filt754,
    ratio129 / filt129,
    ratio754 / filt754,
    100 *
      (1 -
        min(imp_phi$min_ratio[
          grepl("^GSE", imp_phi$context) &
            !grepl("scale_counts", imp_phi$context, fixed = TRUE)
        ])),
    if (max(kd_common$n_lost + kd_common$n_gained) == 0) {
      "no DTE call changes"
    } else {
      sprintf(
        "DTE calls change by at most %d transcripts",
        max(kd_common$n_lost + kd_common$n_gained)
      )
    },
    max(usage_rows$n_lost + usage_rows$n_gained)
  ),
  paste(
    "   Proposal: in data-raw/build_knockdown_example.R keep the genome-wide",
    "   moments (or only the scalar genome-wide prior) next to draws.parquet, so",
    "   that inferential_overdispersion() on the example reports the prior that",
    "   catchSalmon would give on the full runs; or state in the uncertainty",
    "   article that the example's prior comes from a 39-gene panel with an",
    "   expression filter and is larger than the all-targets prior of the full",
    "   runs. Low priority for the calls; relevant for any text that quotes the",
    "   prior."
  ),
  "",
  sprintf(
    paste(
      "2. Placenta (8 PE vs 15 control units): replacing the genome-wide prior",
      "(%.3f) by the median random-panel prior changes DTE calls from %d to %d",
      "on the %d transcripts both runs test (lost %d, gained %d). With separate",
      "filterByExpr runs the call changes are larger, but they are mostly",
      "transcripts that one run no longer tests, not changed decisions."
    ),
    imp_pl$prior[imp_pl$label == "genome"],
    pl_med$n_calls_reference,
    pl_med$n_calls_alternative,
    pl_med$n_tested_both,
    pl_med$n_lost,
    pl_med$n_gained
  ),
  "",
  sprintf(
    paste(
      "3. The prior depends on which transcripts enter the median. Placenta:",
      "%.3f with every target, %.3f with mean NumReads >= 1, %.3f with >= 100.",
      "The package already uses the moments of every target (the catchSalmon",
      "convention; checked to %.1e). Proposal: document this dependence and",
      "inform when the prior falls back to the stored (subset) replicates",
      "because the moments attribute is missing."
    ),
    a$prior,
    pl1,
    pl100,
    max(cs$max_rel_diff_phi)
  ),
  "",
  sprintf(
    paste(
      "4. Prior df: the moment fit of limma (fitFDistRobustly returns the same",
      "df2 as fitFDist) on the placenta raw overdispersion gives df0 = %.2f",
      "(jackknife 95%% interval %.2f-%.2f), close to edgeR's fixed 3. The F",
      "model for Pearson phi over Gibbs draws is approximate, so this is an",
      "empirical check, not a test; no practical effect, since every expressed",
      "transcript has df >= 29. No change proposed."
    ),
    jd$point,
    jd$ci_lo,
    jd$ci_hi
  ),
  "",
  sprintf(
    paste(
      "Context (not a WP-D result): raw-count vs divided-count DTE on placenta,",
      "on the %d transcripts both filters keep, one BH correction: %d vs %d",
      "calls, %d shared (lost %d, gained %d): divided counts call almost every",
      "raw-count call plus more; whether the extra calls are true or false is",
      "the question of the null and replication analyses (WP-A), not of this",
      "table. The earlier statement that the two disagree on about a third of",
      "their calls compared different tested sets (74k vs 33k transcripts, each",
      "with its own BH correction) and is withdrawn."
    ),
    pl_raw$n_tested_both,
    pl_raw$n_calls_reference,
    pl_raw$n_calls_alternative,
    pl_raw$n_shared,
    pl_raw$n_lost,
    pl_raw$n_gained
  )
)
writeLines(proposal, file.path(out, "proposal.txt"))
cat(proposal, sep = "\n")
