# 90_article_tables.R: the small result tables that the uncertainty evidence
# article reads.
#
# Purpose: copy, without re-estimating anything, every number that the
# article draft (and REPORT.md) states into a few small Parquet tables, so
# that the article computes each number inline from a table and no number is
# typed by hand. Each value keeps its 95% CI, the kind of CI and its source
# table. Every key is unique (the script stops otherwise).
#
# Inputs (all under $ISOFORMIC_UE_OUT; produced by the scripts named):
#   prereg.yaml                                   01_preregister.R
#   S0/stratum_sizes, S0/run_summary, store/strata_tx      10_placenta_prepare.R
#   null_placenta/summary/*.parquet               29_null_summarise.R (WP-A)
#   placenta_audit/must_report.parquet, replication/summary_flags.parquet,
#     replication/summary_sign_concordance.parquet,
#     story/flt1_isoforms.parquet, story/story_genes_*.parquet,
#     story/rule_selected_*.parquet               39, 31b, 32 (WP-B)
#   placenta_audit/prevalence/summary_family.parquet      30 (WP-B)
#   placenta_audit/post_hoc/*.parquet             33 (post hoc, A9, A21), 34 (A11),
#                                                 92 (A24)
#   null_mouse/results_mouse_null.parquet, control_batches.parquet,
#     wpc_summary/wpc_rule_RA.parquet, rhalf/replication/summary_*.parquet
#                                                 41, 44, 49 (WP-C)
#   od_prior/summary.parquet, od_prior/impact_calls.parquet   59, 51 (WP-D)
#
# Outputs: $ISOFORMIC_UE_OUT/article_tables/
#   evidence_numbers.parquet  long table: key, estimate, ci_lo, ci_hi,
#                             ci_type, n, source (one row per number)
#   flt1_isoforms.parquet     the FLT1 isoforms of the story (one row each)
#   story_genes.parquet       the other story genes of the plan (ENG, LEP and
#                             the fragile- and robust-rule picks), by gene and
#                             by isoform
#   replication_flags.parquet held-out replication by flag (placenta
#                             half-splits and knockdown read halves)
#   session_info.parquet
# Proposed destination in the package: vignettes/articles/data/uncertainty_evidence/
# (the four result tables only).
#
# Env vars: ISOFORMIC_UE_OUT (required), ISOFORMIC_UE_PKG (ue_common.R).
# Resources (estimate): about 1 min, < 1 GB.
# Run from the repository root, inside a memory cap of about 4 GB (see
# README.md):
#   pixi run --frozen Rscript data-raw/uncertainty_evidence/90_article_tables.R

source("data-raw/uncertainty_evidence/ue_common.R")
ue <- ue_setup("90_article_tables")
# jarl-ignore library_call: attached after ue_setup() on purpose, so that its
# names mask those of the package snapshot, as in the runs of the evaluation
library(data.table)
out_root <- ue_out_root()
out_dir <- file.path(out_root, "article_tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

rd <- function(rel) {
  path <- file.path(out_root, rel)
  if (!file.exists(path)) {
    stop("missing input: ", path, call. = FALSE)
  }
  as.data.table(read_parquet(path))
}
slug <- function(x) {
  x <- gsub("[^A-Za-z0-9.+-]+", "_", x)
  gsub("^_+|_+$", "", x)
}
# Parse the first "[lo, hi]" of a CI string
ci_parse <- function(s) {
  m <- regmatches(s, regexpr("\\[[^]]*\\]", s))
  out <- matrix(NA_real_, length(s), 2L)
  has <- !is.na(s) & grepl("\\[[^]]*\\]", s)
  if (any(has)) {
    v <- strsplit(gsub("\\[|\\]", "", m), ",")
    out[has, ] <- t(vapply(
      v,
      function(z) suppressWarnings(as.numeric(z)),
      numeric(2)
    ))
  }
  out
}
rows <- list()
add <- function(
  key,
  estimate,
  ci_lo = NA_real_,
  ci_hi = NA_real_,
  ci_type = NA_character_,
  n = NA_real_,
  source
) {
  # data.frame(), not data.table(): data.table() reads `key =` as its own
  # argument
  # jarl-ignore undesirable_operator: appends to a result list of the enclosing scope
  rows[[length(rows) + 1L]] <<- data.frame(
    key = as.character(key),
    estimate = as.numeric(estimate),
    ci_lo = as.numeric(ci_lo),
    ci_hi = as.numeric(ci_hi),
    ci_type = as.character(ci_type),
    n = as.numeric(n),
    source = as.character(source),
    stringsAsFactors = FALSE
  )
}
# WP-A tables share the jackknife columns: estimate, ci_lo, ci_hi, ci_type,
# n_splits
add_jk <- function(d, prefix, keycols, source) {
  if (nrow(d) == 0L) {
    stop("no rows for ", prefix, call. = FALSE)
  }
  k <- do.call(
    paste,
    c(list(prefix), lapply(keycols, function(cc) slug(d[[cc]])), sep = ".")
  )
  add(
    k,
    d$estimate,
    d$ci_lo,
    d$ci_hi,
    d$ci_type,
    if ("n_splits" %in% names(d)) d$n_splits else NA_real_,
    source
  )
}

# Design facts (prereg.yaml) ---------------------------------------------------------
pre <- ue$prereg
add("design.n_pe", length(pre$data$preeclampsia), source = "prereg.yaml data")
add("design.n_control", length(pre$data$control), source = "prereg.yaml data")
add(
  "design.n_control_units",
  length(pre$data$control) - length(pre$data$technical_pairs),
  source = "prereg.yaml control_units"
)
add(
  "design.null_splits_per_size",
  pre$splits$null_splits_per_size,
  source = "prereg.yaml splits"
)
add(
  "design.half_splits",
  pre$splits$replication_half_splits,
  source = "prereg.yaml splits"
)
add(
  "design.directed_pairs",
  pre$splits$replication_directed_total,
  source = "prereg.yaml splits"
)
add(
  "design.n_samples",
  length(pre$data$preeclampsia) + length(pre$data$control),
  source = "prereg.yaml data"
)

# S0: universe and strata -------------------------------------------------------------
ss <- rd("S0/stratum_sizes.parquet")
add(
  paste("S0", slug(ss$stratum), slug(ss$level), sep = "."),
  ss$n,
  source = "S0/stratum_sizes.parquet"
)
rs <- rd("S0/run_summary.parquet")
add(
  paste("S0.run", slug(rs$metric), sep = "."),
  rs$value,
  source = "S0/run_summary.parquet"
)
cp <- rd("store/strata_cutpoints.parquet")
add(
  paste("S0.cut", slug(cp$what), slug(cp$cut), sep = "."),
  cp$value,
  source = "store/strata_cutpoints.parquet"
)
st <- rd("store/strata_tx.parquet")
stu <- st[in_u %in% TRUE]
add(
  "S0.tertile_agreement_15_vs_30",
  mean(stu$infrv_tertile == stu$infrv_tertile_15),
  n = nrow(stu),
  source = "store/strata_tx.parquet infrv_tertile vs infrv_tertile_15 (in_u)"
)
mi_files <- file.path(
  pre$data$placenta_quants,
  c(pre$data$preeclampsia, pre$data$control),
  "aux_info",
  "meta_info.json"
)
n_boot <- vapply(
  mi_files,
  function(f) as.numeric(jsonlite::read_json(f)$num_bootstraps),
  numeric(1)
)
if (length(unique(n_boot)) != 1L) {
  stop("placenta runs differ in num_bootstraps", call. = FALSE)
}
add(
  "S0.n_draws",
  n_boot[[1]],
  n = length(n_boot),
  source = "placenta quants aux_info/meta_info.json num_bootstraps (all runs equal)"
)

# WP-A: placenta null -----------------------------------------------------------------
A <- "null_placenta/summary/"
g <- rd(paste0(A, "q3_glm_rate_ratios.parquet"))
g[,
  term_s := fifelse(
    term == "infrv_tertilehigh",
    "top",
    fifelse(term == "infrv_tertilemid", "mid", term)
  )
]
add_jk(
  g,
  "A.adj_rr",
  c("design_set", "method", "model", "term_s"),
  paste0(A, "q3_glm_rate_ratios.parquet")
)
# the raw delete-one-unit jackknife next to the MC-corrected one (amendment
# A5): where the correction removes most of the variance, the raw interval
# is the conservative one, and the article quotes both
add_rawjk <- function(d, prefix, keycols, source) {
  d <- d[!is.na(jk_lo_raw)]
  k <- do.call(
    paste,
    c(list(prefix), lapply(keycols, function(cc) slug(d[[cc]])), sep = ".")
  )
  add(
    k,
    d$estimate,
    d$jk_lo_raw,
    d$jk_hi_raw,
    "delete-one-unit jackknife (raw)",
    d$n_splits,
    source
  )
}
add_rawjk(
  g,
  "A.adj_rr_rawjk",
  c("design_set", "method", "model", "term_s"),
  paste0(A, "q3_glm_rate_ratios.parquet")
)
h2 <- rd(paste0(A, "q3_h2_method_interaction.parquet"))
h2[, cmp := paste0(alternative, "_vs_", reference)]
# model "prereg" is the pre-registered covariate set (A6); "ident" adds
# non-identifiable (sensitivity); the M6_vs_M7 rows are a diagnostic only (A7)
add_jk(
  h2,
  "A.rr_ratio",
  c("design_set", "model", "cmp"),
  paste0(A, "q3_h2_method_interaction.parquet")
)
if ("jk_lo_raw" %in% names(h2)) {
  add_rawjk(
    h2,
    "A.rr_ratio_rawjk",
    c("design_set", "model", "cmp"),
    paste0(A, "q3_h2_method_interaction.parquet")
  )
}
cm <- rd(paste0(A, "q3_cmh.parquet"))
add_jk(cm, "A.cmh", c("design_set", "method"), paste0(A, "q3_cmh.parquet"))
q1 <- rd(paste0(A, "q1_q2_calls_and_inflation.parquet"))
q1 <- q1[
  (stratum == "all" | stratum == "infrv_tertile") &
    quantity %in%
      c(
        "mean_calls_q05",
        "p_any_q05",
        "inflation_p01",
        "mean_tested",
        "mean_calls_p01"
      )
]
q1[, lev := fifelse(stratum == "all", "all", paste0("tertile_", level))]
add_jk(
  q1,
  "A",
  c("quantity", "design_set", "method", "set", "lev"),
  paste0(A, "q1_q2_calls_and_inflation.parquet")
)
q2 <- rd(paste0(A, "q2_tertile_rate_ratio.parquet"))[threshold == "n_p01"]
add_jk(
  q2,
  "A.crude_rr",
  c("design_set", "method"),
  paste0(A, "q2_tertile_rate_ratio.parquet")
)
q4 <- rd(paste0(A, "q4_enrichment_top_tertile.parquet"))[threshold == "n_p01"]
add_jk(
  q4,
  "A.enrich_top",
  c("design_set", "method"),
  paste0(A, "q4_enrichment_top_tertile.parquet")
)
q5 <- rd(paste0(A, "q5_m9_vs_m10_removal.parquet"))[threshold == "n_p01"]
add_jk(
  q5,
  "A.m9_vs_m10",
  c("design_set", "set", "m9_rule", "quantity"),
  paste0(A, "q5_m9_vs_m10_removal.parquet")
)
q7 <- rd(paste0(A, "q7_dtu.parquet"))[stratum == "all"]
q7[, sc := fifelse(as.logical(scaled), "scaled", "raw")]
add_jk(
  q7,
  "A.dtu",
  c("backend", "design", "sc", "quantity"),
  paste0(A, "q7_dtu.parquet")
)
nmd <- rd(paste0(A, "q7_nmd.parquet"))
add_jk(nmd, "A.nmd", c("arm", "quantity"), paste0(A, "q7_nmd.parquet"))
sw <- rd(paste0(A, "q7_switches.parquet"))
add_jk(sw, "A.switch", c("arm", "quantity"), paste0(A, "q7_switches.parquet"))
th <- rd(paste0(A, "q8_thinning_affine.parquet"))
add_jk(
  th,
  "A.thin",
  c("design", "method", "subset", "call_rule", "quantity"),
  paste0(A, "q8_thinning_affine.parquet")
)
th1 <- th[quantity == "fdr_pooled"]
for (cc in c("total_calls", "total_false", "total_true")) {
  add(
    paste(
      "A.thin",
      slug(th1$design),
      slug(th1$method),
      slug(th1$subset),
      slug(th1$call_rule),
      cc,
      sep = "."
    ),
    th1[[cc]],
    n = th1$n_splits,
    source = paste0(A, "q8_thinning_affine.parquet ", cc)
  )
}
q10 <- rd(paste0(A, "q10_tail_sensitivity.parquet"))
for (cc in c(
  "mean_tested",
  "mean_q05",
  "max_q05",
  "share_from_max_split",
  "mean_q05_without_max_split",
  "p_any_q05"
)) {
  add(
    paste(
      "A.tail",
      slug(q10$design_set),
      slug(q10$method),
      slug(q10$set),
      cc,
      sep = "."
    ),
    q10[[cc]],
    n = q10$n_splits,
    source = paste0(A, "q10_tail_sensitivity.parquet ", cc)
  )
}
q11 <- rd(paste0(A, "q11_m9_class_vs_null_calls.parquet"))
add_jk(
  q11,
  "A.m9class",
  c("design_set", "quantity"),
  paste0(A, "q11_m9_class_vs_null_calls.parquet")
)
rb <- rd(paste0(A, "robustness.parquet"))[
  quantity == "inflation_p01" & method == "M1" & set == "own"
]
add_jk(
  rb,
  "A.robust",
  c("analysis", "method", "stratum", "level"),
  paste0(A, "robustness.parquet")
)
im <- rd(paste0(A, "robustness_imbalance_glm.parquet"))[
  method == "M1" & grepl("^abs", term)
]
add(
  paste("A.imbalance_rr.M1", slug(im$term), sep = "."),
  im$rr,
  im$jk_lo_raw,
  im$jk_hi_raw,
  "delete-one-unit jackknife (raw)",
  im$n_splits,
  paste0(A, "robustness_imbalance_glm.parquet")
)
mi <- rd(paste0(A, "robustness_tertile_ratio_by_mapping_imbalance.parquet"))
add_jk(
  mi,
  "A.tertile_rr_by_mapping",
  c("design_set", "method", "mapping_imbalance"),
  paste0(A, "robustness_tertile_ratio_by_mapping_imbalance.parquet")
)
vr <- rd(paste0(A, "verdict_R_A.parquet"))
add(
  paste("A.verdict", slug(vr$criterion), sep = "."),
  as.numeric(as.logical(vr$value)),
  source = paste0(
    A,
    "verdict_R_A.parquet (1 TRUE, 0 FALSE, NA not informative)"
  )
)

# WP-B: placenta audit ----------------------------------------------------------------
mr <- rd("placenta_audit/must_report.parquet")
ci <- ci_parse(mr$ci)
add(
  paste("B", slug(mr$item), slug(mr$statistic), sep = "."),
  mr$value,
  ci[, 1],
  ci[, 2],
  sub(" *\\[.*", "", mr$ci),
  source = paste0("placenta_audit/must_report.parquet <- ", mr$source)
)

# The planned WP-B bootstrap resamples the half-splits, then the genes within
# each directed pair (amendment A15)
wpb_ci_type <- paste(
  "two-level bootstrap (half-splits, then genes within each directed pair;",
  "A15)"
)
pf0 <- rd("placenta_audit/replication/summary_flags.parquet")
for (cc in c(
  "rate_flagged",
  "rate_unflagged",
  "deficit",
  "matched_deficit",
  "adjusted_or"
)) {
  cic <- paste0(cc, "_ci")
  ci <- if (cic %in% names(pf0)) {
    ci_parse(pf0[[cic]])
  } else {
    matrix(NA_real_, nrow(pf0), 2L)
  }
  add(
    paste("B.flags", slug(pf0$flag), cc, sep = "."),
    pf0[[cc]],
    ci[, 1],
    ci[, 2],
    if (cic %in% names(pf0)) {
      wpb_ci_type
    } else {
      NA_character_
    },
    pf0$n_calls,
    "placenta_audit/replication/summary_flags.parquet"
  )
}
add(
  paste("B.flags", slug(pf0$flag), "pct_flagged", sep = "."),
  pf0$pct_flagged,
  n = pf0$n_calls,
  source = "placenta_audit/replication/summary_flags.parquet"
)
# The same flag on the full cohort (F1: transcript calls at q < 0.05 of the
# 8 PE against 15 control units comparison), next to its share among the
# calls of the discovery halves above
sfam <- rd("placenta_audit/prevalence/summary_family.parquet")[
  family == "F1" & statistic == "pct_sign_agree_lt_0.95"
]
ci <- ci_parse(sfam$ci_gene_cluster)
add(
  "B.prevalence.F1.pct_sign_agree_lt_0.95",
  sfam$estimate,
  ci[, 1],
  ci[, 2],
  "gene-cluster bootstrap",
  n = rd("placenta_audit/prevalence/summary_family.parquet")[
    family == "F1" & statistic == "n"
  ]$estimate,
  source = "placenta_audit/prevalence/summary_family.parquet <- 30_audit_prevalence.R"
)
scn <- rd("placenta_audit/replication/summary_sign_concordance.parquet")
ci <- ci_parse(scn$ci)
add(
  paste("B.signconc", slug(scn$flag), "adjusted_or", sep = "."),
  scn$adjusted_or_sign_concordance,
  ci[, 1],
  ci[, 2],
  wpb_ci_type,
  source = "placenta_audit/replication/summary_sign_concordance.parquet"
)
for (cc in c("sign_concordance_flagged", "sign_concordance_unflagged")) {
  add(
    paste("B.signconc", slug(scn$flag), cc, sep = "."),
    scn[[cc]],
    source = "placenta_audit/replication/summary_sign_concordance.parquet"
  )
}
# Post hoc checks (33_post_hoc_checks.R, amendments A9 and A21): long tables
# with key and estimate, no CI
post_hoc_amendment <- c(
  point_vs_draws = "A9",
  power_sensitivity = "A9",
  unstable_same_size_controls = "A21"
)
for (f in names(post_hoc_amendment)) {
  ph <- rd(file.path("placenta_audit", "post_hoc", paste0(f, ".parquet")))
  add(
    paste0("B.posthoc.", ph$key),
    ph$estimate,
    n = ph$n,
    ci_type = paste0(
      "post hoc, not pre-registered (",
      post_hoc_amendment[[f]],
      "); no CI"
    ),
    source = paste0(
      "placenta_audit/post_hoc/",
      f,
      ".parquet <- 33_post_hoc_checks.R"
    )
  )
}

# Post hoc gene-clustered bootstrap, pair fixed effects, per-pair rates and
# read-half fallback (34_cluster_bootstrap.R, amendment A11)
cb <- rd("placenta_audit/post_hoc/cluster_bootstrap.parquet")
add(
  paste0("B.cluster.", cb$key),
  cb$estimate,
  cb$ci_lo,
  cb$ci_hi,
  ifelse(
    is.na(cb$ci_type),
    "post hoc, not pre-registered (A11)",
    paste0("post hoc, not pre-registered (A11); ", cb$ci_type)
  ),
  cb$n,
  "placenta_audit/post_hoc/cluster_bootstrap.parquet <- 34_cluster_bootstrap.R"
)

# Post hoc check that the package flag is the evaluated flag
# (92_package_flag_check.R, amendment A24): counts, no CI
pc <- rd("placenta_audit/post_hoc/package_flag_check.parquet")
add(
  paste0("B.pkgcheck.", pc$key),
  pc$estimate,
  n = pc$n,
  ci_type = "post hoc check of the package code, not pre-registered (A24); no CI",
  source = "placenta_audit/post_hoc/package_flag_check.parquet <- 92_package_flag_check.R"
)

# WP-C: mouse null, read halves -------------------------------------------------------
mo <- rd("null_mouse/results_mouse_null.parquet")
add(
  paste("C", slug(mo$method), slug(mo$statistic), sep = "."),
  mo$estimate,
  mo$ci_lo,
  mo$ci_hi,
  mo$jackknife,
  mo$n_splits_full,
  "null_mouse/results_mouse_null.parquet"
)
ra <- rd("wpc_summary/wpc_rule_RA.parquet")
for (cc in c("m1_adj_rr", "m1_m2_diff", "m1_m6_diff", "m11_adj_rr")) {
  add(
    paste("C.RA", slug(ra$splits), cc, sep = "."),
    ra[[cc]],
    ra[[sub("(_adj_rr|_diff)$", "_ci_lo", cc)]],
    ra[[sub("(_adj_rr|_diff)$", "_ci_hi", cc)]],
    "exact jackknife",
    source = "wpc_summary/wpc_rule_RA.parquet"
  )
}
add(
  paste("C.RA", slug(ra$splits), "passes", sep = "."),
  as.numeric(as.logical(ra$R_A_passes)),
  source = "wpc_summary/wpc_rule_RA.parquet"
)
cb <- rd("null_mouse/control_batches.parquet")
cbs <- cb[,
  .(lo = min(frag_length_mean), hi = max(frag_length_mean), n = .N),
  by = batch
]
add(
  paste("C.batch", slug(cbs$batch), "frag_min", sep = "."),
  cbs$lo,
  n = cbs$n,
  source = "null_mouse/control_batches.parquet"
)
add(
  paste("C.batch", slug(cbs$batch), "frag_max", sep = "."),
  cbs$hi,
  n = cbs$n,
  source = "null_mouse/control_batches.parquet"
)
spl <- rd("null_mouse/splits.parquet")
if ("batch_balanced" %in% names(spl)) {
  full <- if ("design" %in% names(spl)) spl[grepl("full|5v5", design)] else spl
  add(
    "C.n_splits_batch_balanced",
    sum(full$batch_balanced %in% TRUE),
    source = "null_mouse/splits.parquet batch_balanced"
  )
  add(
    "C.n_splits_batch_aligned",
    sum(full$batch_balanced %in% FALSE),
    source = "null_mouse/splits.parquet batch_balanced"
  )
}
hf <- rd("rhalf/replication/summary_flags.parquet")
for (cc in c(
  "n_calls",
  "pct_flagged",
  "replication_overall",
  "n_untestable_v"
)) {
  add(
    paste(
      "R",
      slug(hf$halving),
      slug(hf$series),
      slug(hf$method),
      slug(hf$flag),
      cc,
      sep = "."
    ),
    hf[[cc]],
    source = "rhalf/replication/summary_flags.parquet"
  )
}
for (cc in c(
  "rate_flagged",
  "rate_unflagged",
  "deficit",
  "adjusted_or",
  "rate_flagged_testable",
  "rate_unflagged_testable",
  "deficit_testable"
)) {
  ci <- ci_parse(hf[[paste0(cc, "_ci")]])
  add(
    paste(
      "R",
      slug(hf$halving),
      slug(hf$series),
      slug(hf$method),
      slug(hf$flag),
      cc,
      sep = "."
    ),
    hf[[cc]],
    ci[, 1],
    ci[, 2],
    "bootstrap (44_rhalf_replicate.R)",
    source = "rhalf/replication/summary_flags.parquet"
  )
}
pl <- rd("rhalf/replication/summary_placebo.parquet")
for (cc in c("observed", "placebo_mean", "empirical_p_one_sided")) {
  add(
    paste(
      "R.placebo",
      slug(pl$halving),
      slug(pl$series),
      slug(pl$method),
      slug(pl$statistic),
      cc,
      sep = "."
    ),
    pl[[cc]],
    if (cc == "placebo_mean") pl[["placebo_2.5"]] else NA_real_,
    if (cc == "placebo_mean") pl[["placebo_97.5"]] else NA_real_,
    if (cc == "placebo_mean") "2.5-97.5% of placebos" else NA_character_,
    source = "rhalf/replication/summary_placebo.parquet"
  )
}
pm <- rd("rhalf/replication/summary_perm.parquet")
pms <- pm[,
  .(
    true = mean(replication[true_partition]),
    perm = mean(replication[!true_partition]),
    perm_min = min(replication[!true_partition]),
    perm_max = max(replication[!true_partition]),
    n_perm = sum(!true_partition)
  ),
  by = .(halving, series)
]
add(
  paste(
    "R.perm",
    slug(pms$halving),
    slug(pms$series),
    "replication_true",
    sep = "."
  ),
  pms$true,
  source = "rhalf/replication/summary_perm.parquet"
)
add(
  paste(
    "R.perm",
    slug(pms$halving),
    slug(pms$series),
    "replication_permuted",
    sep = "."
  ),
  pms$perm,
  pms$perm_min,
  pms$perm_max,
  "range over permutations (not a CI)",
  pms$n_perm,
  "rhalf/replication/summary_perm.parquet"
)

# WP-D: overdispersion prior ----------------------------------------------------------
od <- rd("od_prior/summary.parquet")
add(
  paste("D", slug(od$item), sep = "."),
  od$value,
  od$ci_lo,
  od$ci_hi,
  od$ci_method,
  source = paste0("od_prior/summary.parquet <- ", od$source)
)
ic <- rd("od_prior/impact_calls.parquet")
for (cc in c(
  "n_tested_both",
  "n_calls_reference",
  "n_calls_alternative",
  "n_lost",
  "n_gained",
  "n_lost_untested"
)) {
  add(
    paste(
      "D.impact",
      slug(ic$context),
      slug(ic$comparison),
      slug(ic$tested_set),
      cc,
      sep = "."
    ),
    ic[[cc]],
    source = "od_prior/impact_calls.parquet"
  )
}

ev <- rbindlist(rows)
dup <- ev$key[duplicated(ev$key)]
if (length(dup) > 0L) {
  stop(
    "duplicated keys: ",
    paste(utils::head(unique(dup), 20), collapse = ", "),
    call. = FALSE
  )
}
ue_write(ev, file.path(out_dir, "evidence_numbers.parquet"), ue)

# FLT1 isoforms (story) ----------------------------------------------------------------
fl <- rd("placenta_audit/story/flt1_isoforms.parquet")
keep <- c(
  "transcript_name",
  "transcript_type",
  "mean_count",
  "mean_inf_rv",
  "tested_edgeR",
  "naive_log2FC",
  "naive_ci_lo",
  "naive_ci_hi",
  "draw_refit_lfc_lo_2.5",
  "draw_refit_lfc_hi_97.5",
  "swish_log2FC",
  "swish_q",
  "naive_q",
  "confidence",
  "n_called_discovery",
  "replicated_when_called"
)
miss <- setdiff(keep, names(fl))
if (length(miss) > 0L) {
  stop("flt1_isoforms lacks: ", paste(miss, collapse = ", "), call. = FALSE)
}
fl <- fl[, ..keep]
setnames(
  fl,
  c("draw_refit_lfc_lo_2.5", "draw_refit_lfc_hi_97.5"),
  c("draw_refit_lo", "draw_refit_hi")
)
ue_write(fl, file.path(out_dir, "flt1_isoforms.parquet"), ue)

# The other story genes (prereg.yaml story_genes) ---------------------------------------
# ENG and LEP (named in the plan), and the genes picked by the fragile and
# robust rules (32_story.R). One row per gene (level "gene", the summed
# isoforms, plus the rule and its number of candidates) and one row per
# isoform (level "isoform"). The pre-registered fragile rule picked no gene;
# its row has gene_name NA and n_candidates 0.
iso_keep <- c(
  "gene_name",
  "transcript_name",
  "mean_count",
  "mean_inf_rv",
  "tested_edgeR",
  "naive_log2FC",
  "naive_q",
  "draw_refit_lfc_lo_2.5",
  "draw_refit_lfc_hi_97.5",
  "sign_agree",
  "confidence",
  "n_called_discovery",
  "replicated_when_called"
)
gene_keep <- c(
  "gene_name",
  "n_isoforms_U",
  "naive_log2FC",
  "naive_q",
  "draw_refit_lfc_lo_2.5",
  "draw_refit_lfc_hi_97.5",
  "sign_agree",
  "n_called_discovery",
  "replicated_when_called"
)
sg_iso <- rbind(
  rd("placenta_audit/story/story_genes_isoforms.parquet")[
    gene_name %in% c("ENG", "LEP"),
    ..iso_keep
  ],
  rd("placenta_audit/story/rule_selected_isoforms.parquet")[, ..iso_keep]
)
sg_gene <- rbind(
  rd("placenta_audit/story/story_genes_gene_level.parquet")[, ..gene_keep],
  rd("placenta_audit/story/rule_selected_gene_level.parquet")[, ..gene_keep]
)
rules <- rd("placenta_audit/story/rule_selected_genes.parquet")[,
  .(rule, gene_name, n_candidates)
]
sg_gene <- merge(sg_gene, rules, by = "gene_name", all = TRUE)
sg_gene[
  is.na(rule) & gene_name %in% c("ENG", "LEP"),
  rule := "named in the plan (gene level)"
]
sg_gene[
  is.na(rule) & gene_name == "FLT1",
  rule := "named in the plan (isoform level)"
]
sg <- rbind(
  cbind(level = "gene", sg_gene, transcript_name = NA_character_),
  cbind(level = "isoform", sg_iso, n_isoforms_U = NA_integer_),
  fill = TRUE
)
setnames(
  sg,
  c("draw_refit_lfc_lo_2.5", "draw_refit_lfc_hi_97.5"),
  c("draw_refit_lo", "draw_refit_hi")
)
ue_write(sg, file.path(out_dir, "story_genes.parquet"), ue)

# Replication by flag ------------------------------------------------------------------
pf <- rd("placenta_audit/replication/summary_flags.parquet")
ci_or <- ci_parse(pf$adjusted_or_ci)
ci_df <- ci_parse(pf$deficit_ci)
ci_md <- ci_parse(pf$matched_deficit_ci)
rep_pl <- data.table(
  dataset = "placenta half-splits",
  flag = pf$flag,
  n_calls = pf$n_calls,
  pct_flagged = pf$pct_flagged,
  rate_flagged = pf$rate_flagged,
  rate_unflagged = pf$rate_unflagged,
  deficit = pf$deficit,
  deficit_lo = ci_df[, 1],
  deficit_hi = ci_df[, 2],
  matched_deficit = pf$matched_deficit,
  matched_deficit_lo = ci_md[, 1],
  matched_deficit_hi = ci_md[, 2],
  adjusted_or = pf$adjusted_or,
  or_lo = ci_or[, 1],
  or_hi = ci_or[, 2]
)
hm <- hf[method == "M1"]
ci_or <- ci_parse(hm$adjusted_or_ci)
ci_df <- ci_parse(hm$deficit_ci)
ci_md <- ci_parse(hm$matched_deficit_ci)
rep_kd <- data.table(
  dataset = paste("read halves", hm$series, hm$halving),
  flag = hm$flag,
  n_calls = hm$n_calls,
  pct_flagged = hm$pct_flagged,
  rate_flagged = hm$rate_flagged,
  rate_unflagged = hm$rate_unflagged,
  deficit = hm$deficit,
  deficit_lo = ci_df[, 1],
  deficit_hi = ci_df[, 2],
  matched_deficit = hm$matched_deficit,
  matched_deficit_lo = ci_md[, 1],
  matched_deficit_hi = ci_md[, 2],
  adjusted_or = hm$adjusted_or,
  or_lo = ci_or[, 1],
  or_hi = ci_or[, 2]
)
ue_write(
  rbind(rep_pl, rep_kd),
  file.path(out_dir, "replication_flags.parquet"),
  ue
)
ue_write_session(out_dir, ue)

cat("evidence_numbers:", nrow(ev), "rows\n")
cat(
  "flt1_isoforms:",
  nrow(fl),
  "rows; story_genes:",
  nrow(sg),
  "rows; replication_flags:",
  nrow(rep_pl) + nrow(rep_kd),
  "rows\n"
)
print(file.info(list.files(out_dir, full.names = TRUE))[, "size", drop = FALSE])
