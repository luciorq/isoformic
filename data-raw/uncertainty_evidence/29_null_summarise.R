# WP-A step 5: inference and summary tables of the genome-wide placenta null.
#
# Purpose: turn the per-split outputs of 21_null_dte.R, 22_null_dtu.R and
# 23_thin_dte.R into the reported numbers, each with a 95% CI and its
# source table:
#   q1  false calls per split (q < 0.05, q < 0.10, p < 0.01, p < 0.001) and
#       P(>= 1 call at q < 0.05), every method, size, set (own / common) and
#       stratum level, zeros included
#   q2  type I inflation (observed share at p < 0.01 / 0.01, and at
#       p < 0.001 / 0.001), and the top against bottom InfRV tertile ratio
#   q3  formal H1: quasi-Poisson GLM of per-transcript false-call counts at
#       p < 0.01 over the 7v7 leave-one-unit-out splits. Primary model
#       "prereg", as pre-registered (R_A): InfRV tertile + expression decile
#       + isoform class (isoform_bin). Sensitivity models (not
#       pre-registered; amendment A6): "ident" adds non-identifiable, and
#       "zero_unique" adds unique_fraction == 0 instead (7v7 only). The
#       Mantel-Haenszel OR of top tertile against called, stratified by
#       expression decile, is the pre-registered check. H2: method
#       interaction in a stacked model, for the "prereg" and "ident"
#       covariate sets (M1 vs M2 and M1 vs M6, the pre-registered
#       contrasts; M7 vs M6 is a diagnostic of swish on tied draws only,
#       amendment A7)
#   q4  enrichment of M1 null calls in the top tertile
#   q5  share of M1 null calls removed by M9 against M10 (common set)
#   q6  M11 gene-level false calls, raw against divided (in q1/q2 tables)
#   q7  DTU false genes per split (scale FALSE / TRUE), false switches and
#       NMD (productive fraction) shifts, naive against aware
#   q8  thinning arm: empirical FDR at q 0.05 / 0.10, TPR, FDR at matched k
#   q9  p-value histograms (bins pooled over splits)
#   q10 tail sensitivity: share of each mean false-call count from the single
#       worst split, and the mean without it (DTE and DTU)
#   q11 de_confidence (M9) classes against M1 null calls: share of tested
#       and of calls in class high / stable / uncertain-or-ambiguous
#   robustness (draws 16-30, chrX/Y, unrestricted sex, other A/B member,
#   unfolded haplotypes, imbalance covariates with a delete-one-unit
#   jackknife CI, top/low tertile ratio by mapping-rate imbalance) and the
#   verdict on R-A (criterion 3 is NA, not informative, when the M11 CI
#   contains the M1 top-tertile inflation).
#
# Inference:
#   jackknife  delete-one-unit over the 15 control units: the estimate
#              without unit u is the statistic on the splits that do not use
#              u (for 7v7_louo, the 14 splits that leave u out). Raw variance
#              (n-1)/n * sum (theta_u - mean)^2. Monte Carlo component: each
#              theta_u carries split-sampling noise, estimated per unit as
#              (theta(half 1) - theta(half 2))^2 / 4 over two halves of its
#              splits (alternating split order); its expected contribution,
#              (n-1)^2/n * mean, is subtracted. Means use the exact
#              s^2/k. Ratios and rate ratios on the log scale. Primary CI
#              (ci_lo, ci_hi, ci_type): the corrected variance when it is
#              positive and not below the split-bootstrap variance, else the
#              raw jackknife variance (the bootstrap is never the primary
#              variance: it ignores the between-unit variance); clipped to
#              the range of the quantity. Raw and corrected CIs are kept
#              (jk_*_raw, jk_lo/jk_hi). For 5v5 and 3v3 the leave-u-out
#              subsets overlap, so the correction is approximate there.
#   bootstrap  percentile bootstrap over splits (B = 2000; 100 for GLMs),
#              labelled a lower bound: splits share units, so it ignores
#              the between-unit variance.
#   wilson     Wilson intervals on pooled counts (thinning), also a lower
#              bound for the same reason.
#
# Inputs: <out>/null_placenta/{splits,dte,dtu,thin} and <out>/store
# (strata). Env vars: ISOFORMIC_UE_OUT (ue_common.R); ISOFORMIC_UE_BOOT
# (bootstrap replicates, default 2000); ISOFORMIC_UE_BOOT_GLM (bootstrap replicates of
# the GLM and CMH statistics, default 100).
# Resources: measured 41 min on one core, peak 2.3 GB. Run from the
# repository root inside a memory cap of about 8 GB (see README.md), and pin
# BLAS to one thread, or the GLM fits use every core:
#   OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
#     pixi run --frozen Rscript data-raw/uncertainty_evidence/29_null_summarise.R
# Outputs: <out>/null_placenta/summary/*.parquet and verdict.txt.

source("data-raw/uncertainty_evidence/null_common.R")
ue <- ue_setup("29_null_summarise")
p <- wpa_paths()
dir.create(p$summary, recursive = TRUE, showWarnings = FALSE)
seed <- ue$seed_base + 1000L * ue$prereg$seeds$package_index$A + 900L
B <- as.integer(Sys.getenv("ISOFORMIC_UE_BOOT", "2000"))
B_glm <- as.integer(Sys.getenv("ISOFORMIC_UE_BOOT_GLM", "100"))
# jarl-ignore library_call: attached after ue_setup() on purpose, so that its
# names mask those of the package snapshot, as in the runs of the evaluation
library(data.table)
z975 <- stats::qnorm(0.975)
W <- function(x, name) {
  ue_write(
    as.data.frame(x),
    file.path(p$summary, paste0(name, ".parquet")),
    ue,
    seed
  )
  wpa_log("wrote ", name, " (", nrow(x), " rows)")
}

splits <- as.data.table(wpa_read(file.path(p$splits, "splits.parquet")))
splits[,
  units := lapply(seq_len(.N), function(i) {
    c(strsplit(units_a[[i]], ",")[[1L]], strsplit(units_b[[i]], ",")[[1L]])
  })
]
all_units <- sort(unique(unlist(splits$units)))
stopifnot(length(all_units) == 15L)
strata <- as.data.table(wpa_read(file.path(p$store, "strata_tx.parquet")))
U <- strata[in_u %in% TRUE]

# Inference helpers ------------------------------------------------------------------

#' f(ids) -> scalar; ids may repeat (bootstrap). Returns estimate and CIs.
# jarl-ignore cyclomatic_complexity: one pass over all methods of a split, kept as run
infer <- function(
  ids,
  f,
  log_scale = FALSE,
  boot = TRUE,
  nboot = B,
  jk = TRUE,
  boot_seed = seed,
  range = if (log_scale) c(0, Inf) else c(0, Inf)
) {
  tr <- if (log_scale) function(x) log(x) else identity
  bk <- if (log_scale) exp else identity
  theta <- f(ids)
  res <- data.table(estimate = theta, n_splits = length(ids))
  if (jk) {
    sp <- splits[match(ids, split_id)]
    n <- length(all_units)
    loo <- rep(NA_real_, n)
    mc <- rep(NA_real_, n)
    for (j in seq_len(n)) {
      use <- ids[
        !vapply(sp$units, function(x) all_units[[j]] %in% x, logical(1L))
      ]
      if (length(use) < 2L) {
        next
      }
      loo[j] <- tr(f(use))
      vv <- attr(f, "values")
      if (!is.null(vv) && !log_scale) {
        # mean of per-split values: exact split-sampling variance s^2 / k
        mc[j] <- stats::var(vv[use]) / length(use)
      } else {
        h <- seq_along(use) %% 2L == 0L
        mc[j] <- (tr(f(use[h])) - tr(f(use[!h])))^2 / 4
      }
    }
    ok <- is.finite(loo)
    nn <- sum(ok)
    if (nn >= 3L && is.finite(tr(theta))) {
      var_raw <- (nn - 1) / nn * sum((loo[ok] - mean(loo[ok]))^2)
      mcv <- (nn - 1)^2 / nn * mean(mc[ok & is.finite(mc)])
      var_cor <- if (is.finite(mcv)) max(var_raw - mcv, 0) else NA_real_
      res[, `:=`(
        jk_units = nn,
        jk_se_raw = sqrt(var_raw),
        jk_lo_raw = bk(tr(theta) - z975 * sqrt(var_raw)),
        jk_hi_raw = bk(tr(theta) + z975 * sqrt(var_raw)),
        jk_mc_var = mcv,
        jk_se = sqrt(var_cor),
        jk_lo = bk(tr(theta) - z975 * sqrt(var_cor)),
        jk_hi = bk(tr(theta) + z975 * sqrt(var_cor)),
        # the Monte Carlo estimate removes all the variance: report the raw CI
        # TRUE also when the Monte Carlo term could not be estimated
        mc_degenerate = !isTRUE(var_cor > 0) && isTRUE(var_raw > 0)
      )]
    } else {
      res[, `:=`(
        jk_units = nn,
        jk_se_raw = NA_real_,
        jk_lo_raw = NA_real_,
        jk_hi_raw = NA_real_,
        jk_mc_var = NA_real_,
        jk_se = NA_real_,
        jk_lo = NA_real_,
        jk_hi = NA_real_
      )]
    }
  }
  if (boot) {
    # the same resamples for every statistic (and cache hits in GLM fits)
    set.seed(boot_seed)
    bs <- vapply(
      seq_len(nboot),
      function(b) f(sample(ids, replace = TRUE)),
      numeric(1L)
    )
    q <- stats::quantile(bs, c(0.025, 0.975), na.rm = TRUE, names = FALSE)
    res[, `:=`(boot_lo_lowerbound = q[[1L]], boot_hi_lowerbound = q[[2L]])]
    tb <- if (log_scale) suppressWarnings(log(bs)) else bs
    res[, boot_var := stats::var(tb[is.finite(tb)])]
  }
  # Primary CI (the pre-registered delete-one-unit jackknife): the
  # MC-corrected variance when the correction is usable (positive and not
  # below the split-bootstrap variance, which is the split-sampling variance
  # alone and so a lower bound of the total); otherwise the raw jackknife
  # variance. The split bootstrap is never used as the primary variance: it
  # ignores the variance between units. Log scale for ratios; clipped to
  # the range of the quantity.
  if (jk && !is.null(res$jk_se_raw) && is.finite(res$jk_se_raw)) {
    tr <- if (log_scale) log else identity
    bk <- if (log_scale) exp else identity
    v_cor <- res$jk_se^2
    v_boot <- if (!is.null(res$boot_var)) res$boot_var else NA_real_
    use_cor <- isTRUE(v_cor > 0) && (is.na(v_boot) || isTRUE(v_cor >= v_boot))
    v <- if (use_cor) v_cor else res$jk_se_raw^2
    th <- tr(res$estimate)
    if (is.finite(th)) {
      res[, `:=`(
        ci_lo = max(range[[1L]], bk(th - z975 * sqrt(v))),
        ci_hi = min(range[[2L]], bk(th + z975 * sqrt(v))),
        ci_type = if (use_cor) {
          "jackknife MC-corrected"
        } else {
          "jackknife raw (MC correction degenerate or below split-bootstrap variance)"
        }
      )]
    } else {
      res[, `:=`(ci_lo = NA_real_, ci_hi = NA_real_, ci_type = "not estimable")]
    }
  }
  res
}

f_mean <- function(v) {
  f <- function(ids) mean(v[ids])
  attr(f, "values") <- v
  f
}
f_ratio <- function(num, den, scale = 1) {
  function(ids) sum(num[ids]) / sum(den[ids]) / scale
}

read_dir <- function(dir, cols = "*") {
  if (!dir.exists(dir) || !length(list.files(dir, pattern = "parquet$"))) {
    return(NULL)
  }
  as.data.table(wpa_query(sprintf(
    "SELECT %s FROM read_parquet('%s/*.parquet', union_by_name = true)",
    cols,
    dir
  )))
}

# 1. Per-split summaries (zeros filled) ---------------------------------------------
summ <- read_dir(file.path(p$dte, "main", "summary"))
summ[, c("git_sha", "snapshot_sha256", "seed", "script") := NULL]
# Fill zeros: a (split, method) pair that ran gets every (set, stratum, level)
# key the method has anywhere in the same design
fill_zero <- function(d) {
  ran <- unique(d[, .(design, split_id, method)])
  keys <- unique(d[, .(design, method, set, stratum, level)])
  grid <- merge(ran, keys, by = c("design", "method"), allow.cartesian = TRUE)
  out <- merge(
    grid,
    d,
    by = c("design", "split_id", "method", "set", "stratum", "level"),
    all.x = TRUE
  )
  for (cc in c("n_tested", "n_q05", "n_q10", "n_p01", "n_p001")) {
    set(out, which(is.na(out[[cc]])), cc, 0L)
  }
  out
}
summ <- fill_zero(summ)
summ[, variant := "main"]
wpa_log("summary rows: ", nrow(summ), "; splits: ", uniqueN(summ$split_id))

size_of <- c(
  "7v7_louo" = "7v7",
  "7v7_mc" = "7v7",
  "5v5" = "5v5",
  "3v3" = "3v3",
  "7v7_unres" = "7v7_unres"
)
design_sets <- list(
  "7v7" = "7v7_louo",
  "7v7_louo+mc" = c("7v7_louo", "7v7_mc"),
  "5v5" = "5v5",
  "3v3" = "3v3",
  "7v7_unres" = "7v7_unres"
)

# q1 and q2 per (design set, method, set, stratum, level)
set.seed(seed)
q12 <- list()
for (ds in names(design_sets)) {
  d <- summ[design %in% design_sets[[ds]]]
  # the 28-per-unit set is only for M1 and M2
  if (ds == "7v7_louo+mc") {
    d <- d[method %in% c("M1", "M2")]
  }
  # strata beyond all/tertile/main for designs other than 7v7 add little; keep all
  for (key in split(
    d,
    by = c("method", "set", "stratum", "level"),
    drop = TRUE
  )) {
    ids <- key$split_id
    v <- function(col) stats::setNames(as.numeric(key[[col]]), ids)
    one <- list()
    bt <- key$stratum[[1L]] %in% c("all", "infrv_tertile")
    for (thr in c("n_q05", "n_q10", "n_p01", "n_p001")) {
      r <- infer(ids, f_mean(v(thr)), boot = bt)
      r[, `:=`(
        quantity = paste0("mean_calls_", sub("n_", "", thr, fixed = TRUE))
      )]
      one[[length(one) + 1L]] <- r
    }
    anyc <- stats::setNames(as.numeric(key$n_q05 >= 1), ids)
    r <- infer(ids, f_mean(anyc), boot = bt, range = c(0, 1))
    r[, quantity := "p_any_q05"]
    one[[length(one) + 1L]] <- r
    one[[length(one) + 1L]] <- infer(ids, f_mean(v("n_tested")), boot = FALSE)[,
      quantity := "mean_tested"
    ]
    if (sum(key$n_tested) > 0) {
      r <- infer(
        ids,
        f_ratio(v("n_p01"), v("n_tested"), 0.01),
        log_scale = TRUE,
        boot = bt
      )
      r[, quantity := "inflation_p01"]
      one[[length(one) + 1L]] <- r
      r <- infer(
        ids,
        f_ratio(v("n_p001"), v("n_tested"), 0.001),
        log_scale = TRUE,
        boot = bt
      )
      r[, quantity := "inflation_p001"]
      one[[length(one) + 1L]] <- r
    }
    one <- rbindlist(one, fill = TRUE)
    one[, `:=`(
      design_set = ds,
      method = key$method[[1L]],
      set = key$set[[1L]],
      stratum = key$stratum[[1L]],
      level = key$level[[1L]]
    )]
    q12[[length(q12) + 1L]] <- one
  }
  wpa_log("q1/q2 done for ", ds)
}
q12 <- rbindlist(q12, fill = TRUE)
setcolorder(
  q12,
  c("design_set", "method", "set", "stratum", "level", "quantity", "estimate")
)
W(q12, "q1_q2_calls_and_inflation")

# q2: top against bottom tertile, ratio of inflations, per method and design
tert_ratio <- list()
for (ds in c("7v7", "5v5", "3v3")) {
  d <- summ[
    design %in% design_sets[[ds]] & stratum == "infrv_tertile" & set == "own"
  ]
  for (m in unique(d$method)) {
    dm <- d[method == m]
    hi <- dm[level == "high"]
    lo <- dm[level == "low"]
    ids <- intersect(hi$split_id, lo$split_id)
    if (!length(ids)) {
      next
    }
    for (thr in c("n_p01", "n_q05")) {
      nh <- stats::setNames(as.numeric(hi[[thr]]), hi$split_id)
      th <- stats::setNames(as.numeric(hi$n_tested), hi$split_id)
      nl <- stats::setNames(as.numeric(lo[[thr]]), lo$split_id)
      tl <- stats::setNames(as.numeric(lo$n_tested), lo$split_id)
      f <- function(ids) {
        (sum(nh[ids]) / sum(th[ids])) / (sum(nl[ids]) / sum(tl[ids]))
      }
      r <- infer(ids, f, log_scale = TRUE)
      r[, `:=`(
        design_set = ds,
        method = m,
        threshold = thr,
        quantity = "rate_ratio_high_vs_low_tertile",
        calls_high = sum(nh),
        calls_low = sum(nl)
      )]
      tert_ratio[[length(tert_ratio) + 1L]] <- r
    }
  }
}
tert_ratio <- rbindlist(tert_ratio, fill = TRUE)
W(tert_ratio, "q2_tertile_rate_ratio")

# q4 enrichment of null calls in the top tertile (share of calls / share of tested)
enr <- list()
for (ds in c("7v7", "5v5", "3v3")) {
  d <- summ[
    design %in% design_sets[[ds]] & stratum == "infrv_tertile" & set == "own"
  ]
  for (m in unique(d$method)) {
    dm <- d[method == m]
    tot <- dm[,
      .(n_tested = sum(n_tested), n_q05 = sum(n_q05), n_p01 = sum(n_p01)),
      by = split_id
    ]
    hi <- dm[level == "high"]
    ids <- hi$split_id
    for (thr in c("n_q05", "n_p01")) {
      ch <- stats::setNames(as.numeric(hi[[thr]]), ids)
      th <- stats::setNames(as.numeric(hi$n_tested), ids)
      ca <- stats::setNames(as.numeric(tot[[thr]]), tot$split_id)
      ta <- stats::setNames(as.numeric(tot$n_tested), tot$split_id)
      f <- function(ids) {
        (sum(ch[ids]) / sum(ca[ids])) / (sum(th[ids]) / sum(ta[ids]))
      }
      if (sum(ca) == 0) {
        enr[[length(enr) + 1L]] <- data.table(
          design_set = ds,
          method = m,
          threshold = thr,
          estimate = NA_real_,
          calls_total = 0
        )
        next
      }
      r <- infer(ids, f, log_scale = TRUE)
      r[, `:=`(
        design_set = ds,
        method = m,
        threshold = thr,
        calls_total = sum(ca),
        calls_top = sum(ch),
        share_calls_top = sum(ch) / sum(ca),
        share_tested_top = sum(th) / sum(ta)
      )]
      enr[[length(enr) + 1L]] <- r
    }
  }
}
W(rbindlist(enr, fill = TRUE), "q4_enrichment_top_tertile")

# q5 share of M1 null calls removed by M9 against M10 (matched removal), common set
rem <- list()
for (ds in c("7v7", "5v5", "3v3")) {
  d <- summ[design %in% design_sets[[ds]] & stratum == "all"]
  for (st in c("common", "own")) {
    for (v in c("a", "b")) {
      for (thr in c("n_q05", "n_p01")) {
        m1 <- d[method == "M1" & set == st]
        m9 <- d[method == paste0("M9", v) & set == st]
        m10 <- d[method == paste0("M10", v) & set == st]
        ids <- Reduce(intersect, list(m1$split_id, m9$split_id, m10$split_id))
        if (!length(ids)) {
          next
        }
        c1 <- stats::setNames(as.numeric(m1[[thr]]), m1$split_id)
        c9 <- stats::setNames(as.numeric(m9[[thr]]), m9$split_id)
        c10 <- stats::setNames(as.numeric(m10[[thr]]), m10$split_id)
        if (sum(c1[ids]) == 0) {
          rem[[length(rem) + 1L]] <- data.table(
            design_set = ds,
            set = st,
            m9_rule = v,
            threshold = thr,
            quantity = "no M1 calls",
            estimate = NA_real_,
            m1_calls = 0
          )
          next
        }
        for (q in c("removed_M9", "removed_M10", "removed_M9_minus_M10")) {
          f <- switch(
            q,
            removed_M9 = function(ids) 1 - sum(c9[ids]) / sum(c1[ids]),
            removed_M10 = function(ids) 1 - sum(c10[ids]) / sum(c1[ids]),
            removed_M9_minus_M10 = function(ids) {
              (sum(c10[ids]) - sum(c9[ids])) / sum(c1[ids])
            }
          )
          r <- infer(
            ids,
            f,
            range = if (q == "removed_M9_minus_M10") c(-1, 1) else c(0, 1)
          )
          r[, `:=`(
            design_set = ds,
            set = st,
            m9_rule = v,
            threshold = thr,
            quantity = q,
            m1_calls = sum(c1[ids]),
            m9_calls = sum(c9[ids]),
            m10_calls = sum(c10[ids])
          )]
          rem[[length(rem) + 1L]] <- r
        }
      }
    }
  }
}
W(rbindlist(rem, fill = TRUE), "q5_m9_vs_m10_removal")

# 2. Per-transcript counts (H1 / H2 GLMs, CMH) --------------------------------------
# Primary: the 210 7v7 leave-one-unit-out splits, M1, M2, M6, M7. Secondary
# (smaller samples; the jackknife over random splits is approximate because
# the leave-u-out subsets overlap): 5v5 (M1, M2, M6, M7) and 3v3 (M1, M2),
# the pre-registered model and the non-identifiable sensitivity model.
Uid <- U$transcript_id
cov <- data.frame(
  infrv_tertile = factor(U$infrv_tertile, levels = c("low", "mid", "high")),
  expr_decile = factor(U$expr_decile),
  isoform_bin = factor(U$isoform_bin),
  non_identifiable = factor(
    U$non_identifiable %in% TRUE,
    levels = c(FALSE, TRUE)
  ),
  zero_unique = factor(U$zero_unique %in% TRUE, levels = c(FALSE, TRUE))
)
wpa_log(
  "covariate NA: non_identifiable ",
  sum(is.na(U$non_identifiable)),
  ", zero_unique ",
  sum(is.na(U$zero_unique)),
  ", isoform_bin levels ",
  paste(levels(cov$isoform_bin), collapse = "/")
)
fit_glm <- function(y, n, form_rhs) {
  d <- cbind(cov, y = y, n = n)
  d <- d[d$n > 0, ]
  suppressWarnings(stats::glm(
    stats::as.formula(paste("y ~", form_rhs, "+ offset(log(n))")),
    family = stats::quasipoisson(),
    data = d
  ))
}
# "prereg" is the R_A model of prereg.yaml (expression decile and isoform
# class) and the primary one; "ident" and "zero_unique" are sensitivity
# models that were not pre-registered (amendment A6)
rhs_all <- c(
  prereg = "infrv_tertile + expr_decile + isoform_bin",
  ident = "infrv_tertile + expr_decile + isoform_bin + non_identifiable",
  zero_unique = "infrv_tertile + expr_decile + isoform_bin + zero_unique"
)
terms_of <- c(
  "infrv_tertilehigh",
  "infrv_tertilemid",
  "non_identifiableTRUE",
  "zero_uniqueTRUE"
)
glm_plan <- list(
  "7v7" = list(
    design = "7v7_louo",
    methods = c("M1", "M2", "M6", "M7"),
    models = c("prereg", "ident", "zero_unique"),
    h2 = list(c("M1", "M2"), c("M1", "M6"), c("M7", "M6"))
  ),
  "5v5" = list(
    design = "5v5",
    methods = c("M1", "M2", "M6", "M7"),
    models = c("prereg", "ident"),
    h2 = list(c("M1", "M2"), c("M1", "M6"), c("M7", "M6"))
  ),
  "3v3" = list(
    design = "3v3",
    methods = c("M1", "M2"),
    models = c("prereg", "ident"),
    h2 = list(c("M1", "M2"))
  )
)
txdir <- file.path(p$dte, "main", "tx")
glm_rows <- list()
h2 <- list()
cmh_rows <- list()
for (ds in names(glm_plan)) {
  gp <- glm_plan[[ds]]
  ids_d <- splits[design == gp$design, split_id]
  have <- ids_d[file.exists(file.path(txdir, paste0(ids_d, ".parquet")))]
  wpa_log("tx tables of ", gp$design, ": ", length(have))
  files <- paste0(
    "'",
    file.path(txdir, paste0(have, ".parquet")),
    "'",
    collapse = ","
  )
  tx <- as.data.table(wpa_query(sprintf(
    "SELECT split_id, transcript_id, %s FROM read_parquet([%s], union_by_name = true)",
    paste0(gp$methods, "_p", collapse = ", "),
    files
  )))
  ti <- match(tx$transcript_id, Uid)
  si <- match(tx$split_id, have)
  mats <- list()
  for (m in gp$methods) {
    pv <- tx[[paste0(m, "_p")]]
    ok <- !is.na(pv)
    Tm <- matrix(0L, length(Uid), length(have), dimnames = list(Uid, have))
    Cm <- Tm
    Tm[cbind(ti[ok], si[ok])] <- 1L
    Cm[cbind(ti[ok], si[ok])] <- as.integer(pv[ok] < 0.01)
    mats[[m]] <- list(T = Tm, C = Cm)
  }
  rm(tx, ti, si)
  gc()
  w_of <- function(ids) tabulate(match(ids, have), length(have))
  for (m in names(mats)) {
    Tm <- mats[[m]]$T
    Cm <- mats[[m]]$C
    if (sum(Tm) == 0) {
      next
    }
    for (mod in gp$models) {
      memo <- new.env()
      fit_coefs <- function(ids) {
        w <- w_of(ids)
        k <- paste(w, collapse = ",")
        if (is.null(memo[[k]])) {
          memo[[k]] <- stats::coef(fit_glm(
            as.vector(Cm %*% w),
            as.vector(Tm %*% w),
            rhs_all[[mod]]
          ))
        }
        memo[[k]]
      }
      coef_fun <- function(term) {
        function(ids) {
          cf <- fit_coefs(ids)
          if (term %in% names(cf)) exp(cf[[term]]) else NA_real_
        }
      }
      full <- fit_glm(rowSums(Cm), rowSums(Tm), rhs_all[[mod]])
      sm <- summary(full)$coefficients
      for (term in intersect(terms_of, rownames(sm))) {
        r <- infer(have, coef_fun(term), log_scale = TRUE, nboot = B_glm)
        r[, `:=`(
          design_set = ds,
          method = m,
          model = mod,
          term = term,
          role = if (mod == "prereg") {
            "pre-registered"
          } else {
            "sensitivity (not pre-registered; A6)"
          },
          model_lo = exp(sm[term, 1L] - z975 * sm[term, 2L]),
          model_hi = exp(sm[term, 1L] + z975 * sm[term, 2L]),
          quasi_dispersion = summary(full)$dispersion,
          n_transcripts = sum(rowSums(Tm) > 0),
          n_calls = sum(Cm),
          n_tests = sum(Tm)
        )]
        glm_rows[[length(glm_rows) + 1L]] <- r
      }
      wpa_log("GLM ", ds, " ", m, " ", mod)
    }
  }
  # H2: method interaction in a stacked model (log ratio of the two RRs), for
  # the pre-registered covariates and the non-identifiable sensitivity set
  h2_models <- intersect(c("prereg", "ident"), gp$models)
  for (pair in gp$h2) {
    for (mod in h2_models) {
      a <- mats[[pair[[1L]]]]
      b <- mats[[pair[[2L]]]]
      if (is.null(a) || is.null(b) || sum(a$T) == 0 || sum(b$T) == 0) {
        next
      }
      h2_form <- stats::as.formula(paste0(
        "y ~ (",
        rhs_all[[mod]],
        ") * method + offset(log(n))"
      ))
      f <- function(ids) {
        w <- w_of(ids)
        d <- rbind(
          cbind(
            cov,
            y = as.vector(a$C %*% w),
            n = as.vector(a$T %*% w),
            method = "ref"
          ),
          cbind(
            cov,
            y = as.vector(b$C %*% w),
            n = as.vector(b$T %*% w),
            method = "alt"
          )
        )
        d$method <- factor(d$method, levels = c("ref", "alt"))
        d <- d[d$n > 0, ]
        fit <- suppressWarnings(stats::glm(
          h2_form,
          family = stats::quasipoisson(),
          data = d
        ))
        exp(stats::coef(fit)[["infrv_tertilehigh:methodalt"]])
      }
      r <- infer(have, f, log_scale = TRUE, nboot = B_glm)
      r[, `:=`(
        design_set = ds,
        model = mod,
        reference = pair[[1L]],
        alternative = pair[[2L]],
        quantity = "ratio of top-tertile adjusted RR (alternative / reference)",
        role = if (identical(pair[[1L]], "M7")) {
          "diagnostic (M7 is swish on tied copies of the point counts, not a calibrated reference; A7)"
        } else if (mod == "prereg") {
          "pre-registered"
        } else {
          "sensitivity (not pre-registered; A6)"
        }
      )]
      h2[[length(h2) + 1L]] <- r
      wpa_log("H2 ", ds, " ", mod, " ", paste(pair, collapse = " vs "))
    }
  }
  # CMH: top tertile vs other, called vs not, stratified by expression decile
  top <- cov$infrv_tertile %in% "high"
  dec <- cov$expr_decile
  for (m in names(mats)) {
    Tm <- mats[[m]]$T
    Cm <- mats[[m]]$C
    if (sum(Tm) == 0) {
      next
    }
    cmh_or <- function(ids, full = FALSE) {
      w <- w_of(ids)
      y <- as.vector(Cm %*% w)
      n <- as.vector(Tm %*% w)
      tab <- array(0, c(2L, 2L, nlevels(dec)))
      for (k in seq_len(nlevels(dec))) {
        s <- dec == levels(dec)[k]
        tab[1L, 1L, k] <- sum(y[s & top])
        tab[1L, 2L, k] <- sum(n[s & top] - y[s & top])
        tab[2L, 1L, k] <- sum(y[s & !top])
        tab[2L, 2L, k] <- sum(n[s & !top] - y[s & !top])
      }
      keep <- apply(tab, 3L, function(x) {
        all(rowSums(x) > 0) && all(colSums(x) > 0)
      })
      if (sum(keep) < 1L) {
        return(if (full) NULL else NA_real_)
      }
      mh <- stats::mantelhaen.test(tab[,, keep, drop = FALSE])
      if (full) mh else unname(mh$estimate)
    }
    mh <- cmh_or(have, TRUE)
    if (is.null(mh)) {
      next
    }
    r <- infer(have, cmh_or, log_scale = TRUE, nboot = B_glm)
    r[, `:=`(
      design_set = ds,
      method = m,
      quantity = "MH OR top tertile vs rest, called at p<0.01, strata expression decile",
      mh_p = mh$p.value,
      mh_lo_independence = mh$conf.int[[1L]],
      mh_hi_independence = mh$conf.int[[2L]]
    )]
    cmh_rows[[length(cmh_rows) + 1L]] <- r
  }
  rm(mats)
  gc()
}
glm_tab <- rbindlist(glm_rows, fill = TRUE)
W(glm_tab, "q3_glm_rate_ratios")
W(rbindlist(h2, fill = TRUE), "q3_h2_method_interaction")
W(rbindlist(cmh_rows, fill = TRUE), "q3_cmh")


# 3. p-value histograms ----------------------------------------------------------------
hist <- read_dir(
  file.path(p$dte, "main", "hist"),
  "design, method, stratum, level, bin, lower, n, split_id"
)
if (!is.null(hist)) {
  hist[, size := size_of[design]]
  h <- hist[,
    .(n = sum(n), n_splits = uniqueN(split_id)),
    by = .(size, design, method, stratum, level, bin, lower)
  ]
  h[, share := n / sum(n), by = .(design, method, stratum, level)]
  W(h, "q9_pvalue_histograms")
}

# 4. DTU, switches, NMD ------------------------------------------------------------------
dtu_rows <- list()
for (bk in c("diffsplice", "satuRn", "DEXSeq")) {
  d <- read_dir(file.path(p$dtu, bk, "summary"))
  if (is.null(d)) {
    next
  }
  d <- merge(d, splits[, .(split_id, design)], by = "split_id")
  for (key in split(
    d,
    by = c("design", "scaled", "stratum", "level"),
    drop = TRUE
  )) {
    ids <- key$split_id
    for (q in c("n_genes_q05", "n_genes_q10", "n_genes_p01", "n_tx_q05")) {
      r <- infer(ids, f_mean(stats::setNames(as.numeric(key[[q]]), ids)))
      r[, quantity := paste0("mean_", q)]
      dtu_rows[[length(dtu_rows) + 1L]] <- cbind(
        r,
        backend = bk,
        design = key$design[[1L]],
        scaled = key$scaled[[1L]],
        stratum = key$stratum[[1L]],
        level = key$level[[1L]]
      )
    }
    if (sum(key$n_genes_tested) > 0) {
      r <- infer(
        ids,
        f_ratio(
          stats::setNames(as.numeric(key$n_genes_p01), ids),
          stats::setNames(as.numeric(key$n_genes_tested), ids),
          0.01
        ),
        log_scale = TRUE
      )
      r[, quantity := "gene_inflation_p01"]
      dtu_rows[[length(dtu_rows) + 1L]] <- cbind(
        r,
        backend = bk,
        design = key$design[[1L]],
        scaled = key$scaled[[1L]],
        stratum = key$stratum[[1L]],
        level = key$level[[1L]]
      )
      # false genes at q < 0.05 per tested gene: the scaled arm tests fewer genes
      r <- infer(
        ids,
        f_ratio(
          stats::setNames(as.numeric(key$n_genes_q05), ids),
          stats::setNames(as.numeric(key$n_genes_tested), ids)
        ),
        log_scale = TRUE
      )
      r[, quantity := "rate_q05_per_tested_gene"]
      dtu_rows[[length(dtu_rows) + 1L]] <- cbind(
        r,
        backend = bk,
        design = key$design[[1L]],
        scaled = key$scaled[[1L]],
        stratum = key$stratum[[1L]],
        level = key$level[[1L]]
      )
    }
    r <- infer(
      ids,
      f_mean(stats::setNames(as.numeric(key$n_genes_tested), ids)),
      boot = FALSE
    )
    r[, quantity := "mean_n_genes_tested"]
    dtu_rows[[length(dtu_rows) + 1L]] <- cbind(
      r,
      backend = bk,
      design = key$design[[1L]],
      scaled = key$scaled[[1L]],
      stratum = key$stratum[[1L]],
      level = key$level[[1L]]
    )
    r <- infer(
      ids,
      f_mean(stats::setNames(as.numeric(key$n_genes_q05 >= 1), ids)),
      range = c(0, 1)
    )
    r[, quantity := "p_any_gene_q05"]
    dtu_rows[[length(dtu_rows) + 1L]] <- cbind(
      r,
      backend = bk,
      design = key$design[[1L]],
      scaled = key$scaled[[1L]],
      stratum = key$stratum[[1L]],
      level = key$level[[1L]]
    )
    r <- infer(
      ids,
      f_mean(stats::setNames(as.numeric(key$seconds), ids)),
      boot = FALSE,
      jk = FALSE
    )
    r[, quantity := "seconds_per_split"]
    dtu_rows[[length(dtu_rows) + 1L]] <- cbind(
      r,
      backend = bk,
      design = key$design[[1L]],
      scaled = key$scaled[[1L]],
      stratum = key$stratum[[1L]],
      level = key$level[[1L]]
    )
  }
}
if (length(dtu_rows)) {
  W(rbindlist(dtu_rows, fill = TRUE), "q7_dtu")
}
sw <- read_dir(file.path(p$dtu, "switches"))
if (!is.null(sw)) {
  rows <- list()
  for (key in split(sw, by = "arm")) {
    ids <- key$split_id
    for (q in c(
      "n_genes_with_sig_isoform",
      "n_opposite_direction",
      "n_det_not_deg",
      "n_consistent",
      "n_opposite_top_tertile",
      "n_opposite_low_tertile",
      "n_deg"
    )) {
      r <- infer(ids, f_mean(stats::setNames(as.numeric(key[[q]]), ids)))
      rows[[length(rows) + 1L]] <- r[, `:=`(
        arm = key$arm[[1L]],
        quantity = paste0("mean_", q)
      )]
    }
    r <- infer(
      ids,
      f_mean(stats::setNames(as.numeric(key$n_opposite_direction >= 1), ids)),
      range = c(0, 1)
    )
    rows[[length(rows) + 1L]] <- r[, `:=`(
      arm = key$arm[[1L]],
      quantity = "p_any_opposite_direction"
    )]
  }
  W(rbindlist(rows, fill = TRUE), "q7_switches")
}
nmd <- read_dir(file.path(p$dtu, "nmd"))
if (!is.null(nmd)) {
  rows <- list()
  for (key in split(nmd, by = "arm")) {
    ids <- key$split_id
    for (q in c(
      "n_padj05",
      "n_padj10",
      "n_p01",
      "n_padj05_top_tertile",
      "n_genes_tested"
    )) {
      r <- infer(
        ids,
        f_mean(stats::setNames(as.numeric(key[[q]]), ids)),
        boot = q != "n_genes_tested"
      )
      rows[[length(rows) + 1L]] <- r[, `:=`(
        arm = key$arm[[1L]],
        quantity = paste0("mean_", q)
      )]
    }
    r <- infer(
      ids,
      f_ratio(
        stats::setNames(as.numeric(key$n_p01), ids),
        stats::setNames(as.numeric(key$n_genes_tested), ids),
        0.01
      ),
      log_scale = TRUE
    )
    rows[[length(rows) + 1L]] <- r[, `:=`(
      arm = key$arm[[1L]],
      quantity = "inflation_p01"
    )]
  }
  W(rbindlist(rows, fill = TRUE), "q7_nmd")
}

# 5. Thinning ------------------------------------------------------------------------------
for (rule in c("affine", "scaled")) {
  met <- read_dir(file.path(p$thin, rule, "metrics"))
  if (is.null(met)) {
    next
  }
  met[, design := splits$design[match(split_id, splits$split_id)]]
  rows <- list()
  for (key in split(
    met,
    by = c("design", "method", "subset", "rule"),
    drop = TRUE
  )) {
    ids <- key$split_id
    nf <- stats::setNames(as.numeric(key$n_false), ids)
    nc <- stats::setNames(as.numeric(key$n_calls), ids)
    nt <- stats::setNames(as.numeric(key$n_true_calls), ids)
    ntr <- stats::setNames(as.numeric(key$n_true), ids)
    fdp <- ifelse(nc > 0, nf / nc, 0)
    out <- list(
      infer(ids, f_ratio(nf, nc), range = c(0, 1))[, quantity := "fdr_pooled"],
      infer(ids, f_mean(fdp), range = c(0, 1))[, quantity := "mean_fdp"],
      infer(ids, f_ratio(nt, ntr), range = c(0, 1))[, quantity := "tpr_pooled"],
      infer(ids, f_mean(nc), boot = FALSE)[, quantity := "mean_calls"]
    )
    wf <- wilson(sum(nf), sum(nc))
    wt <- wilson(sum(nt), sum(ntr))
    out[[1L]][, `:=`(
      wilson_lo_lowerbound = wf[[1L]],
      wilson_hi_lowerbound = wf[[2L]]
    )]
    out[[3L]][, `:=`(
      wilson_lo_lowerbound = wt[[1L]],
      wilson_hi_lowerbound = wt[[2L]]
    )]
    o <- rbindlist(out, fill = TRUE)
    o[, `:=`(
      thin_rule = rule,
      design = key$design[[1L]],
      method = key$method[[1L]],
      subset = key$subset[[1L]],
      call_rule = key$rule[[1L]],
      total_calls = sum(nc),
      total_false = sum(nf),
      total_true = sum(ntr)
    )]
    rows[[length(rows) + 1L]] <- o
  }
  W(rbindlist(rows, fill = TRUE), paste0("q8_thinning_", rule))
}

# 6. Robustness ----------------------------------------------------------------------------
rob <- list()
add_rob <- function(d, label) {
  if (is.null(d) || !nrow(d)) {
    return(invisible())
  }
  for (key in split(
    d,
    by = c("method", "set", "stratum", "level"),
    drop = TRUE
  )) {
    ids <- key$split_id
    v <- function(col) stats::setNames(as.numeric(key[[col]]), ids)
    r1 <- infer(ids, f_mean(v("n_q05")))[, quantity := "mean_calls_q05"]
    r2 <- if (sum(key$n_tested) > 0) {
      infer(ids, f_ratio(v("n_p01"), v("n_tested"), 0.01), log_scale = TRUE)[,
        quantity := "inflation_p01"
      ]
    } else {
      NULL
    }
    r <- rbindlist(list(r1, r2), fill = TRUE)
    r[, `:=`(
      analysis = label,
      method = key$method[[1L]],
      set = key$set[[1L]],
      stratum = key$stratum[[1L]],
      level = key$level[[1L]]
    )]
    # jarl-ignore undesirable_operator: appends to a result list of the enclosing scope
    rob[[length(rob) + 1L]] <<- r
  }
}
# (a) chrX/chrY excluded, (b) unrestricted sex: from the main summaries
add_rob(
  summ[
    design == "7v7_louo" &
      stratum == "chr_xy" &
      level == "FALSE" &
      method %in% c("M1", "M2", "M6", "M7", "M9a", "M11raw", "M11div")
  ],
  "7v7 without chrX/chrY transcripts"
)
add_rob(
  summ[design == "7v7_unres" & stratum %in% c("all", "infrv_tertile")],
  "7v7 unrestricted sex"
)
add_rob(
  summ[
    design == "7v7_louo" &
      split_id %in%
        splits[design == "7v7_louo" & abs(sex_diff) == 0, split_id] &
      stratum %in% c("all", "infrv_tertile") &
      method %in% c("M1", "M2")
  ],
  "7v7 with equal sex counts"
)
for (vr in c("draws16", "member", "unfolded")) {
  d <- read_dir(file.path(p$dte, vr, "summary"))
  if (is.null(d)) {
    next
  }
  d[, c("git_sha", "snapshot_sha256", "seed", "script") := NULL]
  d <- fill_zero(d)
  # the same splits of the main run, for a paired comparison
  main_same <- summ[
    split_id %in% unique(d$split_id) & method %in% unique(d$method)
  ]
  keep_strata <- c("all", "infrv_tertile", "infrv_tertile_15")
  add_rob(d[stratum %in% keep_strata], paste0("variant ", vr))
  add_rob(
    main_same[stratum %in% keep_strata],
    paste0("main run on the splits of variant ", vr)
  )
}
if (length(rob)) {
  W(rbindlist(rob, fill = TRUE), "robustness")
}

# Imbalance covariates: false-call counts per split on depth, mapping and sex imbalance
m1 <- summ[
  method %in%
    c("M1", "M2") &
    set == "own" &
    stratum == "all" &
    design %in% c("7v7_louo", "7v7_mc", "5v5", "3v3", "7v7_unres")
]
m1 <- merge(
  m1,
  splits[, .(split_id, depth_diff, mapping_diff, sex_diff)],
  by = "split_id"
)
imb_rows <- list()
for (m in c("M1", "M2")) {
  d <- m1[method == m]
  fit <- suppressWarnings(stats::glm(
    n_p01 ~ abs(depth_diff) +
      abs(mapping_diff) +
      abs(sex_diff) +
      design +
      offset(log(n_tested)),
    family = stats::quasipoisson(),
    data = d
  ))
  sm <- summary(fit)$coefficients
  # the model SE treats splits as independent, but splits share units: a
  # delete-one-unit jackknife of the coefficients gives the honest CI
  su <- splits$units[match(d$split_id, splits$split_id)]
  loo <- t(vapply(
    all_units,
    function(u) {
      dj <- d[!vapply(su, function(x) u %in% x, logical(1L))]
      cf <- stats::coef(suppressWarnings(stats::glm(
        n_p01 ~ abs(depth_diff) +
          abs(mapping_diff) +
          abs(sex_diff) +
          design +
          offset(log(n_tested)),
        family = stats::quasipoisson(),
        data = dj
      )))
      cf[rownames(sm)]
    },
    numeric(nrow(sm))
  ))
  nn <- nrow(loo)
  jk_se <- sqrt(
    (nn - 1) /
      nn *
      colSums(sweep(loo, 2L, colMeans(loo, na.rm = TRUE))^2, na.rm = TRUE)
  )
  imb_rows[[m]] <- data.table(
    method = m,
    term = rownames(sm),
    estimate_log = sm[, 1L],
    se = sm[, 2L],
    rr = exp(sm[, 1L]),
    rr_lo = exp(sm[, 1L] - z975 * sm[, 2L]),
    rr_hi = exp(sm[, 1L] + z975 * sm[, 2L]),
    p_model_invalid = sm[, 4L],
    jk_se_raw = jk_se,
    jk_lo_raw = exp(sm[, 1L] - z975 * jk_se),
    jk_hi_raw = exp(sm[, 1L] + z975 * jk_se),
    n_splits = nrow(d),
    dispersion = summary(fit)$dispersion,
    note = "rr_lo/rr_hi and p_model_invalid assume independent splits (they share units): use jk_lo_raw/jk_hi_raw"
  )
}
W(rbindlist(imb_rows), "robustness_imbalance_glm")

# 7. Monte Carlo check: 14 against 28 splits per unit (M1, M2) ---------------------------
mc_rows <- q12[
  design_set %in%
    c("7v7", "7v7_louo+mc") &
    method %in% c("M1", "M2") &
    set == "own" &
    stratum %in% c("all", "infrv_tertile") &
    quantity %in% c("mean_calls_q05", "mean_calls_p01", "inflation_p01")
]
W(mc_rows, "mc_check_14_vs_28")

# 7b. Tail sensitivity: how much of each mean comes from the single worst split ---------
tail_rows <- list()
for (ds in c("7v7", "5v5", "3v3", "7v7_unres")) {
  d <- summ[design %in% design_sets[[ds]] & stratum == "all"]
  for (key in split(d, by = c("method", "set"), drop = TRUE)) {
    x <- key$n_q05
    top <- which.max(x)
    tail_rows[[length(tail_rows) + 1L]] <- data.table(
      design_set = ds,
      method = key$method[[1L]],
      set = key$set[[1L]],
      n_splits = length(x),
      mean_tested = mean(key$n_tested),
      total_q05 = sum(x),
      mean_q05 = mean(x),
      median_q05 = stats::median(x),
      max_q05 = max(x),
      max_split = key$split_id[[top]],
      share_from_max_split = if (sum(x) > 0) max(x) / sum(x) else NA_real_,
      mean_q05_without_max_split = mean(x[-top]),
      p_any_q05 = mean(x >= 1),
      p_any_q05_without_max_split = mean(x[-top] >= 1)
    )
  }
}
for (bk in c("diffsplice", "satuRn", "DEXSeq")) {
  d <- read_dir(file.path(p$dtu, bk, "summary"))
  if (is.null(d)) {
    next
  }
  d <- merge(
    d[stratum == "all"],
    splits[, .(split_id, design)],
    by = "split_id"
  )
  for (key in split(d, by = c("design", "scaled"), drop = TRUE)) {
    x <- key$n_genes_q05
    top <- which.max(x)
    tail_rows[[length(tail_rows) + 1L]] <- data.table(
      design_set = key$design[[1L]],
      method = paste0("DTU ", bk, if (key$scaled[[1L]]) " scaled" else " raw"),
      set = "own",
      n_splits = length(x),
      mean_tested = mean(key$n_genes_tested),
      total_q05 = sum(x),
      mean_q05 = mean(x),
      median_q05 = stats::median(x),
      max_q05 = max(x),
      max_split = key$split_id[[top]],
      share_from_max_split = if (sum(x) > 0) max(x) / sum(x) else NA_real_,
      mean_q05_without_max_split = mean(x[-top]),
      p_any_q05 = mean(x >= 1),
      p_any_q05_without_max_split = mean(x[-top] >= 1)
    )
  }
}
W(rbindlist(tail_rows, fill = TRUE), "q10_tail_sensitivity")

# 7c. de_confidence (M9) and the M1 null calls: is "high" depleted of false calls? -------
# Per split, among M1-tested transcripts: counts by class, and the M1 calls at
# p < 0.01 in each class. The sign-agreement criterion ("stable") rewards a
# consistent direction across draws, which a null call with a large |lfc|
# also has.
circ_rows <- list()
for (ds in c("7v7", "5v5", "3v3")) {
  ids_d <- splits[design %in% design_sets[[ds]], split_id]
  have <- ids_d[file.exists(file.path(txdir, paste0(ids_d, ".parquet")))]
  files <- paste0(
    "'",
    file.path(txdir, paste0(have, ".parquet")),
    "'",
    collapse = ","
  )
  ag <- as.data.table(wpa_query(sprintf(
    paste(
      "SELECT split_id, count(*) AS n_tested, sum((M1_p < 0.01)::INT) AS n_call,",
      "sum((M9_high IS TRUE)::INT) AS n_high, sum((M9_high IS TRUE AND M1_p < 0.01)::INT) AS n_call_high,",
      "sum((M9_stable IS TRUE)::INT) AS n_stable, sum((M9_stable IS TRUE AND M1_p < 0.01)::INT) AS n_call_stable,",
      "sum((M9_uncertain IS TRUE OR M9_ambiguous IS TRUE)::INT) AS n_flag,",
      "sum(((M9_uncertain IS TRUE OR M9_ambiguous IS TRUE) AND M1_p < 0.01)::INT) AS n_call_flag",
      "FROM read_parquet([%s]) WHERE M1_p IS NOT NULL GROUP BY split_id"
    ),
    files
  )))
  ids <- ag$split_id
  v <- function(col) stats::setNames(as.numeric(ag[[col]]), ids)
  qs <- list(
    share_tested_high = list(
      f_ratio(v("n_high"), v("n_tested")),
      FALSE,
      c(0, 1)
    ),
    share_calls_high = list(
      f_ratio(v("n_call_high"), v("n_call")),
      FALSE,
      c(0, 1)
    ),
    share_tested_stable = list(
      f_ratio(v("n_stable"), v("n_tested")),
      FALSE,
      c(0, 1)
    ),
    share_calls_stable = list(
      f_ratio(v("n_call_stable"), v("n_call")),
      FALSE,
      c(0, 1)
    ),
    share_tested_uncertain_or_ambiguous = list(
      f_ratio(v("n_flag"), v("n_tested")),
      FALSE,
      c(0, 1)
    ),
    share_calls_uncertain_or_ambiguous = list(
      f_ratio(v("n_call_flag"), v("n_call")),
      FALSE,
      c(0, 1)
    ),
    p_call_given_high = list(
      f_ratio(v("n_call_high"), v("n_high")),
      FALSE,
      c(0, 1)
    ),
    p_call_overall = list(f_ratio(v("n_call"), v("n_tested")), FALSE, c(0, 1)),
    enrichment_calls_in_high = list(
      function(ids) {
        (sum(v("n_call_high")[ids]) / sum(v("n_call")[ids])) /
          (sum(v("n_high")[ids]) / sum(v("n_tested")[ids]))
      },
      TRUE,
      c(0, Inf)
    ),
    enrichment_calls_in_stable = list(
      function(ids) {
        (sum(v("n_call_stable")[ids]) / sum(v("n_call")[ids])) /
          (sum(v("n_stable")[ids]) / sum(v("n_tested")[ids]))
      },
      TRUE,
      c(0, Inf)
    ),
    share_tested_removed_by_M9a = list(
      function(ids) 1 - sum(v("n_high")[ids]) / sum(v("n_tested")[ids]),
      FALSE,
      c(0, 1)
    )
  )
  for (qn in names(qs)) {
    r <- infer(
      ids,
      qs[[qn]][[1L]],
      log_scale = qs[[qn]][[2L]],
      range = qs[[qn]][[3L]]
    )
    r[, `:=`(
      design_set = ds,
      quantity = qn,
      threshold = "M1 p<0.01",
      n_calls_total = sum(ag$n_call),
      n_tested_total = sum(ag$n_tested)
    )]
    circ_rows[[length(circ_rows) + 1L]] <- r
  }
  rm(ag)
}
W(rbindlist(circ_rows, fill = TRUE), "q11_m9_class_vs_null_calls")

# 7d. Top/low tertile rate ratio by mapping-rate imbalance (confounding check) ----------
mr_rows <- list()
for (ds in c("7v7", "5v5", "3v3")) {
  d <- summ[
    design %in%
      design_sets[[ds]] &
      method == "M1" &
      set == "own" &
      stratum == "infrv_tertile" &
      level %in% c("high", "low")
  ]
  d <- merge(d, splits[, .(split_id, mapping_diff)], by = "split_id")
  med <- stats::median(abs(unique(d[, .(split_id, mapping_diff)])$mapping_diff))
  for (half in c("low", "high")) {
    dh <- d[
      if (half == "low") abs(mapping_diff) <= med else abs(mapping_diff) > med
    ]
    hi <- dh[level == "high"]
    lo <- dh[level == "low"]
    ids <- intersect(hi$split_id, lo$split_id)
    nh <- stats::setNames(as.numeric(hi$n_p01), hi$split_id)
    th <- stats::setNames(as.numeric(hi$n_tested), hi$split_id)
    nl <- stats::setNames(as.numeric(lo$n_p01), lo$split_id)
    tl <- stats::setNames(as.numeric(lo$n_tested), lo$split_id)
    f <- function(ids) {
      (sum(nh[ids]) / sum(th[ids])) / (sum(nl[ids]) / sum(tl[ids]))
    }
    r <- infer(ids, f, log_scale = TRUE)
    r[, `:=`(
      design_set = ds,
      method = "M1",
      mapping_imbalance = half,
      median_abs_mapping_diff = med,
      rate_p01_low_and_high_tertiles = (sum(nh[ids]) + sum(nl[ids])) /
        (sum(th[ids]) + sum(tl[ids])),
      quantity = "rate_ratio_high_vs_low_tertile_p01"
    )]
    mr_rows[[length(mr_rows) + 1L]] <- r
  }
}
W(
  rbindlist(mr_rows, fill = TRUE),
  "robustness_tertile_ratio_by_mapping_imbalance"
)

# 8. Decision rule R-A ------------------------------------------------------------------------
# Pre-registered model (amendment A6) and pre-registered contrasts: "M2 or M6
# reduces it", where "it" is the M1 top-tertile RR, so both reductions are
# against M1 (amendment A7). M6/M7 is reported as a diagnostic only: M7 (swish
# on 30 tied copies of the point counts) makes far more null calls than M1 and
# is not a calibrated reference.
rr_of <- function(mod) {
  glm_tab[
    design_set == "7v7" &
      method == "M1" &
      model == mod &
      term == "infrv_tertilehigh"
  ]
}
rr_m1 <- rr_of("prereg")
rr_m1_ident <- rr_of("ident")
h2t <- rbindlist(h2, fill = TRUE)[design_set == "7v7"]
red_m2 <- h2t[model == "prereg" & reference == "M1" & alternative == "M2"]
red_m6 <- h2t[model == "prereg" & reference == "M1" & alternative == "M6"]
diag_m6_m7 <- h2t[model == "prereg" & reference == "M7" & alternative == "M6"]
infl <- q12[design_set == "7v7" & set == "own" & quantity == "inflation_p01"]
m1_top <- infl[method == "M1" & stratum == "infrv_tertile" & level == "high"]
m11 <- infl[method == "M11raw" & stratum == "all"]
crit1 <- nrow(rr_m1) == 1L && isTRUE(rr_m1$ci_lo > 1)
crit2_m2 <- nrow(red_m2) == 1L && isTRUE(red_m2$ci_hi < 1)
crit2_m6 <- nrow(red_m6) == 1L && isTRUE(red_m6$ci_hi < 1)
crit2 <- crit2_m2 || crit2_m6
# (3) TRUE when the M11 CI lies below the M1 top-tertile estimate, FALSE when
# it lies at or above it, NA (not informative) when it contains it
crit3 <- if (nrow(m11) == 1L && nrow(m1_top) == 1L && is.finite(m11$ci_hi)) {
  if (m11$ci_hi < m1_top$estimate) {
    TRUE
  } else if (m11$ci_lo >= m1_top$estimate) {
    FALSE
  } else {
    NA
  }
} else {
  NA
}
# all criteria must hold: FALSE as soon as one is FALSE, NA if undecided
r_a <- if (isFALSE(crit1) || isFALSE(crit2) || isFALSE(crit3)) {
  FALSE
} else if (anyNA(c(crit1, crit2, crit3))) {
  NA
} else {
  TRUE
}
fmt <- function(x, d = 3) {
  if (length(x) && !is.na(x)) format(signif(x, d)) else "NA"
}
lines <- c(
  "Decision rule R-A (prereg.yaml): ignoring uncertainty inflates false isoform calls.",
  sprintf(
    "(1) M1 adjusted RR top vs low InfRV tertile (pre-registered model: quasi-Poisson on expression decile and isoform class, p<0.01 counts, 7v7 louo): %s, 95%% CI %s-%s (%s; raw %s-%s). CI excludes 1 and > 1: %s",
    fmt(rr_m1$estimate),
    fmt(rr_m1$ci_lo),
    fmt(rr_m1$ci_hi),
    rr_m1$ci_type,
    fmt(rr_m1$jk_lo_raw),
    fmt(rr_m1$jk_hi_raw),
    crit1
  ),
  sprintf(
    "(2a) M2/M1 ratio of top-tertile RR: %s, CI %s-%s (%s); (2b) M6/M1: %s, CI %s-%s (%s). A reduction with CI excluding 1: %s (M2/M1 %s, M6/M1 %s)",
    fmt(red_m2$estimate),
    fmt(red_m2$ci_lo),
    fmt(red_m2$ci_hi),
    red_m2$ci_type,
    fmt(red_m6$estimate),
    fmt(red_m6$ci_lo),
    fmt(red_m6$ci_hi),
    red_m6$ci_type,
    crit2,
    crit2_m2,
    crit2_m6
  ),
  sprintf(
    "(3) M11 raw gene-level inflation at p<0.01: %s (CI %s-%s) against M1 top-tertile inflation %s: no comparable excess (M11 upper CI below M1 top estimate; NA = not informative, the CI contains it): %s",
    fmt(m11$estimate),
    fmt(m11$ci_lo),
    fmt(m11$ci_hi),
    fmt(m1_top$estimate),
    crit3
  ),
  "Criterion (2) is moot when (1) fails: a reduction of an excess that is not there is not evidence for R-A.",
  sprintf("R-A holds: %s", r_a),
  sprintf(
    "Sensitivity (not pre-registered, A6): M1 RR with non-identifiable added: %s, CI %s-%s (%s).",
    fmt(rr_m1_ident$estimate),
    fmt(rr_m1_ident$ci_lo),
    fmt(rr_m1_ident$ci_hi),
    rr_m1_ident$ci_type
  ),
  sprintf(
    "Diagnostic only (A7): M6/M7 ratio of top-tertile RR %s, CI %s-%s (%s); M7 is swish on tied copies of the point counts and not a calibrated reference.",
    fmt(diag_m6_m7$estimate),
    fmt(diag_m6_m7$ci_lo),
    fmt(diag_m6_m7$ci_hi),
    diag_m6_m7$ci_type
  )
)
writeLines(lines, file.path(p$summary, "verdict.txt"))
W(
  data.table(
    criterion = c(
      "rr_top_excludes_1",
      "aware_reduces",
      "gene_level_no_excess",
      "R_A"
    ),
    value = c(crit1, crit2, crit3, r_a),
    text = lines[c(2:4, 6)]
  ),
  "verdict_R_A"
)
cat(lines, sep = "\n")
ue_write_session(p$summary, ue)
wpa_log("done")
