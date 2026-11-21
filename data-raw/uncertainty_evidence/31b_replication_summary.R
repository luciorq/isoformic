# WP-B (placenta audit), part 2b: statistics of the held-out replication.
#
# Sourced by 31_halfsplit_replication.R (stage "summarise"), which defines
# units, runs, splits, directed, U, txu, the helpers and say(). Not meant to
# be run alone.
#
# Pooled over the 60 directed pairs. Inference: two-level bootstrap
# (ISOFORMIC_UE_BOOT resamples, default 1,000; 500 for the secondary flags):
# the 30 half-splits are resampled with replacement (a split brings both its
# directions), then the genes of each directed pair with replacement (all the
# calls of a gene move together). The splits overlap (every pair of units
# shares a side in many splits), so the coverage of these intervals is
# approximate. The split-to-split 2.5-97.5% range of a statistic is reported
# as a range, not a CI.
#
# Definitions
#   replicated: held-out log2FC with the discovery sign and one-sided p <
#     0.05 (calls not testable in V count as not replicated; none for DTE,
#     because V tests every discovery call).
#   deficit(flag) = replication rate of unflagged calls - that of flagged
#     calls (positive: flagged calls replicate less).
#   matched deficit: the same within cells of the joint deciles of discovery
#     mean count x |discovery z| (p is a function of |z| in edgeR QL, so the
#     third margin of the plan, p, adds no cell), weighted by the flagged
#     calls of each cell.
#   adjusted OR: logistic regression replicated ~ flag + log2(mean count +
#     0.5) + |z| + log2(isoforms per gene) + transcript type class; OR < 1
#     means flagged calls replicate less.
#   delta-AUC: AUC for predicting non-replication of the logistic model
#     count + |z| + flag minus that of count + |z| (in sample).
#   PPV = P(not replicated | flagged); sensitivity = P(flagged | not
#     replicated); cost = replicated calls that the flag removes.
#   C4: per directed pair, as many calls as the primary flag, the lowest
#     discovery mean counts. This is the control R_B was judged on as
#     implemented; it is not the count-matched flag of prereg.yaml (same
#     mean-count distribution), a departure recorded after the results
#     (amendment A12). The matched deficit and the placebo below are the
#     closest analyses to the planned control.
#   C4b (stricter, not preregistered): per directed pair, as many calls as
#     the primary flag with the lowest replication predicted from count +
#     |z| alone (in-sample logistic fit on all pooled calls). A flag that
#     carries information beyond count and |z| should beat it; R_B is
#     still judged on C4, as implemented (amendment A12).
#   sign concordance outcome: the same adjusted OR with "held-out log2FC
#     has the discovery sign" as the outcome, which depends less on the
#     power of the validation half than replication at p < 0.05.
#   placebo: the primary flag permuted within (directed pair x cell), 200
#     times; its deficit and OR give the null distribution of a flag that
#     carries only count and |z|.
#   C1: label permutations of the discovery side (stage perm). A
#     permutation within fetal-sex blocks swaps n_relabelled samples; with 4
#     PE samples in D, 4 - n_relabelled / 2 true PE samples stay in the
#     permuted PE group, so most permutations keep part of the real
#     contrast (or reverse it). summary_c1_by_relabelled.parquet gives the
#     results per stratum with the residual contrast (share of true PE in
#     the permuted PE group minus that in the permuted control group); the
#     pooled C1 is contaminated with real signal and is not a clean null.
#   H3: the groups were selected because a member call did not replicate,
#     so they are reported as the share of those calls that summing the
#     group rescues, next to the replication of every group that contains
#     a discovery call (not compared with unselected single calls only).
#   R_B passes when the adjusted OR CI of the primary flag lies below 1 AND
#     the CI of deficit(primary) - deficit(C4) lies above 0.

seed_sum <- wpb_seed(ue, 400L)
n_boot2 <- max(100L, n_boot %/% 2L)
read_set <- function(suffix, dir = pair_dir) {
  f <- list.files(
    dir,
    pattern = paste0(suffix, "\\.parquet$"),
    full.names = TRUE
  )
  if (!length(f)) {
    return(NULL)
  }
  do.call(rbind, lapply(f, wpb_read))
}
P <- read_set("_dte")
info <- read_set("_info")
G <- read_set("_gene")
DTU <- read_set("_dtu")
GR <- read_set("_groups")
PM <- read_set("perm[0-9]+", perm_dir)
P$split <- as.integer(substr(P$pair_id, 2, 3))
say(
  "pooled: %d DTE calls over %d directed pairs (%d with fallback); %d with NA class",
  nrow(P),
  length(unique(P$pair_id)),
  sum(info$dte_fallback),
  sum(is.na(P$confidence))
)
ue_write(info, file.path(rep_dir, "pair_info.parquet"), ue, seed_sum)

prep <- function(d) {
  d <- d[!is.na(d$confidence), ]
  d$flag_primary <- d$confidence != "high"
  d$flag_unstable <- d$sign_agreement < 0.95
  d$flag_infrv_top_tertile <- d$infrv_top_tertile
  d$flag_ambiguous <- d$ambiguous
  d$flag_zero_unique <- d$zero_unique
  d$flag_non_identifiable <- d$non_identifiable
  d$flag_primary15 <- d$confidence15 != "high"
  d$log2_count <- log2(d$d_mean_count + 0.5)
  d$log2_iso <- log2(pmax(d$n_isoforms, 1))
  d$type_class <- factor(d$type_class)
  d
}
P <- prep(P)
P$cell <- interaction(
  dplyr::ntile(P$log2_count, 10),
  dplyr::ntile(P$d_absz, 10),
  drop = TRUE
)
c4_flag <- function(d) {
  out <- logical(nrow(d))
  for (ix in split(seq_len(nrow(d)), d$pair_id)) {
    k <- sum(d$flag_primary[ix])
    if (k > 0) out[ix[order(d$d_mean_count[ix])[seq_len(k)]]] <- TRUE
  }
  out
}
P$flag_c4 <- c4_flag(P)
c4b_flag <- function(d) {
  m <- suppressWarnings(stats::glm(
    replicated ~ log2_count + d_absz,
    data = d,
    family = stats::binomial()
  ))
  pr <- stats::fitted(m)
  out <- logical(nrow(d))
  for (ix in split(seq_len(nrow(d)), d$pair_id)) {
    k <- sum(d$flag_primary[ix])
    if (k > 0) out[ix[order(pr[ix])[seq_len(k)]]] <- TRUE
  }
  out
}
P$flag_c4b <- c4b_flag(P)
ue_write(P, file.path(rep_dir, "pooled_dte.parquet"), ue, seed_sum)

# Two-level bootstrap indices ----------------------------------------------------------
rows_by_pair_gene <- lapply(split(seq_len(nrow(P)), P$pair_id), function(ix) {
  split(ix, P$gene_id[ix])
})
pair_split <- as.integer(substr(names(rows_by_pair_gene), 2, 3))
boot_rows <- function() {
  sp <- sample(unique(P$split), replace = TRUE)
  unlist(
    lapply(sp, function(s) {
      lapply(which(pair_split == s), function(k) {
        gl <- rows_by_pair_gene[[k]]
        unlist(gl[sample.int(length(gl), replace = TRUE)], use.names = FALSE)
      })
    }),
    use.names = FALSE
  )
}
set.seed(seed_sum)
B_idx <- lapply(seq_len(n_boot), function(b) boot_rows())

# Statistics ------------------------------------------------------------------------------
cov_terms <- c("log2_count", "d_absz", "log2_iso", "type_class")
glm_or <- function(d, flags) {
  d <- d[stats::complete.cases(d[, flags, drop = FALSE]), ]
  if (length(unique(d$type_class)) < 2L) {
    f <- stats::reformulate(
      c(flags, setdiff(cov_terms, "type_class")),
      "replicated"
    )
  } else {
    f <- stats::reformulate(c(flags, cov_terms), "replicated")
  }
  co <- stats::coef(suppressWarnings(stats::glm(
    f,
    data = d,
    family = stats::binomial()
  )))
  unname(exp(co[paste0(flags, "TRUE")]))
}
deficit <- function(d, flag) {
  f <- d[[flag]]
  mean(d$replicated[f %in% FALSE]) - mean(d$replicated[f %in% TRUE])
}
matched_deficit <- function(d, flag) {
  f <- d[[flag]]
  ok <- !is.na(f)
  d <- d[ok, ]
  f <- f[ok]
  num <- 0
  w <- 0
  for (ix in split(seq_len(nrow(d)), d$cell)) {
    a <- ix[f[ix]]
    b <- ix[!f[ix]]
    if (length(a) && length(b)) {
      num <- num + length(a) * (mean(d$replicated[b]) - mean(d$replicated[a]))
      w <- w + length(a)
    }
  }
  num / w
}
auc_gain <- function(d, flag) {
  d <- d[!is.na(d[[flag]]), ]
  y <- !d$replicated
  m0 <- suppressWarnings(stats::glm(
    y ~ log2_count + d_absz,
    data = d,
    family = stats::binomial()
  ))
  m1 <- suppressWarnings(stats::glm(
    stats::reformulate(c("log2_count", "d_absz", flag), "y"),
    data = d,
    family = stats::binomial()
  ))
  c(auc0 = wpb_auc(stats::fitted(m0), y), auc1 = wpb_auc(stats::fitted(m1), y))
}
flag_stats <- function(d, flag) {
  a <- auc_gain(d, flag)
  c(
    deficit = deficit(d, flag),
    matched_deficit = matched_deficit(d, flag),
    adjusted_or = glm_or(d, flag),
    auc_base = a[["auc0"]],
    delta_auc = a[["auc1"]] - a[["auc0"]],
    deficit_minus_c4 = if (flag == "flag_primary") {
      deficit(d, flag) - deficit(d, "flag_c4")
    } else {
      NA_real_
    }
  )
}
par_boot <- function(fun, idx = B_idx) {
  do.call(
    rbind,
    parallel::mclapply(
      idx,
      function(ix) {
        tryCatch(fun(P[ix, ]), error = function(e) NA_real_)
      },
      mc.cores = workers
    )
  )
}

flags <- c(
  "flag_primary",
  "flag_unstable",
  "flag_infrv_top_tertile",
  "flag_ambiguous",
  "flag_zero_unique",
  "flag_non_identifiable",
  "flag_primary15",
  "flag_c4",
  "flag_c4b"
)
rows <- list()
boot_store <- list()
# Resume: a flag whose bootstrap is already stored is not recomputed. The
# indices B_idx are regenerated from the same seed in the same order, so a
# resumed run gives the same numbers as an uninterrupted one.
bf_path <- file.path(rep_dir, "bootstrap_flags.parquet")
sf_path <- file.path(rep_dir, "summary_flags.parquet")
prev_bt <- if (file.exists(bf_path) && file.exists(sf_path)) {
  wpb_read(bf_path)
} else {
  NULL
}
prev_rows <- if (!is.null(prev_bt)) wpb_read(sf_path) else NULL
for (fl in flags) {
  if (all(is.na(P[[fl]]))) {
    next
  }
  if (!is.null(prev_bt) && fl %in% prev_bt$flag && fl %in% prev_rows$flag) {
    boot_store[[fl]] <- prev_bt[prev_bt$flag == fl, ]
    rows[[fl]] <- prev_rows[prev_rows$flag == fl, ]
    say(
      "%s: bootstrap read from %s (%d resamples)",
      fl,
      bf_path,
      nrow(boot_store[[fl]])
    )
    next
  }
  d <- P[!is.na(P[[fl]]), ]
  f <- d[[fl]]
  t0 <- Sys.time()
  est <- flag_stats(P, fl)
  bt <- par_boot(
    function(x) flag_stats(x, fl),
    if (fl == "flag_primary") B_idx else B_idx[seq_len(n_boot2)]
  )
  boot_store[[fl]] <- data.frame(flag = fl, bt)
  per_pair <- vapply(
    split(d, d$pair_id),
    function(x) deficit(x, fl),
    numeric(1)
  )
  per_pair_or <- vapply(
    split(d, d$pair_id),
    function(x) {
      tryCatch(glm_or(x, fl), error = function(e) NA_real_)
    },
    numeric(1)
  )
  rows[[fl]] <- data.frame(
    flag = fl,
    n_calls = nrow(d),
    n_flagged = sum(f),
    pct_flagged = 100 * mean(f),
    rate_flagged = mean(d$replicated[f]),
    rate_unflagged = mean(d$replicated[!f]),
    deficit = est[["deficit"]],
    deficit_ci = wpb_ci(bt[, "deficit"]),
    matched_deficit = est[["matched_deficit"]],
    matched_deficit_ci = wpb_ci(bt[, "matched_deficit"]),
    adjusted_or = est[["adjusted_or"]],
    adjusted_or_ci = wpb_ci(bt[, "adjusted_or"]),
    auc_count_z = est[["auc_base"]],
    delta_auc = est[["delta_auc"]],
    delta_auc_ci = wpb_ci(bt[, "delta_auc"]),
    deficit_split_range_not_ci = sprintf(
      "[%.3g, %.3g]",
      wpb_q(per_pair, 0.025),
      wpb_q(per_pair, 0.975)
    ),
    adjusted_or_split_range_not_ci = sprintf(
      "[%.3g, %.3g]",
      wpb_q(per_pair_or, 0.025),
      wpb_q(per_pair_or, 0.975)
    ),
    ppv_nonreplication = mean(!d$replicated[f]),
    sensitivity_nonreplication = sum(f & !d$replicated) / sum(!d$replicated),
    replicated_removed_n = sum(f & d$replicated),
    replicated_removed_pct = 100 * sum(f & d$replicated) / sum(d$replicated),
    nonreplicated_removed_n = sum(f & !d$replicated),
    n_boot = nrow(bt)
  )
  ue_write(do.call(rbind, rows), sf_path, ue, seed_sum)
  ue_write(do.call(rbind, boot_store), bf_path, ue, seed_sum)
  say(
    "%s: %.1f%% flagged, OR %.3f %s, deficit %.4f %s, dAUC %.4f (%.1f min)",
    fl,
    100 * mean(f),
    est[["adjusted_or"]],
    wpb_ci(bt[, "adjusted_or"]),
    est[["deficit"]],
    wpb_ci(bt[, "deficit"]),
    est[["delta_auc"]],
    as.numeric(difftime(Sys.time(), t0, units = "mins"))
  )
}
flag_tab <- do.call(rbind, rows)
ue_write(flag_tab, file.path(rep_dir, "summary_flags.parquet"), ue, seed_sum)
ue_write(
  do.call(rbind, boot_store),
  file.path(rep_dir, "bootstrap_flags.parquet"),
  ue,
  seed_sum
)

# Joint model --------------------------------------------------------------------------------
joint_flags <- c(
  "flag_unstable",
  "flag_infrv_top_tertile",
  "flag_ambiguous",
  if (!all(is.na(P$flag_zero_unique))) "flag_zero_unique",
  if (!all(is.na(P$flag_non_identifiable))) "flag_non_identifiable"
)
j_est <- glm_or(P, joint_flags)
j_bt <- par_boot(function(x) glm_or(x, joint_flags), B_idx[seq_len(n_boot2)])
joint <- data.frame(
  flag = joint_flags,
  adjusted_or_joint = j_est,
  ci = apply(j_bt, 2, wpb_ci),
  n_complete = sum(stats::complete.cases(P[, joint_flags]))
)
ue_write(joint, file.path(rep_dir, "summary_joint_model.parquet"), ue, seed_sum)

# C4 ------------------------------------------------------------------------------------------
bt_p <- boot_store$flag_primary
c4 <- data.frame(
  comparison = "deficit(primary flag) - deficit(count-matched flag C4)",
  estimate = deficit(P, "flag_primary") - deficit(P, "flag_c4"),
  ci = wpb_ci(bt_p$deficit_minus_c4),
  deficit_primary = deficit(P, "flag_primary"),
  deficit_c4 = deficit(P, "flag_c4"),
  rate_c4_flagged = mean(P$replicated[P$flag_c4]),
  rate_primary_flagged = mean(P$replicated[P$flag_primary])
)
ue_write(c4, file.path(rep_dir, "summary_c4.parquet"), ue, seed_sum)

# C4b: same-size flag of the lowest predicted replication from count + |z|
# (the flag is fixed on the full data, as C4; the bootstrap resamples calls)
set.seed(seed_sum + 3L)
c4b_bt <- unlist(parallel::mclapply(
  B_idx,
  function(ix) {
    d <- P[ix, ]
    deficit(d, "flag_primary") - deficit(d, "flag_c4b")
  },
  mc.cores = workers
))
c4b <- data.frame(
  comparison = "deficit(primary flag) - deficit(same-size flag of lowest predicted replication from count + |z|, C4b)",
  estimate = deficit(P, "flag_primary") - deficit(P, "flag_c4b"),
  ci = wpb_ci(c4b_bt),
  deficit_primary = deficit(P, "flag_primary"),
  deficit_c4b = deficit(P, "flag_c4b"),
  rate_c4b_flagged = mean(P$replicated[P$flag_c4b]),
  rate_primary_flagged = mean(P$replicated[P$flag_primary]),
  overlap_primary_c4b = mean(P$flag_primary == P$flag_c4b),
  note = "C4b is not preregistered; its predicted replication is fitted in sample on the outcome (favourable to C4b, 2 parameters on 167k calls)"
)
ue_write(c4b, file.path(rep_dir, "summary_c4b.parquet"), ue, seed_sum)
print(c4b)

# Sign concordance as the outcome (less dependent on validation power)
glm_or_sc <- function(d, flag) {
  d <- d[!is.na(d[[flag]]), ]
  f <- stats::reformulate(c(flag, cov_terms), "sign_concordant")
  unname(exp(stats::coef(suppressWarnings(stats::glm(
    f,
    data = d,
    family = stats::binomial()
  )))[paste0(flag, "TRUE")]))
}
sc_rows <- lapply(
  c("flag_primary", "flag_unstable", "flag_non_identifiable"),
  function(fl) {
    d <- P[!is.na(P[[fl]]), ]
    bt <- unlist(parallel::mclapply(
      B_idx[seq_len(n_boot2)],
      function(ix) {
        tryCatch(glm_or_sc(P[ix, ], fl), error = function(e) NA_real_)
      },
      mc.cores = workers
    ))
    data.frame(
      flag = fl,
      sign_concordance_flagged = mean(d$sign_concordant[d[[fl]]]),
      sign_concordance_unflagged = mean(d$sign_concordant[!d[[fl]]]),
      adjusted_or_sign_concordance = glm_or_sc(P, fl),
      ci = wpb_ci(bt),
      n_boot = length(bt)
    )
  }
)
sc_tab <- do.call(rbind, sc_rows)
ue_write(
  sc_tab,
  file.path(rep_dir, "summary_sign_concordance.parquet"),
  ue,
  seed_sum
)
print(sc_tab)

# Placebo ------------------------------------------------------------------------------------------
set.seed(seed_sum + 1L)
cells <- split(seq_len(nrow(P)), interaction(P$pair_id, P$cell, drop = TRUE))
placebo <- do.call(
  rbind,
  parallel::mclapply(
    seq_len(200L),
    function(b) {
      set.seed(seed_sum + 1000L + b)
      f <- P$flag_primary
      for (ix in cells) {
        if (length(ix) > 1L) f[ix] <- f[ix][sample.int(length(ix))]
      }
      d <- P
      d$flag_placebo <- f
      c(
        deficit = deficit(d, "flag_placebo"),
        adjusted_or = glm_or(d, "flag_placebo")
      )
    },
    mc.cores = workers
  )
)
obs <- c(
  deficit = deficit(P, "flag_primary"),
  adjusted_or = glm_or(P, "flag_primary")
)
pl <- data.frame(
  statistic = c("deficit", "adjusted_or"),
  observed = obs,
  placebo_mean = colMeans(placebo, na.rm = TRUE),
  placebo_2.5 = apply(placebo, 2, wpb_q, 0.025),
  placebo_97.5 = apply(placebo, 2, wpb_q, 0.975),
  empirical_p_one_sided = c(
    (1 + sum(placebo[, 1] >= obs[1], na.rm = TRUE)) / 201,
    (1 + sum(placebo[, 2] <= obs[2], na.rm = TRUE)) / 201
  ),
  n_placebo = nrow(placebo)
)
ue_write(pl, file.path(rep_dir, "summary_placebo.parquet"), ue, seed_sum)
ue_write(
  as.data.frame(placebo),
  file.path(rep_dir, "placebo_draws.parquet"),
  ue,
  seed_sum
)

# C1: label permutation baseline -------------------------------------------------------------------
if (!is.null(PM) && nrow(PM)) {
  PM <- prep(PM)
  per_perm <- do.call(
    rbind,
    lapply(split(PM, PM$perm), function(x) {
      data.frame(
        perm = x$perm[1],
        pair_id = x$pair_id[1],
        n_calls = nrow(x),
        n_relabelled = x$n_relabelled[1],
        fallback = x$fallback[1],
        replication_rate = mean(x$replicated),
        pct_flagged = 100 * mean(x$flag_primary),
        deficit = deficit(x, "flag_primary")
      )
    })
  )
  ue_write(
    per_perm,
    file.path(rep_dir, "c1_per_permutation.parquet"),
    ue,
    seed_sum
  )
  strat <- function(d, s) {
    x <- d[d$fixed_infrv_tertile %in% s, ]
    c(
      n = nrow(x),
      rate = mean(x$replicated),
      pct_flagged = 100 * mean(x$flag_primary),
      deficit = deficit(x, "flag_primary")
    )
  }
  st <- do.call(
    rbind,
    lapply(c("low", "mid", "high", "all"), function(s) {
      lev <- if (s == "all") c("low", "mid", "high") else s
      o <- strat(P, lev)
      p <- strat(PM, lev)
      data.frame(
        fixed_infrv_tertile = s,
        n_observed = o[["n"]],
        rate_observed = o[["rate"]],
        pct_flagged_observed = o[["pct_flagged"]],
        deficit_observed = o[["deficit"]],
        n_perm = p[["n"]],
        rate_perm = p[["rate"]],
        pct_flagged_perm = p[["pct_flagged"]],
        deficit_perm = p[["deficit"]],
        excess_rate_over_c1 = o[["rate"]] - p[["rate"]],
        excess_deficit_over_c1 = o[["deficit"]] - p[["deficit"]]
      )
    })
  )
  ue_write(st, file.path(rep_dir, "summary_c1_strata.parquet"), ue, seed_sum)
  c1 <- data.frame(
    n_permutations = nrow(per_perm),
    n_perm_relabelled_0 = sum(per_perm$n_relabelled == 0),
    replication_rate_observed = mean(P$replicated),
    replication_rate_perm = mean(PM$replicated),
    replication_rate_perm_range_not_ci = sprintf(
      "[%.3g, %.3g]",
      wpb_q(per_perm$replication_rate, 0.025),
      wpb_q(per_perm$replication_rate, 0.975)
    ),
    deficit_observed = deficit(P, "flag_primary"),
    deficit_perm_pooled = deficit(PM, "flag_primary"),
    deficit_perm_range_not_ci = sprintf(
      "[%.3g, %.3g]",
      wpb_q(per_perm$deficit, 0.025),
      wpb_q(per_perm$deficit, 0.975)
    ),
    excess_deficit_over_c1 = deficit(P, "flag_primary") -
      deficit(PM, "flag_primary"),
    excess_deficit_ci = wpb_ci(bt_p$deficit - deficit(PM, "flag_primary")),
    note = "CI: two-level bootstrap of the observed deficit minus the pooled permutation deficit (its own variability is shown by the range)"
  )
  ue_write(c1, file.path(rep_dir, "summary_c1.parquet"), ue, seed_sum)
  # Per number of relabelled samples, with the residual real contrast
  n_d <- stats::setNames(
    lengths(strsplit(info$D, ",")),
    info$pair_id
  )
  PM$k_true_pe <- 4L - PM$n_relabelled %/% 2L
  PM$residual_contrast <- PM$k_true_pe /
    4 -
    (4L - PM$k_true_pe) / (n_d[PM$pair_id] - 4L)
  by_rel <- do.call(
    rbind,
    lapply(split(PM, PM$n_relabelled), function(x) {
      data.frame(
        n_relabelled = x$n_relabelled[1],
        true_pe_in_permuted_pe_group = x$k_true_pe[1],
        mean_residual_contrast = mean(tapply(
          x$residual_contrast,
          x$perm,
          `[`,
          1
        )),
        n_permutations = length(unique(x$perm)),
        n_calls = nrow(x),
        replication_rate = mean(x$replicated),
        pct_flagged = 100 * mean(x$flag_primary),
        deficit = deficit(x, "flag_primary")
      )
    })
  )
  by_rel$note <- "residual contrast = share of true PE in the permuted PE group - share in the permuted control group of D (1 = the real contrast, negative = reversed); no stratum is a clean null"
  ue_write(
    by_rel,
    file.path(rep_dir, "summary_c1_by_relabelled.parquet"),
    ue,
    seed_sum
  )
  print(by_rel)
}

# C3: gene level ----------------------------------------------------------------------------------
G$log2_count <- log2(G$d_mean_count + 0.5)
G$log2_maxrv <- log2(G$max_inf_rv)
G$split <- as.integer(substr(G$pair_id, 2, 3))
c3fit <- function(d) {
  unname(stats::coef(suppressWarnings(stats::glm(
    replicated ~ log2_maxrv + log2_count + d_absz,
    data = d,
    family = stats::binomial()
  )))["log2_maxrv"])
}
g_rows <- lapply(split(seq_len(nrow(G)), G$pair_id), function(ix) {
  split(ix, G$gene_id[ix])
})
g_split <- as.integer(substr(names(g_rows), 2, 3))
set.seed(seed_sum + 2L)
gB <- lapply(seq_len(n_boot), function(b) {
  sp <- sample(unique(G$split), replace = TRUE)
  unlist(
    lapply(sp, function(s) {
      lapply(which(g_split == s), function(k) {
        gl <- g_rows[[k]]
        unlist(gl[sample.int(length(gl), replace = TRUE)], use.names = FALSE)
      })
    }),
    use.names = FALSE
  )
})
c3_bt <- unlist(parallel::mclapply(
  gB,
  function(ix) tryCatch(c3fit(G[ix, ]), error = function(e) NA_real_),
  mc.cores = workers
))
tert_g <- cut(
  G$max_inf_rv,
  stats::quantile(G$max_inf_rv, c(0, 1 / 3, 2 / 3, 1), na.rm = TRUE),
  include.lowest = TRUE,
  labels = c("low", "mid", "high")
)
c3 <- data.frame(
  n_gene_calls = nrow(G),
  replication_rate_gene = mean(G$replicated),
  rate_low_max_infrv_tertile = mean(G$replicated[tert_g == "low"]),
  rate_high_max_infrv_tertile = mean(G$replicated[tert_g == "high"]),
  slope_log_odds_per_doubling_max_infrv = c3fit(G),
  ci = wpb_ci(c3_bt),
  replication_rate_transcript = mean(P$replicated),
  note = "logistic replicated ~ log2(max isoform InfRV of the gene, discovery draws) + log2 count + |z|"
)
ue_write(c3, file.path(rep_dir, "summary_c3.parquet"), ue, seed_sum)

# H3: groups of flagged non-replicating calls ----------------------------------------------------
if (!is.null(GR) && nrow(GR)) {
  key <- paste(P$pair_id, P$group_id)
  target <- P$flag_primary & !P$replicated & !is.na(P$group_id)
  tg <- unique(key[target])
  GR$key <- paste(GR$pair_id, GR$group_id)
  gsel <- GR[GR$key %in% tg, ]
  single <- P[!P$flag_primary & is.na(P$group_id), ]
  ta <- unique(key[!is.na(P$group_id)])
  gall <- GR[GR$key %in% ta, ]
  h3 <- data.frame(
    n_flagged_nonreplicating_calls = sum(P$flag_primary & !P$replicated),
    n_of_them_in_a_group = sum(target),
    n_groups_tested = nrow(gsel),
    group_replication_rate = mean(gsel$replicated),
    group_rate_ci_wilson_no_clustering = if (nrow(gsel)) {
      sprintf(
        "[%.3f, %.3f]",
        stats::prop.test(sum(gsel$replicated), nrow(gsel))$conf.int[1],
        stats::prop.test(sum(gsel$replicated), nrow(gsel))$conf.int[2]
      )
    } else {
      NA_character_
    },
    pct_of_flagged_nonreplicating_calls_rescued = 100 *
      sum(target & key %in% gsel$key[gsel$replicated]) /
      sum(P$flag_primary & !P$replicated),
    pct_selected_groups_discovery_p_lt_0.1 = 100 *
      mean(gsel$d_p < 0.1, na.rm = TRUE),
    n_groups_with_any_call = nrow(gall),
    all_groups_with_a_call_replication_rate = mean(gall$replicated),
    unflagged_single_rate = mean(single$replicated),
    n_unflagged_single = nrow(single),
    note = "group = connected component of ambiguous pairs (discovery draws); summed group re-tested on V. The selected groups were chosen because a member call did not replicate (selection on the outcome): read group_replication_rate as the share of those groups that summing rescues, not against unflagged single calls"
  )
  ue_write(h3, file.path(rep_dir, "summary_h3.parquet"), ue, seed_sum)
  ue_write(gsel, file.path(rep_dir, "h3_groups.parquet"), ue, seed_sum)
}

# Secondary outcomes, DTU, ablation ---------------------------------------------------------------
sec <- do.call(
  rbind,
  lapply(split(P, P$pair_id), function(x) {
    data.frame(
      pair_id = x$pair_id[1],
      n_calls = nrow(x),
      fallback = x$fallback[1],
      replication_rate = mean(x$replicated),
      sign_concordance = mean(x$sign_concordant),
      storey_pi1 = wpb_pi1(x$v_p1),
      spearman_logFC = suppressWarnings(stats::cor(
        x$d_logFC,
        x$v_logFC,
        method = "spearman",
        use = "complete.obs"
      )),
      pct_flagged = 100 * mean(x$flag_primary),
      rate_flagged = mean(x$replicated[x$flag_primary]),
      rate_unflagged = mean(x$replicated[!x$flag_primary])
    )
  })
)
ue_write(sec, file.path(rep_dir, "summary_per_pair.parquet"), ue, seed_sum)
sec_tab <- data.frame(
  metric = c(
    "replication_rate",
    "sign_concordance",
    "storey_pi1",
    "spearman_logFC"
  ),
  pooled_or_mean = c(
    mean(P$replicated),
    mean(P$sign_concordant),
    wpb_pi1(P$v_p1),
    mean(sec$spearman_logFC, na.rm = TRUE)
  ),
  split_range_not_ci = vapply(
    c("replication_rate", "sign_concordance", "storey_pi1", "spearman_logFC"),
    function(m) {
      sprintf("[%.3g, %.3g]", wpb_q(sec[[m]], 0.025), wpb_q(sec[[m]], 0.975))
    },
    character(1)
  ),
  flagged = c(
    mean(P$replicated[P$flag_primary]),
    mean(P$sign_concordant[P$flag_primary]),
    wpb_pi1(P$v_p1[P$flag_primary]),
    NA
  ),
  unflagged = c(
    mean(P$replicated[!P$flag_primary]),
    mean(P$sign_concordant[!P$flag_primary]),
    wpb_pi1(P$v_p1[!P$flag_primary]),
    NA
  )
)
ue_write(sec_tab, file.path(rep_dir, "summary_secondary.parquet"), ue, seed_sum)
if (!is.null(DTU) && nrow(DTU)) {
  DTU <- DTU[!is.na(DTU$confidence), ]
  fl <- DTU$confidence != "high"
  dtu <- data.frame(
    n_calls = nrow(DTU),
    n_pairs_fallback = length(unique(DTU$pair_id[DTU$fallback])),
    pct_tested_in_V = 100 * mean(DTU$v_tested),
    replication_rate = mean(DTU$replicated),
    pct_flagged = 100 * mean(fl),
    rate_flagged = mean(DTU$replicated[fl]),
    rate_unflagged = mean(DTU$replicated[!fl]),
    deficit = mean(DTU$replicated[!fl]) - mean(DTU$replicated[fl])
  )
  ue_write(dtu, file.path(rep_dir, "summary_dtu.parquet"), ue, seed_sum)
}
abl <- data.frame(
  class_agreement_30_vs_15 = mean(P$confidence == P$confidence15, na.rm = TRUE),
  primary_flag_agreement_30_vs_15 = mean(
    P$flag_primary == P$flag_primary15,
    na.rm = TRUE
  ),
  pct_flagged_30 = 100 * mean(P$flag_primary),
  pct_flagged_15 = 100 * mean(P$flag_primary15, na.rm = TRUE),
  infrv_tertile_flag_agreement_30_vs_15 = mean(
    P$infrv_top_tertile == P$infrv_top_tertile15,
    na.rm = TRUE
  ),
  unstable_flag_agreement_30_vs_15 = mean(
    (P$sign_agreement < 0.95) == (P$sign_agreement15 < 0.95),
    na.rm = TRUE
  )
)
ue_write(abl, file.path(rep_dir, "summary_ablation.parquet"), ue, seed_sum)

# Verdict R-B --------------------------------------------------------------------------------------
prim <- flag_tab[flag_tab$flag == "flag_primary", ]
or_ci <- wpb_q(bt_p$adjusted_or, c(0.025, 0.975))
c4_ci <- wpb_q(bt_p$deficit_minus_c4, c(0.025, 0.975))
verdict <- data.frame(
  rule = "R_B",
  adjusted_or = prim$adjusted_or,
  adjusted_or_ci = prim$adjusted_or_ci,
  or_ci_below_1 = or_ci[2] < 1,
  or_ci_excludes_1 = or_ci[2] < 1 | or_ci[1] > 1,
  deficit_minus_c4 = c4$estimate,
  deficit_minus_c4_ci = c4$ci,
  beats_count_matched = c4_ci[1] > 0,
  passed = (or_ci[2] < 1) && (c4_ci[1] > 0),
  definition = "passed = adjusted OR CI below 1 (flagged calls replicate less) AND CI of deficit(primary) - deficit(C4) above 0"
)
ue_write(
  verdict,
  file.path(rep_dir, "summary_verdict_RB.parquet"),
  ue,
  seed_sum
)
print(flag_tab[, c(
  "flag",
  "pct_flagged",
  "rate_flagged",
  "rate_unflagged",
  "deficit",
  "deficit_ci",
  "adjusted_or",
  "adjusted_or_ci",
  "delta_auc",
  "delta_auc_ci"
)])
print(verdict)
say("replication summaries written")
