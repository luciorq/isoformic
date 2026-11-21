# WP-B (placenta audit), part 1b: summaries of the prevalence stage.
#
# Sourced by 30_audit_prevalence.R (stage "summarise"), which loads the
# frozen package and the store objects (cohort, groups, U, gene_of, tx,
# counts_U, arr, ...). Not meant to be run alone. Label of every figure:
# sensitivity of conclusions to quantification uncertainty
# (self-consistency), not an error rate.
#
# Definitions
#   naive conclusion: F1 transcript q < 0.05; F2 switch gene; F3 transcript
#     q < 0.05 (raw usage test); F4 gene q < 0.05; G gene q < 0.05.
#   sign_agree: share of the 30 per-draw refits whose sign (log2FC, usage
#     log2FC, delta; for F2 the up member up and the down member down) is the
#     naive sign. sig_frac: share of refits with the naive sign and q < 0.05
#     (F2: the gene is still a switch).
#   lost by divided / swish / scaled: not a call (q < 0.05) of the
#     alternative, including features the alternative's filter removed.
#   class high: de_confidence() class of the package (F1, F3: the
#     transcript; F2: every naive up and down member high); F4: |delta_z| >=
#     1.96 (productive_fraction() over the draws).
#   p-matched gene reference: for each conclusion the gene-level test (any
#     gene tested in G, called or not) with the nearest naive p on the log10
#     scale (ties at random), and its own sign_agree / sig_frac / swish loss.
#   changed: sig_frac < 0.5, or lost by any alternative of the family.
#   q-matched gene reference: the same, matched on the naive BH q instead of
#     p. The families have different BH cut-offs (F3 about p 0.0016, gene
#     level about 0.007), so a p-matched reference can sit at a stricter q
#     than the conclusions; both references are reported.
#   lost reason: a lost conclusion is either not tested by the alternative
#     (its filter, e.g. filterByExpr on the divided counts, removed it) or
#     tested and not called (summary_lost_reason.parquet).
#
# Writes to placenta_audit/prevalence/:
#   conclusions_<family>.parquet, reference_gene_matched_<family>.parquet,
#   summary_family.parquet (MUST REPORT 1), summary_strata.parquet,
#   summary_changed_q001.parquet, sex_sensitivity.parquet,
#   null_link.parquet, null_link_by_decile.parquet, summary_lost_reason.parquet,
#   reference_gene_qmatched_<family>.parquet, f5_summary.parquet

seed_sum <- wpb_seed(ue, 4L)
set.seed(seed_sum)
n_boot_gene <- 2000L
n_boot_draw <- 500L
rd <- function(name) wpb_read(pq(name))
f1_raw <- rd("naive_f1_raw")
f1_div <- rd("naive_f1_div")
f1_swish <- rd("naive_f1_swish")
g_raw <- rd("naive_g_raw")
g_swish <- rd("naive_g_swish")
f3_raw <- rd("naive_f3_raw")
f3_scaled <- rd("naive_f3_scaled")
f4_raw <- rd("naive_f4_raw")
f4_div <- rd("naive_f4_div")
f2_all <- rd("naive_f2_switches")
conf <- rd("package_confidence")
conf <- conf[conf$feature_type == "transcript", ]
pfz <- rd("package_productive_fraction")
pd <- list(
  f1 = rd("perdraw_f1"),
  g = rd("perdraw_g"),
  f3 = rd("perdraw_f3"),
  f4 = rd("perdraw_f4")
)

# Per-draw matrix (rows ids, columns draws) of one statistic
draw_mat <- function(tab, id_col, stat_col, ids) {
  m <- matrix(NA_real_, length(ids), n_draw)
  i <- match(tab[[id_col]], ids)
  ok <- !is.na(i)
  m[cbind(i[ok], tab$draw[ok])] <- tab[[stat_col]][ok]
  m
}
agree_sig <- function(naive_sign, lfc_m, q_m) {
  agree <- sign(lfc_m) == naive_sign
  agree[is.na(agree)] <- FALSE
  list(agree = agree, sig = agree & !is.na(q_m) & q_m < 0.05)
}
called <- function(tab, id_col, ids) {
  q <- tab$qvalue[match(ids, tab[[id_col]])]
  !is.na(q) & q < 0.05
}
lost_reason <- function(tab, id_col, ids) {
  hit <- match(ids, tab[[id_col]])
  ifelse(
    is.na(hit),
    "not_tested",
    ifelse((tab$qvalue[hit] < 0.05) %in% TRUE, "called", "tested_not_called")
  )
}
tx_dec <- stats::setNames(tx$expr_decile, U)
sg <- wpb_read(file.path(ue_out_root(), "store", "strata_gene.parquet"))
sg <- sg[sg$gene_id %in% unique(gene_of), ]
g_dec <- stats::setNames(dplyr::ntile(sg$gene_cpm_pe, 10), sg$gene_id)
g_sig <- stats::setNames(g_raw$qvalue < 0.05, g_raw$feature_id)
cls <- stats::setNames(conf$confidence, conf$feature_id)

fam_tables <- list()
fam_draws <- list()

# G: every gene tested (reference pool); conclusions are q < 0.05
gl <- draw_mat(pd$g, "feature_id", "logFC", g_raw$feature_id)
gq <- draw_mat(pd$g, "feature_id", "qvalue", g_raw$feature_id)
rs_g_all <- agree_sig(sign(g_raw$logFC), gl, gq)
g_pool <- data.frame(
  family = "G_gene_DGE",
  feature_id = g_raw$feature_id,
  gene_id = g_raw$feature_id,
  naive_effect = g_raw$logFC,
  naive_p = g_raw$pvalue,
  naive_q = g_raw$qvalue,
  lost_divided = NA,
  lost_swish = !called(g_swish, "feature_id", g_raw$feature_id),
  sign_agree = rowMeans(rs_g_all$agree),
  sig_frac = rowMeans(rs_g_all$sig),
  confidence = NA_character_,
  mean_inf_rv = NA_real_,
  lfc_z = NA_real_,
  expr_decile = unname(g_dec[g_raw$feature_id]),
  gene_dge = g_raw$qvalue < 0.05,
  swish_status = lost_reason(g_swish, "feature_id", g_raw$feature_id)
)
gsel <- g_raw$qvalue < 0.05
fam_tables$G <- g_pool[gsel, ]
fam_draws$G <- list(
  agree = rs_g_all$agree[gsel, , drop = FALSE],
  sig = rs_g_all$sig[gsel, , drop = FALSE]
)

# F1 DTE
f1 <- f1_raw[f1_raw$qvalue < 0.05, ]
rs1 <- agree_sig(
  sign(f1$logFC),
  draw_mat(pd$f1, "feature_id", "logFC", f1$feature_id),
  draw_mat(pd$f1, "feature_id", "qvalue", f1$feature_id)
)
cf <- conf[match(f1$feature_id, conf$feature_id), ]
fam_tables$F1 <- data.frame(
  family = "F1_DTE",
  feature_id = f1$feature_id,
  gene_id = gene_of[f1$feature_id],
  naive_effect = f1$logFC,
  naive_p = f1$pvalue,
  naive_q = f1$qvalue,
  lost_divided = !called(f1_div, "feature_id", f1$feature_id),
  lost_swish = !called(f1_swish, "feature_id", f1$feature_id),
  sign_agree = rowMeans(rs1$agree),
  sig_frac = rowMeans(rs1$sig),
  confidence = cf$confidence,
  mean_inf_rv = cf$mean_inf_rv,
  lfc_z = cf$lfc_median / cf$lfc_sd,
  expr_decile = unname(tx_dec[f1$feature_id]),
  gene_dge = g_sig[gene_of[f1$feature_id]] %in% TRUE,
  divided_status = lost_reason(f1_div, "feature_id", f1$feature_id),
  swish_status = lost_reason(f1_swish, "feature_id", f1$feature_id)
)
fam_draws$F1 <- rs1
gains <- list(
  F1 = c(
    divided = sum(f1_div$qvalue < 0.05 & !f1_div$feature_id %in% f1$feature_id),
    swish = sum(
      f1_swish$qvalue < 0.05 & !f1_swish$feature_id %in% f1$feature_id,
      na.rm = TRUE
    )
  )
)

# F2 switches: per draw, the strongest naive up member stays up and the
# strongest down member stays down (agree); the gene is still a switch (sig)
sw <- f2_all$gene_id[f2_all$method == "raw"]
det_sig <- f1_raw[f1_raw$qvalue < 0.05 & abs(f1_raw$logFC) >= 1, ]
det_sig$gene_id <- gene_of[det_sig$feature_id]
pick_tx <- function(g, up) {
  d <- det_sig[det_sig$gene_id == g & (det_sig$logFC > 0) == up, ]
  d$feature_id[which.max(abs(d$logFC))]
}
if (length(sw)) {
  sw_up <- vapply(sw, pick_tx, character(1), up = TRUE)
  sw_dn <- vapply(sw, pick_tx, character(1), up = FALSE)
  pl <- draw_mat(pd$f1, "feature_id", "logFC", unique(c(sw_up, sw_dn)))
  rownames(pl) <- unique(c(sw_up, sw_dn))
  pqm <- draw_mat(pd$f1, "feature_id", "qvalue", rownames(pl))
  rownames(pqm) <- rownames(pl)
  f2_agree <- matrix(
    pl[sw_up, , drop = FALSE] > 0 & pl[sw_dn, , drop = FALSE] < 0,
    nrow = length(sw)
  )
  f2_agree[is.na(f2_agree)] <- FALSE
  # still a switch in draw k: rule applied to the per-draw F1 results
  f2_sig <- vapply(
    seq_len(n_draw),
    function(k) {
      x <- pd$f1[pd$f1$draw == k, ]
      sw %in%
        wpb_switch_genes(data.frame(
          transcript_id = x$feature_id,
          gene_id = gene_of[x$feature_id],
          logFC = x$logFC,
          qvalue = x$qvalue
        ))
    },
    logical(length(sw))
  )
  f2_sig <- matrix(f2_sig, nrow = length(sw))
  all_high <- vapply(
    sw,
    function(g) {
      all(cls[det_sig$feature_id[det_sig$gene_id == g]] %in% "high")
    },
    logical(1)
  )
  sw_p <- vapply(
    sw,
    function(g) min(det_sig$pvalue[det_sig$gene_id == g]),
    numeric(1)
  )
  sw_q <- vapply(
    sw,
    function(g) min(det_sig$qvalue[det_sig$gene_id == g]),
    numeric(1)
  )
  fam_tables$F2 <- data.frame(
    family = "F2_switch",
    feature_id = sw,
    gene_id = sw,
    naive_effect = NA_real_,
    naive_p = sw_p,
    naive_q = sw_q,
    lost_divided = !sw %in% f2_all$gene_id[f2_all$method == "divided"],
    lost_swish = !sw %in% f2_all$gene_id[f2_all$method == "swish"],
    sign_agree = rowMeans(f2_agree),
    sig_frac = rowMeans(f2_sig),
    confidence = ifelse(all_high, "high", "not_high"),
    mean_inf_rv = NA_real_,
    lfc_z = NA_real_,
    expr_decile = unname(g_dec[sw]),
    gene_dge = g_sig[sw] %in% TRUE,
    up_member = sw_up,
    down_member = sw_dn
  )
  fam_draws$F2 <- list(agree = f2_agree, sig = f2_sig)
}
gains$F2 <- c(
  divided = length(setdiff(f2_all$gene_id[f2_all$method == "divided"], sw)),
  swish = length(setdiff(f2_all$gene_id[f2_all$method == "swish"], sw))
)

# F3 DTU
f3 <- f3_raw[f3_raw$qvalue < 0.05, ]
rs3 <- agree_sig(
  sign(f3$usage_log2FC),
  draw_mat(pd$f3, "transcript_id", "usage_log2FC", f3$transcript_id),
  draw_mat(pd$f3, "transcript_id", "qvalue", f3$transcript_id)
)
cf3 <- conf[match(f3$transcript_id, conf$feature_id), ]
fam_tables$F3 <- data.frame(
  family = "F3_DTU",
  feature_id = f3$transcript_id,
  gene_id = f3$gene_id,
  naive_effect = f3$usage_log2FC,
  naive_p = f3$pvalue,
  naive_q = f3$qvalue,
  lost_divided = !called(f3_scaled, "transcript_id", f3$transcript_id),
  lost_swish = NA,
  sign_agree = rowMeans(rs3$agree),
  sig_frac = rowMeans(rs3$sig),
  confidence = cf3$confidence,
  mean_inf_rv = cf3$mean_inf_rv,
  lfc_z = cf3$lfc_median / cf3$lfc_sd,
  expr_decile = unname(tx_dec[f3$transcript_id]),
  gene_dge = g_sig[f3$gene_id] %in% TRUE,
  divided_status = lost_reason(f3_scaled, "transcript_id", f3$transcript_id)
)
fam_draws$F3 <- rs3
gains$F3 <- c(
  divided = sum(
    f3_scaled$qvalue < 0.05 & !f3_scaled$transcript_id %in% f3$transcript_id
  ),
  swish = NA
)

# F4 productive fraction
f4 <- f4_raw[!is.na(f4_raw$qvalue) & f4_raw$qvalue < 0.05, ]
rs4 <- agree_sig(
  sign(f4$delta),
  draw_mat(pd$f4, "gene_id", "delta", f4$gene_id),
  draw_mat(pd$f4, "gene_id", "qvalue", f4$gene_id)
)
dz <- pfz$delta_z[match(f4$gene_id, pfz$gene_id)]
fam_tables$F4 <- data.frame(
  family = "F4_productive_fraction",
  feature_id = f4$gene_id,
  gene_id = f4$gene_id,
  naive_effect = f4$delta,
  naive_p = f4$pvalue,
  naive_q = f4$qvalue,
  lost_divided = !called(f4_div, "gene_id", f4$gene_id),
  lost_swish = NA,
  sign_agree = rowMeans(rs4$agree),
  sig_frac = rowMeans(rs4$sig),
  confidence = ifelse(abs(dz) >= 1.96, "high", "not_high"),
  mean_inf_rv = NA_real_,
  lfc_z = dz,
  expr_decile = unname(g_dec[f4$gene_id]),
  gene_dge = g_sig[f4$gene_id] %in% TRUE,
  divided_status = lost_reason(f4_div, "gene_id", f4$gene_id)
)
fam_draws$F4 <- rs4
gains$F4 <- c(
  divided = sum(
    f4_div$qvalue < 0.05 & !f4_div$gene_id %in% f4$gene_id,
    na.rm = TRUE
  ),
  swish = NA
)
gains$G <- c(
  divided = NA,
  swish = sum(
    g_swish$qvalue < 0.05 &
      !g_swish$feature_id %in%
        fam_tables$G$feature_id,
    na.rm = TRUE
  )
)

# p-matched gene reference (pool: every gene tested at gene level)
match_ref <- function(p, pool = g_pool$naive_p) {
  lg <- log10(pmax(pool, 1e-300))
  vapply(
    log10(pmax(p, 1e-300)),
    function(v) {
      d <- abs(lg - v)
      cand <- which(d == min(d))
      cand[sample.int(length(cand), 1L)]
    },
    integer(1)
  )
}

fams <- intersect(c("F1", "F2", "F3", "F4", "G"), names(fam_tables))
for (fam in fams) {
  ue_write(fam_tables[[fam]], pq(paste0("conclusions_", fam)), ue, seed_sum)
}

stats_of <- function(tab, agree = NULL, sig = NULL) {
  sa <- if (is.null(agree)) tab$sign_agree else rowMeans(agree)
  sf <- if (is.null(sig)) tab$sig_frac else rowMeans(sig)
  na_all <- function(x, v) if (all(is.na(x))) NA_real_ else v
  c(
    n = nrow(tab),
    pct_lost_divided = na_all(tab$lost_divided, 100 * mean(tab$lost_divided)),
    pct_lost_swish = na_all(tab$lost_swish, 100 * mean(tab$lost_swish)),
    pct_sign_agree_lt_0.95 = 100 * mean(sa < 0.95),
    pct_sign_agree_lt_0.5 = 100 * mean(sa < 0.5),
    mean_sig_frac = mean(sf),
    pct_sig_frac_lt_0.5 = 100 * mean(sf < 0.5),
    pct_high = na_all(tab$confidence, 100 * mean(tab$confidence %in% "high"))
  )
}
draw_stats <- c(
  "pct_sign_agree_lt_0.95",
  "pct_sign_agree_lt_0.5",
  "mean_sig_frac",
  "pct_sig_frac_lt_0.5"
)

summ <- list()
changed_tab <- list()
for (fam in fams) {
  tab <- fam_tables[[fam]]
  if (!nrow(tab)) {
    next
  }
  dr <- fam_draws[[fam]]
  est <- stats_of(tab)
  set.seed(seed_sum + match(fam, fams))
  ref_idx <- match_ref(tab$naive_p)
  ref <- g_pool[ref_idx, ]
  ref_dr <- list(
    agree = rs_g_all$agree[ref_idx, , drop = FALSE],
    sig = rs_g_all$sig[ref_idx, , drop = FALSE]
  )
  ref$confidence <- NA_character_
  ref$lost_divided <- NA
  ref_est <- stats_of(ref)
  ue_write(
    data.frame(
      family = fam,
      conclusion = tab$feature_id,
      conclusion_p = tab$naive_p,
      ref_gene = ref$feature_id,
      ref_p = ref$naive_p,
      ref_q = ref$naive_q,
      ref_sign_agree = ref$sign_agree,
      ref_sig_frac = ref$sig_frac,
      ref_lost_swish = ref$lost_swish
    ),
    pq(paste0("reference_gene_matched_", fam)),
    ue,
    seed_sum
  )
  genes <- unique(tab$gene_id)
  rows_by_gene <- split(seq_len(nrow(tab)), tab$gene_id)
  bg <- t(vapply(
    seq_len(n_boot_gene),
    function(b) {
      ix <- unlist(
        rows_by_gene[sample(genes, length(genes), replace = TRUE)],
        use.names = FALSE
      )
      c(stats_of(tab[ix, ]), stats_of(ref[ix, ]))
    },
    numeric(2L * length(est))
  ))
  bd <- t(vapply(
    seq_len(n_boot_draw),
    function(b) {
      k <- sample.int(n_draw, n_draw, replace = TRUE)
      c(
        stats_of(tab, dr$agree[, k, drop = FALSE], dr$sig[, k, drop = FALSE]),
        stats_of(
          ref,
          ref_dr$agree[, k, drop = FALSE],
          ref_dr$sig[, k, drop = FALSE]
        )
      )
    },
    numeric(2L * length(est))
  ))
  # q-matched reference: own seeds, after the p-matched bootstraps, so the
  # p-matched numbers do not depend on it
  set.seed(seed_sum + 100L + match(fam, fams))
  qref_idx <- match_ref(tab$naive_q, g_pool$naive_q)
  qref <- g_pool[qref_idx, ]
  qref_dr <- list(
    agree = rs_g_all$agree[qref_idx, , drop = FALSE],
    sig = rs_g_all$sig[qref_idx, , drop = FALSE]
  )
  qref$confidence <- NA_character_
  qref$lost_divided <- NA
  qref_est <- stats_of(qref)
  ue_write(
    data.frame(
      family = fam,
      conclusion = tab$feature_id,
      conclusion_q = tab$naive_q,
      ref_gene = qref$feature_id,
      ref_p = qref$naive_p,
      ref_q = qref$naive_q,
      ref_sign_agree = qref$sign_agree,
      ref_sig_frac = qref$sig_frac,
      ref_lost_swish = qref$lost_swish
    ),
    pq(paste0("reference_gene_qmatched_", fam)),
    ue,
    seed_sum
  )
  bgq <- t(vapply(
    seq_len(n_boot_gene),
    function(b) {
      ix <- unlist(
        rows_by_gene[sample(genes, length(genes), replace = TRUE)],
        use.names = FALSE
      )
      stats_of(qref[ix, ])
    },
    numeric(length(est))
  ))
  bdq <- t(vapply(
    seq_len(n_boot_draw),
    function(b) {
      k <- sample.int(n_draw, n_draw, replace = TRUE)
      stats_of(
        qref,
        qref_dr$agree[, k, drop = FALSE],
        qref_dr$sig[, k, drop = FALSE]
      )
    },
    numeric(length(est))
  ))
  nm <- names(est)
  for (i in seq_along(nm)) {
    summ[[length(summ) + 1L]] <- data.frame(
      family = fam,
      statistic = nm[i],
      estimate = unname(est[i]),
      ci_gene_cluster = wpb_ci(bg[, i]),
      ci_draw = if (nm[i] %in% draw_stats) wpb_ci(bd[, i]) else NA_character_,
      ref_gene_matched = unname(ref_est[i]),
      ref_ci_gene_cluster = wpb_ci(bg[, length(est) + i]),
      ref_ci_draw = if (nm[i] %in% draw_stats) {
        wpb_ci(bd[, length(est) + i])
      } else {
        NA_character_
      },
      ref_gene_q_matched = unname(qref_est[i]),
      ref_q_ci_gene_cluster = wpb_ci(bgq[, i]),
      ref_q_ci_draw = if (nm[i] %in% draw_stats) {
        wpb_ci(bdq[, i])
      } else {
        NA_character_
      },
      max_naive_p_of_conclusions = max(tab$naive_p),
      gained_divided = unname(gains[[fam]]["divided"]),
      gained_swish = unname(gains[[fam]]["swish"]),
      label = wpb_label
    )
  }
  changed <- tab$sig_frac < 0.5 |
    tab$lost_divided %in% TRUE |
    tab$lost_swish %in% TRUE
  changed_tab[[fam]] <- data.frame(
    family = fam,
    n_conclusions = nrow(tab),
    n_changed = sum(changed),
    pct_changed = 100 * mean(changed),
    n_changed_naive_q_lt_0.001 = sum(
      changed & tab$naive_q < 0.001,
      na.rm = TRUE
    ),
    pct_of_changed_with_naive_q_lt_0.001 = 100 *
      mean(tab$naive_q[changed] < 0.001, na.rm = TRUE),
    definition = "changed = sig_frac < 0.5 or lost by divided counts / scaled test or swish"
  )
}
summary_family <- do.call(rbind, summ)
ue_write(summary_family, pq("summary_family"), ue, seed_sum)
ue_write(do.call(rbind, changed_tab), pq("summary_changed_q001"), ue, seed_sum)
print(summary_family[, c(
  "family",
  "statistic",
  "estimate",
  "ci_gene_cluster",
  "ci_draw",
  "ref_gene_matched",
  "ref_gene_q_matched"
)])

# Why conclusions are lost: not tested by the alternative or tested and not called
lr_one <- function(fam, alternative, status) {
  if (is.null(status)) {
    return(NULL)
  }
  n <- nrow(fam_tables[[fam]])
  data.frame(
    family = fam,
    alternative = alternative,
    n_conclusions = n,
    n_lost = sum(status != "called"),
    n_lost_not_tested = sum(status == "not_tested"),
    n_lost_tested_not_called = sum(status == "tested_not_called"),
    pct_lost = 100 * mean(status != "called"),
    pct_lost_not_tested = 100 * mean(status == "not_tested"),
    pct_lost_tested_not_called = 100 * mean(status == "tested_not_called"),
    note = "not_tested: removed by the filter of the alternative (filterByExpr on divided counts, min_gene_count on divided F4 counts, or the scaled usage filter)"
  )
}
lost_tab <- rbind(
  lr_one("F1", "divided counts", fam_tables$F1$divided_status),
  lr_one("F1", "swish", fam_tables$F1$swish_status),
  lr_one("F3", "scaled test_usage", fam_tables$F3$divided_status),
  lr_one("F4", "divided counts", fam_tables$F4$divided_status),
  lr_one("G", "swish", fam_tables$G$swish_status)
)
ue_write(lost_tab, pq("summary_lost_reason"), ue, seed_sum)
print(lost_tab[, 1:9])

# Strata: expression decile and gene DGE status --------------------------------------
one_stratum <- function(t, fam, stratum, level) {
  data.frame(
    family = fam,
    stratum = stratum,
    level = level,
    n = nrow(t),
    pct_sign_agree_lt_0.95 = 100 * mean(t$sign_agree < 0.95),
    mean_sig_frac = mean(t$sig_frac),
    pct_sig_frac_lt_0.5 = 100 * mean(t$sig_frac < 0.5),
    pct_lost_divided = if (all(is.na(t$lost_divided))) {
      NA
    } else {
      100 * mean(t$lost_divided)
    },
    pct_lost_swish = if (all(is.na(t$lost_swish))) {
      NA
    } else {
      100 * mean(t$lost_swish)
    },
    pct_high = if (all(is.na(t$confidence))) {
      NA
    } else {
      100 * mean(t$confidence %in% "high")
    }
  )
}
strata_tab <- do.call(
  rbind,
  lapply(fams, function(fam) {
    tab <- fam_tables[[fam]]
    if (!nrow(tab)) {
      return(NULL)
    }
    rbind(
      do.call(
        rbind,
        lapply(split(tab, tab$expr_decile), function(t) {
          one_stratum(
            t,
            fam,
            "expression_decile",
            as.character(t$expr_decile[1])
          )
        })
      ),
      do.call(
        rbind,
        lapply(split(tab, tab$gene_dge), function(t) {
          one_stratum(t, fam, "gene_DGE_q_lt_0.05", as.character(t$gene_dge[1]))
        })
      )
    )
  })
)
ue_write(strata_tab, pq("summary_strata"), ue, seed_sum)

# Sex covariate sensitivity --------------------------------------------------------------
f1_sex <- rd("naive_f1_sex")
g_sex <- rd("naive_g_sex")
ov <- function(a, b) {
  length(intersect(
    a$feature_id[a$qvalue < 0.05],
    b$feature_id[b$qvalue < 0.05]
  ))
}
sx <- data.frame(
  analysis = c("F1 transcript", "G gene"),
  n_calls_condition = c(sum(f1_raw$qvalue < 0.05), sum(g_raw$qvalue < 0.05)),
  n_calls_sex_condition = c(
    sum(f1_sex$qvalue < 0.05),
    sum(g_sex$qvalue < 0.05)
  ),
  n_shared = c(ov(f1_raw, f1_sex), ov(g_raw, g_sex)),
  sex_by_group = paste(
    apply(table(sex, groups), 1, paste, collapse = "/"),
    collapse = "; "
  )
)
ue_write(sx, pq("sex_sensitivity"), ue, seed_sum)

# Non-circular link: flag rate of real calls against null calls ----------------------------
nls <- rd("null_link_splits")
nlc <- rd("null_link_calls")
real_p001 <- f1_raw$feature_id[f1_raw$pvalue < 0.001]
flag_rate <- function(ids) 100 * mean(!cls[ids] %in% "high")
# base rate: every transcript tested in F1 (assessed)
tested_cls <- cls[f1_raw$feature_id]
wpa_files <- list.files(
  file.path(ue_out_root(), "null_placenta"),
  pattern = "\\.parquet$",
  recursive = TRUE,
  full.names = TRUE
)
null_n <- sum(nls$n_calls_p001)
null_flag <- sum(nls$n_flag_p001)
wl <- wilson_ci <- function(x, n) {
  if (!n) {
    return(c(NA, NA))
  }
  z <- 1.959964
  p <- x / n
  den <- 1 + z^2 / n
  c0 <- (p + z^2 / (2 * n)) / den
  h <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / den
  100 * c(c0 - h, c0 + h)
}
nl <- data.frame(
  set = c(
    "real F1 calls q < 0.05",
    "real F1 calls p < 0.001",
    "null calls p < 0.001 (20 control-only 7 v 8 splits)",
    "null calls q < 0.05 (20 control-only splits)",
    "all transcripts tested in F1 (base rate)"
  ),
  n_calls = c(
    nrow(fam_tables$F1),
    length(real_p001),
    null_n,
    sum(nls$n_calls_q05),
    sum(!is.na(tested_cls))
  ),
  pct_flagged = c(
    flag_rate(fam_tables$F1$feature_id),
    flag_rate(real_p001),
    100 * null_flag / max(1, null_n),
    100 * sum(nls$n_flag_q05) / max(1, sum(nls$n_calls_q05)),
    100 * mean(!tested_cls[!is.na(tested_cls)] %in% "high")
  )
)
ci <- rbind(
  wl(sum(!cls[fam_tables$F1$feature_id] %in% "high"), nrow(fam_tables$F1)),
  wl(sum(!cls[real_p001] %in% "high"), length(real_p001)),
  wl(null_flag, null_n),
  wl(sum(nls$n_flag_q05), sum(nls$n_calls_q05)),
  wl(sum(!tested_cls[!is.na(tested_cls)] %in% "high"), sum(!is.na(tested_cls)))
)
nl$wilson_lo <- ci[, 1]
nl$wilson_hi <- ci[, 2]
nl$note <- c(
  rep("class of the package on the 23 cohort samples", 2),
  rep(
    "class from the draws of the 15 control samples of each split (calls pooled over splits; Wilson CI ignores clustering)",
    2
  ),
  "class of the package on the 23 cohort samples"
)
nl$wpa_parquet_files_found <- length(wpa_files)
print(nl)
ue_write(nl, pq("null_link"), ue, seed_sum)

# The same within expression deciles (deciles of the cohort mean count of the
# transcripts tested in F1): real calls are more highly expressed than null
# calls, and class high depends on expression. The null rate standardised to
# the decile mix of the real calls and the Mantel-Haenszel OR (flagged, null
# against real, stratified by decile) remove that difference. Caveats: real
# effects are larger and therefore more stable over the draws; the null
# classes come from the 15 control samples of a split, the real classes from
# the 23 cohort samples.
mc <- stats::setNames(conf$mean_count, conf$feature_id)
br <- stats::quantile(mc[f1_raw$feature_id], 0:10 / 10, na.rm = TRUE)
br[1] <- -Inf
br[11] <- Inf
dec_of <- function(ids) as.integer(cut(mc[ids], br, labels = FALSE))
real_ids <- fam_tables$F1$feature_id
real_dec <- dec_of(real_ids)
real_flag <- !cls[real_ids] %in% "high"
null_dec <- dec_of(nlc$transcript_id)
null_flag <- !nlc$confidence %in% "high"
nld <- do.call(
  rbind,
  lapply(1:10, function(k) {
    data.frame(
      decile = k,
      n_real = sum(real_dec == k, na.rm = TRUE),
      pct_flagged_real = 100 * mean(real_flag[real_dec %in% k]),
      n_null = sum(null_dec == k, na.rm = TRUE),
      pct_flagged_null = if (any(null_dec %in% k)) {
        100 * mean(null_flag[null_dec %in% k])
      } else {
        NA_real_
      }
    )
  })
)
w <- nld$n_real / sum(nld$n_real)
ok <- !is.na(nld$pct_flagged_null) & nld$n_null > 0 & nld$n_real > 0
mh_tab <- table(
  factor(
    c(rep("real", sum(!is.na(real_dec))), rep("null", sum(!is.na(null_dec)))),
    c("null", "real")
  ),
  factor(
    c(real_flag[!is.na(real_dec)], null_flag[!is.na(null_dec)]),
    c(TRUE, FALSE)
  ),
  factor(c(real_dec[!is.na(real_dec)], null_dec[!is.na(null_dec)]), 1:10)
)
mh_tab <- mh_tab[,, ok, drop = FALSE]
mh <- stats::mantelhaen.test(mh_tab)
nld_sum <- data.frame(
  pct_flagged_real = 100 * mean(real_flag),
  pct_flagged_null_crude = 100 * mean(null_flag),
  pct_flagged_null_standardised_to_real_deciles = sum(
    w[ok] * nld$pct_flagged_null[ok]
  ) /
    sum(w[ok]),
  n_deciles_null_flag_rate_above_real = sum(
    nld$pct_flagged_null[ok] > nld$pct_flagged_real[ok]
  ),
  n_deciles_compared = sum(ok),
  mh_or_flagged_null_vs_real = unname(mh$estimate),
  mh_or_ci = sprintf("[%.4g, %.4g]", mh$conf.int[1], mh$conf.int[2]),
  mh_p = mh$p.value,
  note = "MH OR > 1: null (false) calls are flagged more often than real calls of the same expression decile; CI ignores clustering of null calls within splits"
)
ue_write(nld, pq("null_link_by_decile"), ue, seed_sum)
ue_write(nld_sum, pq("null_link_decile_summary"), ue, seed_sum)
print(nld)
print(nld_sum)

# F5 summary (secondary) ------------------------------------------------------------------
e <- rd("f5_enrichment")
if (nrow(e)) {
  e$key <- paste(e$experiment, e$pathway)
  sig_l <- e$key[e$rank_by == "log2FC" & e$padj < 0.05]
  sig_z <- e$key[e$rank_by == "lfc_z" & e$padj < 0.05]
  f5 <- data.frame(
    n_tests_log2FC = sum(e$rank_by == "log2FC"),
    n_tests_lfc_z = sum(e$rank_by == "lfc_z"),
    collections = paste(sort(unique(e$experiment)), collapse = ","),
    n_sig_log2FC = length(sig_l),
    n_sig_lfc_z = length(sig_z),
    n_shared = length(intersect(sig_l, sig_z)),
    pct_log2FC_sig_kept_by_lfc_z = 100 * mean(sig_l %in% sig_z),
    note = "fgsea padj < 0.05; a test is one gene set in one transcript-type collection (experiment); no per-draw fgsea"
  )
  ue_write(f5, pq("f5_summary"), ue, seed_sum)
  print(f5)
}
say("summaries written")
