# WP-B (placenta audit): shared helpers (sourced, not run).
#
# Purpose: the loader of the S0 placenta store (10_placenta_prepare.R) into
# plain matrices, and the matrix-level statistics that the per-draw refits
# (30_audit_prevalence.R) and the per-split loops
# (31_halfsplit_replication.R) run thousands of times, so that no
# IsoformicExperiment is built or subset inside a loop. Every helper either
# calls the package code of the frozen snapshot directly (internal functions
# are available through `devtools::load_all(export_all = TRUE)` in
# ue_setup()) or restates a package rule line by line;
# 30_audit_prevalence.R (stage `validate`) checks each restatement against
# the package function on the full cohort and writes the comparison to
# placenta_audit/prevalence/validation.parquet.
#
# Inputs: the S0 store, $ISOFORMIC_UE_OUT/store/ (qc, control_units,
#   strata_tx, strata_gene, point_counts, moments, draws/<run>.parquet,
#   identifiability_pe/<run>.parquet).
# Env vars: those of ue_common.R; ISOFORMIC_UE_PILOT ("true": restrict U to
#   the transcripts of 1,500 random genes plus the story genes and write to
#   placenta_audit_pilot/, for code checks).
# Use: source("data-raw/uncertainty_evidence/ue_common.R") first (ue_setup()
#   loads the snapshot), then this file.

wpb_pilot <- function() {
  identical(tolower(Sys.getenv("ISOFORMIC_UE_PILOT")), "true")
}

wpb_dir <- function(...) {
  base <- if (wpb_pilot()) "placenta_audit_pilot" else "placenta_audit"
  file.path(ue_out_root(), base, ...)
}

wpb_label <- paste(
  "sensitivity of conclusions to quantification uncertainty",
  "(self-consistency), not an error rate"
)

# Seeds of WP-B: base + 1000 * 2 + unit (prereg rule, package index B = 2)
wpb_seed <- function(ue, unit) {
  as.integer(ue$seed_base + 1000L * ue$prereg$seeds$package_index$B + unit)
}

wpb_logger <- function(path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  function(...) {
    msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " ", sprintf(...))
    cat(msg, "\n", file = path, append = TRUE)
    message(msg)
  }
}

# Parquet reader (DuckDB, the package reader), provenance columns dropped
wpb_read <- function(path) {
  x <- as.data.frame(read_parquet(path))
  x[,
    setdiff(colnames(x), c("git_sha", "snapshot_sha256", "seed", "script")),
    drop = FALSE
  ]
}

wpb_query <- function(sql, tables = list()) {
  con <- DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  DBI::dbExecute(con, sprintf("SET threads = %d", ue_threads()))
  for (nm in names(tables)) {
    duckdb::duckdb_register(con, nm, tables[[nm]])
  }
  DBI::dbGetQuery(con, sql)
}

# Store loader --------------------------------------------------------------------

#' The cohort of WP-B from the S0 store: the 8 preeclampsia (PE) samples and
#' the 15 primary control units (the deeper member of each A/B pair), the
#' universe U (strata_tx$in_u), point counts and draws of U, library sizes
#' (sum of the folded counts of all units, as in S0), and optionally the
#' moments of every folded unit (for the overdispersion, as catchSalmon).
wpb_load_store <- function(draws = TRUE, moments = TRUE, say = message) {
  store <- file.path(ue_out_root(), "store")
  qc <- wpb_read(file.path(store, "qc.parquet"))
  units <- qc[
    qc$condition == "preeclampsia" | qc$control_unit_primary %in% TRUE,
  ]
  units <- units[order(units$condition != "preeclampsia", units$sample_id), ]
  runs <- units$sample_id
  stopifnot(length(runs) == 23L, sum(units$condition == "preeclampsia") == 8L)
  strata <- wpb_read(file.path(store, "strata_tx.parquet"))
  tx <- strata[strata$in_u %in% TRUE, ]
  if (wpb_pilot()) {
    set.seed(20261117L)
    genes <- unique(tx$gene_id)
    keep_g <- unique(c(
      sample(genes, min(1500L, length(genes))),
      tx$gene_id[tx$gene_name %in% c("FLT1", "ENG", "LEP")]
    ))
    tx <- tx[tx$gene_id %in% keep_g, ]
  }
  tx <- tx[order(tx$gene_id, tx$transcript_id), ]
  U <- tx$transcript_id
  rownames(tx) <- U
  say(
    "store: %d units (%d PE, %d control), |U| = %d transcripts of %d genes",
    length(runs),
    sum(units$condition == "preeclampsia"),
    sum(units$condition == "control"),
    length(U),
    length(unique(tx$gene_id))
  )
  u_tab <- data.frame(transcript_id = U)
  r_tab <- data.frame(sample_id = runs)
  pc <- wpb_query(
    sprintf(
      "SELECT p.transcript_id, p.sample_id, p.count
       FROM read_parquet('%s') p
       JOIN u_tab USING (transcript_id) JOIN r_tab USING (sample_id)",
      file.path(store, "point_counts.parquet")
    ),
    list(u_tab = u_tab, r_tab = r_tab)
  )
  counts <- matrix(0, length(U), length(runs), dimnames = list(U, runs))
  counts[cbind(
    match(pc$transcript_id, U),
    match(pc$sample_id, runs)
  )] <- pc$count
  rm(pc)
  lib <- wpb_query(
    sprintf(
      "SELECT sample_id, sum(count) AS lib FROM read_parquet('%s')
       JOIN r_tab USING (sample_id) GROUP BY sample_id",
      file.path(store, "point_counts.parquet")
    ),
    list(r_tab = r_tab)
  )
  lib_size <- stats::setNames(lib$lib[match(runs, lib$sample_id)], runs)
  out <- list(
    units = units,
    runs = runs,
    groups = factor(units$condition, c("control", "preeclampsia")),
    sex = factor(units$fetal_sex),
    tx = tx,
    U = U,
    gene_of = stats::setNames(tx$gene_id, U),
    counts = counts,
    lib_size = lib_size,
    strata = strata
  )
  if (moments) {
    mm <- wpb_query(
      sprintf(
        "SELECT m.transcript_id, m.sample_id, m.mean, m.var
         FROM read_parquet('%s') m JOIN r_tab USING (sample_id)",
        file.path(store, "moments.parquet")
      ),
      list(r_tab = r_tab)
    )
    ids <- sort(unique(mm$transcript_id))
    i <- cbind(match(mm$transcript_id, ids), match(mm$sample_id, runs))
    mu <- matrix(0, length(ids), length(runs), dimnames = list(ids, runs))
    v <- mu
    mu[i] <- mm$mean
    v[i] <- mm$var
    rm(mm, i)
    out$mom <- list(mean = mu, var = v, n = 30L)
    say("moments: %d folded units x %d samples", length(ids), length(runs))
  }
  if (draws) {
    arr <- array(
      0,
      c(length(U), length(runs), 30L),
      dimnames = list(U, runs, NULL)
    )
    for (j in seq_along(runs)) {
      d <- wpb_query(
        sprintf(
          "SELECT d.transcript_id, d.draw, d.count FROM read_parquet('%s') d
           JOIN u_tab USING (transcript_id)",
          file.path(store, "draws", paste0(runs[j], ".parquet"))
        ),
        list(u_tab = u_tab)
      )
      stopifnot(nrow(d) == 30L * length(U))
      arr[cbind(match(d$transcript_id, U), j, d$draw)] <- d$count
    }
    rm(d)
    out$arr <- arr
    say(
      "draws: %d x %d x %d array (%.2f GB)",
      dim(arr)[1],
      dim(arr)[2],
      dim(arr)[3],
      prod(dim(arr)) * 8 / 1024^3
    )
  }
  invisible(gc(verbose = FALSE))
  out
}

# Per-sample identifiability tables of the PE samples (S0), summed over a set
# of samples as equivalence_identifiability() does: not identifiable when it
# fails in at least one sample with fragments.
wpb_identifiability <- function(samples) {
  dir <- file.path(ue_out_root(), "store", "identifiability_pe")
  tabs <- lapply(samples, function(s) {
    f <- file.path(dir, paste0(s, ".parquet"))
    if (!file.exists(f)) {
      return(NULL)
    }
    wpb_read(f)[, c("transcript_id", "n_samples", "n_identifiable")]
  })
  all <- do.call(rbind, tabs)
  if (is.null(all)) {
    return(NULL)
  }
  ns <- tapply(all$n_samples, all$transcript_id, sum)
  ni <- tapply(all$n_identifiable, all$transcript_id, sum)
  data.frame(
    transcript_id = names(ns),
    n_samples = as.integer(ns),
    non_identifiable = as.integer(ni) < as.integer(ns)
  )
}

# edgeR quasi-likelihood ---------------------------------------------------------

#' edgeR v4 QL (DGEList, normLibSizes, glmQLFit, glmQLFTest), the route of
#' Baldoni et al. 2024. Rows: `keep` (logical or IDs) when given, else
#' edgeR::filterByExpr(group = groups) on the counts tested. Returns logFC,
#' p, BH q over the rows tested and a signed z (from the two-sided p).
wpb_edger <- function(
  counts,
  groups,
  keep = NULL,
  covariates = NULL,
  coef = NULL
) {
  if (is.null(keep)) {
    keep <- edgeR::filterByExpr(counts, group = groups)
  }
  if (is.character(keep)) {
    keep <- rownames(counts) %in% keep
  }
  y <- edgeR::DGEList(counts[keep, , drop = FALSE])
  y <- edgeR::normLibSizes(y)
  design <- wpb_design(groups, covariates)
  fit <- edgeR::glmQLFit(y, design)
  res <- edgeR::glmQLFTest(fit, coef = coef %||% ncol(design))
  tab <- res$table
  p <- tab$PValue
  data.frame(
    feature_id = rownames(y),
    logFC = tab$logFC,
    logCPM = tab$logCPM,
    pvalue = p,
    qvalue = stats::p.adjust(p, "BH"),
    z = sign(tab$logFC) * stats::qnorm(pmax(p, 1e-300) / 2, lower.tail = FALSE)
  )
}

# Gene sums of a transcript matrix
wpb_gene_sums <- function(m, gene_of) {
  rowsum(m, gene_of, reorder = TRUE)
}

# Design with an optional covariate data frame; the tested term is last
wpb_design <- function(groups, covariates = NULL) {
  if (is.null(covariates)) {
    stats::model.matrix(~groups)
  } else {
    df <- data.frame(covariates, groups = groups)
    stats::model.matrix(
      stats::reformulate(c(colnames(covariates), "groups")),
      data = df
    )
  }
}

# Overdispersion (package estimator, frozen; equal to edgeR::catchSalmon, S0
# validation) on a subset of samples, from the moments of all folded units.
wpb_overdispersion <- function(mom, samples) {
  od <- ue_overdispersion(
    mom$mean[, samples, drop = FALSE],
    mom$var[, samples, drop = FALSE],
    mom$n
  )
  out <- stats::setNames(od$overdispersion, od$transcript_id)
  attr(out, "prior") <- attr(od, "prior")
  out
}

# Divided counts; units without an estimate (no draws variance) keep 1
wpb_divide <- function(counts, od) {
  d <- od[rownames(counts)]
  d[is.na(d)] <- 1
  counts / d
}

# Isoform switch rule of isoform_switches() (package): a gene with a
# significant transcript up and one down (switch_type "opposite_direction";
# det_sig = BH q < 0.05 and |log2FC| >= 1, de_cutoffs() defaults).
wpb_switch_genes <- function(det, q = 0.05, lfc = 1) {
  sig <- !is.na(det$qvalue) & det$qvalue < q & abs(det$logFC) >= lfc
  up <- tapply(sig & det$logFC > 0, det$gene_id, any)
  down <- tapply(sig & det$logFC < 0, det$gene_id, any)
  names(up)[up & down]
}

# Usage test: the diffSplice backend of test_usage() (package functions
# filter_usage(), usage_diffsplice(), adjust_usage()); `od` divides the
# counts (scale = TRUE). `kept` fixes the rows (per-draw refits), and
# `gene_p = FALSE` skips the Simes gene p-values (BH over transcripts only,
# identical to adjust_usage()$qvalue).
wpb_usage <- function(
  counts,
  gene_of,
  groups,
  od = NULL,
  kept = NULL,
  gene_p = TRUE,
  min_count = 10,
  min_prop = 0.1
) {
  if (!is.null(od)) {
    counts <- wpb_divide(counts, od)
  }
  if (is.null(kept)) {
    kept <- rownames(counts)[filter_usage(
      counts,
      gene_of[rownames(counts)],
      groups,
      min_count,
      min_prop
    )]
  }
  counts <- counts[kept, , drop = FALSE]
  g <- unname(gene_of[kept])
  contrasts <- resolve_contrasts(NULL, levels(groups))
  res <- as.data.frame(usage_diffsplice(counts, g, groups, contrasts))
  res$gene_id <- g[match(res$transcript_id, kept)]
  if (gene_p) {
    res <- as.data.frame(adjust_usage(res))
  } else {
    res$qvalue <- stats::p.adjust(res$pvalue, "BH")
  }
  res <- res[, intersect(
    c(
      "transcript_id",
      "gene_id",
      "usage_log2FC",
      "pvalue",
      "qvalue",
      "gene_pvalue",
      "gene_qvalue"
    ),
    colnames(res)
  )]
  attr(res, "kept") <- kept
  res
}

# Productive fraction test (package two_group_proportion_test(), as in
# test_productive_fraction(); is_prod = productive or altered_productive,
# coding genes only, genes in `keep_genes`).
wpb_prod_test <- function(counts, gene_of, is_prod, groups, keep_genes) {
  total <- rowsum(counts, gene_of[rownames(counts)])
  productive <- rowsum(counts * is_prod, gene_of[rownames(counts)])
  k <- rownames(total) %in% keep_genes
  total <- total[k, , drop = FALSE]
  productive <- productive[k, , drop = FALSE]
  test <- two_group_proportion_test(productive, total, groups)
  data.frame(
    gene_id = rownames(total),
    mean_gene_count = rowMeans(total),
    delta = test$delta,
    statistic = test$statistic,
    pvalue = test$pvalue,
    qvalue = test$padj
  )
}

# Genes of test_productive_fraction(min_gene_count = m): coding genes whose
# mean gene count reaches m in every group
wpb_prod_genes <- function(counts, gene_of, groups, coding_genes, m = 10) {
  tot <- rowsum(counts, gene_of[rownames(counts)])
  pass <- rep(TRUE, nrow(tot))
  for (lv in levels(groups)) {
    pass <- pass & rowMeans(tot[, groups == lv, drop = FALSE]) >= m
  }
  intersect(rownames(tot)[pass], coding_genes)
}

# Confidence classes (compute_de_confidence(), transcripts, no collapse) ----

# Moments over the draws of the samples `s` (array tx x sample x draw)
wpb_moments <- function(arr, s = seq_len(dim(arr)[2L])) {
  a <- arr[, s, , drop = FALSE]
  n <- dim(a)[3L]
  mu <- rowMeans(a, dims = 2L)
  v <- (rowSums(a^2, dims = 2L) - n * mu^2) / (n - 1L)
  v <- pmax(v, 0)
  list(mean = mu, var = v, n = n)
}

# Per-draw Welch log2 fold-change on library-size-scaled draws, as
# compute_de_confidence(): scale = mean(lib) / lib, log2(x * scale + 5).
wpb_draw_welch <- function(arr, s, in_b, lib, pc = 5) {
  d <- dim(arr)[3L]
  scale <- mean(lib) / lib
  lfc <- matrix(NA_real_, dim(arr)[1L], d)
  for (j in seq_len(d)) {
    x <- matrix(arr[, s, j], nrow = dim(arr)[1L])
    x <- log2(sweep(x, 2L, scale, "*") + pc)
    lfc[, j] <- welch_matrix(x, in_b)$lfc
  }
  lfc
}

# Mean correlation of the draws of `ids` with their siblings (same gene, rows
# with a mean count of at least `min_count` over the samples `s`), averaged
# over the samples where both are expressed and vary
# (inferential_correlation() restricted to the pairs that touch `ids`).
wpb_sibling_cor <- function(arr, s, ids, gene_of, assessed, min_count = 1) {
  rn <- dimnames(arr)[[1L]]
  empty <- data.frame(
    gene_id = character(),
    transcript_id_1 = character(),
    transcript_id_2 = character(),
    mean_cor = numeric()
  )
  ids <- intersect(ids, rn[assessed])
  genes <- unique(gene_of[match(ids, rn)])
  cand <- which(assessed & gene_of %in% genes)
  if (!length(ids) || length(cand) < 2L) {
    return(empty)
  }
  g <- gene_of[cand]
  pairs <- do.call(
    rbind,
    lapply(split(cand, g), function(ix) {
      if (length(ix) < 2L) {
        return(NULL)
      }
      cmb <- utils::combn(ix, 2L)
      data.frame(a = cmb[1L, ], b = cmb[2L, ])
    })
  )
  if (is.null(pairs)) {
    return(empty)
  }
  pairs <- pairs[rn[pairs$a] %in% ids | rn[pairs$b] %in% ids, , drop = FALSE]
  if (!nrow(pairs)) {
    return(empty)
  }
  d <- dim(arr)[3L]
  r_sum <- numeric(nrow(pairs))
  r_n <- integer(nrow(pairs))
  rows <- sort(unique(c(pairs$a, pairs$b)))
  pa <- match(pairs$a, rows)
  pb <- match(pairs$b, rows)
  for (j in s) {
    x <- matrix(arr[rows, j, ], nrow = length(rows))
    mu <- rowMeans(x)
    cen <- x - mu
    sds <- sqrt(rowSums(cen^2) / (d - 1L))
    z <- cen / sds
    z[sds == 0 | mu < min_count, ] <- NA_real_
    for (start in seq.int(1L, length(pa), by = 50000L)) {
      ix <- start:min(start + 49999L, length(pa))
      r <- rowSums(z[pa[ix], , drop = FALSE] * z[pb[ix], , drop = FALSE]) /
        (d - 1L)
      ok <- !is.na(r)
      r_sum[ix[ok]] <- r_sum[ix[ok]] + r[ok]
      r_n[ix[ok]] <- r_n[ix[ok]] + 1L
    }
  }
  keep <- r_n > 0L
  data.frame(
    gene_id = gene_of[pairs$a][keep],
    transcript_id_1 = rn[pairs$a][keep],
    transcript_id_2 = rn[pairs$b][keep],
    mean_cor = r_sum[keep] / r_n[keep]
  )
}

#' Classes of the transcripts `ids` from the draws of the samples `s`, the
#' rules of add_uncertainty() + de_confidence() with their defaults:
#' stable = sign agreement of the per-draw Welch log2FC with its median >=
#' 0.95; uncertain = mean InfRV over `s` >= 0.5; ambiguous = a sibling with
#' mean draw correlation <= -0.5; class high / moderate / low for 0 / 1 / 2+
#' failures, NA when the mean count over `s` is below 1. Identifiability is
#' not part of the class (add_uncertainty(identifiability = FALSE), the
#' default). `draws` restricts the draws (ablation).
wpb_confidence <- function(
  arr,
  s,
  in_b,
  lib,
  ids,
  gene_of,
  draws = NULL,
  min_count = 1,
  inf_rv_cutoff = 0.5,
  cor_cutoff = -0.5,
  min_agreement = 0.95
) {
  rn <- dimnames(arr)[[1L]]
  sub <- if (is.null(draws)) {
    arr[, s, , drop = FALSE]
  } else {
    arr[, s, draws, drop = FALSE]
  }
  ss <- seq_along(s)
  mom <- wpb_moments(sub)
  mean_count <- rowMeans(mom$mean)
  rv_m <- inferential_rv(mom)
  rv <- rowMeans(rv_m)
  assessed <- mean_count >= min_count
  pairs <- wpb_sibling_cor(sub, ss, ids, gene_of, assessed, min_count)
  pairs$ambiguous <- pairs$mean_cor <= cor_cutoff
  amb_ids <- unique(c(
    pairs$transcript_id_1[pairs$ambiguous],
    pairs$transcript_id_2[pairs$ambiguous]
  ))
  long <- rbind(
    data.frame(id = pairs$transcript_id_1, cor = pairs$mean_cor),
    data.frame(id = pairs$transcript_id_2, cor = pairs$mean_cor)
  )
  min_cor <- if (nrow(long)) tapply(long$cor, long$id, min) else numeric()
  idx <- match(ids, rn)
  lfc <- wpb_draw_welch(sub[idx, , , drop = FALSE], ss, in_b, lib)
  med <- apply(lfc, 1L, stats::median)
  agree <- rowMeans(sign(lfc) == sign(med))
  stable <- agree >= min_agreement
  uncertain <- rv[idx] >= inf_rv_cutoff
  ambiguous <- ids %in% amb_ids
  n_fail <- as.integer(!stable) + as.integer(uncertain) + as.integer(ambiguous)
  cls <- c("high", "moderate", "low")[pmin(n_fail, 2L) + 1L]
  cls[!assessed[idx]] <- NA_character_
  list(
    table = data.frame(
      transcript_id = ids,
      mean_count = mean_count[idx],
      mean_inf_rv = rv[idx],
      lfc_median = med,
      lfc_sd = apply(lfc, 1L, stats::sd),
      sign_agreement = agree,
      stable = stable,
      uncertain = uncertain,
      ambiguous = ambiguous,
      min_sibling_cor = as.vector(min_cor[ids]),
      confidence = cls
    ),
    pairs = pairs,
    mean_inf_rv_all = stats::setNames(rv, rn),
    mean_count_all = stats::setNames(mean_count, rn)
  )
}

# Connected components of ambiguous pairs (package assign_ambiguity_groups())
wpb_groups <- function(pairs) {
  if (!nrow(pairs) || !any(pairs$ambiguous)) {
    return(data.frame(
      transcript_id = character(),
      ambiguity_group = character()
    ))
  }
  as.data.frame(assign_ambiguity_groups(pairs))
}

# Bootstrap and small helpers ------------------------------------------------------

wpb_ci <- function(x, probs = c(0.025, 0.975)) {
  if (all(is.na(x))) {
    return(NA_character_)
  }
  q <- stats::quantile(x, probs, na.rm = TRUE, names = FALSE)
  sprintf("[%.4g, %.4g]", q[1], q[2])
}

wpb_q <- function(x, p) {
  if (all(is.na(x))) {
    return(NA_real_)
  }
  stats::quantile(x, p, na.rm = TRUE, names = FALSE)
}

# Area under the ROC curve (Mann-Whitney), y logical, score numeric
wpb_auc <- function(score, y) {
  ok <- !is.na(score) & !is.na(y)
  score <- score[ok]
  y <- y[ok]
  # doubles: n1 * n0 overflows integers above ~46,000 calls per class
  n1 <- as.numeric(sum(y))
  n0 <- as.numeric(sum(!y))
  if (!n1 || !n0) {
    return(NA_real_)
  }
  r <- rank(score)
  (sum(r[y]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

# Storey pi1 = 1 - pi0 at lambda = 0.5
wpb_pi1 <- function(p) {
  p <- p[!is.na(p)]
  if (!length(p)) {
    return(NA_real_)
  }
  1 - min(1, mean(p > 0.5) / 0.5)
}

# One-sided p of the validation test in the direction of the discovery sign
wpb_one_sided <- function(lfc_v, p_v, sign_d) {
  ifelse(
    is.na(p_v) | is.na(lfc_v),
    NA_real_,
    ifelse(sign(lfc_v) == sign_d, p_v / 2, 1 - p_v / 2)
  )
}
