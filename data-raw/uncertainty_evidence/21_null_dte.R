# WP-A step 2: differential transcript expression on the placenta null splits.
#
# Purpose: run the DTE methods on control-versus-control splits
# (20_null_splits.R), where every call is false, and store per split the
# p-values of every tested transcript plus per-stratum call counts, so that
# 29_null_summarise.R can estimate false-call rates, where they fall, and
# whether uncertainty-aware methods remove them. Matrices only: no
# IsoformicExperiment subsetting in the loop. Per split, the tested set of
# the raw-count methods is U (S0 universe, from the PE runs) intersected with
# edgeR::filterByExpr() without design or group (so the filter does not see
# the labels); divided counts are filtered the same way on the divided
# counts (Baldoni et al. 2024).
#
# Methods (WP-A numbering; prereg.yaml `methods` uses another numbering:
# prereg M1 = M1, M2 = M2, M3 = M3, M4 = M4, M5 = M5, M6 = M6, prereg M7
# (de_confidence) = M9 here, prereg M11 = M11; M7, M8 and M10 here are not
# in the prereg list):
#   M1  edgeR QL (v4), raw counts
#   M2  edgeR QL, counts divided by the overdispersion of the split's own
#       samples (package overdispersion_from_moments() on the S0 moments of
#       every folded transcript, so the prior is that of catchSalmon)
#   M3, M4  limma-voom raw, divided (7v7 designs only)
#   M5  DESeq2 Wald (20 7v7 splits)
#   M6  fishpond swish on the 30 draws (7v7 and 5v5; nperms 100; seed per split)
#   M7  swish with every draw replaced by the point counts (same splits as M6;
#       30 copies: the pilot found swish p-values with 2 copies differ by up
#       to 0.099 from those with 30, because swish breaks ties with random
#       jitter per replicate)
#   M8  Wilcoxon rank-sum on TMM log CPM
#   M9  de_confidence (package rules, matrix form) over M1: (a) M1 calls with
#       class high; (b) M1 calls minus the uncertain or ambiguous ones
#   M10 M1 calls after removing as many lowest-expressed tested transcripts
#       as M9 removes ((a) and (b) matched separately)
#   M11 gene level: edgeR QL on gene sums, raw and divided by the gene
#       overdispersion from the summed draws
#
# Inputs: <out>/store (10_placenta_prepare.R), <out>/null_placenta/splits
# (20_null_splits.R). The first run builds a cache of matrices in
# <out>/null_placenta/cache (Parquet tables plus an .rds bundle of the same
# matrices for fast loading by the workers).
#
# Env vars:
#   ISOFORMIC_UE_OUT      output root (ue_common.R)
#   ISOFORMIC_UE_SPLITS   "<from>:<to>" row range of the selected splits
#                         (default all); or comma-separated split IDs
#   ISOFORMIC_UE_DESIGNS  comma-separated designs to run (default
#                         7v7_louo,7v7_mc,5v5,3v3,7v7_unres)
#   ISOFORMIC_UE_WORKERS  forked workers (default 4; each about 1-3 GB)
#   ISOFORMIC_UE_VARIANT  main (default), member (sensitivity members),
#                         draws16 (methods on draws 16-30), unfolded
#                         (haplotype copies not folded: primary rows of the
#                         unfolded store)
#   ISOFORMIC_UE_METHODS  comma-separated subset (default: per design, see
#                         methods_for())
#   ISOFORMIC_UE_SWISH_SPLITS  max swish splits per design (default all)
#   ISOFORMIC_UE_FIRST_LOUO  n: only the first n 7v7_louo splits in (rank in
#                         unit, left-out unit) order (robustness variants)
#   ISOFORMIC_UE_PILOT    1: run 3 splits per size through every method
#                         sequentially, record wall time and peak RSS per
#                         method, validate M9 against the package and M7 with
#                         2 against 30 copies; writes dte/pilot_*.parquet
#
# Resources (measured, peak RSS of the largest process): main variant 81 min
# with 4 workers, 4.8 GB (a rerun from the cache 10 min, 6.2 GB); the
# draws16, member and unfolded variants 3 to 15 min, 2.3 to 4.4 GB; pilot
# 8 min, 6.1 GB. The main run used a memory cap of 16 GB with 4 forked
# workers (see README.md); lower ISOFORMIC_UE_WORKERS under a smaller cap.
# Run from the repository root, for example:
#   ISOFORMIC_UE_PILOT=1 pixi run --frozen Rscript data-raw/uncertainty_evidence/21_null_dte.R
#   ISOFORMIC_UE_DESIGNS=7v7_louo ISOFORMIC_UE_WORKERS=4 \
#     pixi run --frozen Rscript data-raw/uncertainty_evidence/21_null_dte.R
# Outputs: <out>/null_placenta/dte/<variant>/tx/<split>.parquet (per
# transcript p-values and M9 flags), gene/<split>.parquet,
# summary/<split>.parquet (calls per method, set and stratum), hist/,
# timing/<split>.parquet.

source("data-raw/uncertainty_evidence/null_common.R")
ue <- ue_setup("21_null_dte")
p <- wpa_paths()
for (d in c(p$cache, p$dte, p$logs)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}
seed_pkg <- ue$seed_base + 1000L * ue$prereg$seeds$package_index$A
n_draws <- 30L
variant <- Sys.getenv("ISOFORMIC_UE_VARIANT", "main")
pilot <- identical(Sys.getenv("ISOFORMIC_UE_PILOT"), "1")
workers <- as.integer(Sys.getenv("ISOFORMIC_UE_WORKERS", "4"))

# 1. Cache of matrices -------------------------------------------------------------
cache_rds <- file.path(p$cache, "dte_cache.rds")
build_cache <- function() {
  wpa_log("building the cache")
  cu <- wpa_read(file.path(p$store, "control_units.parquet"))
  ctrl <- cu$sample_id
  strata <- wpa_read(file.path(p$store, "strata_tx.parquet"))
  U <- strata$transcript_id[strata$in_u]
  all_units <- strata$transcript_id
  gene_all <- stats::setNames(strata$gene_id, strata$transcript_id)
  inlist <- paste0("'", ctrl, "'", collapse = ",")
  pc <- wpa_query(sprintf(
    "SELECT transcript_id, sample_id, count, effective_length FROM read_parquet('%s') WHERE sample_id IN (%s)",
    file.path(p$store, "point_counts.parquet"),
    inlist
  ))
  counts <- long_to_matrix(
    pc$transcript_id,
    pc$sample_id,
    pc$count,
    all_units,
    ctrl
  )
  efflen <- long_to_matrix(
    pc$transcript_id,
    pc$sample_id,
    pc$effective_length,
    U,
    ctrl
  )
  rm(pc)
  mo <- wpa_query(sprintf(
    "SELECT transcript_id, sample_id, mean, var, mean_15, var_15 FROM read_parquet('%s') WHERE sample_id IN (%s)",
    file.path(p$store, "moments.parquet"),
    inlist
  ))
  mom_mean <- long_to_matrix(
    mo$transcript_id,
    mo$sample_id,
    mo$mean,
    all_units,
    ctrl
  )
  mom_var <- long_to_matrix(
    mo$transcript_id,
    mo$sample_id,
    mo$var,
    all_units,
    ctrl
  )
  m15 <- long_to_matrix(
    mo$transcript_id,
    mo$sample_id,
    mo$mean_15,
    all_units,
    ctrl
  )
  v15 <- long_to_matrix(
    mo$transcript_id,
    mo$sample_id,
    mo$var_15,
    all_units,
    ctrl
  )
  rm(mo)
  # Moments of draws 16-30 from those of 1-30 and 1-15 (sums of squares)
  m16 <- 2 * mom_mean - m15
  s2_30 <- 29 * mom_var + 30 * mom_mean^2
  s2_15 <- 14 * v15 + 15 * m15^2
  v16 <- pmax((s2_30 - s2_15) - 15 * m16^2, 0) / 14
  m16[m16 < 0] <- 0
  rm(m15, v15, s2_30, s2_15)
  # Draws of U, and gene sums of the draws of every stored unit
  arr <- array(
    0,
    c(length(U), length(ctrl), n_draws),
    dimnames = list(U, ctrl, NULL)
  )
  gene_ids <- sort(unique(gene_all))
  gmean <- matrix(
    0,
    length(gene_ids),
    length(ctrl),
    dimnames = list(gene_ids, ctrl)
  )
  gvar <- gmean
  gmean16 <- gmean
  gvar16 <- gmean
  for (j in seq_along(ctrl)) {
    f <- file.path(p$store, "draws", paste0(ctrl[[j]], ".parquet"))
    d <- wpa_query(sprintf(
      "SELECT transcript_id, draw, count FROM read_parquet('%s')",
      f
    ))
    i <- match(d$transcript_id, U)
    ok <- !is.na(i)
    arr[cbind(i[ok], j, d$draw[ok])] <- d$count[ok]
    g <- match(gene_all[d$transcript_id], gene_ids)
    gs <- matrix(0, length(gene_ids), n_draws)
    agg <- rowsum(
      as.numeric(d$count),
      g + (d$draw - 1L) * length(gene_ids),
      reorder = FALSE
    )
    gs[as.integer(rownames(agg))] <- agg[, 1L]
    gm <- replicate_moments(gs, n_draws)
    gmean[, j] <- gm$mean
    gvar[, j] <- gm$var
    gm <- replicate_moments(gs[, 16:30], 15L)
    gmean16[, j] <- gm$mean
    gvar16[, j] <- gm$var
    rm(d, gs, agg)
    wpa_log("draws read: ", ctrl[[j]])
  }
  # Within-gene pairs of U and their per-sample draw correlation (NA when a
  # member has no spread or a mean below 1 in that sample), on draws 1-30 and
  # 16-30
  gU <- gene_all[U]
  ord <- order(gU)
  gs <- gU[ord]
  rl <- rle(gs)
  end_of <- rep(cumsum(rl$lengths), rl$lengths)
  remaining <- end_of - seq_along(gs)
  first <- rep(seq_along(gs), remaining)
  second <- first + sequence(remaining)
  pairs <- data.frame(i1 = ord[first], i2 = ord[second])
  pair_cor <- function(slice) {
    out <- matrix(
      NA_real_,
      nrow(pairs),
      length(ctrl),
      dimnames = list(NULL, ctrl)
    )
    for (j in seq_along(ctrl)) {
      x <- arr[, j, slice]
      dd <- ncol(x)
      mu <- rowMeans(x)
      cen <- x - mu
      sds <- sqrt(rowSums(cen^2) / (dd - 1L))
      z <- cen / sds
      z[sds == 0 | mu < 1, ] <- NA_real_
      for (s in seq.int(1L, nrow(pairs), by = 500000L)) {
        idx <- s:min(s + 499999L, nrow(pairs))
        out[idx, j] <- rowSums(
          z[pairs$i1[idx], , drop = FALSE] * z[pairs$i2[idx], , drop = FALSE]
        ) /
          (dd - 1L)
      }
    }
    out
  }
  pc30 <- pair_cor(1:30)
  pc16 <- pair_cor(16:30)
  wpa_log("pairs: ", nrow(pairs))
  # Unfolded store: primary rows only (copies not added), for the sensitivity
  pu <- wpa_query(sprintf(
    "SELECT transcript_id, sample_id, count FROM read_parquet('%s') WHERE sample_id IN (%s)",
    file.path(p$store, "point_counts_unfolded.parquet"),
    inlist
  ))
  ids_u <- sort(unique(pu$transcript_id))
  counts_unf <- long_to_matrix(
    pu$transcript_id,
    pu$sample_id,
    pu$count,
    ids_u,
    ctrl
  )
  rm(pu)
  mu <- wpa_query(sprintf(
    "SELECT transcript_id, sample_id, mean, var FROM read_parquet('%s') WHERE sample_id IN (%s)",
    file.path(p$store, "moments_unfolded.parquet"),
    inlist
  ))
  umean <- long_to_matrix(mu$transcript_id, mu$sample_id, mu$mean, ids_u, ctrl)
  uvar <- long_to_matrix(mu$transcript_id, mu$sample_id, mu$var, ids_u, ctrl)
  rm(mu)
  cache <- list(
    ctrl = ctrl,
    U = U,
    all_units = all_units,
    gene_all = gene_all,
    counts = counts,
    efflen = efflen,
    mom_mean = mom_mean,
    mom_var = mom_var,
    mom_mean16 = m16,
    mom_var16 = v16,
    arr = arr,
    gene_ids = gene_ids,
    gmean = gmean,
    gvar = gvar,
    gmean16 = gmean16,
    gvar16 = gvar16,
    pairs = pairs,
    pc30 = pc30,
    pc16 = pc16,
    counts_unf = counts_unf,
    umean = umean,
    uvar = uvar
  )
  # Parquet copies of the derived tables (canonical; the .rds is a cache)
  ue_write(
    data.frame(transcript_id = U),
    file.path(p$cache, "universe.parquet"),
    ue
  )
  ue_write(
    data.frame(
      gene_id = rep(gene_ids, length(ctrl)),
      sample_id = rep(ctrl, each = length(gene_ids)),
      mean = as.vector(gmean),
      var = as.vector(gvar),
      mean_16_30 = as.vector(gmean16),
      var_16_30 = as.vector(gvar16)
    ),
    file.path(p$cache, "gene_moments.parquet"),
    ue
  )
  ue_write(
    data.frame(
      transcript_id_1 = U[pairs$i1],
      transcript_id_2 = U[pairs$i2],
      gene_id = gene_all[U[pairs$i1]],
      pc30
    ),
    file.path(p$cache, "pair_cor_draws1_30.parquet"),
    ue
  )
  ue_write(
    data.frame(
      transcript_id_1 = U[pairs$i1],
      transcript_id_2 = U[pairs$i2],
      pc16
    ),
    file.path(p$cache, "pair_cor_draws16_30.parquet"),
    ue
  )
  saveRDS(cache, cache_rds, compress = FALSE)
  cache
}
C <- if (file.exists(cache_rds)) readRDS(cache_rds) else build_cache()
wpa_log(
  "cache: |U| = ",
  length(C$U),
  ", controls ",
  length(C$ctrl),
  ", pairs ",
  nrow(C$pairs)
)

strata <- wpa_read(file.path(p$store, "strata_tx.parquet"))
strata <- strata[match(C$U, strata$transcript_id), ]
gstrata <- wpa_read(file.path(p$store, "strata_gene.parquet"))
splits <- wpa_read(file.path(p$splits, "splits.parquet"))
split_samples <- wpa_read(file.path(p$splits, "split_samples.parquet"))

# Strata (fixed, from the PE runs) -----------------------------------------------
st_vars <- list(
  all = rep("all", length(C$U)),
  infrv_tertile = strata$infrv_tertile,
  infrv_tertile_15 = strata$infrv_tertile_15,
  uncertain = as.character(strata$uncertain),
  non_identifiable = as.character(strata$non_identifiable),
  zero_unique = as.character(strata$zero_unique),
  isoform_bin = strata$isoform_bin,
  expr_decile = as.character(strata$expr_decile),
  chr_xy = as.character(strata$seqname %in% c("chrX", "chrY"))
)

# 2. Methods per design --------------------------------------------------------------
deseq_ids <- {
  l <- splits[splits$design == "7v7_louo", ]
  l <- l[order(l$left_out, l$rank_in_unit), ]
  c(l$split_id[l$rank_in_unit == 1L], l$split_id[l$rank_in_unit == 2L][1:5])
}
methods_for <- function(design, split_id) {
  m <- c("M1", "M2", "M8", "M9", "M11")
  if (startsWith(design, "7v7")) {
    m <- c(m, "M3", "M4")
  }
  # the second 14 splits per unit only serve the Monte Carlo check of M1, M2
  if (design == "7v7_mc") {
    m <- c("M1", "M2")
  }
  if (design %in% c("7v7_louo", "5v5")) {
    m <- c(m, "M6", "M7")
  }
  if (split_id %in% deseq_ids) {
    m <- c(m, "M5")
  }
  if (variant == "member") {
    m <- c("M1", "M2")
  }
  if (variant == "unfolded") {
    m <- c("M1", "M2")
  }
  if (variant == "draws16") {
    m <- intersect(m, c("M1", "M2", "M6", "M9"))
  }
  only <- Sys.getenv("ISOFORMIC_UE_METHODS")
  if (nzchar(only)) {
    m <- intersect(m, strsplit(only, ",", fixed = TRUE)[[1L]])
  }
  m
}

design_free_keep <- function(counts, lib) {
  suppressMessages(edgeR::filterByExpr(
    edgeR::DGEList(counts, lib.size = lib),
    group = NULL
  ))
}

clear_peak <- function() {
  try(writeLines("5", "/proc/self/clear_refs"), silent = TRUE)
}

# 3. One split ---------------------------------------------------------------------
# jarl-ignore cyclomatic_complexity: one pass over all methods of a split, kept as run
run_split <- function(sid, m7_copies = 30L, timing_detail = FALSE, CC = C) {
  t_split <- Sys.time()
  row <- splits[splits$split_id == sid, ]
  ms <- if (variant == "member") "sensitivity" else "primary"
  sm <- split_samples[
    split_samples$split_id == sid & split_samples$member_set == ms,
  ]
  s <- c(sm$sample_id[sm$arm == "a"], sm$sample_id[sm$arm == "b"])
  in_b <- c(rep(FALSE, sum(sm$arm == "a")), rep(TRUE, sum(sm$arm == "b")))
  seed_split <- seed_pkg + 100L + match(sid, splits$split_id)
  meths <- methods_for(row$design, sid)
  timing <- list()
  tick <- function(m, t0) {
    # jarl-ignore undesirable_operator: appends to a result list of the enclosing scope
    timing[[length(timing) + 1L]] <<- data.frame(
      split_id = sid,
      method = m,
      seconds = as.numeric(difftime(Sys.time(), t0, units = "secs")),
      peak_rss_gb = ue_peak_rss_gb()
    )
    if (timing_detail) clear_peak()
  }
  if (timing_detail) {
    clear_peak()
  }
  t0 <- Sys.time()
  if (variant == "unfolded") {
    cnt_all <- CC$counts_unf[, s, drop = FALSE]
    mm <- CC$umean[, s, drop = FALSE]
    mv <- CC$uvar[, s, drop = FALSE]
    U <- intersect(CC$U, rownames(cnt_all))
  } else {
    cnt_all <- CC$counts[, s, drop = FALSE]
    mm <- if (variant == "draws16") {
      CC$mom_mean16[, s, drop = FALSE]
    } else {
      CC$mom_mean[, s, drop = FALSE]
    }
    mv <- if (variant == "draws16") {
      CC$mom_var16[, s, drop = FALSE]
    } else {
      CC$mom_var[, s, drop = FALSE]
    }
    U <- CC$U
  }
  nd <- if (variant == "draws16") 15L else n_draws
  lib <- colSums(cnt_all)
  od_tab <- overdispersion_from_moments(list(mean = mm, var = mv, n = nd))
  od <- stats::setNames(od_tab$overdispersion, od_tab$transcript_id)
  prior <- attr(od_tab, "prior")
  cntU <- cnt_all[U, , drop = FALSE]
  keep1 <- U[design_free_keep(cntU, lib)]
  scaled_all <- cnt_all / od[rownames(cnt_all)]
  lib2 <- colSums(scaled_all)
  scU <- scaled_all[U, , drop = FALSE]
  keep2 <- U[design_free_keep(scU, lib2)]
  tick("setup", t0)
  res <- list()
  if ("M1" %in% meths) {
    t0 <- Sys.time()
    res$M1 <- edger_ql(cntU[keep1, , drop = FALSE], in_b)
    tick("M1", t0)
  }
  if ("M2" %in% meths) {
    t0 <- Sys.time()
    res$M2 <- edger_ql(scU[keep2, , drop = FALSE], in_b)
    tick("M2", t0)
  }
  if ("M3" %in% meths) {
    t0 <- Sys.time()
    res$M3 <- limma_voom(cntU[keep1, , drop = FALSE], in_b)
    tick("M3", t0)
  }
  if ("M4" %in% meths) {
    t0 <- Sys.time()
    res$M4 <- limma_voom(scU[keep2, , drop = FALSE], in_b)
    tick("M4", t0)
  }
  if ("M5" %in% meths) {
    t0 <- Sys.time()
    res$M5 <- deseq2_wald(cntU[keep1, , drop = FALSE], in_b)
    tick("M5", t0)
  }
  if ("M8" %in% meths) {
    t0 <- Sys.time()
    res$M8 <- wilcox_logcpm(cntU[keep1, , drop = FALSE], in_b)
    tick("M8", t0)
  }
  arr <- NULL
  if (any(c("M6", "M7", "M9") %in% meths)) {
    slice <- if (variant == "draws16") 16:30 else seq_len(n_draws)
    arr <- CC$arr[keep1, s, slice, drop = FALSE]
  }
  if ("M6" %in% meths) {
    t0 <- Sys.time()
    res$M6 <- swish_test(
      cntU[keep1, , drop = FALSE],
      arr,
      CC$efflen[keep1, s, drop = FALSE],
      in_b,
      seed_split
    )
    tick("M6", t0)
  }
  if ("M7" %in% meths) {
    t0 <- Sys.time()
    pt <- array(
      rep(cntU[keep1, , drop = FALSE], m7_copies),
      c(length(keep1), length(s), m7_copies)
    )
    res$M7 <- swish_test(
      cntU[keep1, , drop = FALSE],
      pt,
      CC$efflen[keep1, s, drop = FALSE],
      in_b,
      seed_split
    )
    rm(pt)
    tick("M7", t0)
  }
  flags <- NULL
  if ("M9" %in% meths) {
    t0 <- Sys.time()
    flags <- m9_flags(keep1, s, in_b, lib, arr, mm, mv, CC)
    tick("M9", t0)
  }
  gene <- NULL
  if ("M11" %in% meths) {
    t0 <- Sys.time()
    gc_all <- rowsum(cnt_all, CC$gene_all[rownames(cnt_all)], reorder = TRUE)
    glib <- colSums(gc_all)
    kg <- design_free_keep(gc_all, glib)
    g1 <- edger_ql(gc_all[kg, , drop = FALSE], in_b)
    gmm <- if (variant == "draws16") CC$gmean16 else CC$gmean
    gmv <- if (variant == "draws16") CC$gvar16 else CC$gvar
    god_tab <- overdispersion_from_moments(list(
      mean = gmm[, s, drop = FALSE],
      var = gmv[, s, drop = FALSE],
      n = nd
    ))
    god <- stats::setNames(god_tab$overdispersion, god_tab$transcript_id)
    gfac <- god[rownames(gc_all)]
    gfac[is.na(gfac)] <- attr(god_tab, "prior")
    gsc <- gc_all / gfac
    kg2 <- design_free_keep(gsc, colSums(gsc))
    g2 <- edger_ql(gsc[kg2, , drop = FALSE], in_b)
    gene <- merge(
      data.frame(gene_id = g1$id, M11raw_p = g1$pvalue, M11raw_lfc = g1$log2FC),
      data.frame(
        gene_id = g2$id,
        M11div_p = g2$pvalue,
        M11div_lfc = g2$log2FC,
        gene_overdispersion = gfac[g2$id]
      ),
      by = "gene_id",
      all = TRUE
    )
    gene$split_id <- sid
    tick("M11", t0)
  }
  # Per transcript table over U rows tested by any method
  tested <- sort(unique(c(keep1, keep2)))
  tab <- data.frame(split_id = sid, transcript_id = tested)
  tab$overdispersion <- od[tested]
  lc <- log2(t(t(cntU[tested, , drop = FALSE]) / lib) * 1e6 + 0.5)
  tab$mean_logcpm_split <- rowMeans(lc)
  for (m in names(res)) {
    r <- res[[m]]
    i <- match(tested, r$id)
    tab[[paste0(m, "_p")]] <- r$pvalue[i]
    tab[[paste0(m, "_lfc")]] <- r$log2FC[i]
    if (!is.null(r$swish_qvalue)) {
      tab[[paste0(m, "_swishq")]] <- r$swish_qvalue[i]
    }
  }
  if (!is.null(flags)) {
    i <- match(tested, flags$transcript_id)
    for (f in c(
      "assessed",
      "mean_inf_rv",
      "uncertain",
      "ambiguous",
      "sign_agreement",
      "stable",
      "high"
    )) {
      tab[[paste0("M9_", f)]] <- flags[[f]][i]
    }
  }
  attr(tab, "prior") <- prior
  list(
    tab = tab,
    gene = gene,
    timing = do.call(rbind, timing),
    meta = data.frame(
      split_id = sid,
      design = row$design,
      variant = variant,
      member_set = ms,
      n_keep1 = length(keep1),
      n_keep2 = length(keep2),
      prior_overdispersion = prior,
      methods = paste(meths, collapse = ","),
      seed = seed_split,
      seconds = as.numeric(difftime(Sys.time(), t_split, units = "secs")),
      peak_rss_gb = ue_peak_rss_gb()
    )
  )
}

#' de_confidence rules of the package, matrix form, on the split's samples:
#' uncertain = mean InfRV >= 0.5; ambiguous = in a within-gene pair of U
#' with mean draw correlation <= -0.5 over the split samples (both members
#' assessed); stable = sign agreement of the per-draw log2FC >= 0.95 (counts
#' scaled to the mean library size, pseudocount 5); high = assessed, stable,
#' not uncertain, not ambiguous. Identifiability is off (package default).
m9_flags <- function(ids, s, in_b, lib, arr, mm, mv, CC = C) {
  defaults <- uncertainty_defaults()
  # Assessed and ambiguity over all of U (the object holds U)
  uidx <- match(ids, CC$U)
  mean_count_U <- rowMeans(mm[CC$U, , drop = FALSE])
  assessed_U <- mean_count_U >= defaults$min_count
  pc <- if (variant == "draws16") {
    CC$pc16[, s, drop = FALSE]
  } else {
    CC$pc30[, s, drop = FALSE]
  }
  okp <- assessed_U[CC$pairs$i1] & assessed_U[CC$pairs$i2]
  mc <- rowMeans(pc, na.rm = TRUE)
  amb_pair <- okp & is.finite(mc) & mc <= defaults$cor_cutoff
  amb_U <- logical(length(CC$U))
  amb_U[c(CC$pairs$i1[amb_pair], CC$pairs$i2[amb_pair])] <- TRUE
  rv <- rowMeans(inferential_rv(list(
    mean = mm[ids, , drop = FALSE],
    var = mv[ids, , drop = FALSE]
  )))
  scale <- mean(lib) / lib
  lfc <- matrix(NA_real_, length(ids), dim(arr)[3L])
  for (k in seq_len(dim(arr)[3L])) {
    x <- log2(sweep(arr[,, k], 2L, scale, "*") + defaults$pseudocount)
    lfc[, k] <- rowMeans(x[, in_b, drop = FALSE]) -
      rowMeans(x[, !in_b, drop = FALSE])
  }
  med <- apply(lfc, 1L, stats::median)
  agreement <- rowMeans(sign(lfc) == sign(med))
  assessed <- assessed_U[uidx]
  uncertain <- rv >= defaults$inf_rv_cutoff
  ambiguous <- amb_U[uidx]
  stable <- agreement >= defaults$min_agreement
  data.frame(
    transcript_id = ids,
    assessed = assessed,
    mean_inf_rv = rv,
    uncertain = uncertain,
    ambiguous = ambiguous,
    sign_agreement = agreement,
    stable = stable,
    high = assessed & stable & !uncertain & !ambiguous
  )
}

# 4. Per-split summary --------------------------------------------------------------
bh <- function(p) stats::p.adjust(p, "BH")
method_calls <- function(tab) {
  # p and q of every method on its own tested set; M9/M10 derived from M1
  out <- list()
  for (m in c("M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8")) {
    pc <- paste0(m, "_p")
    if (!pc %in% colnames(tab)) {
      next
    }
    p <- tab[[pc]]
    q <- rep(NA_real_, length(p))
    ok <- !is.na(p)
    q[ok] <- bh(p[ok])
    if (paste0(m, "_swishq") %in% colnames(tab)) {
      q <- tab[[paste0(m, "_swishq")]]
    }
    out[[m]] <- list(p = p, q = q)
  }
  if (all(c("M1_p", "M9_high") %in% colnames(tab))) {
    p1 <- out$M1$p
    q1 <- out$M1$q
    t1 <- !is.na(p1)
    keep_a <- tab$M9_high %in% TRUE
    keep_b <- !(tab$M9_uncertain %in% TRUE | tab$M9_ambiguous %in% TRUE)
    out$M9a <- list(
      p = ifelse(keep_a, p1, 1),
      q = ifelse(keep_a, q1, 1),
      tested = t1
    )
    out$M9b <- list(
      p = ifelse(keep_b, p1, 1),
      q = ifelse(keep_b, q1, 1),
      tested = t1
    )
    # M10: drop as many lowest-expressed tested transcripts as M9 drops
    expr <- tab$mean_logcpm_split
    expr[!t1] <- Inf
    o <- order(expr)
    for (v in c("a", "b")) {
      n_drop <- sum(t1 & !(if (v == "a") keep_a else keep_b))
      drop <- o[seq_len(n_drop)]
      kk <- rep(TRUE, length(p1))
      kk[drop] <- FALSE
      out[[paste0("M10", v)]] <- list(
        p = ifelse(kk, p1, 1),
        q = ifelse(kk, q1, 1),
        tested = t1,
        n_drop = n_drop
      )
    }
  }
  out
}

# jarl-ignore cyclomatic_complexity: one pass over all methods of a split, kept as run
summarise_split <- function(tab, gene, meta) {
  calls <- method_calls(tab)
  idx <- match(tab$transcript_id, C$U)
  base <- c("M1", "M2", "M6", "M7", "M8")
  in_common <- Reduce(
    `&`,
    lapply(intersect(base, names(calls)), function(m) !is.na(calls[[m]]$p))
  )
  rows <- list()
  hist_rows <- list()
  for (m in names(calls)) {
    p <- calls[[m]]$p
    q <- calls[[m]]$q
    tested <- calls[[m]]$tested %||% !is.na(p)
    for (set in c("own", "common")) {
      sel <- if (set == "own") tested else tested & in_common
      qq <- q
      if (set == "common") {
        # BH within the common set (swish: BH of its permutation p-values)
        qq <- rep(NA_real_, length(p))
        pc <- p
        if (m %in% c("M9a", "M9b", "M10a", "M10b")) {
          pc <- calls$M1$p
        }
        qq[sel] <- bh(pc[sel])
        if (m %in% c("M9a", "M9b", "M10a", "M10b")) {
          qq[sel & p == 1 & pc < 1] <- 1
        }
      }
      for (v in names(st_vars)) {
        lev <- st_vars[[v]][idx]
        lev[is.na(lev)] <- "NA"
        for (l in unique(lev[sel])) {
          w <- sel & lev == l
          rows[[length(rows) + 1L]] <- data.frame(
            method = m,
            set = set,
            stratum = v,
            level = l,
            n_tested = sum(w),
            n_q05 = sum(qq[w] < 0.05, na.rm = TRUE),
            n_q10 = sum(qq[w] < 0.10, na.rm = TRUE),
            n_p01 = sum(p[w] < 0.01, na.rm = TRUE),
            n_p001 = sum(p[w] < 0.001, na.rm = TRUE)
          )
        }
        if (set == "own" && v %in% c("all", "infrv_tertile")) {
          for (l in unique(lev[sel])) {
            w <- sel & lev == l
            h <- tabulate(pmin(floor(p[w] * 20) + 1L, 20L), 20L)
            hist_rows[[length(hist_rows) + 1L]] <- data.frame(
              method = m,
              stratum = v,
              level = l,
              bin = seq_len(20L),
              lower = (seq_len(20L) - 1) / 20,
              n = h
            )
          }
        }
      }
    }
  }
  if (!is.null(gene)) {
    for (gm in c("M11raw", "M11div")) {
      pg <- gene[[paste0(gm, "_p")]]
      ok <- !is.na(pg)
      qg <- rep(NA_real_, length(pg))
      qg[ok] <- bh(pg[ok])
      both <- !is.na(gene$M11raw_p) & !is.na(gene$M11div_p)
      for (set in c("own", "common")) {
        sel <- if (set == "own") ok else both
        qq <- qg
        if (set == "common") {
          qq <- rep(NA_real_, length(pg))
          qq[sel] <- bh(pg[sel])
        }
        gl <- gstrata$max_infrv_tertile[match(gene$gene_id, gstrata$gene_id)]
        gl[is.na(gl)] <- "not_in_U_multi"
        for (v in c("all", "gene_max_infrv_tertile")) {
          lev <- if (v == "all") rep("all", length(pg)) else gl
          for (l in unique(lev[sel])) {
            w <- sel & lev == l
            rows[[length(rows) + 1L]] <- data.frame(
              method = gm,
              set = set,
              stratum = v,
              level = l,
              n_tested = sum(w),
              n_q05 = sum(qq[w] < 0.05, na.rm = TRUE),
              n_q10 = sum(qq[w] < 0.10, na.rm = TRUE),
              n_p01 = sum(pg[w] < 0.01, na.rm = TRUE),
              n_p001 = sum(pg[w] < 0.001, na.rm = TRUE)
            )
          }
        }
      }
    }
  }
  out <- do.call(rbind, rows)
  out$split_id <- meta$split_id
  out$design <- meta$design
  out$variant <- meta$variant
  hist <- do.call(rbind, hist_rows)
  hist$split_id <- meta$split_id
  hist$design <- meta$design
  hist$variant <- meta$variant
  list(summary = out, hist = hist)
}

write_split <- function(r) {
  sid <- r$meta$split_id
  vdir <- file.path(p$dte, variant)
  seed <- r$meta$seed
  ue_write(r$tab, file.path(vdir, "tx", paste0(sid, ".parquet")), ue, seed)
  if (!is.null(r$gene)) {
    ue_write(r$gene, file.path(vdir, "gene", paste0(sid, ".parquet")), ue, seed)
  }
  sm <- summarise_split(r$tab, r$gene, r$meta)
  ue_write(
    sm$summary,
    file.path(vdir, "summary", paste0(sid, ".parquet")),
    ue,
    seed
  )
  ue_write(sm$hist, file.path(vdir, "hist", paste0(sid, ".parquet")), ue, seed)
  ue_write(
    r$timing,
    file.path(vdir, "timing", paste0(sid, ".parquet")),
    ue,
    seed
  )
  ue_write(r$meta, file.path(vdir, "meta", paste0(sid, ".parquet")), ue, seed)
  invisible(sid)
}

# 5. Pilot ---------------------------------------------------------------------------
# 23_thin_dte.R sources this file with options(wpa.library_only = TRUE) to
# reuse the cache and run_split(); nothing below runs then.
library_only <- isTRUE(getOption("wpa.library_only"))
if (pilot && !library_only) {
  pick <- c(
    splits$split_id[
      splits$design == "7v7_louo" & splits$split_id %in% deseq_ids
    ][1:3],
    splits$split_id[splits$design == "5v5"][1:3],
    splits$split_id[splits$design == "3v3"][1:3]
  )
  all_t <- list()
  all_m <- list()
  for (sid in pick) {
    wpa_log("pilot split ", sid)
    d <- splits$design[splits$split_id == sid]
    # In the pilot every split runs M5 and swish where allowed
    deseq_ids <- unique(c(deseq_ids, sid))
    r <- run_split(sid, m7_copies = 30L, timing_detail = TRUE)
    write_split(r)
    all_t[[sid]] <- cbind(r$timing, design = d)
    all_m[[sid]] <- r$meta
  }
  ue_write(do.call(rbind, all_t), file.path(p$dte, "pilot_timing.parquet"), ue)
  ue_write(do.call(rbind, all_m), file.path(p$dte, "pilot_meta.parquet"), ue)
  print(do.call(rbind, all_t))
  # (a) M7 with 2 copies of the point counts against 30 copies
  sid <- pick[[1L]]
  r2 <- run_split(sid, m7_copies = 2L)
  t30 <- wpa_read(file.path(p$dte, variant, "tx", paste0(sid, ".parquet")))
  v_m7 <- data.frame(
    check = "M7 swish p-values: 2 copies vs 30 copies of the point counts",
    split_id = sid,
    n = sum(!is.na(t30$M7_p)),
    max_abs_diff = max(
      abs(
        t30$M7_p - r2$tab$M7_p[match(t30$transcript_id, r2$tab$transcript_id)]
      ),
      na.rm = TRUE
    )
  )
  # (b) M9 flags against the package (IsoformicExperiment, add_uncertainty)
  row <- splits[splits$split_id == sid, ]
  sm <- split_samples[
    split_samples$split_id == sid & split_samples$member_set == "primary",
  ]
  s <- c(sm$sample_id[sm$arm == "a"], sm$sample_id[sm$arm == "b"])
  lib <- colSums(C$counts[, s])
  reps <- lapply(seq_len(n_draws), function(k) {
    x <- C$arr[, s, k]
    dimnames(x) <- list(C$U, s)
    x
  })
  annot_dir <- Sys.getenv(
    "ISOFORMIC_UE_ANNOT",
    path.expand(
      "~/.cache/R/isoformic/gencode_v50_chr_patch_hapl_scaff_annotation_gff3"
    )
  )
  iso <- IsoformicExperiment(
    annot_path = annot_dir,
    assays = list(counts = C$counts[C$U, s]),
    col_data = data.frame(
      sample_id = s,
      cond = rep(c("a", "b"), c(sum(sm$arm == "a"), sum(sm$arm == "b"))),
      library_size = lib
    )
  )
  iso@metadata$inf_reps <- reps
  iso <- add_uncertainty(
    iso,
    group = "cond",
    levels = c("a", "b"),
    run_qc = FALSE,
    equivalence = FALSE
  )
  pk <- as.data.frame(de_confidence(iso))
  ours <- t30[!is.na(t30$M9_high), ]
  pk <- pk[match(ours$transcript_id, pk$feature_id), ]
  v_m9 <- data.frame(
    check = c(
      "M9 uncertain vs package de_confidence",
      "M9 ambiguous vs package",
      "M9 sign_agreement max abs diff",
      "M9 class high vs package"
    ),
    split_id = sid,
    n = nrow(ours),
    max_abs_diff = c(
      sum(ours$M9_uncertain != (pk$mean_inf_rv >= 0.5), na.rm = TRUE),
      sum(ours$M9_ambiguous != pk$ambiguous, na.rm = TRUE),
      max(abs(ours$M9_sign_agreement - pk$sign_agreement), na.rm = TRUE),
      sum(ours$M9_high != (pk$confidence %in% "high"), na.rm = TRUE)
    )
  )
  val <- rbind(v_m7, v_m9)
  ue_write(val, file.path(p$dte, "pilot_validation.parquet"), ue)
  print(val)
  ue_write_session(file.path(p$dte), ue)
  quit(save = "no")
}

# 6. Full run ------------------------------------------------------------------------
if (!library_only) {
  designs <- strsplit(
    Sys.getenv("ISOFORMIC_UE_DESIGNS", "7v7_louo,7v7_mc,5v5,3v3,7v7_unres"),
    ",",
    fixed = TRUE
  )[[1L]]
  sel <- splits[splits$design %in% designs, ]
  rng <- Sys.getenv("ISOFORMIC_UE_SPLITS")
  if (nzchar(rng)) {
    if (grepl("^[0-9]+:[0-9]+$", rng)) {
      ab <- as.integer(strsplit(rng, ":", fixed = TRUE)[[1L]])
      sel <- sel[seq.int(ab[[1L]], min(ab[[2L]], nrow(sel))), ]
    } else {
      sel <- sel[sel$split_id %in% strsplit(rng, ",", fixed = TRUE)[[1L]], ]
    }
  }
  # ISOFORMIC_UE_FIRST_LOUO=<n>: the first n 7v7_louo splits in (rank in unit,
  # left-out unit) order, so that a subset is balanced over the units (n = 20:
  # rank 1 of every unit and rank 2 of five; the same splits as M5)
  first_louo <- as.integer(Sys.getenv("ISOFORMIC_UE_FIRST_LOUO", "0"))
  if (first_louo > 0L) {
    l <- splits[splits$design == "7v7_louo", ]
    l <- l[order(l$rank_in_unit, l$left_out), ]
    sel <- sel[sel$split_id %in% l$split_id[seq_len(first_louo)], ]
  }
  if (variant == "member") {
    sel <- sel[sel$has_sensitivity_member, ]
  }
  done <- file.exists(file.path(
    p$dte,
    variant,
    "meta",
    paste0(sel$split_id, ".parquet")
  ))
  todo <- sel$split_id[!done]
  wpa_log(
    "variant ",
    variant,
    ": ",
    length(todo),
    " splits to run (",
    sum(done),
    " done) with ",
    workers,
    " workers"
  )
  res <- parallel::mclapply(
    todo,
    function(sid) {
      out <- try({
        r <- run_split(sid)
        write_split(r)
        r$meta
      })
      if (inherits(out, "try-error")) {
        cat("FAILED ", sid, ": ", out, "\n")
      }
      out
    },
    mc.cores = workers,
    mc.preschedule = FALSE
  )
  failed <- todo[vapply(res, inherits, logical(1L), "try-error")]
  wpa_log(
    "finished; failed: ",
    length(failed),
    if (length(failed)) paste0(" (", paste(failed, collapse = ","), ")")
  )
  ue_write_session(file.path(p$dte, variant), ue)
}
