# WP-D, step 2c: jackknife intervals of the overdispersion prior
#
# Purpose: the sample bootstrap of 50_od_prior.R is biased for this
# statistic. The prior decreases with the number of distinct samples (null
# splits of placenta control units: about 1.61 with 6 samples, 1.52 with 14,
# 1.48 with 25), and a resample with replacement holds only about 63% distinct
# samples, so the bootstrap distribution sits above the point estimate and
# its percentile interval can exclude it. This script gives delete-one-block
# jackknife standard errors and normal intervals (point +/- 1.96 SE) for the
# sample-level statistics: the placenta prior (23 blocks; the two technical
# A/B pairs are one block each), the folded and primary-chromosome
# variants, the 15/30 draw ratio, the limma::fitFDistRobustly df and scale,
# and for each knockdown series the genome-wide prior and the panel /
# genome ratio (6 runs, delete one run). It also reports the bootstrap
# standard error and bias from bootstrap_replicates.parquet. Jackknife
# intervals of a median-type statistic are approximate: tagged "estimate".
#
# Inputs: the Salmon directories read by 50_od_prior.R (same env vars), the
#   S0 fold map, $ISOFORMIC_UE_OUT/od_prior/{priors,placenta_samples,
#   bootstrap_replicates,od_transcripts_GSE*}.parquet
# Env vars: ISOFORMIC_V50_DIR, ISOFORMIC_UE_OUT, ISOFORMIC_UE_RLIB,
#   ISOFORMIC_OD_OUT.
# Outputs: od_prior/jackknife.parquet (one row per statistic),
#   od_prior/jackknife_replicates.parquet
# Resources (estimate): about 10 min, under 4 GB.
# Run from the repository root, inside a memory cap of about 6 GB (see
# README.md):
#   pixi run --frozen Rscript data-raw/uncertainty_evidence/53_od_prior_jackknife.R

source("data-raw/uncertainty_evidence/ue_common.R")
v50 <- ue_v50_dir()
ue_use_rlib()
ue <- ue_setup("53_od_prior_jackknife")
out <- Sys.getenv("ISOFORMIC_OD_OUT", file.path(ue_out_root(), "od_prior"))
prior_df <- 3
strip_bar <- function(x) sub("\\|.*$", "", x)

read_moments_light <- function(dirs, n_draws = NULL) {
  first <- readLines(
    gzfile(file.path(dirs[[1]], "aux_info", "bootstrap", "names.tsv.gz")),
    n = 1L,
    warn = FALSE
  )
  first <- strip_bar(strsplit(first, "\t", fixed = TRUE)[[1]][1])
  reps <- read_inferential_replicates(
    dirs,
    transcripts = first,
    n_draws = n_draws
  )
  attr(reps, "moments", exact = TRUE)
}
terms <- function(m) {
  ex <- m$mean > 0
  p <- m$var * (m$n - 1) / m$mean
  p[!ex] <- 0
  storage.mode(ex) <- "double"
  list(P = p, E = ex, n = m$n)
}
phi_df <- function(t, w) {
  df <- drop(t$E %*% w) * (t$n - 1)
  list(phi = drop(t$P %*% w) / df, df = df)
}
prior_pd <- function(x, rows = NULL) {
  phi <- x$phi
  df <- x$df
  if (!is.null(rows)) {
    phi <- phi[rows]
    df <- df[rows]
  }
  has <- df > 0
  max(
    1,
    stats::median(phi[has]) /
      stats::qf(0.5, df1 = stats::median(df[has]), df2 = prior_df)
  )
}
fd <- function(x) {
  ok <- x$df > 0 & is.finite(x$phi) & x$phi > 0
  f <- limma::fitFDistRobustly(x$phi[ok], df1 = x$df[ok])
  c(df2 = f$df2, scale = f$scale)
}
jk_rows <- list()
rep_rows <- list()
jk <- function(dataset, statistic, point, reps, blocks, note = "") {
  n <- length(reps)
  se <- sqrt((n - 1) / n * sum((reps - mean(reps))^2))
  # jarl-ignore undesirable_operator: appends to a result list of the enclosing scope
  jk_rows[[length(jk_rows) + 1L]] <<- data.frame(
    dataset = dataset,
    statistic = statistic,
    point = point,
    n_blocks = n,
    jackknife_se = se,
    jackknife_bias = (n - 1) * (mean(reps) - point),
    ci_lo = point - 1.959964 * se,
    ci_hi = point + 1.959964 * se,
    note = note
  )
  # jarl-ignore undesirable_operator: appends to a result list of the enclosing scope
  rep_rows[[length(rep_rows) + 1L]] <<- data.frame(
    dataset = dataset,
    statistic = statistic,
    left_out = blocks,
    value = reps
  )
}

# Placenta -----------------------------------------------------------------
ss <- read_parquet(file.path(out, "placenta_samples.parquet"))
dirs <- stats::setNames(
  file.path(v50, "human", "placenta", "quants", ss$sample_id),
  ss$sample_id
)
m30 <- read_moments_light(dirs)
ids <- rownames(m30$mean)
t30 <- terms(m30)
fm <- read_parquet(file.path(ue_out_root(), "store", "fold_map.parquet"))
fid <- fm$folded_id[match(ids, fm$transcript_id)]
fid[is.na(fid)] <- ids[is.na(fid)]
prim <- fm$primary[match(ids, fm$transcript_id)] %in% TRUE
rm(m30, fm)
gc(verbose = FALSE)
m15 <- read_moments_light(dirs, n_draws = 15L)
stopifnot(m15$n == 15L, identical(rownames(m15$mean), ids))
t15 <- terms(m15)
rm(m15)
gc(verbose = FALSE)
blocks <- unique(ss$block)
w_all <- rep(1, nrow(ss))
pd30 <- phi_df(t30, w_all)
pt_a <- prior_pd(pd30)
pt_b <- prior_pd(pd30, which(prim))
pt_15 <- prior_pd(phi_df(t15, w_all)) / pt_a
pt_fd <- fd(pd30)
pri <- read_parquet(file.path(out, "priors.parquet"))
stopifnot(
  abs(
    pt_a - pri$prior[pri$dataset == "placenta" & pri$variant == "a_all_targets"]
  ) <
    1e-10,
  abs(
    pt_b -
      pri$prior[
        pri$dataset == "placenta" & pri$variant == "b_primary_chromosomes"
      ]
  ) <
    1e-10
)
reps <- t(vapply(
  blocks,
  function(b) {
    w <- as.numeric(ss$block != b)
    x <- phi_df(t30, w)
    a <- prior_pd(x)
    c(
      a = a,
      b = prior_pd(x, which(prim)),
      r15 = prior_pd(phi_df(t15, w)) / a,
      fd(x)
    )
  },
  numeric(5)
))
jk("placenta", "prior a_all_targets", pt_a, reps[, "a"], blocks)
jk("placenta", "prior b_primary_chromosomes", pt_b, reps[, "b"], blocks)
jk("placenta", "ratio prior first 15 / 30 draws", pt_15, reps[, "r15"], blocks)
jk(
  "placenta",
  "fitFDistRobustly df2 (prior df)",
  pt_fd[["df2"]],
  reps[, "df2"],
  blocks
)
jk(
  "placenta",
  "fitFDistRobustly scale",
  pt_fd[["scale"]],
  reps[, "scale"],
  blocks
)
rm(t15)
# Folded variant (S0 map): moments are not additive over transcripts, so the
# draws are folded again (as in 50_od_prior.R) before the jackknife
fgrp <- factor(fid, levels = unique(fid))
mf <- list(mean = matrix(0, nlevels(fgrp), nrow(ss)), n = 30L)
mf$var <- mf$mean
for (j in seq_along(dirs)) {
  x <- read_bootstrap_dir(dirs[[j]])
  cnt <- x$counts[match(ids, strip_bar(rownames(x$counts))), , drop = FALSE]
  rm(x)
  folded <- rowsum(cnt, fgrp, reorder = FALSE)
  mm <- replicate_moments(folded, 30L)
  mf$mean[, j] <- mm$mean
  mf$var[, j] <- mm$var
  rm(cnt, folded, mm)
}
tf <- terms(mf)
rm(mf)
pt_c <- prior_pd(phi_df(tf, w_all))
stopifnot(
  abs(
    pt_c -
      pri$prior[pri$dataset == "placenta" & pri$variant == "c_haplotype_folded"]
  ) <
    1e-10
)
reps_c <- vapply(
  blocks,
  function(b) prior_pd(phi_df(tf, as.numeric(ss$block != b))),
  numeric(1)
)
jk("placenta", "prior c_haplotype_folded (S0 map)", pt_c, reps_c, blocks)
rm(tf, t30)
gc(verbose = FALSE)

# Knockdown series ------------------------------------------------------------
kd_runs <- list(
  GSE129540 = c(
    control_1 = "SRR8870780",
    control_2 = "SRR8870781",
    control_3 = "SRR8870782",
    knockdown_1 = "SRR8870783",
    knockdown_2 = "SRR8870784",
    knockdown_3 = "SRR8870785"
  ),
  GSE75489 = c(
    control_1 = "SRR2966415",
    control_2 = "SRR2966416",
    control_3 = "SRR2966417",
    knockdown_1 = "SRR2966412",
    knockdown_2 = "SRR2966413",
    knockdown_3 = "SRR2966414"
  )
)
for (s in names(kd_runs)) {
  d <- stats::setNames(
    file.path(v50, "human", "quants_sub", paste0(s, "_", kd_runs[[s]])),
    names(kd_runs[[s]])
  )
  m <- read_moments_light(d)
  t <- terms(m)
  otx <- read_parquet(file.path(out, paste0("od_transcripts_", s, ".parquet")))
  prow <- match(otx$transcript_id[otx$in_panel], rownames(m$mean))
  rm(m)
  x <- phi_df(t, rep(1, 6))
  g <- prior_pd(x)
  p <- prior_pd(x, prow)
  stopifnot(
    abs(g - pri$prior[pri$dataset == s & pri$variant == "genome_all_targets"]) <
      1e-10
  )
  rr <- t(vapply(
    seq_len(6),
    function(i) {
      w <- rep(1, 6)
      w[i] <- 0
      y <- phi_df(t, w)
      gg <- prior_pd(y)
      c(g = gg, r = prior_pd(y, prow) / gg)
    },
    numeric(2)
  ))
  jk(s, "genome-wide prior", g, rr[, "g"], names(d))
  jk(
    s,
    "ratio panel (unrounded moments) / genome-wide prior",
    p / g,
    rr[, "r"],
    names(d),
    note = "the bundled-draw panel prior differs from the unrounded one by < 1e-4"
  )
  rm(t, x)
  gc(verbose = FALSE)
}

# Bootstrap bias and SE (from 50_od_prior.R) for comparison -----------------
bt <- read_parquet(file.path(out, "bootstrap_replicates.parquet"))
bpoint <- c(
  "placenta|a_all_targets" = pt_a,
  "placenta|c_haplotype_folded" = pt_c,
  "placenta|b_primary_chromosomes" = pt_b
)
bsum <- do.call(
  rbind,
  lapply(names(bpoint), function(k) {
    dv <- strsplit(k, "|", fixed = TRUE)[[1]]
    v <- bt$value[
      bt$dataset == dv[1] & bt$variant == dv[2] & bt$kind == "sample_block"
    ]
    data.frame(
      dataset = dv[1],
      variant = dv[2],
      point = bpoint[[k]],
      boot_mean = mean(v),
      boot_bias = mean(v) - bpoint[[k]],
      boot_se = stats::sd(v),
      frac_boot_above_point = mean(v > bpoint[[k]])
    )
  })
)
ue_write(bsum, file.path(out, "bootstrap_bias.parquet"), ue)
jkt <- dplyr::bind_rows(jk_rows)
ue_write(jkt, file.path(out, "jackknife.parquet"), ue)
ue_write(
  dplyr::bind_rows(rep_rows),
  file.path(out, "jackknife_replicates.parquet"),
  ue
)
options(width = 200)
print(jkt[, 1:8])
print(bsum)
