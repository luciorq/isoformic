# WP-D, step 2c: gene-clustered bootstrap of the transcript-level intervals
#
# Purpose: 50_od_prior.R reports a "transcript bootstrap" of the
# genome-wide prior (transcripts resampled independently). The isoforms of a
# gene share reads, so their overdispersion is correlated, and that interval
# is too narrow (verifier of WP-D: SD about 0.003 on placenta against about
# 0.016 when genes are resampled). This script replaces it with a
# gene-clustered bootstrap: genes are resampled with replacement and every
# transcript of a drawn gene enters with the gene's multiplicity. Same
# statistic as the snapshot code (prior = max(1, median(phi) /
# qf(0.5, median(df), 3)) over transcripts with df > 0). It covers:
#   - placenta, all targets (25 runs);
#   - GSE129540 and GSE75489 genome-wide priors (quants_sub, 6 runs each);
#   - the draw-count ratios of the 3 g100 runs (100 / first 30 draws and
#     100 draws / the independent 30-draw quants_sub runs), with the moments
#     re-read here through the snapshot reader (50_od_prior.R did not keep the
#     per-transcript values of these runs).
# Genes are the gene_name of the GENCODE 50 annotation (from
# od_transcripts_placenta.parquet): haplotype/patch copies of a gene share
# its name and so fall in the same cluster. The transcript bootstrap of
# 50_od_prior.R stays in bootstrap_replicates.parquet and is reported only
# as a lower bound of the uncertainty.
#
# Inputs: $ISOFORMIC_UE_OUT/od_prior/od_transcripts_{placenta,GSE129540,
#   GSE75489}.parquet and priors.parquet (50_od_prior.R);
#   $ISOFORMIC_V50_DIR/human/g100/<run> and human/quants_sub/GSE129540_<run>.
# Env vars: ISOFORMIC_V50_DIR, ISOFORMIC_UE_OUT, ISOFORMIC_UE_RLIB,
#   ISOFORMIC_OD_OUT (as in 50_od_prior.R), ISOFORMIC_OD_NBOOT (resamples,
#   default 2000).
# Seeds: base + 1000 x 5 (prior package) + unit; units 70 (placenta),
#   71 and 72 (GSE129540, GSE75489), 73 (draw count): 20266187, 20266188,
#   20266189, 20266190 with the pre-registered base 20261117.
# Outputs: od_prior/gene_bootstrap.parquet (one row per statistic: point,
#   percentile interval, SD, number of genes) and
#   od_prior/gene_bootstrap_replicates.parquet.
# Memory: one 100-draw Salmon run in the reader is 655k x 100 x 8 B = 0.52
#   GB; the per-transcript vectors are a few MB. Measured: 13.5 min, peak
#   2.4 GB (one thread is enough).
# Run from the repository root, inside a memory cap of about 6 GB (see
# README.md):
#   pixi run --frozen Rscript data-raw/uncertainty_evidence/54_od_prior_gene_bootstrap.R

source("data-raw/uncertainty_evidence/ue_common.R")
v50 <- ue_v50_dir()
ue_use_rlib()
ue <- ue_setup("54_od_prior_gene_bootstrap")
out <- Sys.getenv("ISOFORMIC_OD_OUT", file.path(ue_out_root(), "od_prior"))
n_boot <- as.integer(Sys.getenv("ISOFORMIC_OD_NBOOT", "2000"))
prior_df <- 3
seed0 <- ue$seed_base + 1000L * ue$prereg$seeds$package_index$prior
seed_of <- function(unit) as.integer(seed0 + unit)
strip_bar <- function(x) sub("\\|.*$", "", x)
t0 <- Sys.time()
log_step <- function(...) {
  message(
    "[",
    format(round(difftime(Sys.time(), t0, units = "mins"), 1)),
    ", ",
    round(ue_peak_rss_gb(), 2),
    " GB peak] ",
    ...
  )
}

# Prior with transcript weights (multiplicities of the resampled genes)
prior_w <- function(phi, df, w) {
  has <- df > 0 & w > 0
  med_phi <- stats::median(rep.int(phi[has], w[has]))
  med_df <- stats::median(rep.int(df[has], w[has]))
  max(1, med_phi / stats::qf(0.5, df1 = med_df, df2 = prior_df))
}

# Gene-clustered bootstrap of one or more statistics that share the
# resampled genes; `stat` takes the transcript weights and returns a vector
gene_boot <- function(gene, n, seed, stat) {
  g <- match(gene, unique(gene))
  ng <- max(g)
  set.seed(seed)
  res <- lapply(seq_len(n), function(b) {
    wg <- tabulate(sample.int(ng, ng, replace = TRUE), ng)
    stat(wg[g])
  })
  do.call(rbind, res)
}

priors <- read_parquet(file.path(out, "priors.parquet"))
get_prior <- function(dataset, variant) {
  p <- priors$prior[priors$dataset == dataset & priors$variant == variant]
  stopifnot(length(p) == 1L)
  p
}
pl <- read_parquet(file.path(out, "od_transcripts_placenta.parquet"))
gene_of <- stats::setNames(pl$gene_name, pl$transcript_id)
stopifnot(!anyNA(gene_of))

rows <- list()
reps <- list()
record <- function(dataset, statistic, point, ref, values, seed, ngenes) {
  stopifnot(abs(point - ref) <= 1e-9 * abs(ref))
  # jarl-ignore undesirable_operator: appends to a result list of the enclosing scope
  rows[[length(rows) + 1L]] <<- data.frame(
    dataset = dataset,
    statistic = statistic,
    point = point,
    point_in_priors_parquet = ref,
    ci_lo = stats::quantile(values, 0.025, names = FALSE),
    ci_hi = stats::quantile(values, 0.975, names = FALSE),
    boot_sd = stats::sd(values),
    boot_mean = mean(values),
    n_boot = length(values),
    n_genes = ngenes,
    boot_seed = seed,
    method = "gene-clustered bootstrap (gene_name), percentile interval"
  )
  # jarl-ignore undesirable_operator: appends to a result list of the enclosing scope
  reps[[length(reps) + 1L]] <<- data.frame(
    dataset = dataset,
    statistic = statistic,
    replicate = seq_along(values),
    value = values,
    boot_seed = seed
  )
}

# 1. Placenta and the knockdown series ---------------------------------------
sets <- list(
  list(
    dataset = "placenta",
    file = "od_transcripts_placenta.parquet",
    variant = "a_all_targets",
    unit = 70L
  ),
  list(
    dataset = "GSE129540",
    file = "od_transcripts_GSE129540.parquet",
    variant = "genome_all_targets",
    unit = 71L
  ),
  list(
    dataset = "GSE75489",
    file = "od_transcripts_GSE75489.parquet",
    variant = "genome_all_targets",
    unit = 72L
  )
)
for (st in sets) {
  log_step(st$dataset)
  otx <- read_parquet(file.path(out, st$file))
  gene <- unname(gene_of[otx$transcript_id])
  stopifnot(!anyNA(gene))
  phi <- otx$phi_raw
  df <- otx$df
  point <- prior_w(phi, df, rep(1L, length(phi)))
  sd <- seed_of(st$unit)
  b <- gene_boot(gene, n_boot, sd, function(w) prior_w(phi, df, w))
  record(
    st$dataset,
    paste("prior", st$variant),
    point,
    get_prior(st$dataset, st$variant),
    b[, 1],
    sd,
    length(unique(gene[df > 0]))
  )
}

# 2. Draw count: 100 vs 30 draws ----------------------------------------------
log_step("draw count")
read_moments_light <- function(dirs, n_draws = NULL) {
  first <- readLines(
    gzfile(file.path(dirs[[1]], "aux_info", "bootstrap", "names.tsv.gz")),
    n = 1L,
    warn = FALSE
  )
  first <- strip_bar(strsplit(first, "\t", fixed = TRUE)[[1]][1])
  r <- read_inferential_replicates(
    dirs,
    transcripts = first,
    n_draws = n_draws
  )
  attr(r, "moments", exact = TRUE)
}
phi_of <- function(m) {
  ex <- m$mean > 0
  pear <- m$var * (m$n - 1) / m$mean
  pear[!ex] <- 0
  df <- rowSums(ex) * (m$n - 1)
  list(phi = rowSums(pear) / df, df = df, ids = rownames(m$mean))
}
g100_runs <- c("SRR8870780", "SRR8870781", "SRR8870783")
g100_dirs <- stats::setNames(
  file.path(v50, "human", "g100", g100_runs),
  g100_runs
)
sub_dirs <- stats::setNames(
  file.path(v50, "human", "quants_sub", paste0("GSE129540_", g100_runs)),
  g100_runs
)
dc <- list()
m <- read_moments_light(g100_dirs)
stopifnot(m$n == 100L)
dc$g100_100_draws <- phi_of(m)
m <- read_moments_light(g100_dirs, n_draws = 30L)
stopifnot(m$n == 30L)
dc$g100_first_30_draws <- phi_of(m)
m <- read_moments_light(sub_dirs)
stopifnot(m$n == 30L)
dc$quants_sub_30_draws <- phi_of(m)
rm(m)
gc(verbose = FALSE)
ids <- dc[[1]]$ids
for (p in dc) {
  stopifnot(identical(p$ids, ids))
}
gene <- unname(gene_of[strip_bar(ids)])
stopifnot(!anyNA(gene))
one <- rep(1L, length(ids))
pts <- vapply(dc, function(p) prior_w(p$phi, p$df, one), numeric(1))
for (nm in names(dc)) {
  ref <- get_prior("GSE129540_g100_runs", nm)
  stopifnot(abs(pts[[nm]] - ref) <= 1e-9 * ref)
}
log_step("draw count bootstrap")
sdc <- seed_of(73L)
b <- gene_boot(gene, n_boot, sdc, function(w) {
  vapply(dc, function(p) prior_w(p$phi, p$df, w), numeric(1))
})
ng <- length(unique(gene[dc[[1]]$df > 0]))
record(
  "GSE129540_g100_runs",
  "ratio prior 100 / first 30 draws",
  pts[[1]] / pts[[2]],
  get_prior("GSE129540_g100_runs", "g100_100_draws") /
    get_prior("GSE129540_g100_runs", "g100_first_30_draws"),
  b[, 1] / b[, 2],
  sdc,
  ng
)
record(
  "GSE129540_g100_runs",
  "ratio prior 100 draws / quants_sub 30 draws",
  pts[[1]] / pts[[3]],
  get_prior("GSE129540_g100_runs", "g100_100_draws") /
    get_prior("GSE129540_g100_runs", "quants_sub_30_draws"),
  b[, 1] / b[, 3],
  sdc,
  ng
)

res <- dplyr::bind_rows(rows)
ue_write(res, file.path(out, "gene_bootstrap.parquet"), ue)
ue_write(
  dplyr::bind_rows(reps),
  file.path(out, "gene_bootstrap_replicates.parquet"),
  ue
)
log_step("done")
options(width = 200)
print(res[, c(
  "dataset",
  "statistic",
  "point",
  "ci_lo",
  "ci_hi",
  "boot_sd",
  "n_genes",
  "boot_seed"
)])
