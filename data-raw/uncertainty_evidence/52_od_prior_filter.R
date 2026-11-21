# WP-D, step 2b: why a 39-gene panel prior is larger than the genome-wide
# prior: expression filter or gene selection?
#
# Purpose: 50_od_prior.R found that the prior of the bundled 39-gene
# knockdown panels, and of random 39-gene panels on placenta, is 2-3 times
# the genome-wide prior. A panel keeps only transcripts with mean NumReads
# >= 1 (the rule of data-raw/build_knockdown_example.R), while the
# genome-wide prior uses every target with at least one nonzero draw mean.
# This script computes the prior of the whole genome under the panel's
# expression rule (and stricter count floors), which separates the filter
# effect from the gene-selection and panel-size effects. It also breaks down
# the placenta impact of 51_od_prior_impact.R (genome prior vs the median
# random-panel prior) by the df of each transcript and by whether edgeR's
# filterByExpr keeps it, because transcripts without expression (df = 0)
# get phi = prior by construction.
#
# Inputs: $ISOFORMIC_UE_OUT/od_prior/od_transcripts_<dataset>.parquet
#   (phi_raw, df, double precision; 50_od_prior.R), impact_dte_placenta.parquet
#   and impact_placenta_priors.parquet (51_od_prior_impact.R),
#   moments_placenta.parquet, placenta_samples.parquet, the quant.sf files.
# Env vars: ISOFORMIC_V50_DIR, ISOFORMIC_UE_OUT, ISOFORMIC_UE_RLIB,
#   ISOFORMIC_OD_OUT (as in 50_od_prior.R).
# Outputs: od_prior/filter_priors.parquet, od_prior/impact_phi_by_df.parquet
# Resources (estimate): a few minutes, under 4 GB.
# Run from the repository root, inside a memory cap of about 6 GB (see
# README.md):
#   pixi run --frozen Rscript data-raw/uncertainty_evidence/52_od_prior_filter.R

source("data-raw/uncertainty_evidence/ue_common.R")
v50 <- ue_v50_dir()
ue_use_rlib()
ue <- ue_setup("52_od_prior_filter")
out <- Sys.getenv("ISOFORMIC_OD_OUT", file.path(ue_out_root(), "od_prior"))
prior_df <- 3
strip_bar <- function(x) sub("\\|.*$", "", x)

prior_of <- function(phi, df) {
  has <- df > 0
  c(
    prior = max(
      1,
      stats::median(phi[has]) /
        stats::qf(0.5, df1 = stats::median(df[has]), df2 = prior_df)
    ),
    median_phi = stats::median(phi[has]),
    median_df = stats::median(df[has]),
    n_transcripts = sum(has)
  )
}
moderate_phi <- function(phi, df, prior) {
  has <- df > 0
  o <- rep(prior, length(phi))
  o[has] <- pmax(
    (prior_df * prior + df[has] * phi[has]) / (prior_df + df[has]),
    1
  )
  o
}
read_counts <- function(dirs, ids) {
  m <- matrix(0, length(ids), length(dirs), dimnames = list(ids, names(dirs)))
  for (i in seq_along(dirs)) {
    q <- data.table::fread(
      file.path(dirs[[i]], "quant.sf"),
      select = c("Name", "NumReads")
    )
    m[, i] <- q$NumReads[match(ids, strip_bar(q$Name))]
  }
  m
}

# 1. Genome-wide prior under expression floors ------------------------------
ss <- read_parquet(file.path(out, "placenta_samples.parquet"))
kd_runs <- list(
  GSE129540 = c(
    "SRR8870780",
    "SRR8870781",
    "SRR8870782",
    "SRR8870783",
    "SRR8870784",
    "SRR8870785"
  ),
  GSE75489 = c(
    "SRR2966415",
    "SRR2966416",
    "SRR2966417",
    "SRR2966412",
    "SRR2966413",
    "SRR2966414"
  )
)
dsets <- list(
  placenta = stats::setNames(
    file.path(v50, "human", "placenta", "quants", ss$sample_id),
    ss$sample_id
  ),
  GSE129540 = stats::setNames(
    file.path(
      v50,
      "human",
      "quants_sub",
      paste0("GSE129540_", kd_runs$GSE129540)
    ),
    kd_runs$GSE129540
  ),
  GSE75489 = stats::setNames(
    file.path(
      v50,
      "human",
      "quants_sub",
      paste0("GSE75489_", kd_runs$GSE75489)
    ),
    kd_runs$GSE75489
  )
)
rows <- list()
for (d in names(dsets)) {
  otx <- read_parquet(file.path(out, paste0("od_transcripts_", d, ".parquet")))
  cnt <- read_counts(dsets[[d]], otx$transcript_id)
  mc <- rowMeans(cnt)
  for (floor in c(0, 1, 10, 100)) {
    keep <- if (floor == 0) rep(TRUE, length(mc)) else mc >= floor
    st <- prior_of(otx$phi_raw[keep], otx$df[keep])
    rows[[length(rows) + 1L]] <- data.frame(
      dataset = d,
      rule = if (floor == 0) {
        "all targets (genome-wide prior)"
      } else {
        sprintf("mean NumReads >= %g", floor)
      },
      count_floor = floor,
      prior = st[["prior"]],
      median_phi = st[["median_phi"]],
      median_df = st[["median_df"]],
      n_transcripts = st[["n_transcripts"]]
    )
  }
  # Expressed in every sample (full df), the transcripts that dominate a
  # panel of well expressed genes
  full <- otx$df == max(otx$df)
  st <- prior_of(otx$phi_raw[full], otx$df[full])
  rows[[length(rows) + 1L]] <- data.frame(
    dataset = d,
    rule = "nonzero draw mean in every sample (full df)",
    count_floor = NA_real_,
    prior = st[["prior"]],
    median_phi = st[["median_phi"]],
    median_df = st[["median_df"]],
    n_transcripts = st[["n_transcripts"]]
  )
  rm(otx, cnt)
}
fp <- dplyr::bind_rows(rows)
ue_write(fp, file.path(out, "filter_priors.parquet"), ue)
print(fp)

# 2. Placenta impact by df class (23 samples of 51_od_prior_impact.R) -------
use <- ss$sample_id[ss$condition == "preeclampsia" | ss$control_unit]
otx <- read_parquet(file.path(out, "od_transcripts_placenta.parquet"))
ids <- otx$transcript_id
m <- read_parquet(file.path(out, "moments_placenta.parquet"))
m <- m[m$sample_id %in% use, ]
num <- tapply(
  m$var * (unique(m$n_draws) - 1) / m$mean,
  factor(m$transcript_id, ids),
  sum
)
nexp <- tapply(m$mean > 0, factor(m$transcript_id, ids), sum)
rm(m)
num[is.na(num)] <- 0
nexp[is.na(nexp)] <- 0
df <- as.numeric(nexp) * 29
phi <- as.numeric(num) / df
pri <- read_parquet(file.path(out, "impact_placenta_priors.parquet"))
p_gen <- pri$prior[pri$label == "genome"]
stopifnot(abs(prior_of(phi, df)[["prior"]] - p_gen) < 1e-4 * p_gen)
dte <- read_parquet(file.path(out, "impact_dte_placenta.parquet"))
tested <- unique(dte$transcript_id[dte$analysis == "genome_prior"])
phi_g <- moderate_phi(phi, df, p_gen)
df_class <- cut(
  df,
  c(-Inf, 0, 29 * 5, 29 * 15, 29 * 22, Inf),
  labels = c(
    "df 0 (no nonzero draw mean)",
    "1-5 samples",
    "6-15 samples",
    "16-22 samples",
    "all 23 samples"
  )
)
by <- list()
for (lab in pri$label[pri$label != "genome"]) {
  p_alt <- pri$prior[pri$label == lab]
  r <- moderate_phi(phi, df, p_alt) / phi_g
  for (cl in c(levels(df_class), "edgeR filterByExpr kept (tested)")) {
    sel <- if (cl == "edgeR filterByExpr kept (tested)") {
      ids %in% tested
    } else {
      df_class == cl
    }
    by[[length(by) + 1L]] <- data.frame(
      alternative = lab,
      prior_genome = p_gen,
      prior_alternative = p_alt,
      subset = cl,
      n_transcripts = sum(sel),
      n_changed_gt_10pct = sum(abs(r[sel] - 1) > 0.10),
      n_changed_gt_1pct = sum(abs(r[sel] - 1) > 0.01),
      median_ratio = stats::median(r[sel]),
      max_ratio = max(r[sel])
    )
  }
}
byd <- dplyr::bind_rows(by)
ue_write(byd, file.path(out, "impact_phi_by_df.parquet"), ue)
options(width = 200)
print(byd)
