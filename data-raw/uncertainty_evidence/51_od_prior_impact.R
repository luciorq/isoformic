# WP-D, step 2: impact of replacing the panel overdispersion prior with the
# genome-wide prior
#
# Purpose: the bundled knockdown draws are a 39-gene panel, so the package
# estimates their overdispersion prior from the panel alone. This script
# replaces it with the genome-wide prior of the same series (from
# 50_od_prior.R) and measures what changes:
#   - per transcript moderated overdispersion (median and max ratio, number
#     of transcripts changed by more than 10%), two ways: (i) prior swapped,
#     same raw overdispersion of the bundled draws; (ii) the genome-wide
#     moments of the Salmon runs attached to the replicates (the `moments`
#     attribute), which is how the package itself would use them;
#   - the counts of scale_counts();
#   - DTE calls (edgeR QL on divided counts, BH q < 0.05) and
#     test_usage(scale = TRUE) calls (diffSplice, transcript q < 0.05) on the
#     bundled panel: overlap, lost and gained.
# On placenta (8 preeclampsia vs the 15 control units), DTE on divided
# counts with the genome-wide prior against the prior of a random 39-gene
# panel (the panels of 50_od_prior.R recomputed on these 23 samples: the
# median panel and the 2.5% and 97.5% panels), and raw counts as reference.
# Tested sets: every DTE comparison is reported twice. "separate": each
# count matrix gets its own filterByExpr and its own BH correction (what a
# user running each method alone would get; lost and gained calls then mix
# real changes with transcripts that one run did not test, counted in
# n_lost_untested / n_gained_untested). "common": both methods are refitted
# on the intersection of the two filterByExpr sets, so they test the same
# transcripts with the same BH denominator; this is the comparison of the
# methods themselves.
# Nothing under R/ or inst/ is changed; the proposal is in 59_od_summarise.R.
#
# Inputs: $ISOFORMIC_UE_OUT/od_prior/ from 50_od_prior.R (priors.parquet,
#   panels_prior.parquet, moments_<dataset>.parquet, placenta_samples.parquet,
#   od_transcripts_placenta.parquet); the knockdown example of the snapshot;
#   the placenta quant.sf files.
# Env vars: ISOFORMIC_V50_DIR, ISOFORMIC_UE_OUT, ISOFORMIC_UE_RLIB,
#   ISOFORMIC_UE_THREADS, ISOFORMIC_OD_OUT (as in 50_od_prior.R).
# Resources (estimate): a few minutes, under 4 GB.
# Run from the repository root, inside a memory cap of about 6 GB (see
# README.md):
#   pixi run --frozen Rscript data-raw/uncertainty_evidence/51_od_prior_impact.R

source("data-raw/uncertainty_evidence/ue_common.R")
v50 <- ue_v50_dir()
ue_use_rlib()
ue <- ue_setup("51_od_prior_impact")
out <- Sys.getenv("ISOFORMIC_OD_OUT", file.path(ue_out_root(), "od_prior"))
prior_df <- 3
q_cut <- 0.05

priors <- read_parquet(file.path(out, "priors.parquet"))
get_prior <- function(dataset, variant) {
  p <- priors$prior[priors$dataset == dataset & priors$variant == variant]
  stopifnot(length(p) == 1L)
  p
}
strip_bar <- function(x) sub("\\|.*$", "", x)

moderate_phi <- function(phi, df, prior) {
  has <- df > 0
  out <- rep(prior, length(phi))
  out[has] <- pmax(
    (prior_df * prior + df[has] * phi[has]) /
      (prior_df + df[has]),
    1
  )
  out
}

# Moments of a dataset from the long Parquet table (zeros where absent)
moments_from_long <- function(path, ids, samples) {
  m <- read_parquet(path)
  m <- m[m$sample_id %in% samples, ]
  mean <- matrix(0, length(ids), length(samples), dimnames = list(ids, samples))
  var <- mean
  idx <- cbind(match(m$transcript_id, ids), match(m$sample_id, samples))
  keep <- !is.na(idx[, 1])
  mean[idx[keep, , drop = FALSE]] <- m$mean[keep]
  var[idx[keep, , drop = FALSE]] <- m$var[keep]
  list(mean = mean, var = var, n = unique(m$n_draws))
}

# `rows`: transcripts to test; NULL applies filterByExpr to `counts`
edger_dte <- function(
  counts,
  groups,
  lib_size = NULL,
  norm = TRUE,
  rows = NULL
) {
  y <- edgeR::DGEList(counts, group = groups, lib.size = lib_size)
  keep <- if (is.null(rows)) {
    edgeR::filterByExpr(y, group = groups, lib.size = lib_size)
  } else {
    rownames(counts) %in% rows
  }
  y <- y[keep, , keep.lib.sizes = TRUE]
  if (norm) {
    y <- edgeR::normLibSizes(y)
  }
  fit <- edgeR::glmQLFit(y, stats::model.matrix(~groups))
  tt <- edgeR::topTags(edgeR::glmQLFTest(fit, coef = 2), n = Inf)$table
  data.frame(
    transcript_id = rownames(tt),
    logFC = tt$logFC,
    pvalue = tt$PValue,
    qvalue = tt$FDR
  )
}

# ref_tested / alt_tested: the ids each run tested; lost (gained) calls that
# the other run did not test are counted apart
overlap_row <- function(
  context,
  comparison,
  ref_ids,
  alt_ids,
  ref_tested,
  alt_tested,
  tested_set = "separate"
) {
  lost <- setdiff(ref_ids, alt_ids)
  gained <- setdiff(alt_ids, ref_ids)
  data.frame(
    context = context,
    comparison = comparison,
    tested_set = tested_set,
    n_tested_reference = length(ref_tested),
    n_tested_alternative = length(alt_tested),
    n_tested_both = length(intersect(ref_tested, alt_tested)),
    n_calls_reference = length(ref_ids),
    n_calls_alternative = length(alt_ids),
    n_shared = length(intersect(ref_ids, alt_ids)),
    n_lost = length(lost),
    n_gained = length(gained),
    n_lost_untested = sum(!lost %in% alt_tested),
    n_gained_untested = sum(!gained %in% ref_tested),
    lost = paste(sort(lost), collapse = ";"),
    gained = paste(sort(gained), collapse = ";")
  )
}

dte_calls <- function(d) d$transcript_id[d$qvalue < q_cut]

# The same comparison on the intersection of the two filterByExpr sets, one
# BH correction over it; returns the overlap row and both result tables
common_compare <- function(
  context,
  comparison,
  ref_counts,
  alt_counts,
  groups,
  lib_size = NULL,
  norm = TRUE
) {
  kept <- function(x) {
    y <- edgeR::DGEList(x, group = groups, lib.size = lib_size)
    rownames(x)[edgeR::filterByExpr(y, group = groups, lib.size = lib_size)]
  }
  rows <- intersect(kept(ref_counts), kept(alt_counts))
  d_r <- edger_dte(ref_counts, groups, lib_size, norm, rows = rows)
  d_a <- edger_dte(alt_counts, groups, lib_size, norm, rows = rows)
  stopifnot(setequal(d_r$transcript_id, d_a$transcript_id))
  list(
    row = overlap_row(
      context,
      comparison,
      dte_calls(d_r),
      dte_calls(d_a),
      d_r$transcript_id,
      d_a$transcript_id,
      tested_set = "common"
    ),
    ref = d_r,
    alt = d_a
  )
}

phi_change_row <- function(context, ref, alt, prior_ref, prior_alt) {
  r <- alt / ref
  data.frame(
    context = context,
    prior_reference = prior_ref,
    prior_alternative = prior_alt,
    n_transcripts = length(r),
    median_ratio = stats::median(r),
    min_ratio = min(r),
    max_ratio = max(r),
    max_abs_log2_ratio = max(abs(log2(r))),
    n_changed_gt_10pct = sum(abs(r - 1) > 0.10),
    n_changed_gt_1pct = sum(abs(r - 1) > 0.01)
  )
}

extdata <- file.path(ue$pkg_path, "inst", "extdata")
phi_rows <- list()
call_rows <- list()
per_tx <- list()

# 1. Knockdown panels -------------------------------------------------------
for (s in c("GSE129540", "GSE75489")) {
  message(s)
  dir <- file.path(extdata, "knockdown", s)
  col <- read_parquet(file.path(dir, "col_data.parquet"))
  build <- function() {
    iso <- IsoformicExperiment(
      experiment_name = s,
      annot_path = file.path(extdata, "knockdown_annotation"),
      assays = list(
        abundance = read_parquet(file.path(dir, "abundance.parquet")),
        counts = read_parquet(file.path(dir, "counts.parquet"))
      ),
      col_data = col
    )
    iso@metadata$inf_reps <- as_inferential_replicates(
      read_parquet(file.path(dir, "draws.parquet"))
    )
    iso
  }
  iso_p <- build()
  od_p <- inferential_overdispersion(iso_p)
  prior_panel <- attr(od_p, "prior")
  prior_genome <- get_prior(s, "genome_all_targets")

  # (i) prior swapped on the raw overdispersion of the bundled draws
  pm <- draw_moments(draws_view(iso_p))
  ex <- pm$mean > 0
  pear <- ifelse(ex, pm$var * (pm$n - 1) / pm$mean, 0)
  df <- rowSums(ex) * (pm$n - 1)
  phi_raw <- rowSums(pear) / df
  phi_swap <- moderate_phi(phi_raw, df, prior_genome)
  phi_panel <- od_p$overdispersion[match(rownames(pm$mean), od_p$transcript_id)]
  stopifnot(!anyNA(phi_panel))
  phi_rows[[length(phi_rows) + 1L]] <- phi_change_row(
    paste0(s, ": prior swapped, bundled draws"),
    phi_panel,
    phi_swap,
    prior_panel,
    prior_genome
  )

  # (ii) genome-wide moments attached to the replicates
  iso_g <- build()
  reps <- iso_g@metadata$inf_reps
  all_ids <- unique(
    read_parquet(file.path(
      out,
      paste0("od_transcripts_", s, ".parquet")
    ))$transcript_id
  )
  mom <- moments_from_long(
    file.path(out, paste0("moments_", s, ".parquet")),
    all_ids,
    col$sample_id
  )
  attr(reps, "moments") <- mom
  iso_g@metadata$inf_reps <- reps
  od_g <- inferential_overdispersion(iso_g)
  # The long moment tables are single precision: the prior recomputed from
  # them matches the double-precision value of 50_od_prior.R to ~1e-7
  float_rel <- abs(attr(od_g, "prior") - prior_genome) / prior_genome
  message(
    s,
    ": prior from float moments, relative difference ",
    signif(float_rel, 3)
  )
  stopifnot(float_rel < 1e-4)
  phi_g <- od_g$overdispersion[match(od_p$transcript_id, od_g$transcript_id)]
  phi_rows[[length(phi_rows) + 1L]] <- phi_change_row(
    paste0(s, ": genome-wide moments attached"),
    od_p$overdispersion,
    phi_g,
    prior_panel,
    prior_genome
  )
  per_tx[[s]] <- data.frame(
    dataset = s,
    transcript_id = od_p$transcript_id,
    df = od_p$df,
    phi_panel_prior = od_p$overdispersion,
    phi_genome_prior_swapped = phi_swap[match(
      od_p$transcript_id,
      rownames(pm$mean)
    )],
    phi_genome_moments = phi_g
  )

  # scale_counts()
  sc_p <- attr(scale_counts(iso_p), "assays")$counts_scaled
  sc_g <- attr(scale_counts(iso_g), "assays")$counts_scaled
  sc_p <- as.matrix(sc_p)
  sc_g <- as.matrix(sc_g)[rownames(sc_p), colnames(sc_p)]
  nz <- sc_p > 0
  r <- sc_g[nz] / sc_p[nz]
  phi_rows[[length(phi_rows) + 1L]] <- data.frame(
    context = paste0(
      s,
      ": scale_counts() counts_scaled, genome moments / panel"
    ),
    prior_reference = prior_panel,
    prior_alternative = prior_genome,
    n_transcripts = nrow(sc_p),
    median_ratio = stats::median(r),
    min_ratio = min(r),
    max_ratio = max(r),
    max_abs_log2_ratio = max(abs(log2(r))),
    n_changed_gt_10pct = sum(apply(
      abs(sc_g / sc_p - 1) > 0.10,
      1,
      any,
      na.rm = TRUE
    )),
    n_changed_gt_1pct = sum(apply(
      abs(sc_g / sc_p - 1) > 0.01,
      1,
      any,
      na.rm = TRUE
    ))
  )

  # DTE: edgeR QL on divided counts (library sizes of the full run)
  counts <- as.matrix(attr(iso_p, "assays")$counts)
  groups <- factor(col$condition, levels = c("control", "knockdown"))
  lib <- col$library_size
  div_p <- counts /
    od_p$overdispersion[match(rownames(counts), od_p$transcript_id)]
  div_g <- counts /
    od_g$overdispersion[match(rownames(counts), od_g$transcript_id)]
  div_s <- counts / phi_swap[match(rownames(counts), rownames(pm$mean))]
  d_p <- edger_dte(div_p, groups, lib, norm = FALSE)
  d_g <- edger_dte(div_g, groups, lib, norm = FALSE)
  d_s <- edger_dte(div_s, groups, lib, norm = FALSE)
  calls <- dte_calls
  call_rows[[length(call_rows) + 1L]] <- overlap_row(
    paste0(s, ": DTE edgeR QL, divided counts"),
    "panel prior -> genome moments attached",
    calls(d_p),
    calls(d_g),
    d_p$transcript_id,
    d_g$transcript_id
  )
  call_rows[[length(call_rows) + 1L]] <- overlap_row(
    paste0(s, ": DTE edgeR QL, divided counts"),
    "panel prior -> genome prior swapped",
    calls(d_p),
    calls(d_s),
    d_p$transcript_id,
    d_s$transcript_id
  )
  cc_g <- common_compare(
    paste0(s, ": DTE edgeR QL, divided counts"),
    "panel prior -> genome moments attached",
    div_p,
    div_g,
    groups,
    lib,
    norm = FALSE
  )
  cc_s <- common_compare(
    paste0(s, ": DTE edgeR QL, divided counts"),
    "panel prior -> genome prior swapped",
    div_p,
    div_s,
    groups,
    lib,
    norm = FALSE
  )
  call_rows[[length(call_rows) + 1L]] <- cc_g$row
  call_rows[[length(call_rows) + 1L]] <- cc_s$row
  # Usage: test_usage(scale = TRUE)
  u_p <- test_usage(
    iso_p,
    group = "condition",
    levels = c("control", "knockdown"),
    scale = TRUE
  )
  u_g <- test_usage(
    iso_g,
    group = "condition",
    levels = c("control", "knockdown"),
    scale = TRUE
  )
  ucalls <- function(u) u$transcript_id[!is.na(u$qvalue) & u$qvalue < q_cut]
  call_rows[[length(call_rows) + 1L]] <- overlap_row(
    paste0(s, ": test_usage(scale = TRUE), diffSplice"),
    "panel prior -> genome moments attached",
    ucalls(u_p),
    ucalls(u_g),
    u_p$transcript_id,
    u_g$transcript_id
  )
  ugcalls <- function(u) {
    unique(u$gene_id[!is.na(u$gene_qvalue) & u$gene_qvalue < q_cut])
  }
  call_rows[[length(call_rows) + 1L]] <- overlap_row(
    paste0(s, ": test_usage(scale = TRUE), diffSplice, genes"),
    "panel prior -> genome moments attached",
    ugcalls(u_p),
    ugcalls(u_g),
    unique(u_p$gene_id),
    unique(u_g$gene_id)
  )
  ue_write(
    dplyr::bind_rows(
      cbind(dataset = s, analysis = "dte_panel_prior", d_p),
      cbind(dataset = s, analysis = "dte_genome_moments", d_g),
      cbind(dataset = s, analysis = "dte_genome_prior_swapped", d_s),
      cbind(dataset = s, analysis = "common_dte_panel_prior", cc_g$ref),
      cbind(dataset = s, analysis = "common_dte_genome_moments", cc_g$alt),
      cbind(dataset = s, analysis = "common_dte_genome_prior_swapped", cc_s$alt)
    ),
    file.path(out, paste0("impact_dte_", s, ".parquet")),
    ue
  )
  ue_write(
    dplyr::bind_rows(
      cbind(dataset = s, analysis = "usage_panel_prior", as.data.frame(u_p)),
      cbind(dataset = s, analysis = "usage_genome_moments", as.data.frame(u_g))
    ),
    file.path(out, paste0("impact_usage_", s, ".parquet")),
    ue
  )
}

# 2. Placenta: 8 PE vs 15 control units -------------------------------------
message("placenta")
ss <- read_parquet(file.path(out, "placenta_samples.parquet"))
use <- ss$sample_id[ss$condition == "preeclampsia" | ss$control_unit]
grp <- factor(
  ss$condition[match(use, ss$sample_id)],
  levels = c("control", "preeclampsia")
)
otx <- read_parquet(file.path(out, "od_transcripts_placenta.parquet"))
ids <- otx$transcript_id
mom <- moments_from_long(file.path(out, "moments_placenta.parquet"), ids, use)
ex <- mom$mean > 0
pear <- ifelse(ex, mom$var * (mom$n - 1) / mom$mean, 0)
df <- rowSums(ex) * (mom$n - 1)
phi_raw <- rowSums(pear) / df
rm(pear, ex)
prior_of <- function(rows) {
  has <- rows[df[rows] > 0]
  max(
    1,
    stats::median(phi_raw[has]) /
      stats::qf(0.5, df1 = stats::median(df[has]), df2 = prior_df)
  )
}
prior_genome <- prior_of(seq_along(ids))
od_check <- attr(ue_overdispersion(mom$mean, mom$var, mom$n), "prior")
stopifnot(abs(od_check - prior_genome) < 1e-9 * prior_genome)
counts <- matrix(0, length(ids), length(use), dimnames = list(ids, use))
for (r in use) {
  q <- data.table::fread(
    file.path(v50, "human", "placenta", "quants", r, "quant.sf"),
    select = c("Name", "NumReads")
  )
  counts[, r] <- q$NumReads[match(ids, strip_bar(q$Name))]
}
# The random panels of 50_od_prior.R, recomputed on these 23 samples
panels <- read_parquet(file.path(out, "panels_prior.parquet"))
panels <- panels[
  panels$dataset == "placenta" & panels$panel_type == "random_39",
]
keep_tx <- rowMeans(counts) >= 1
rows_by_gene <- split(which(keep_tx), otx$gene_name[keep_tx])
panels$prior_23 <- vapply(
  strsplit(panels$genes, ";", fixed = TRUE),
  function(g) prior_of(unlist(rows_by_gene[g], use.names = FALSE)),
  numeric(1)
)
ord <- order(panels$prior_23)
pick <- function(p) ord[max(1L, ceiling(p * length(ord)))]
chosen <- data.frame(
  label = c("panel_p2.5", "panel_median", "panel_p97.5"),
  panel = panels$panel[c(pick(0.025), pick(0.5), pick(0.975))],
  prior = panels$prior_23[c(pick(0.025), pick(0.5), pick(0.975))]
)
ue_write(
  panels[, c("panel", "genes", "prior", "prior_23")],
  file.path(out, "impact_placenta_panels_23.parquet"),
  ue
)
d_raw <- edger_dte(counts, grp)
phi_gw <- moderate_phi(phi_raw, df, prior_genome)
d_gw <- edger_dte(counts / phi_gw, grp)
calls <- dte_calls
call_rows[[length(call_rows) + 1L]] <- overlap_row(
  "placenta 8 PE vs 15 control units: DTE edgeR QL",
  "raw counts -> divided counts, genome-wide prior",
  calls(d_raw),
  calls(d_gw),
  d_raw$transcript_id,
  d_gw$transcript_id
)
cc_raw <- common_compare(
  "placenta 8 PE vs 15 control units: DTE edgeR QL",
  "raw counts -> divided counts, genome-wide prior",
  counts,
  counts / phi_gw,
  grp
)
call_rows[[length(call_rows) + 1L]] <- cc_raw$row
dte_all <- list(
  cbind(analysis = "raw_counts", prior = NA_real_, d_raw),
  cbind(analysis = "genome_prior", prior = prior_genome, d_gw),
  cbind(
    analysis = "common_raw_vs_genome:raw_counts",
    prior = NA_real_,
    cc_raw$ref
  ),
  cbind(
    analysis = "common_raw_vs_genome:genome_prior",
    prior = prior_genome,
    cc_raw$alt
  )
)
for (i in seq_len(nrow(chosen))) {
  phi_i <- moderate_phi(phi_raw, df, chosen$prior[i])
  phi_rows[[length(phi_rows) + 1L]] <- phi_change_row(
    paste0("placenta 23 samples: genome prior -> ", chosen$label[i]),
    phi_gw,
    phi_i,
    prior_genome,
    chosen$prior[i]
  )
  d_i <- edger_dte(counts / phi_i, grp)
  call_rows[[length(call_rows) + 1L]] <- overlap_row(
    "placenta 8 PE vs 15 control units: DTE edgeR QL, divided counts",
    paste0("genome prior -> ", chosen$label[i]),
    calls(d_gw),
    calls(d_i),
    d_gw$transcript_id,
    d_i$transcript_id
  )
  cc_i <- common_compare(
    "placenta 8 PE vs 15 control units: DTE edgeR QL, divided counts",
    paste0("genome prior -> ", chosen$label[i]),
    counts / phi_gw,
    counts / phi_i,
    grp
  )
  call_rows[[length(call_rows) + 1L]] <- cc_i$row
  dte_all[[length(dte_all) + 1L]] <- cbind(
    analysis = chosen$label[i],
    prior = chosen$prior[i],
    d_i
  )
  dte_all[[length(dte_all) + 1L]] <- cbind(
    analysis = paste0("common_genome_vs_", chosen$label[i], ":genome_prior"),
    prior = prior_genome,
    cc_i$ref
  )
  dte_all[[length(dte_all) + 1L]] <- cbind(
    analysis = paste0(
      "common_genome_vs_",
      chosen$label[i],
      ":",
      chosen$label[i]
    ),
    prior = chosen$prior[i],
    cc_i$alt
  )
}
ue_write(
  dplyr::bind_rows(dte_all),
  file.path(out, "impact_dte_placenta.parquet"),
  ue
)
ue_write(
  rbind(data.frame(label = "genome", panel = NA, prior = prior_genome), chosen),
  file.path(out, "impact_placenta_priors.parquet"),
  ue
)

ue_write(dplyr::bind_rows(phi_rows), file.path(out, "impact_phi.parquet"), ue)
ue_write(
  dplyr::bind_rows(call_rows),
  file.path(out, "impact_calls.parquet"),
  ue
)
ue_write(
  dplyr::bind_rows(per_tx),
  file.path(out, "impact_transcripts_knockdown.parquet"),
  ue
)
print(dplyr::bind_rows(phi_rows)[, 1:10])
print(dplyr::bind_rows(call_rows)[, 1:15])
