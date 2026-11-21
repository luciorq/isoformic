# WP-D, step 1: genome-wide estimate of the RTA overdispersion prior
#
# Purpose: estimate the prior of the read-to-transcript ambiguity
# overdispersion (Baldoni et al. 2024; edgeR::catchSalmon) from every Salmon
# target, and compare it with the prior that the bundled knockdown example
# gets from its 39-gene panel. The code under test is
# overdispersion_from_moments() of the frozen package snapshot
# (R/uncertainty.R): prior = max(1, median(phi) / qf(0.5, median(df), 3)),
# phi moderated with prior_df = 3.
#
# What it computes (all tables in $ISOFORMIC_UE_OUT/od_prior/):
#   1. Placenta (25 runs, 30 Gibbs draws): prior of (a) all targets,
#      (b) primary-chromosome targets, (c) haplotype/patch copies folded per
#      draw onto their primary transcript with the pre-registered fold map of
#      the S0 store (prereg amendment A1; a 31-mer fold map built here is a
#      sensitivity), (d) preeclampsia only, controls
#      only, the 15 control units, (e) null splits of the control units
#      (3v3, 5v5, 7v7; 100 seeded splits per size, plus the 15 distinct
#      14-unit subsets that every 7v7 split uses). 95% CIs from a sample
#      bootstrap (the 2 technical A/B pairs resampled as blocks; 2,000
#      resamples) and a transcript bootstrap (2,000 resamples).
#   2. Knockdown series GSE129540 and GSE75489: the prior of the bundled
#      panel draws (inst/extdata/knockdown of the snapshot, through the
#      snapshot functions, exactly as the uncertainty article builds the
#      object) and the genome-wide prior of the same 6 runs per series
#      (human/quants_sub, 12 M read pairs, 30 draws), with the panel/genome
#      ratio and a paired sample bootstrap.
#   3. Random 39-gene panels (1,000, seeded) on placenta and on each
#      knockdown series, and panels matched to the composition of the
#      knockdown panel (19 genes matched on isoform count to the 19 known
#      ESRP targets + the 20 placenta genes with most significant
#      transcripts) on placenta.
#   4. Draw count: 100 vs the first 30 draws of the g100 runs (SRR8870780,
#      SRR8870781, SRR8870783), and vs their independent 30-draw runs of
#      quants_sub; placenta first 15 vs all 30 draws.
#   5. The prior df and scale estimated by limma::fitFDistRobustly on the raw
#      overdispersion with its df (placenta), with a sample-bootstrap CI.
#   6. Mouse BNST (19 runs, GENCODE vM39) genome-wide prior.
#   7. Cross-check against edgeR::catchSalmon on the same directories.
#
# Inputs:
#   $ISOFORMIC_V50_DIR/human/placenta/quants/<run>     25 Salmon directories
#   $ISOFORMIC_V50_DIR/human/quants_sub/<series>_<run> 12 Salmon directories
#   $ISOFORMIC_V50_DIR/human/g100/<run>                3 Salmon directories
#   $ISOFORMIC_V50_DIR/human/ref/gencode.v50.chr_patch_hapl_scaff.annotation.gtf.gz
#   $ISOFORMIC_V50_DIR/human/ref/transcripts.fa.gz
#   $ISOFORMIC_V50_DIR/human/ref/salmon_index_v50/duplicate_clusters.tsv
#   $ISOFORMIC_V50_DIR/cases/mouse_bnst/quants/<run>   19 Salmon directories
#   $ISOFORMIC_UE_OUT/prereg.yaml (01_preregister.R): snapshot, seeds, samples
#   $ISOFORMIC_UE_OUT/store/fold_map.parquet, control_units.parquet
#     (10_placenta_prepare.R, S0): the fold map of variant (c), the 15 units
#
# Env vars:
#   ISOFORMIC_V50_DIR     work directory with the Salmon runs (required)
#   ISOFORMIC_UE_OUT      output root (required)
#   ISOFORMIC_UE_RLIB     extra R library with edgeR/limma, default
#                         $ISOFORMIC_V50_DIR/rlib (prepended when it exists)
#   ISOFORMIC_UE_THREADS  threads (default 4)
#   ISOFORMIC_OD_NBOOT    bootstrap resamples (default 2000)
#   ISOFORMIC_OD_NPANEL   random panels per dataset (default 1000)
#   ISOFORMIC_OD_NBOOT_DF bootstrap resamples of the fitFDist fit (default 200)
#   ISOFORMIC_OD_OUT      output directory, default $ISOFORMIC_UE_OUT/od_prior
#   ISOFORMIC_OD_MOUSE_DIR mouse quants, default $ISOFORMIC_V50_DIR/cases/mouse_bnst/quants
#   ISOFORMIC_OD_SPLITS   optional Parquet of null splits written by WP-A
#                         (columns split_id, size, arm, sample_id); default:
#                         splits drawn here (seed in the tables)
#
# Resources: measured 96 min and a peak of 4.8 GB (default threads).
# Run from the repository root, inside a memory cap of about 8 GB (see
# README.md):
#   pixi run --frozen Rscript data-raw/uncertainty_evidence/50_od_prior.R
# Then 51_od_prior_impact.R and 59_od_summarise.R.

source("data-raw/uncertainty_evidence/ue_common.R")
v50 <- ue_v50_dir()
ue_use_rlib()
ue <- ue_setup("50_od_prior")
out <- Sys.getenv("ISOFORMIC_OD_OUT", file.path(ue_out_root(), "od_prior"))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
n_boot <- as.integer(Sys.getenv("ISOFORMIC_OD_NBOOT", "2000"))
n_panel <- as.integer(Sys.getenv("ISOFORMIC_OD_NPANEL", "1000"))
n_boot_df <- as.integer(Sys.getenv("ISOFORMIC_OD_NBOOT_DF", "200"))
# Seeds: base + 1000 * package_index ("prior" = 5) + unit index
seed0 <- ue$seed_base + 1000L * ue$prereg$seeds$package_index$prior
seed_of <- function(unit) as.integer(seed0 + unit)
prior_df <- 3
t_start <- Sys.time()
log_step <- function(...) {
  message(
    format(Sys.time(), "%H:%M:%S"),
    " [",
    round(ue_peak_rss_gb(), 2),
    " GB peak] ",
    ...
  )
}

# Helpers -------------------------------------------------------------------

strip_bar <- function(x) sub("\\|.*$", "", x)

# Per sample Pearson sums and expressed indicators from the moments: the
# terms of overdispersion_from_moments(), so that subsets and bootstrap
# weights of samples are matrix products
pearson_terms <- function(mean, var, n) {
  expressed <- mean > 0
  pear <- var * (n - 1) / mean
  pear[!expressed] <- 0
  storage.mode(expressed) <- "double"
  list(P = pear, E = expressed, n = n)
}

# Prior from Pearson terms with sample weights w (multiplicities) and an
# optional row subset; same formula as the snapshot code
prior_from_terms <- function(pt, w = NULL, rows = NULL) {
  P <- pt$P
  E <- pt$E
  if (!is.null(rows)) {
    P <- P[rows, , drop = FALSE]
    E <- E[rows, , drop = FALSE]
  }
  if (is.null(w)) {
    w <- rep(1, ncol(P))
  }
  num <- drop(P %*% w)
  df <- drop(E %*% w) * (pt$n - 1)
  prior_from_phi(num / df, df)
}

prior_from_phi <- function(phi, df) {
  has <- df > 0
  med_phi <- stats::median(phi[has])
  med_df <- stats::median(df[has])
  raw <- med_phi / stats::qf(0.5, df1 = med_df, df2 = prior_df)
  c(
    prior = max(1, raw),
    prior_unfloored = raw,
    median_phi = med_phi,
    median_df = med_df,
    n_transcripts = sum(has)
  )
}

# Moderated phi for a given prior (the snapshot formula)
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

pct_ci <- function(x, probs = c(0.025, 0.975)) {
  unname(stats::quantile(x, probs, na.rm = TRUE, names = FALSE))
}

# Read the moments of every target through the snapshot reader
read_moments <- function(dirs, n_draws = NULL) {
  reps <- read_inferential_replicates(
    dirs,
    transcripts = NULL,
    n_draws = n_draws
  )
  # Only the moments are needed; keep one transcript of the replicates
  mom <- attr(reps, "moments", exact = TRUE)
  rm(reps)
  gc(verbose = FALSE)
  mom
}

# The replicate matrices are not needed: request one transcript only
read_moments_light <- function(dirs, n_draws = NULL) {
  first <- readLines(
    gzfile(file.path(
      dirs[[1]],
      "aux_info",
      "bootstrap",
      "names.tsv.gz"
    )),
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

num_reads <- function(dirs) {
  m <- NULL
  for (i in seq_along(dirs)) {
    q <- data.table::fread(
      file.path(dirs[[i]], "quant.sf"),
      select = c("Name", "NumReads")
    )
    if (is.null(m)) {
      m <- matrix(
        0,
        nrow(q),
        length(dirs),
        dimnames = list(strip_bar(q$Name), names(dirs))
      )
    }
    m[, i] <- q$NumReads[match(rownames(m), strip_bar(q$Name))]
  }
  m
}

meta_num <- function(dir, key) {
  jsonlite::read_json(file.path(dir, "aux_info", "meta_info.json"))[[key]]
}

# Rows of a long table of the nonzero moments (recomputable input)
moments_long <- function(mom, dataset) {
  idx <- which(mom$mean > 0, arr.ind = TRUE)
  data.frame(
    dataset = dataset,
    transcript_id = rownames(mom$mean)[idx[, 1]],
    sample_id = colnames(mom$mean)[idx[, 2]],
    n_draws = mom$n,
    mean = mom$mean[idx],
    var = mom$var[idx]
  )
}

priors <- list()
add_prior <- function(
  dataset,
  variant,
  samples,
  stats,
  official = NA_real_,
  note = ""
) {
  # jarl-ignore undesirable_operator: appends to a result list of the enclosing scope
  priors[[length(priors) + 1L]] <<- data.frame(
    dataset = dataset,
    variant = variant,
    n_samples = length(samples),
    samples = paste(samples, collapse = ";"),
    prior = unname(stats["prior"]),
    prior_unfloored = unname(stats["prior_unfloored"]),
    median_phi = unname(stats["median_phi"]),
    median_df = unname(stats["median_df"]),
    n_transcripts_df = unname(stats["n_transcripts"]),
    prior_snapshot_function = official,
    note = note
  )
}
boots <- list()
add_boot <- function(dataset, variant, kind, values, seed) {
  # jarl-ignore undesirable_operator: appends to a result list of the enclosing scope
  boots[[length(boots) + 1L]] <<- data.frame(
    dataset = dataset,
    variant = variant,
    kind = kind,
    replicate = seq_along(values),
    value = values,
    boot_seed = seed
  )
}

# Sample bootstrap: blocks of samples resampled with replacement; returns
# the weight matrix (samples x resamples)
block_weights <- function(samples, blocks, n, seed) {
  set.seed(seed)
  ub <- unique(blocks)
  w <- matrix(0, length(samples), n)
  for (b in seq_len(n)) {
    pick <- table(factor(sample(ub, length(ub), replace = TRUE), levels = ub))
    w[, b] <- as.numeric(pick[blocks])
  }
  w
}

transcript_boot <- function(phi, df, n, seed) {
  set.seed(seed)
  has <- which(df > 0)
  phi <- phi[has]
  df <- df[has]
  vapply(
    seq_len(n),
    function(b) {
      i <- sample.int(length(phi), replace = TRUE)
      prior_from_phi(phi[i], df[i])[["prior"]]
    },
    numeric(1)
  )
}

# 0. Annotation, sample sheets ------------------------------------------------
log_step("annotation")
ref <- file.path(v50, "human", "ref")
gtf_cmd <- sprintf(
  "zcat %s | awk -F'\\t' '$3==\"transcript\"{print $1\"\\t\"$9}'",
  shQuote(file.path(ref, "gencode.v50.chr_patch_hapl_scaff.annotation.gtf.gz"))
)
gtf <- data.table::fread(
  cmd = gtf_cmd,
  sep = "\t",
  header = FALSE,
  col.names = c("seqid", "attr"),
  quote = ""
)
attr_of <- function(x, key) {
  m <- regmatches(x, regexpr(paste0(key, " \"[^\"]+\""), x))
  out <- rep(NA_character_, length(x))
  has <- grepl(paste0(key, " \""), x)
  out[has] <- sub(paste0(key, " \"([^\"]+)\""), "\\1", m)
  out
}
tx <- data.frame(
  transcript_id = attr_of(gtf$attr, "transcript_id"),
  gene_id = attr_of(gtf$attr, "gene_id"),
  gene_name = attr_of(gtf$attr, "gene_name"),
  transcript_name = attr_of(gtf$attr, "transcript_name"),
  transcript_type = attr_of(gtf$attr, "transcript_type"),
  seqid = gtf$seqid
)
rm(gtf)
primary_contigs <- c(paste0("chr", c(1:22, "X", "Y", "M")))
tx$primary <- tx$seqid %in% primary_contigs

pre <- ue$prereg
pe_samples <- unlist(pre$data$preeclampsia)
ctrl_samples <- unlist(pre$data$control)
placenta_samples <- c(pe_samples, ctrl_samples)
pairs <- lapply(pre$data$technical_pairs, unlist)
pl_dirs <- stats::setNames(
  file.path(v50, "human", "placenta", "quants", placenta_samples),
  placenta_samples
)
stopifnot(dir.exists(pl_dirs))
# Blocks: the members of an A/B pair share a block
pl_blocks <- stats::setNames(placenta_samples, placenta_samples)
for (p in names(pairs)) {
  pl_blocks[pairs[[p]]] <- p
}
# Control units: the shallower member (fewer mapped fragments) of each pair
# is left out (prereg)
mapped <- vapply(pl_dirs, meta_num, numeric(1), key = "num_mapped")
drop_members <- vapply(
  pairs,
  function(p) p[which.min(mapped[p])],
  character(1)
)
ctrl_units <- setdiff(ctrl_samples, drop_members)
stopifnot(length(ctrl_units) == 15L)
# Same units as the S0 store (the deeper member of each pair)
cu_s0 <- read_parquet(file.path(
  ue_out_root(),
  "store",
  "control_units.parquet"
))
stopifnot(setequal(ctrl_units, cu_s0$sample_id[cu_s0$unit_primary]))
sample_sheet <- data.frame(
  sample_id = placenta_samples,
  condition = ifelse(
    placenta_samples %in% pe_samples,
    "preeclampsia",
    "control"
  ),
  block = unname(pl_blocks[placenta_samples]),
  num_mapped = unname(mapped[placenta_samples]),
  control_unit = placenta_samples %in% ctrl_units,
  num_draws = vapply(pl_dirs, meta_num, numeric(1), key = "num_bootstraps")
)
ue_write(sample_sheet, file.path(out, "placenta_samples.parquet"), ue)

# 1. Placenta moments ---------------------------------------------------------
log_step("placenta moments (snapshot reader)")
mom <- read_moments_light(pl_dirs)
ids <- rownames(mom$mean)
stopifnot(length(ids) == 654828L, mom$n == 30L)
counts_pl <- num_reads(pl_dirs)[ids, ]
ue_write_float(
  moments_long(mom, "placenta"),
  file.path(out, "moments_placenta.parquet"),
  float_cols = c("mean", "var"),
  ue = ue
)

# Folding map: haplotype/patch copies onto their primary transcript --------
log_step("folding map")
fold_path <- file.path(out, "fold_map.parquet")
tx_t <- tx[match(ids, tx$transcript_id), ]
stopifnot(!anyNA(tx_t$transcript_id))
dup <- data.table::fread(file.path(
  ref,
  "salmon_index_v50",
  "duplicate_clusters.tsv"
))
dup$retained <- strip_bar(dup$RetainedRef)
dup$duplicate <- strip_bar(dup$DuplicateRef)
fold <- data.frame(
  transcript_id = ids,
  gene_name = tx_t$gene_name,
  primary = tx_t$primary,
  fold_id = ids,
  rule = ifelse(tx_t$primary, "primary", "unmatched")
)
# (i) Salmon kept a non-primary copy of a duplicate cluster that contains a
# primary transcript: the target carries the reads of the primary sequence
dup_primary <- dup[dup$duplicate %in% tx$transcript_id[tx$primary], ]
dup_primary <- dup_primary[order(dup_primary$duplicate), ]
dup_primary <- dup_primary[!duplicated(dup_primary$retained), ]
i_dup <- match(dup_primary$retained, fold$transcript_id)
i_dup_alt <- i_dup[!is.na(i_dup) & !fold$primary[i_dup]]
fold$fold_id[i_dup_alt] <- dup_primary$duplicate[match(
  fold$transcript_id[i_dup_alt],
  dup_primary$retained
)]
fold$rule[i_dup_alt] <- "duplicate_cluster"
# (ii) Other non-primary targets: the primary target of the same gene_name
# with the highest share of 31-mers of the copy (40 sampled positions), at
# least 0.9; ties broken by the closest length
alt <- which(!fold$primary & fold$rule == "unmatched")
fa <- Biostrings::readDNAStringSet(file.path(ref, "transcripts.fa.gz"))
names(fa) <- strip_bar(names(fa))
cand_all <- fold$transcript_id[fold$primary]
cand_gene <- fold$gene_name[fold$primary]
need <- unique(c(
  fold$transcript_id[alt],
  cand_all[cand_gene %in% fold$gene_name[alt]]
))
seqs <- as.character(fa[need])
rm(fa)
gc(verbose = FALSE)
k <- 31L
match_score <- rep(NA_real_, nrow(fold))
by_gene <- split(cand_all, cand_gene)
for (i in alt) {
  cands <- by_gene[[fold$gene_name[i]]]
  if (is.null(cands)) {
    next
  }
  s <- seqs[[fold$transcript_id[i]]]
  len <- nchar(s)
  if (len < k) {
    next
  }
  pos <- unique(round(seq(
    1,
    len - k + 1L,
    length.out = min(40L, len - k + 1L)
  )))
  kmers <- substring(s, pos, pos + k - 1L)
  cs <- seqs[cands]
  score <- vapply(
    cs,
    function(cseq) mean(stringi::stri_detect_fixed(cseq, kmers)),
    numeric(1)
  )
  best <- which(score == max(score))
  if (length(best) > 1L) {
    best <- best[which.min(abs(nchar(cs[best]) - len))]
  }
  match_score[i] <- score[best]
  if (score[best] >= 0.9) {
    fold$fold_id[i] <- cands[best]
    fold$rule[i] <- "kmer_match"
  }
}
fold$kmer_score <- match_score
rm(seqs)
ue_write(fold, fold_path, ue)
log_step(
  "fold rules: ",
  paste(names(table(fold$rule)), table(fold$rule), collapse = ", ")
)

# Custom pass over the draws: folded moments, first-15-draw moments and a
# check that the unfolded moments equal those of the snapshot reader
log_step("placenta folded and 15-draw moments")
fold_group <- factor(fold$fold_id, levels = unique(fold$fold_id))
mom_fold <- list(
  mean = matrix(
    0,
    nlevels(fold_group),
    length(pl_dirs),
    dimnames = list(levels(fold_group), names(pl_dirs))
  ),
  n = 30L
)
mom_fold$var <- mom_fold$mean
mom15 <- list(mean = mom$mean * 0, n = 15L)
mom15$var <- mom15$mean
# The fold map of the shared placenta store (10_placenta_prepare.R: gene_name
# plus exon structure), when it exists, as a second folding
s0_map_path <- file.path(ue_out_root(), "store", "fold_map.parquet")
has_s0_map <- file.exists(s0_map_path)
stopifnot(has_s0_map)
if (has_s0_map) {
  s0_map <- read_parquet(s0_map_path)
  s0_fold <- s0_map$folded_id[match(ids, s0_map$transcript_id)]
  s0_fold[is.na(s0_fold)] <- ids[is.na(s0_fold)]
  fold$fold_id_s0 <- s0_fold
  fold_group_s0 <- factor(s0_fold, levels = unique(s0_fold))
  mom_fold_s0 <- list(
    mean = matrix(
      0,
      nlevels(fold_group_s0),
      length(pl_dirs),
      dimnames = list(levels(fold_group_s0), names(pl_dirs))
    ),
    n = 30L
  )
  mom_fold_s0$var <- mom_fold_s0$mean
}
max_check <- 0
for (j in seq_along(pl_dirs)) {
  x <- read_bootstrap_dir(pl_dirs[[j]])
  rows <- match(ids, strip_bar(rownames(x$counts)))
  cnt <- x$counts[rows, , drop = FALSE]
  rm(x)
  m30 <- replicate_moments(cnt, 30L)
  max_check <- max(
    max_check,
    abs(m30$mean - mom$mean[, j]),
    abs(m30$var - mom$var[, j])
  )
  m15 <- replicate_moments(cnt, 15L)
  mom15$mean[, j] <- m15$mean
  mom15$var[, j] <- m15$var
  folded <- rowsum(cnt, fold_group, reorder = TRUE)
  mf <- replicate_moments(folded, 30L)
  mom_fold$mean[, j] <- mf$mean[match(levels(fold_group), rownames(folded))]
  mom_fold$var[, j] <- mf$var[match(levels(fold_group), rownames(folded))]
  if (has_s0_map) {
    folded <- rowsum(cnt, fold_group_s0, reorder = TRUE)
    mf <- replicate_moments(folded, 30L)
    lv <- match(levels(fold_group_s0), rownames(folded))
    mom_fold_s0$mean[, j] <- mf$mean[lv]
    mom_fold_s0$var[, j] <- mf$var[lv]
  }
  rm(cnt, folded, m30, m15, mf)
  gc(verbose = FALSE)
}
log_step("custom vs snapshot reader max abs difference: ", signif(max_check, 3))

# Variants (a)-(d) --------------------------------------------------------------
log_step("placenta variants")
pt <- pearson_terms(mom$mean, mom$var, mom$n)
# pt_fold: the pre-registered S0 fold map; pt_fold_k: the 31-mer fold map
pt_fold <- pearson_terms(mom_fold_s0$mean, mom_fold_s0$var, mom_fold_s0$n)
pt_fold_k <- pearson_terms(mom_fold$mean, mom_fold$var, mom_fold$n)
pt15 <- pearson_terms(mom15$mean, mom15$var, mom15$n)
official <- function(m, rows = NULL, cols = NULL) {
  rows <- rows %||% seq_len(nrow(m$mean))
  cols <- cols %||% seq_len(ncol(m$mean))
  attr(
    ue_overdispersion(
      m$mean[rows, cols, drop = FALSE],
      m$var[rows, cols, drop = FALSE],
      m$n
    ),
    "prior"
  )
}
w_of <- function(samples, all = placenta_samples) {
  as.numeric(all %in% samples)
}
od_a <- ue_overdispersion(mom$mean, mom$var, mom$n)
st_a <- prior_from_terms(pt)
add_prior(
  "placenta",
  "a_all_targets",
  placenta_samples,
  st_a,
  attr(od_a, "prior")
)
prim_rows <- which(fold$primary)
add_prior(
  "placenta",
  "b_primary_chromosomes",
  placenta_samples,
  prior_from_terms(pt, rows = prim_rows),
  official(mom, rows = prim_rows),
  note = "targets on chr1-22, X, Y, M (non-primary copies dropped)"
)
od_c <- ue_overdispersion(mom_fold_s0$mean, mom_fold_s0$var, mom_fold_s0$n)
add_prior(
  "placenta",
  "c_haplotype_folded",
  placenta_samples,
  prior_from_terms(pt_fold),
  attr(od_c, "prior"),
  note = paste(
    "non-primary copies summed per draw onto their primary target with the",
    "S0 fold map (prereg amendment A1); unmapped and ambiguous copies kept"
  )
)
add_prior(
  "placenta",
  "c_haplotype_folded_kmer",
  placenta_samples,
  prior_from_terms(pt_fold_k),
  official(mom_fold),
  note = "sensitivity: fold map of this script (duplicate clusters, 31-mer match >= 0.9)"
)
fold$agree_with_s0 <- fold$fold_id == fold$fold_id_s0
ue_write(fold, fold_path, ue)
rm(mom_fold_s0)
add_prior(
  "placenta",
  "d_preeclampsia_only",
  pe_samples,
  prior_from_terms(pt, w = w_of(pe_samples)),
  official(mom, cols = match(pe_samples, placenta_samples))
)
add_prior(
  "placenta",
  "d_controls_only",
  ctrl_samples,
  prior_from_terms(pt, w = w_of(ctrl_samples)),
  official(mom, cols = match(ctrl_samples, placenta_samples))
)
add_prior(
  "placenta",
  "d_control_units",
  ctrl_units,
  prior_from_terms(pt, w = w_of(ctrl_units)),
  official(mom, cols = match(ctrl_units, placenta_samples))
)
add_prior(
  "placenta",
  "draws_first_15",
  placenta_samples,
  prior_from_terms(pt15),
  official(mom15),
  note = "first 15 of the 30 Gibbs draws of every run"
)

# Per transcript table (placenta, all targets and folded)
phi_raw_a <- drop(pt$P %*% rep(1, ncol(pt$P))) /
  (drop(pt$E %*% rep(1, ncol(pt$E))) * (pt$n - 1))
ue_write(
  data.frame(
    dataset = "placenta",
    transcript_id = ids,
    phi_raw = phi_raw_a,
    df = od_a$df,
    phi_moderated = od_a$overdispersion,
    fold_id = fold$fold_id,
    primary = fold$primary,
    gene_name = fold$gene_name
  ),
  file.path(out, "od_transcripts_placenta.parquet"),
  ue
)
ue_write(
  data.frame(
    dataset = "placenta_folded",
    transcript_id = od_c$transcript_id,
    df = od_c$df,
    phi_moderated = od_c$overdispersion
  ),
  file.path(out, "od_transcripts_placenta_folded.parquet"),
  ue
)

# Free the moment matrices: the bootstraps work on the Pearson terms only
# (the 8 GB cap killed an earlier run here)
rm(mom, mom15, mom_fold, pt_fold_k)
pt_prim <- list(P = pt$P[prim_rows, ], E = pt$E[prim_rows, ], n = pt$n)
gc(verbose = FALSE)

# Bootstraps (1): sample blocks and transcripts
log_step("placenta bootstraps")
s_boot <- seed_of(1)
W <- block_weights(placenta_samples, unname(pl_blocks), n_boot, s_boot)
b_a <- apply(W, 2, function(w) prior_from_terms(pt, w = w)[["prior"]])
b_c <- apply(W, 2, function(w) prior_from_terms(pt_fold, w = w)[["prior"]])
b_b <- apply(W, 2, function(w) {
  prior_from_terms(pt_prim, w = w)[["prior"]]
})
rm(pt_prim)
b_15 <- apply(W, 2, function(w) prior_from_terms(pt15, w = w)[["prior"]])
add_boot("placenta", "a_all_targets", "sample_block", b_a, s_boot)
add_boot("placenta", "c_haplotype_folded", "sample_block", b_c, s_boot)
add_boot("placenta", "b_primary_chromosomes", "sample_block", b_b, s_boot)
add_boot("placenta", "ratio_15_over_30", "sample_block", b_15 / b_a, s_boot)
s_tb <- seed_of(2)
add_boot(
  "placenta",
  "a_all_targets",
  "transcript",
  transcript_boot(phi_raw_a, od_a$df, n_boot, s_tb),
  s_tb
)
rm(pt_fold, pt15)
gc(verbose = FALSE)

# (e) Null splits of the control units ----------------------------------------
log_step("placenta splits")
split_file <- Sys.getenv("ISOFORMIC_OD_SPLITS", "")
split_rows <- list()
if (nzchar(split_file) && file.exists(split_file)) {
  sp <- read_parquet(split_file)
  for (sid in unique(sp$split_id)) {
    s <- sp[sp$split_id == sid, ]
    split_rows[[length(split_rows) + 1L]] <- data.frame(
      source = "wp_a",
      size = s$size[1],
      split_id = sid,
      units = paste(sort(s$sample_id), collapse = ";")
    )
  }
} else {
  s_split <- seed_of(3)
  set.seed(s_split)
  for (size in unlist(pre$splits$null_sizes)) {
    for (r in seq_len(pre$splits$null_splits_per_size)) {
      pick <- sample(ctrl_units, 2L * size)
      split_rows[[length(split_rows) + 1L]] <- data.frame(
        source = paste0("drawn_seed_", s_split),
        size = size,
        split_id = sprintf("s%dv%d_%03d", size, size, r),
        units = paste(sort(pick), collapse = ";"),
        arm_a = paste(sort(pick[seq_len(size)]), collapse = ";")
      )
    }
  }
}
# Every 7v7 split of 15 units uses 14 of them: 15 distinct subsets
for (u in ctrl_units) {
  split_rows[[length(split_rows) + 1L]] <- data.frame(
    source = "leave_one_unit_out",
    size = 7L,
    split_id = paste0("loo_", u),
    units = paste(sort(setdiff(ctrl_units, u)), collapse = ";"),
    arm_a = NA_character_
  )
}
splits <- dplyr::bind_rows(split_rows)
split_stats <- t(vapply(
  strsplit(splits$units, ";", fixed = TRUE),
  function(s) prior_from_terms(pt, w = w_of(s)),
  numeric(5)
))
splits <- cbind(splits, split_stats)
ue_write(splits, file.path(out, "splits_prior.parquet"), ue, seed = seed_of(3))

# 2. Knockdown series -----------------------------------------------------------
log_step("knockdown series")
kd_series <- list(
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
extdata <- file.path(ue$pkg_path, "inst", "extdata")
kd <- list()
for (s in names(kd_series)) {
  col <- read_parquet(file.path(extdata, "knockdown", s, "col_data.parquet"))
  stopifnot(identical(unname(kd_series[[s]][col$sample_id]), col$run))
  dirs <- stats::setNames(
    file.path(v50, "human", "quants_sub", paste0(s, "_", kd_series[[s]])),
    names(kd_series[[s]])
  )
  m <- read_moments_light(dirs)
  stopifnot(m$n == 30L)
  ue_write_float(
    moments_long(m, s),
    file.path(out, paste0("moments_", s, ".parquet")),
    float_cols = c("mean", "var"),
    ue = ue
  )
  # The panel object exactly as the uncertainty article builds it
  dir <- file.path(extdata, "knockdown", s)
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
  od_panel <- inferential_overdispersion(iso)
  panel_tx <- draws_ids(draws_view(iso))
  kd[[s]] <- list(
    dirs = dirs,
    mom = m,
    pt = pearson_terms(m$mean, m$var, m$n),
    counts = num_reads(dirs)[rownames(m$mean), ],
    panel_tx = panel_tx,
    od_panel = od_panel,
    panel_prior = attr(od_panel, "prior")
  )
  od_g <- ue_overdispersion(m$mean, m$var, m$n)
  kd[[s]]$od_genome <- od_g
  add_prior(
    s,
    "genome_all_targets",
    names(dirs),
    prior_from_terms(kd[[s]]$pt),
    attr(od_g, "prior"),
    note = "quants_sub, 12 M read pairs, 30 draws"
  )
  # The panel prior from the bundled draws (rounded to 2 decimals)
  pm <- draw_moments(draws_view(iso))
  add_prior(
    s,
    "panel_bundled_draws",
    names(dirs),
    prior_from_terms(pearson_terms(pm$mean, pm$var, pm$n)),
    kd[[s]]$panel_prior,
    note = paste0(
      length(panel_tx),
      " panel transcripts, bundled draws through inferential_overdispersion()"
    )
  )
  prow <- match(panel_tx, rownames(m$mean))
  stopifnot(!anyNA(prow))
  add_prior(
    s,
    "panel_unrounded_moments",
    names(dirs),
    prior_from_terms(kd[[s]]$pt, rows = prow),
    official(m, rows = prow),
    note = "same panel transcripts, moments of the Salmon draws"
  )
  # Paired sample bootstrap (6 samples, no blocks) and transcript bootstrap
  sb <- seed_of(10 + match(s, names(kd_series)))
  Wk <- block_weights(names(dirs), names(dirs), n_boot, sb)
  bg <- apply(Wk, 2, function(w) prior_from_terms(kd[[s]]$pt, w = w)[["prior"]])
  bp <- apply(
    Wk,
    2,
    function(w) prior_from_terms(kd[[s]]$pt, w = w, rows = prow)[["prior"]]
  )
  add_boot(s, "genome_all_targets", "sample", bg, sb)
  add_boot(s, "panel_unrounded_moments", "sample", bp, sb)
  add_boot(s, "ratio_panel_over_genome", "sample", bp / bg, sb)
  phi_k <- drop(kd[[s]]$pt$P %*% rep(1, 6)) /
    (drop(kd[[s]]$pt$E %*% rep(1, 6)) * 29)
  phi_k[!is.finite(phi_k)] <- NA
  tb <- seed_of(20 + match(s, names(kd_series)))
  add_boot(
    s,
    "genome_all_targets",
    "transcript",
    transcript_boot(phi_k, od_g$df, n_boot, tb),
    tb
  )
  # Panel bootstrap over the panel genes (resample the 39 gene names)
  pg <- tx$gene_name[match(panel_tx, tx$transcript_id)]
  set.seed(tb + 100L)
  ug <- unique(pg)
  bpg <- vapply(
    seq_len(n_boot),
    function(b) {
      gsel <- sample(ug, length(ug), replace = TRUE)
      rows <- unlist(lapply(gsel, function(g) prow[pg == g]))
      prior_from_phi(phi_k[rows], od_g$df[rows])[["prior"]]
    },
    numeric(1)
  )
  add_boot(s, "panel_unrounded_moments", "panel_gene", bpg, tb + 100L)
  ue_write(
    data.frame(
      dataset = s,
      transcript_id = od_g$transcript_id,
      phi_raw = phi_k,
      df = od_g$df,
      phi_moderated_genome = od_g$overdispersion,
      in_panel = od_g$transcript_id %in% panel_tx,
      phi_moderated_panel = od_panel$overdispersion[match(
        od_g$transcript_id,
        od_panel$transcript_id
      )]
    ),
    file.path(out, paste0("od_transcripts_", s, ".parquet")),
    ue
  )
  rm(iso, pm)
}

# 3. Random and matched 39-gene panels ------------------------------------------
log_step("panels")
gene_of_target <- tx$gene_name[match(ids, tx$transcript_id)]
n_iso_primary <- table(tx$gene_name[tx$primary])
panel_genes_kd <- unique(tx$gene_name[match(
  kd$GSE129540$panel_tx,
  tx$transcript_id
)])
panel_genes_kd <- unique(c(
  panel_genes_kd,
  tx$gene_name[match(kd$GSE75489$panel_tx, tx$transcript_id)]
))
known_targets <- c(
  "ESRP1",
  "ESRP2",
  "CD44",
  "CTNND1",
  "SCRIB",
  "ENAH",
  "ARHGEF11",
  "MAP3K7",
  "FGFR2",
  "SLK",
  "EXOC7",
  "NUMB",
  "ITGA6",
  "FLNB",
  "KIF13A",
  "PLOD2",
  "CLSTN1",
  "TCF7L2",
  "MBNL1"
)
iso_bin <- function(g) {
  n <- as.numeric(n_iso_primary[g])
  n[is.na(n)] <- 1
  cut(
    n,
    c(0, 1, 4, 9, 19, 49, Inf),
    labels = c("1", "2-4", "5-9", "10-19", "20-49", "50+")
  )
}
# Panel prior of a set of gene names on a dataset: transcripts of those
# genes with mean NumReads >= 1 over the samples (the bundled panel rule)
panel_prior_fun <- function(pt_d, counts_d) {
  keep <- rowMeans(counts_d) >= 1
  num <- drop(pt_d$P %*% rep(1, ncol(pt_d$P)))
  df <- drop(pt_d$E %*% rep(1, ncol(pt_d$E))) * (pt_d$n - 1)
  phi <- num / df
  genes <- gene_of_target
  rows_by_gene <- split(which(keep), genes[keep])
  list(
    genes = names(rows_by_gene),
    prior = function(g) {
      rows <- unlist(rows_by_gene[g], use.names = FALSE)
      c(prior_from_phi(phi[rows], df[rows]), n_rows = length(rows))
    }
  )
}
panel_rows <- list()
run_panels <- function(dataset, pp, n, seed, type, sampler) {
  set.seed(seed)
  res <- lapply(seq_len(n), function(b) {
    g <- sampler()
    st <- pp$prior(g)
    data.frame(
      dataset = dataset,
      panel_type = type,
      panel = b,
      prior = st[["prior"]],
      prior_unfloored = st[["prior_unfloored"]],
      median_phi = st[["median_phi"]],
      median_df = st[["median_df"]],
      n_transcripts = st[["n_rows"]],
      genes = paste(g, collapse = ";"),
      panel_seed = seed
    )
  })
  # jarl-ignore undesirable_operator: appends to a result list of the enclosing scope
  panel_rows[[length(panel_rows) + 1L]] <<- dplyr::bind_rows(res)
}
stopifnot(
  identical(rownames(kd$GSE129540$mom$mean), ids),
  identical(rownames(kd$GSE75489$mom$mean), ids)
)
pp_pl <- panel_prior_fun(pt, counts_pl)
run_panels("placenta", pp_pl, n_panel, seed_of(30), "random_39", function() {
  sample(pp_pl$genes, 39L)
})
# The knockdown panel genes on placenta (one panel)
st_kdp <- pp_pl$prior(intersect(panel_genes_kd, pp_pl$genes))
panel_rows[[length(panel_rows) + 1L]] <- data.frame(
  dataset = "placenta",
  panel_type = "knockdown_panel_genes",
  panel = 1L,
  prior = st_kdp[["prior"]],
  prior_unfloored = st_kdp[["prior_unfloored"]],
  median_phi = st_kdp[["median_phi"]],
  median_df = st_kdp[["median_df"]],
  n_transcripts = st_kdp[["n_rows"]],
  genes = paste(intersect(panel_genes_kd, pp_pl$genes), collapse = ";"),
  panel_seed = NA_integer_
)
# Matched panels on placenta: top-20 genes by significant transcripts
# (edgeR QL, 8 PE vs 17 controls, q < 0.05, mean count >= 20) plus 19
# genes drawn in the isoform-count bins of the 19 known ESRP targets
log_step("placenta edgeR for the matched panels")
grp <- factor(
  ifelse(placenta_samples %in% pe_samples, "pe", "control"),
  levels = c("control", "pe")
)
y <- edgeR::DGEList(counts_pl, group = grp)
y <- y[edgeR::filterByExpr(y, group = grp), , keep.lib.sizes = FALSE]
y <- edgeR::normLibSizes(y)
fit <- edgeR::glmQLFit(y, stats::model.matrix(~grp))
tt <- edgeR::topTags(edgeR::glmQLFTest(fit, coef = 2), n = Inf)$table
tt$transcript_id <- rownames(tt)
tt$gene_name <- gene_of_target[match(tt$transcript_id, ids)]
tt$mean_count <- rowMeans(counts_pl[tt$transcript_id, ])
sig <- tt[tt$FDR < 0.05 & tt$mean_count >= 20, ]
rank_genes <- stats::aggregate(
  PValue ~ gene_name,
  data = tt[tt$mean_count >= 20, ],
  FUN = min
)
rank_genes$n_sig <- as.numeric(table(factor(
  sig$gene_name,
  levels = rank_genes$gene_name
)))
rank_genes <- rank_genes[order(-rank_genes$n_sig, rank_genes$PValue), ]
top20 <- head(rank_genes$gene_name, 20L)
ue_write(
  data.frame(
    rank = seq_along(top20),
    gene_name = top20,
    n_sig = head(rank_genes$n_sig, 20L),
    min_p = head(rank_genes$PValue, 20L),
    n_transcripts_tested = nrow(tt),
    n_sig_total = nrow(sig)
  ),
  file.path(out, "placenta_top_de_genes.parquet"),
  ue
)
target_bins <- iso_bin(known_targets)
pool <- setdiff(pp_pl$genes, top20)
pool_bin <- iso_bin(pool)
run_panels("placenta", pp_pl, n_panel, seed_of(31), "matched_39", function() {
  c(
    top20,
    unlist(lapply(levels(target_bins), function(b) {
      nb <- sum(target_bins == b)
      if (nb == 0L) character(0) else sample(pool[pool_bin == b], nb)
    }))
  )
})
for (s in names(kd)) {
  pp_k <- panel_prior_fun(kd[[s]]$pt, kd[[s]]$counts)
  run_panels(
    s,
    pp_k,
    n_panel,
    seed_of(32 + match(s, names(kd))),
    "random_39",
    function() sample(pp_k$genes, 39L)
  )
}
panels <- dplyr::bind_rows(panel_rows)
ue_write(panels, file.path(out, "panels_prior.parquet"), ue)
# Only the genome-wide overdispersion of the knockdown runs is used below
for (s in names(kd)) {
  kd[[s]]$mom <- NULL
  kd[[s]]$pt <- NULL
  kd[[s]]$counts <- NULL
}
rm(counts_pl, y, fit)
gc(verbose = FALSE)

# 4. Draw count -------------------------------------------------------------------
log_step("draw count")
g100_runs <- c("SRR8870780", "SRR8870781", "SRR8870783")
g100_dirs <- stats::setNames(
  file.path(v50, "human", "g100", g100_runs),
  g100_runs
)
sub_dirs <- stats::setNames(
  file.path(v50, "human", "quants_sub", paste0("GSE129540_", g100_runs)),
  g100_runs
)
m100 <- read_moments_light(g100_dirs)
m100_30 <- read_moments_light(g100_dirs, n_draws = 30L)
msub <- read_moments_light(sub_dirs)
stopifnot(m100$n == 100L, m100_30$n == 30L, msub$n == 30L)
dc <- list(
  g100_100_draws = m100,
  g100_first_30_draws = m100_30,
  quants_sub_30_draws = msub
)
dc_pt <- lapply(dc, function(m) pearson_terms(m$mean, m$var, m$n))
for (nm in names(dc)) {
  add_prior(
    "GSE129540_g100_runs",
    nm,
    g100_runs,
    prior_from_terms(dc_pt[[nm]]),
    official(dc[[nm]])
  )
}
phi_of <- function(p) {
  df <- drop(p$E %*% rep(1, ncol(p$E))) * (p$n - 1)
  list(phi = drop(p$P %*% rep(1, ncol(p$P))) / df, df = df)
}
dcp <- lapply(dc_pt, phi_of)
sdc <- seed_of(40)
set.seed(sdc)
nt <- length(dcp[[1]]$phi)
dc_boot <- t(vapply(
  seq_len(n_boot),
  function(b) {
    i <- sample.int(nt, replace = TRUE)
    vapply(
      dcp,
      function(p) prior_from_phi(p$phi[i], p$df[i])[["prior"]],
      numeric(1)
    )
  },
  numeric(3)
))
add_boot(
  "GSE129540_g100_runs",
  "ratio_100_over_first30",
  "transcript",
  dc_boot[, 1] / dc_boot[, 2],
  sdc
)
add_boot(
  "GSE129540_g100_runs",
  "ratio_100_over_sub30",
  "transcript",
  dc_boot[, 1] / dc_boot[, 3],
  sdc
)
both <- dcp$g100_100_draws$df > 0 & dcp$g100_first_30_draws$df > 0
ue_write(
  data.frame(
    comparison = c("100_over_first30", "100_over_sub30"),
    median_transcript_phi_ratio = c(
      stats::median(dcp[[1]]$phi[both] / dcp[[2]]$phi[both]),
      stats::median((dcp[[1]]$phi / dcp[[3]]$phi)[
        dcp[[1]]$df > 0 & dcp[[3]]$df > 0
      ])
    ),
    n_transcripts = c(sum(both), sum(dcp[[1]]$df > 0 & dcp[[3]]$df > 0))
  ),
  file.path(out, "draw_count_transcripts.parquet"),
  ue
)
rm(m100, m100_30, msub, dc, dc_pt, dcp)
gc(verbose = FALSE)

# 5. Prior df: limma::fitFDistRobustly on the raw overdispersion -----------------
log_step("fitFDist")
fdist <- function(phi, df) {
  ok <- df > 0 & is.finite(phi) & phi > 0
  f <- limma::fitFDistRobustly(phi[ok], df1 = df[ok])
  f0 <- limma::fitFDist(phi[ok], df1 = df[ok])
  full_df <- ok & df == max(df[ok])
  f_full <- limma::fitFDistRobustly(phi[full_df], df1 = max(df[ok]))
  c(
    robust_scale = f$scale,
    robust_df2 = f$df2,
    nonrobust_scale = f0$scale,
    nonrobust_df2 = f0$df2,
    full_df_scale = f_full$scale,
    full_df_df2 = f_full$df2,
    n = sum(ok),
    n_full_df = sum(full_df)
  )
}
fd_point <- fdist(phi_raw_a, od_a$df)
sfd <- seed_of(50)
Wd <- block_weights(placenta_samples, unname(pl_blocks), n_boot_df, sfd)
# Serial: forked workers of this 5 GB process exceeded the 8 GB cap
fd_boot <- lapply(seq_len(n_boot_df), function(b) {
  w <- Wd[, b]
  df <- drop(pt$E %*% w) * (pt$n - 1)
  fdist(drop(pt$P %*% w) / df, df)
})
fd_boot <- do.call(rbind, fd_boot)
ue_write(
  data.frame(
    statistic = names(fd_point),
    value = unname(fd_point),
    boot_lo = apply(fd_boot, 2, pct_ci)[1, ],
    boot_hi = apply(fd_boot, 2, pct_ci)[2, ],
    n_boot = n_boot_df,
    note = "fitFDistRobustly(phi_raw, df1 = df); df2 = estimated prior df"
  ),
  file.path(out, "prior_df_fit.parquet"),
  ue,
  seed = sfd
)
ue_write(
  cbind(data.frame(replicate = seq_len(nrow(fd_boot))), as.data.frame(fd_boot)),
  file.path(out, "prior_df_boot.parquet"),
  ue,
  seed = sfd
)

# 6. edgeR::catchSalmon cross-check ---------------------------------------------
log_step("catchSalmon")
cs_rows <- list()
cs_check <- function(dataset, dirs, od) {
  cs <- edgeR::catchSalmon(unname(dirs), verbose = FALSE)
  cs_od <- cs$annotation$Overdispersion[match(
    od$transcript_id,
    strip_bar(rownames(cs$annotation))
  )]
  rel <- abs(od$overdispersion - cs_od) / cs_od
  # jarl-ignore undesirable_operator: appends to a result list of the enclosing scope
  cs_rows[[length(cs_rows) + 1L]] <<- data.frame(
    dataset = dataset,
    n_transcripts = length(rel),
    max_rel_diff_phi = max(rel),
    prior_package = attr(od, "prior"),
    prior_catchsalmon = cs$overdispersion.prior,
    rel_diff_prior = abs(attr(od, "prior") - cs$overdispersion.prior) /
      cs$overdispersion.prior
  )
}
cs_check("placenta", pl_dirs, od_a)
for (s in names(kd)) {
  cs_check(s, kd[[s]]$dirs, kd[[s]]$od_genome)
}
ue_write(
  dplyr::bind_rows(cs_rows),
  file.path(out, "catchsalmon_check.parquet"),
  ue
)

# 7. Mouse BNST -------------------------------------------------------------------
log_step("mouse")
mouse_dir <- Sys.getenv(
  "ISOFORMIC_OD_MOUSE_DIR",
  file.path(v50, "cases", "mouse_bnst", "quants")
)
mouse_sheet <- utils::read.delim("data-raw/case_studies/mouse_bnst_samples.tsv")
m_dirs <- stats::setNames(
  file.path(mouse_dir, mouse_sheet$run),
  mouse_sheet$run
)
m_dirs <- m_dirs[dir.exists(m_dirs)]
if (length(m_dirs) >= 2L) {
  mm <- read_moments_light(m_dirs)
  mpt <- pearson_terms(mm$mean, mm$var, mm$n)
  add_prior(
    "mouse_bnst",
    "genome_all_targets",
    names(m_dirs),
    prior_from_terms(mpt),
    official(mm),
    note = "GENCODE vM39, Salmon 2.11.0, 30 Gibbs draws"
  )
  sm <- seed_of(60)
  Wm <- block_weights(names(m_dirs), names(m_dirs), n_boot, sm)
  add_boot(
    "mouse_bnst",
    "genome_all_targets",
    "sample",
    apply(Wm, 2, function(w) prior_from_terms(mpt, w = w)[["prior"]]),
    sm
  )
  rm(mm, mpt)
}

# Write the summaries -----------------------------------------------------------
pri <- dplyr::bind_rows(priors)
bt <- dplyr::bind_rows(boots)
ue_write(pri, file.path(out, "priors.parquet"), ue, seed = seed0)
ue_write(bt, file.path(out, "bootstrap_replicates.parquet"), ue, seed = seed0)
ue_write(
  data.frame(
    item = c(
      "custom_vs_snapshot_moments_max_abs_diff",
      "runtime_min",
      "peak_rss_gb",
      "n_boot",
      "n_panel",
      "n_boot_df"
    ),
    value = c(
      max_check,
      as.numeric(difftime(Sys.time(), t_start, units = "mins")),
      ue_peak_rss_gb(),
      n_boot,
      n_panel,
      n_boot_df
    )
  ),
  file.path(out, "run_info_50.parquet"),
  ue
)
ue_write_session(out, ue)
print(pri[, c(
  "dataset",
  "variant",
  "n_samples",
  "prior",
  "prior_unfloored",
  "median_phi",
  "median_df",
  "prior_snapshot_function"
)])
log_step(
  "done in ",
  round(difftime(Sys.time(), t_start, units = "mins"), 1),
  " min"
)
