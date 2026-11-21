# 40_mouse_prepare.R: inputs of the mouse BNST null (WP-C1).
#
# Purpose: prepare a second, independent null for the uncertainty evidence:
# mouse bed nucleus of the stria terminalis, PRJNA984605, male animals only
# (10 control males are split 5 vs 5 by 41_mouse_null.R; the 9 stressed
# males are never tested and define the universe U_mouse and the strata,
# exactly as the 8 preeclampsia samples do for the placenta null, see
# prereg.yaml `universe` and `strata`).
#
# What it does
#   1. Checks the sample sheet (10 control and 9 stressed males) and the
#      Salmon runs (samp_type, number of draws, version, mapping).
#   2. tx2gene, gene names, biotypes and isoforms per gene from the headers of
#      gencode.vM39.transcripts.fa.gz (no M39 GTF on disk). Patch and
#      haplotype handling is limited to what the headers give: they hold no
#      chromosome, so nothing can be folded; the headers are checked for
#      PAR_Y copies, and the transcripts Salmon dropped as duplicate
#      sequences are listed (the log records both).
#   3. Reads quant.sf and the 30 Gibbs draws of the 19 runs (snapshot reader
#      read_bootstrap_dir()), writes counts, per-sample draw moments, and the
#      draws of the 10 control males for transcripts with >= 10 reads summed
#      over the controls (any transcript that a 5 vs 5 or 4 vs 5 split can
#      test has more, since edgeR::filterByExpr needs a total of 15).
#   4. Universe and strata from the 9 stressed males only: U_mouse = CPM >= 1
#      (library size = sum of NumReads) in >= 3 of the 9; InfRV (fishpond
#      form) averaged over the 9, tertiles on U, uncertain = InfRV >= 0.5,
#      robustness InfRV from draws 1-15; expression deciles of mean CPM on U;
#      isoforms per gene in the index: 1, 2-4, 5-9, 10+; RTA overdispersion of
#      the 9; equivalence_support() and equivalence_identifiability() on the
#      9 (snapshot functions); gene-level max InfRV and any non-identifiable.
#   5. Validation of the matrix-level code: overdispersion of the 10 controls
#      against edgeR::catchSalmon(), InfRV against fishpond::computeInfRV().
#   6. Pairs of sibling isoforms of the 10 controls: per-sample draw
#      correlation (the statistic of the snapshot's inferential_correlation(),
#      validated against it), so that 41 can average any subset of animals.
#
# Inputs
#   MOUSE_QUANTS   default $ISOFORMIC_V50_DIR/cases/mouse_bnst/quants
#   MOUSE_FASTA    default
#                  $ISOFORMIC_V50_DIR/cases/mouse_bnst/ref/gencode.vM39.transcripts.fa.gz
#   (ISOFORMIC_V50_DIR is required unless both are set)
#   MOUSE_VM39_ANNOT  optional, the isoformic-prepared vM39 transcript table
#                  (default in ~/.cache/R/isoformic/), only to log the contigs
#   data-raw/case_studies/mouse_bnst_samples.tsv
#   prereg.yaml and the package snapshot (via ue_common.R)
# Env vars: ISOFORMIC_UE_OUT, ISOFORMIC_UE_THREADS, ISOFORMIC_UE_PKG (see
#   ue_common.R); MOUSE_IDENT=false skips identifiability (slow step: one
#   sample at a time, checkpointed in prepare/identifiability/; rerun with
#   MOUSE_IDENT=true to add it, the other outputs are rewritten identically).
# Outputs: $ISOFORMIC_UE_OUT/null_mouse/prepare/*.parquet and prepare.log
#
# Resources: estimated about 1 h and 8 GB (one thread); the last recorded
# run took 4.5 min with a peak of 4.5 GB (identifiability is checkpointed,
# so a rerun is shorter than a first run).
# Run from the repository root, inside a memory cap of about 8 GB (see
# README.md):
#   pixi run --frozen Rscript data-raw/uncertainty_evidence/40_mouse_prepare.R

source("data-raw/uncertainty_evidence/ue_common.R")
ue <- ue_setup("40_mouse_prepare")
seed <- ue$seed_base + 1000L * ue$prereg$seeds$package_index$C + 0L
set.seed(seed)

quants <- Sys.getenv("MOUSE_QUANTS")
if (!nzchar(quants)) {
  quants <- file.path(ue_v50_dir(), "cases/mouse_bnst/quants")
}
fasta <- Sys.getenv("MOUSE_FASTA")
if (!nzchar(fasta)) {
  fasta <- file.path(
    ue_v50_dir(),
    "cases/mouse_bnst/ref/gencode.vM39.transcripts.fa.gz"
  )
}
do_ident <- !identical(tolower(Sys.getenv("MOUSE_IDENT", "true")), "false")
out <- file.path(ue_out_root(), "null_mouse", "prepare")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(out, "prepare.log")
cat("", file = log_file)
say <- function(...) {
  msg <- paste0(format(Sys.time(), "%H:%M:%S"), " ", sprintf(...))
  cat(msg, "\n", file = log_file, append = TRUE)
  message(msg)
}
say("snapshot %s (sha256 %s), seed %d", ue$git_sha, ue$snapshot_sha256, seed)

# 1. Samples ---------------------------------------------------------------------
sheet <- utils::read.delim("data-raw/case_studies/mouse_bnst_samples.tsv")
males <- sheet[sheet$sex == "male", ]
ctrl <- males$run[males$condition == "control"]
stre <- males$run[males$condition == "stressed"]
stopifnot(length(ctrl) == 10L, length(stre) == 9L)
runs <- c(ctrl, stre)
stopifnot(dir.exists(file.path(quants, runs)))
say("males: %d control, %d stressed (verified)", length(ctrl), length(stre))
meta <- do.call(
  rbind,
  lapply(runs, function(r) {
    m <- jsonlite::read_json(file.path(quants, r, "aux_info", "meta_info.json"))
    data.frame(
      run = r,
      condition = males$condition[match(r, males$run)],
      role = if (r %in% ctrl) "tested" else "strata",
      salmon_version = m$salmon_version,
      samp_type = m$samp_type,
      n_draws = m$num_bootstraps,
      num_processed = m$num_processed,
      num_mapped = m$num_mapped,
      percent_mapped = m$percent_mapped
    )
  })
)
stopifnot(meta$samp_type == "gibbs", meta$n_draws == 30L)
say(
  "all 19 runs: gibbs, 30 draws, salmon %s",
  paste(unique(meta$salmon_version), collapse = ",")
)
ue_write(meta, file.path(out, "samples.parquet"), ue, seed)

# 2. Annotation from the FASTA headers ------------------------------------------------
hdr <- system2(
  "bash",
  c("-c", shQuote(sprintf("zcat %s | grep '^>'", fasta))),
  stdout = TRUE
)
f <- strsplit(sub("^>", "", hdr), "|", fixed = TRUE)
stopifnot(lengths(f) == 8L)
ann <- data.frame(
  transcript_id = vapply(f, `[`, "", 1L),
  gene_id = vapply(f, `[`, "", 2L),
  transcript_name = vapply(f, `[`, "", 5L),
  gene_name = vapply(f, `[`, "", 6L),
  tx_length = as.integer(vapply(f, `[`, "", 7L)),
  transcript_type = vapply(f, `[`, "", 8L)
)
n_par <- sum(grepl("PAR_Y", hdr, fixed = TRUE))
say(
  "FASTA headers: %d transcripts, %d genes, %d PAR_Y entries; headers carry no chromosome, so no patch/haplotype folding is possible",
  nrow(ann),
  length(unique(ann$gene_id)),
  n_par
)
stopifnot(!anyDuplicated(ann$transcript_id))
# Contig check (log only): the isoformic-prepared GENCODE vM39 annotation in
# the user cache, when present, gives the contig of each transcript; GENCODE
# transcripts.fa covers the reference chromosomes only, which this verifies.
vm39 <- Sys.getenv(
  "MOUSE_VM39_ANNOT",
  path.expand(
    "~/.cache/R/isoformic/gencode_vm39_chr_patch_hapl_scaff_annotation_gff3/annot_data_transcripts.parquet"
  )
)
if (file.exists(vm39)) {
  a39 <- as.data.frame(read_parquet(vm39))[, c("transcript_id", "seqid")]
  primary <- paste0("chr", c(1:19, "X", "Y", "M"))
  sq <- a39$seqid[match(ann$transcript_id, a39$transcript_id)]
  ann$seqid <- sq
  other <- sq[!is.na(sq) & !sq %in% primary]
  say(
    "contigs of the FASTA transcripts from %s: %d on chr1-19/X/Y/M, %d on %d other contigs (%s), %d not in that annotation; other contigs that are unplaced or unlocalised scaffolds (GL*, JH*) belong to the reference assembly and are not haplotype copies, so nothing is folded",
    vm39,
    sum(sq %in% primary),
    length(other),
    length(unique(other)),
    paste(sort(unique(sub("[0-9]+\\.[0-9]+$", "*", other))), collapse = ", "),
    sum(is.na(sq))
  )
  rm(a39)
} else {
  ann$seqid <- NA_character_
  say("no prepared vM39 annotation at %s: contigs not checked", vm39)
}

# 3. Counts and draws -----------------------------------------------------------------
qsf <- lapply(runs, function(r) {
  data.table::fread(
    file.path(quants, r, "quant.sf"),
    select = c("Name", "EffectiveLength", "NumReads")
  )
})
ids <- qsf[[1]]$Name
stopifnot(vapply(qsf, function(q) identical(q$Name, ids), TRUE))
dropped <- setdiff(ann$transcript_id, ids)
stopifnot(ids %in% ann$transcript_id)
say(
  "index: %d transcripts; %d FASTA transcripts absent from the index (duplicate sequences removed by Salmon)",
  length(ids),
  length(dropped)
)
counts <- vapply(qsf, function(q) q$NumReads, numeric(length(ids)))
dimnames(counts) <- list(ids, runs)
efflen <- vapply(qsf, function(q) q$EffectiveLength, numeric(length(ids)))
dimnames(efflen) <- list(ids, runs)
rm(qsf)

ann$in_index <- ann$transcript_id %in% ids
ann$n_isoforms_annot <- as.integer(stats::ave(
  seq_len(nrow(ann)),
  ann$gene_id,
  FUN = length
))
ann$n_isoforms_index <- as.integer(stats::ave(
  as.integer(ann$in_index),
  ann$gene_id,
  FUN = sum
))
iso_class <- function(n) {
  as.character(cut(n, c(0, 1, 4, 9, Inf), labels = c("1", "2-4", "5-9", "10+")))
}
ann$isoform_class <- iso_class(ann$n_isoforms_index)
ann$isoform_class[!ann$in_index] <- NA_character_
ue_write(ann, file.path(out, "tx2gene.parquet"), ue, seed)
ue_write(
  data.frame(
    transcript_id = dropped,
    gene_id = ann$gene_id[match(dropped, ann$transcript_id)]
  ),
  file.path(out, "dropped_duplicates.parquet"),
  ue,
  seed
)

n_draws <- 30L
mom_mean <- mom_var <- matrix(
  0,
  length(ids),
  length(runs),
  dimnames = list(ids, runs)
)
ctrl_total <- rowSums(counts[, ctrl])
keep_draws <- ctrl_total >= 10
draws_ctrl <- array(
  0,
  c(sum(keep_draws), length(ctrl), n_draws),
  dimnames = list(ids[keep_draws], ctrl, NULL)
)
draws_strata15 <- matrix(
  0,
  length(ids),
  length(stre),
  dimnames = list(ids, stre)
)
draws_strata15_var <- draws_strata15
t0 <- Sys.time()
for (r in runs) {
  b <- read_bootstrap_dir(file.path(quants, r))
  stopifnot(identical(rownames(b$counts), ids), ncol(b$counts) == n_draws)
  x <- b$counts
  mom_mean[, r] <- rowMeans(x)
  mom_var[, r] <- matrixStats::rowVars(x)
  if (r %in% ctrl) {
    draws_ctrl[, r, ] <- x[keep_draws, ]
  }
  if (r %in% stre) {
    draws_strata15[, r] <- rowMeans(x[, 1:15])
    draws_strata15_var[, r] <- matrixStats::rowVars(x[, 1:15])
  }
  rm(b, x)
}
say(
  "draws read in %.1f min; %d transcripts with >= 10 control reads keep their control draws",
  as.numeric(difftime(Sys.time(), t0, units = "mins")),
  sum(keep_draws)
)

ue_write(
  data.frame(transcript_id = ids, counts, check.names = FALSE),
  file.path(out, "counts.parquet"),
  ue,
  seed
)
ue_write(
  data.frame(transcript_id = ids, efflen, check.names = FALSE),
  file.path(out, "effective_length.parquet"),
  ue,
  seed
)
mom_long <- data.frame(
  transcript_id = rep(ids, length(runs)),
  run = rep(runs, each = length(ids)),
  draw_mean = as.vector(mom_mean),
  draw_var = as.vector(mom_var)
)
ue_write(mom_long, file.path(out, "draw_moments.parquet"), ue, seed)
rm(mom_long)
dl <- do.call(
  rbind,
  lapply(ctrl, function(r) {
    m <- draws_ctrl[, r, ]
    colnames(m) <- sprintf("d%02d", seq_len(n_draws))
    data.frame(
      transcript_id = rownames(draws_ctrl),
      run = r,
      m,
      check.names = FALSE
    )
  })
)
ue_write(dl, file.path(out, "draws_control.parquet"), ue, seed)
rm(dl)
gc()

# 4. Universe and strata from the 9 stressed males ---------------------------------------
lib <- colSums(counts[, stre])
cpm <- sweep(counts[, stre], 2, lib / 1e6, "/")
in_U <- rowSums(cpm >= 1) >= 3L
infrv_s <- ue_infrv(mom_mean[, stre], mom_var[, stre])
infrv <- rowMeans(infrv_s)
infrv15 <- rowMeans(ue_infrv(draws_strata15, draws_strata15_var))
od_s <- ue_overdispersion(mom_mean[, stre], mom_var[, stre], n_draws)
mean_cpm <- rowMeans(cpm)
tert <- function(x, u) {
  q <- stats::quantile(x[u], c(1 / 3, 2 / 3), names = FALSE)
  out <- as.character(cut(x, c(-Inf, q, Inf), labels = c("low", "mid", "high")))
  out[!u] <- NA_character_
  out
}
dec <- function(x, u) {
  q <- stats::quantile(x[u], seq(0.1, 0.9, 0.1), names = FALSE)
  out <- as.integer(cut(x, c(-Inf, q, Inf), labels = FALSE))
  out[!u] <- NA_integer_
  out
}
strata <- data.frame(
  transcript_id = ids,
  gene_id = ann$gene_id[match(ids, ann$transcript_id)],
  in_U = in_U,
  mean_cpm_strata = mean_cpm,
  expr_decile = dec(mean_cpm, in_U),
  infrv_strata = infrv,
  infrv_tertile = tert(infrv, in_U),
  uncertain = infrv >= ue$prereg$strata$infrv_uncertain_cutoff,
  infrv_strata_d15 = infrv15,
  infrv_tertile_d15 = tert(infrv15, in_U),
  od_strata = od_s$overdispersion,
  isoform_class = ann$isoform_class[match(ids, ann$transcript_id)],
  ctrl_total_count = ctrl_total
)
say(
  "U_mouse: %d transcripts (CPM >= 1 in >= 3 of 9 stressed); InfRV tertile cutoffs %.3f, %.3f; uncertain (InfRV >= 0.5) %.1f%% of U",
  sum(in_U),
  stats::quantile(infrv[in_U], 1 / 3),
  stats::quantile(infrv[in_U], 2 / 3),
  100 * mean(strata$uncertain[in_U])
)
say("RTA overdispersion prior of the 9 stressed: %.4f", attr(od_s, "prior"))

sup <- equivalence_support(file.path(quants, stre))
strata$unique_fraction <- sup$unique_fraction[match(ids, sup$transcript_id)]
strata$zero_unique <- strata$unique_fraction %in% 0
strata$identifiable <- NA
strata$non_identifiable <- NA
ue_write(strata, file.path(out, "strata.parquet"), ue, seed)

# 5. Validation of the matrix-level code ---------------------------------------------
od_c <- ue_overdispersion(mom_mean[, ctrl], mom_var[, ctrl], n_draws)
cs <- edgeR::catchSalmon(file.path(quants, ctrl), verbose = FALSE)
cs_od <- cs$annotation$Overdispersion[match(ids, rownames(cs$annotation))]
rel <- abs(od_c$overdispersion - cs_od) / cs_od
say(
  "overdispersion of the 10 controls vs edgeR::catchSalmon: prior %.6f vs %.6f, max relative difference %.3g over %d transcripts",
  attr(od_c, "prior"),
  cs$overdispersion.prior,
  max(rel, na.rm = TRUE),
  sum(!is.na(rel))
)
v_ids <- ids[keep_draws][seq(1, sum(keep_draws), length.out = 5000)]
se <- SummarizedExperiment::SummarizedExperiment(
  assays = c(
    list(counts = counts[v_ids, ctrl]),
    stats::setNames(
      lapply(seq_len(n_draws), function(k) draws_ctrl[v_ids, , k]),
      paste0("infRep", seq_len(n_draws))
    )
  )
)
se <- fishpond::computeInfRV(se, pc = 5, shift = 0.01)
my_rv <- rowMeans(ue_infrv(mom_mean[v_ids, ctrl], mom_var[v_ids, ctrl]))
rv_diff <- max(abs(SummarizedExperiment::mcols(se)$meanInfRV - my_rv))
say(
  "InfRV vs fishpond::computeInfRV on 5,000 transcripts x 10 controls: max abs difference %.3g",
  rv_diff
)
ue_write(
  data.frame(
    check = c(
      "overdispersion_vs_catchSalmon_max_rel_diff",
      "overdispersion_prior_isoformic",
      "overdispersion_prior_catchSalmon",
      "infrv_vs_computeInfRV_max_abs_diff"
    ),
    value = c(
      max(rel, na.rm = TRUE),
      attr(od_c, "prior"),
      cs$overdispersion.prior,
      rv_diff
    )
  ),
  file.path(out, "validation.parquet"),
  ue,
  seed
)
ue_write(
  data.frame(
    transcript_id = ids,
    od_ctrl = od_c$overdispersion,
    od_catchsalmon = cs_od
  ),
  file.path(out, "overdispersion_control.parquet"),
  ue,
  seed
)
rm(cs, se)

# 6. Per-sample draw correlations of sibling isoforms (10 controls) ------------------------
# Same statistic as inferential_correlation(): per sample, correlation of the
# standardised draws, a transcript masked where its draws do not vary or its
# sample mean is below min_count (1). Pairs within genes, among transcripts
# with control draws stored.
ids_d <- rownames(draws_ctrl)
gene_d <- ann$gene_id[match(ids_d, ann$transcript_id)]
ord <- order(gene_d)
g <- gene_d[ord]
rl <- rle(g)
end_of <- rep(cumsum(rl$lengths), rl$lengths)
remaining <- end_of - seq_along(g)
first <- rep(seq_along(g), remaining)
second <- first + sequence(remaining)
say(
  "sibling pairs among %d transcripts with control draws: %d",
  length(ids_d),
  length(first)
)
rmat <- matrix(
  NA_real_,
  length(first),
  length(ctrl),
  dimnames = list(NULL, ctrl)
)
for (s in seq_along(ctrl)) {
  x <- draws_ctrl[ord, s, ]
  mu <- rowMeans(x)
  cen <- x - mu
  sds <- sqrt(rowSums(cen^2) / (n_draws - 1L))
  z <- cen / sds
  z[sds == 0 | mu < 1, ] <- NA_real_
  for (st in seq.int(1L, length(first), by = 200000L)) {
    idx <- st:min(st + 199999L, length(first))
    rmat[idx, s] <- rowSums(
      z[first[idx], , drop = FALSE] * z[second[idx], , drop = FALSE]
    ) /
      (n_draws - 1L)
  }
}
pairs <- data.frame(
  gene_id = g[first],
  transcript_id_1 = ids_d[ord][first],
  transcript_id_2 = ids_d[ord][second],
  rmat,
  check.names = FALSE
)
pairs <- pairs[rowSums(!is.na(rmat)) > 0, ]
ue_write(pairs, file.path(out, "pair_correlation_control.parquet"), ue, seed)
# Validation on the genes of 500 random pairs against the snapshot function
vg <- unique(sample(pairs$gene_id, 500))
vrows <- which(gene_d %in% vg)
mc <- rowMeans(mom_mean[ids_d[vrows], ctrl])
ref <- inferential_correlation(
  draws_ctrl[vrows, , , drop = FALSE],
  sample_means = mom_mean[ids_d[vrows], ctrl],
  gene_id = gene_d[vrows],
  assessed = mc >= 1,
  min_count = 1
)
mine <- pairs[pairs$gene_id %in% vg, ]
mine$mean_cor <- rowMeans(as.matrix(mine[, ctrl]), na.rm = TRUE)
key_r <- paste(ref$transcript_id_1, ref$transcript_id_2)
key_m1 <- paste(mine$transcript_id_1, mine$transcript_id_2)
key_m2 <- paste(mine$transcript_id_2, mine$transcript_id_1)
pos <- ifelse(
  is.na(match(key_r, key_m1)),
  match(key_r, key_m2),
  match(key_r, key_m1)
)
cor_diff <- max(abs(ref$mean_cor - mine$mean_cor[pos]), na.rm = TRUE)
say(
  "pair correlation vs inferential_correlation(): %d pairs compared, %d unmatched, max abs difference %.3g",
  nrow(ref),
  sum(is.na(pos)),
  cor_diff
)
ue_write(
  data.frame(
    check = c(
      "pair_cor_vs_inferential_correlation_max_abs_diff",
      "pair_cor_unmatched"
    ),
    value = c(cor_diff, sum(is.na(pos)))
  ),
  file.path(out, "validation_pairs.parquet"),
  ue,
  seed
)
# Identifiability last (slow); strata.parquet is rewritten with it
# One sample at a time, checkpointed in prepare/identifiability/ (a run of
# all 9 samples at once was killed for memory); the per-sample tables are
# summed as equivalence_identifiability() does over several samples, the
# same as 10_placenta_prepare.R.
if (do_ident) {
  t0 <- Sys.time()
  id_dir <- file.path(out, "identifiability")
  dir.create(id_dir, showWarnings = FALSE)
  id_all <- do.call(
    rbind,
    lapply(stre, function(r) {
      f <- file.path(id_dir, paste0(r, ".parquet"))
      if (!file.exists(f)) {
        t1 <- Sys.time()
        one <- equivalence_identifiability(file.path(quants, r))
        ue_write(
          as.data.frame(one)[, c(
            "transcript_id",
            "n_samples",
            "n_identifiable"
          )],
          f,
          ue,
          seed
        )
        say(
          "identifiability %s: %.1f min; peak RSS %.2f GB",
          r,
          as.numeric(difftime(Sys.time(), t1, units = "mins")),
          ue_peak_rss_gb()
        )
      }
      as.data.frame(read_parquet(f))[, c(
        "transcript_id",
        "n_samples",
        "n_identifiable"
      )]
    })
  )
  n_s <- tapply(id_all$n_samples, id_all$transcript_id, sum)
  n_i <- tapply(id_all$n_identifiable, id_all$transcript_id, sum)
  m <- match(ids, names(n_s))
  strata$identifiable_n_samples <- as.integer(n_s[m])
  strata$identifiable_fraction <- as.numeric(n_i[m] / n_s[m])
  strata$identifiable <- as.logical(n_i[m] == n_s[m])
  strata$non_identifiable <- strata$identifiable %in% FALSE
  say(
    "identifiability on 9 samples: %.1f min; non-identifiable in U: %d of %d",
    as.numeric(difftime(Sys.time(), t0, units = "mins")),
    sum(strata$non_identifiable[in_U]),
    sum(in_U)
  )
} else {
  strata$identifiable <- NA
  strata$non_identifiable <- NA
  say("identifiability skipped (MOUSE_IDENT=false)")
}
gU <- strata[strata$in_U, ]
g_max <- tapply(gU$infrv_strata, gU$gene_id, max)
g_nid <- tapply(gU$non_identifiable, gU$gene_id, function(v) {
  any(v, na.rm = TRUE)
})
strata$gene_max_infrv <- unname(g_max[strata$gene_id])
strata$gene_any_non_identifiable <- unname(g_nid[strata$gene_id])
ue_write(strata, file.path(out, "strata.parquet"), ue, seed)
say("peak RSS %.2f GB", ue_peak_rss_gb())
ue_write_session(out, ue)
say("done")
