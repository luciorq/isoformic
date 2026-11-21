# 41_mouse_null.R: genome-wide null of the mouse BNST control males (WP-C1).
#
# Purpose: replicate the placenta null (WP-A) in an independent dataset. The
# 10 control males have no true difference between any two subsets, so every
# call of a 5 vs 5 split is false. All 126 distinct 5 vs 5 splits are tested
# with the pre-registered methods M1, M2, M6, M7, M8, M9, M11 (M10, satuRn
# and DEXSeq, on 10 splits only, as prereg `deseq2_dexseq_satuRn_splits`;
# the pilot measured DEXSeq at 2.9 min per split, under the 20 min stop rule).
# The universe U_mouse and the strata come from the 9 stressed males, which
# are never tested (40_mouse_prepare.R).
#
# Methods use the WP-A implementations (lib_null_methods.R sources
# null_common.R): M1/M2 test U_mouse rows kept by a design-free
# edgeR::filterByExpr() (labels unseen), M2 divides by the RTA
# overdispersion of the design's own animals, M6 is swish_test() (length
# corrected, 100 permutations) on the rows M1 tests, M11 is edgeR QL on gene
# sums of all transcripts; M7 keeps M1 p-values of class-high transcripts
# (confidence classes from the design's own animals). False-call rates are
# reported per tested feature (the WP-A definition) and per feature of U.
# M10 = satuRn and DEXSeq gene-level usage on 10 random full splits.
#
# Inference: delete-one-animal jackknife. All 126 5 vs 5 splits use all 10
# animals, so no 5 vs 5 split lacks an animal (the plan's "absent from 63 of
# the 126 splits" does not hold); the statistic is therefore defined as a
# functional of the animals available: the mean over the balanced splits of
# the available animals (5 vs 5 for 10, 4 vs 5 for 9; 126 distinct splits in
# both cases). For each deleted animal, 28 of its 126 4 vs 5 splits were
# sampled first (prereg "28-per-unit"; fits/<block>/jk<a>.rds), and the
# other 98 are fitted into fits/<block>_rest/jk<a>.rds (MOUSE_JK_EXACT), so
# that the jackknife is exact: the verification of WP-C showed that the
# inner-bootstrap Monte Carlo correction of 28-split units removed only part
# of their Monte Carlo variance (M1 adjusted rate ratio: 28-split unit
# estimates 4 times as dispersed as the exact ones). A method without all
# 126 splits per unit gets the uncorrected (conservative) 28-split interval;
# the `jackknife` column says which. M1full (M1 on 126 splits per unit, the
# m1jk block) is kept as an independent check of M1. The overdispersion
# (M2, M9) and the confidence classes (M7) are recomputed from the animals
# of each design.
#
# Library batch: the control animals form two library batches (Salmon
# fragment-length means 119-150 vs 172-179 bp; the short-fragment batch is
# the accession block SRR24948951-955). Most false calls come from the 26
# full splits aligned with the batch, so every statistic is also reported
# over the batch-balanced (3-2) and the batch-aligned splits only
# (*_batch_balanced, *_batch_aligned); null_mouse/control_batches.parquet
# holds the fragment lengths and the batch of every control animal.
#
# Inputs: $ISOFORMIC_UE_OUT/null_mouse/prepare/ (40_mouse_prepare.R),
#   prereg.yaml and the snapshot (ue_common.R), lib_null_methods.R.
# Env vars (besides those of ue_common.R):
#   MOUSE_STAGE    fit, analyse or all (default all)
#   MOUSE_METHODS  comma list of method blocks to fit, default
#                  "dte,usage,swish,m10" (dte = M1, M2, M7, M11; usage = M8, M9;
#                  swish = M6; m10 = satuRn and DEXSeq on 10 splits)
#   MOUSE_PILOT    if set to n > 0: fit only the first n full splits into
#                  null_mouse/pilot/ and report timings (no analysis)
#   MOUSE_WORKERS  parallel workers (fork), default 4
#   MOUSE_SWISH_WORKERS  workers of the swish block (default MOUSE_WORKERS)
#   MOUSE_M1_FULLJK  "true" (default) runs M1 on all 126 splits per deleted unit
#   MOUSE_JK_EXACT comma list of blocks whose remaining 98 splits per deleted
#                  animal are fitted (default "dte,swish,usage")
#   MOUSE_JK_BOOT  inner bootstrap replicates of the Monte Carlo correction of
#                  28-split units (default 0 = no correction, conservative)
#   MOUSE_QUANTS   Salmon quants of the mouse runs (fragment lengths), default
#                  $ISOFORMIC_V50_DIR/cases/mouse_bnst/quants
# Outputs: $ISOFORMIC_UE_OUT/null_mouse/fits/<block>/<design>.rds (checkpoints:
#   sparse per-split results: features with p < 0.05 or q < 0.1),
#   splits.parquet (splits, depth and mapping-rate imbalance, batch balance),
#   control_batches.parquet,
#   per_split_counts.parquet, per_split_sig.parquet, gene_strata.parquet,
#   m10_splits.parquet, results_mouse_null.parquet (+ .tsv), 41_mouse_null.log
#
# Memory: draws are held for U_mouse only (66,259 x 10 x 30 doubles, 0.16 GB);
# the pilot peaked at 2.7 GB (main process, 2 workers).
# Time (pilot timings): about 2.5 h with 3 workers for the fit stage, and
# about 4 h more for the exact-jackknife splits (estimate); the analyse stage
# took 26 to 36 min. Peak RSS of the main process 1.8 GB.
# Run from the repository root, each stage inside a memory cap of about 8 GB
# (see README.md):
#   MOUSE_STAGE=fit MOUSE_WORKERS=3 MOUSE_SWISH_WORKERS=2 \
#     pixi run --frozen Rscript data-raw/uncertainty_evidence/41_mouse_null.R
#   MOUSE_STAGE=analyse \
#     pixi run --frozen Rscript data-raw/uncertainty_evidence/41_mouse_null.R

source("data-raw/uncertainty_evidence/ue_common.R")
ue <- ue_setup("41_mouse_null")
source("data-raw/uncertainty_evidence/lib_null_methods.R")
seed0 <- ue$seed_base + 1000L * ue$prereg$seeds$package_index$C
stage <- Sys.getenv("MOUSE_STAGE", "all")
blocks <- strsplit(Sys.getenv("MOUSE_METHODS", "dte,usage,swish,m10"), ",")[[1]]
pilot <- as.integer(Sys.getenv("MOUSE_PILOT", "0"))
workers <- as.integer(Sys.getenv("MOUSE_WORKERS", "4"))
swish_workers <- as.integer(Sys.getenv("MOUSE_SWISH_WORKERS", workers))
m1_fulljk <- !identical(tolower(Sys.getenv("MOUSE_M1_FULLJK", "true")), "false")
jk_exact <- strsplit(Sys.getenv("MOUSE_JK_EXACT", "dte,swish,usage"), ",")[[1]]
root <- file.path(ue_out_root(), "null_mouse")
prep <- file.path(root, "prepare")
outd <- if (pilot > 0) file.path(root, "pilot") else root
dir.create(file.path(outd, "fits"), recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(outd, "41_mouse_null.log")
say <- function(...) {
  msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " ", sprintf(...))
  cat(msg, "\n", file = log_file, append = TRUE)
  message(msg)
}
rd <- function(f) as.data.frame(read_parquet(file.path(prep, f)))
say(
  "start stage=%s blocks=%s pilot=%d workers=%d sha=%s",
  stage,
  paste(blocks, collapse = ","),
  pilot,
  workers,
  ue$git_sha
)

# Data ---------------------------------------------------------------------------------
samples <- rd("samples.parquet")
ctrl <- samples$run[samples$role == "tested"]
stre <- samples$run[samples$role == "strata"]
stopifnot(length(ctrl) == 10L)
ann <- rd("tx2gene.parquet")
strata <- rd("strata.parquet")
cnt <- rd("counts.parquet")
counts <- as.matrix(cnt[, ctrl])
rownames(counts) <- cnt$transcript_id
counts_stre <- as.matrix(cnt[, stre])
rownames(counts_stre) <- cnt$transcript_id
rm(cnt)
gene_all <- ann$gene_id[match(rownames(counts), ann$transcript_id)]
# Moments and draws through DuckDB queries, one control animal at a time, so
# that no long data frame of all runs is held (memory: the machine is shared)
dq <- function(sql) {
  con <- DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  DBI::dbExecute(con, "SET threads = 2")
  DBI::dbGetQuery(con, sql)
}
pq <- function(f) DBI::dbQuoteString(DBI::ANSI(), file.path(prep, f))
ids_all <- rownames(counts)
mom_mean <- mom_var <- matrix(
  NA_real_,
  length(ids_all),
  length(ctrl),
  dimnames = list(ids_all, ctrl)
)
dcols <- sprintf("d%02d", 1:30)
# Draws are held only for the transcripts of U_mouse (M1, M2, M6 and M7 test
# U only, as WP-A); InfRV and the `assessed` rule of the confidence classes
# come from the moments of every transcript (nm_conf_design(mean =, var =)).
U <- strata$transcript_id[strata$in_U]
ids_d <- dq(sprintf(
  "SELECT DISTINCT transcript_id FROM read_parquet(%s) ORDER BY transcript_id",
  pq("draws_control.parquet")
))$transcript_id
n_draws_stored <- length(ids_d)
ids_d <- ids_d[ids_d %in% U]
arr <- array(
  0,
  c(length(ids_d), length(ctrl), 30L),
  dimnames = list(ids_d, ctrl, NULL)
)
for (r in ctrl) {
  m <- dq(sprintf(
    "SELECT transcript_id, draw_mean, draw_var FROM read_parquet(%s) WHERE run = '%s'",
    pq("draw_moments.parquet"),
    r
  ))
  i <- match(m$transcript_id, ids_all)
  mom_mean[i, r] <- m$draw_mean
  mom_var[i, r] <- m$draw_var
  sub <- dq(sprintf(
    "SELECT transcript_id, %s FROM read_parquet(%s) WHERE run = '%s'",
    paste(dcols, collapse = ", "),
    pq("draws_control.parquet"),
    r
  ))
  sub <- sub[sub$transcript_id %in% ids_d, ]
  arr[match(sub$transcript_id, ids_d), r, ] <- as.matrix(sub[, dcols])
  rm(m, sub)
}
stopifnot(!anyNA(mom_mean), !anyNA(mom_var))
pairs <- dq(sprintf(
  "SELECT transcript_id_1, transcript_id_2, %s FROM read_parquet(%s)",
  paste0('"', ctrl, '"', collapse = ", "),
  pq("pair_correlation_control.parquet")
))
# Check the moments against the stored draws
chk <- max(abs(rowMeans(arr[, ctrl[1], ]) - mom_mean[ids_d, ctrl[1]]))
stopifnot(chk < 1e-6)
invisible(gc())

st_U <- strata[strata$in_U, ]
gene_U <- unique(st_U$gene_id)
say(
  "U_mouse %d transcripts, %d genes; draws held for %d transcripts of U (of %d stored)",
  length(U),
  length(gene_U),
  length(ids_d),
  n_draws_stored
)

# Designs and splits ----------------------------------------------------------------------
# full: the 10 animals, 126 splits = 5-subsets containing animal 1 (group b =
# the complement); jk_a: animal a deleted, 126 4-subsets of the remaining 9
# as group b.
mk_splits <- function(animals, size_b) {
  cmb <- utils::combn(length(animals), size_b)
  if (2L * size_b == length(animals)) {
    cmb <- cmb[, cmb[1, ] != 1L, drop = FALSE] # one of each complementary pair
  }
  lapply(seq_len(ncol(cmb)), function(j) animals[cmb[, j]])
}
designs <- list(full = list(animals = ctrl, splits = mk_splits(ctrl, 5L)))
for (a in seq_along(ctrl)) {
  an <- ctrl[-a]
  designs[[paste0("jk", a)]] <- list(animals = an, splits = mk_splits(an, 4L))
}
stopifnot(vapply(designs, function(d) length(d$splits), 1L) == 126L)
set.seed(seed0 + 200L)
for (a in seq_along(ctrl)) {
  designs[[paste0("jk", a)]]$sample28 <- sort(sample.int(126L, 28L))
}
depth <- stats::setNames(samples$num_processed, samples$run)
maprate <- stats::setNames(samples$percent_mapped, samples$run)
# Library batch of the control animals. The verification of WP-C found that
# nearly all false calls come from splits that put the accession block
# SRR24948951-955 against the other five controls. The Salmon fragment
# length distributions separate the same two groups (mean 119-150 vs
# 172-179 bp), a property of the library preparation, so the batch is
# defined from the data: the two groups on either side of the largest gap
# in the sorted fragment-length means of the 10 controls. A split is
# batch-balanced when the counts of batch-1 animals in its two groups
# differ by at most 1 (3-2 for 5 vs 5), batch-aligned otherwise. Animal 1
# (SRR24948951) is in group a of every full split (mk_splits()), so the
# signed depth difference of a split is not symmetric; use its absolute value.
quants_dir <- Sys.getenv("MOUSE_QUANTS")
if (!nzchar(quants_dir)) {
  quants_dir <- file.path(ue_v50_dir(), "cases", "mouse_bnst", "quants")
}
frag <- vapply(
  ctrl,
  function(r) {
    jsonlite::read_json(file.path(
      quants_dir,
      r,
      "aux_info",
      "meta_info.json"
    ))$frag_length_mean
  },
  numeric(1)
)
fs <- sort(frag)
cut_frag <- mean(fs[which.max(diff(fs)) + 0:1])
batch1 <- names(frag)[frag < cut_frag]
say(
  "fragment length means of the controls: %s; largest gap at %.1f bp; batch 1 (short fragments): %s",
  paste(sprintf("%s %.1f", names(frag), frag), collapse = ", "),
  cut_frag,
  paste(batch1, collapse = ", ")
)
split_tab <- do.call(
  rbind,
  lapply(names(designs), function(dn) {
    d <- designs[[dn]]
    do.call(
      rbind,
      lapply(seq_along(d$splits), function(j) {
        b <- d$splits[[j]]
        a <- setdiff(d$animals, b)
        data.frame(
          design = dn,
          split = j,
          in_sample28 = if (dn == "full") TRUE else j %in% d$sample28,
          group_b = paste(b, collapse = ";"),
          log2_depth_diff = mean(log2(depth[b])) - mean(log2(depth[a])),
          maprate_diff = mean(maprate[b]) - mean(maprate[a]),
          batch1_in_a = sum(a %in% batch1),
          batch1_in_b = sum(b %in% batch1),
          batch_balanced = abs(sum(a %in% batch1) - sum(b %in% batch1)) <= 1L
        )
      })
    )
  })
)
if (pilot == 0) {
  ue_write(split_tab, file.path(outd, "splits.parquet"), ue, seed0 + 200L)
  ue_write(
    data.frame(
      run = ctrl,
      frag_length_mean = unname(frag),
      batch = ifelse(ctrl %in% batch1, "batch1_short_fragments", "batch2"),
      depth = unname(depth[ctrl]),
      percent_mapped = unname(maprate[ctrl])
    ),
    file.path(outd, "control_batches.parquet"),
    ue,
    seed0
  )
}

# Per design quantities -------------------------------------------------------------------
od_of <- function(animals) {
  od <- ue_overdispersion(mom_mean[, animals], mom_var[, animals], 30L)
  stats::setNames(od$overdispersion, od$transcript_id)
}

# Sparse result of one method on one split: features with p < 0.05 or
# q < 0.1, plus the tested set as indices into the method's universe (U for
# transcript methods, the genes of U for gene methods). Every filter is
# label-free (design-free filterByExpr, filter_usage() and labelKeep() use
# only the group sizes), so the tested set is the same in all splits of a
# design: run_block() keeps the indices of the first split only and a
# checksum of every split, which the analysis checks.
sparse <- function(res, method, uid) {
  keep <- (!is.na(res$pvalue) & res$pvalue < 0.05) |
    (!is.na(res$qvalue) & res$qvalue < 0.10)
  ti <- match(res$feature_id[!is.na(res$pvalue)], uid)
  ti <- sort(ti[!is.na(ti)])
  list(
    method = method,
    n_tested = sum(!is.na(res$pvalue)),
    n_tested_uni = length(ti),
    tested_sum = sum(as.numeric(ti)),
    tested_idx = ti,
    sig = res[keep, , drop = FALSE]
  )
}
gene_U <- sort(gene_U)
efflen <- {
  e <- dq(sprintf(
    "SELECT transcript_id, %s FROM read_parquet(%s)",
    paste0('"', ctrl, '"', collapse = ", "),
    pq("effective_length.parquet")
  ))
  m <- as.matrix(e[, ctrl])
  rownames(m) <- e$transcript_id
  m[rownames(counts), , drop = FALSE]
}

fit_split_dte <- function(d, j, od, cd) {
  b <- d$splits[[j]]
  an <- d$animals
  g <- factor(ifelse(an %in% b, "b", "a"))
  cs <- counts[, an]
  t0 <- proc.time()[[3]]
  m1 <- nm_edger_dte(cs, g, rows = U)
  t1 <- proc.time()[[3]]
  m2 <- nm_edger_dte(cs, g, od, rows = U)
  t2 <- proc.time()[[3]]
  m11 <- nm_gene_dte(cs, gene_all, g)
  t3 <- proc.time()[[3]]
  conf <- nm_conf_split(cd, an %in% b)
  m1$confidence <- conf$confidence[match(m1$feature_id, conf$transcript_id)]
  t4 <- proc.time()[[3]]
  # M7 keeps the M1 call only when the class is high; its p-value is M1's
  # p-value for class-high transcripts and 1 otherwise (for the p < 0.01 rate)
  m7 <- m1
  m7$pvalue[!(m7$confidence %in% "high")] <- 1
  m7$qvalue[!(m7$confidence %in% "high")] <- 1
  list(
    M1 = sparse(m1, "M1", U),
    M2 = sparse(m2, "M2", U),
    M7 = sparse(m7, "M7", U),
    M11 = sparse(m11, "M11", gene_U),
    n_high_tested = sum(m1$confidence %in% "high"),
    n_conf_na_tested = sum(is.na(m1$confidence)),
    secs = c(M1 = t1 - t0, M2 = t2 - t1, M11 = t3 - t2, M7 = t4 - t3)
  )
}

fit_split_m1 <- function(d, j) {
  b <- d$splits[[j]]
  an <- d$animals
  g <- factor(ifelse(an %in% b, "b", "a"))
  list(M1 = sparse(nm_edger_dte(counts[, an], g, rows = U), "M1", U))
}

fit_split_usage <- function(d, j, od) {
  b <- d$splits[[j]]
  an <- d$animals
  g <- factor(ifelse(an %in% b, "b", "a"))
  cs <- counts[, an]
  t0 <- proc.time()[[3]]
  m8 <- nm_usage(cs, gene_all, g)
  t1 <- proc.time()[[3]]
  m9 <- nm_usage(cs, gene_all, g, od = od)
  t2 <- proc.time()[[3]]
  gs <- function(x, m) {
    gg <- x$genes
    names(gg) <- c("feature_id", "pvalue", "qvalue")
    sparse(gg, m, gene_U)
  }
  list(
    M8 = gs(m8, "M8"),
    M9 = gs(m9, "M9"),
    M8_tx = sparse(m8$transcripts, "M8_tx", U),
    M9_tx = sparse(m9$transcripts, "M9_tx", U),
    secs = c(M8 = t1 - t0, M9 = t2 - t1)
  )
}

# M6 on the rows M1 tests in this design (as 21_null_dte.R: swish on keep1)
fit_split_swish <- function(d, j, rows, unit_seed) {
  b <- d$splits[[j]]
  an <- d$animals
  g <- factor(ifelse(an %in% b, "b", "a"))
  t0 <- proc.time()[[3]]
  m6 <- nm_swish(
    arr[rows, an, , drop = FALSE],
    counts[rows, an],
    efflen[rows, an],
    g,
    seed = unit_seed
  )
  list(M6 = sparse(m6, "M6", U), secs = c(M6 = proc.time()[[3]] - t0))
}

run_block <- function(block, dn, idx, dir_block = block) {
  f <- file.path(outd, "fits", dir_block, paste0(dn, ".rds"))
  if (file.exists(f)) {
    return(invisible(readRDS(f)))
  }
  dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)
  d <- designs[[dn]]
  t0 <- Sys.time()
  dnum <- if (dn == "full") 0L else as.integer(sub("jk", "", dn, fixed = TRUE))
  res <- switch(
    block,
    dte = {
      od <- od_of(d$animals)
      lib <- colSums(counts[, d$animals])
      cd <- nm_conf_design(
        arr,
        match(d$animals, ctrl),
        lib,
        ids_d,
        pairs,
        mean = mom_mean[, d$animals],
        var = mom_var[, d$animals]
      )
      parallel::mclapply(
        idx,
        function(j) fit_split_dte(d, j, od, cd),
        mc.cores = workers,
        mc.preschedule = FALSE
      )
    },
    m1jk = parallel::mclapply(
      idx,
      function(j) fit_split_m1(d, j),
      mc.cores = workers,
      mc.preschedule = FALSE
    ),
    usage = {
      od <- od_of(d$animals)
      parallel::mclapply(
        idx,
        function(j) fit_split_usage(d, j, od),
        mc.cores = workers,
        mc.preschedule = FALSE
      )
    },
    swish = {
      cu <- counts[rownames(counts) %in% U, d$animals]
      rows <- rownames(cu)[nm_design_free_keep(
        cu,
        colSums(counts[, d$animals])
      )]
      stopifnot(rows %in% ids_d)
      parallel::mclapply(
        idx,
        function(j) {
          fit_split_swish(d, j, rows, seed0 + 300L + 126L * dnum + j)
        },
        mc.cores = swish_workers,
        mc.preschedule = FALSE
      )
    }
  )
  # Tested indices kept for the first split of the block only (see sparse())
  for (k in seq_along(res)[-1L]) {
    for (m in names(res[[k]])) {
      if (is.list(res[[k]][[m]]) && !is.null(res[[k]][[m]]$tested_idx)) {
        res[[k]][[m]]$tested_idx <- NULL
      }
    }
  }
  bad <- vapply(res, function(x) inherits(x, "try-error") || is.null(x), TRUE)
  if (any(bad)) {
    stop(
      block,
      " ",
      dn,
      ": ",
      sum(bad),
      " splits failed: ",
      paste(utils::head(unlist(res[bad])), collapse = " | "),
      call. = FALSE
    )
  }
  names(res) <- as.character(idx)
  say(
    "%s %s: %d splits in %.1f min",
    block,
    dn,
    length(idx),
    as.numeric(difftime(Sys.time(), t0, units = "mins"))
  )
  saveRDS(res, f)
  invisible(res)
}

# Fit -----------------------------------------------------------------------------------
if (stage %in% c("fit", "all")) {
  if (pilot > 0) {
    # Check the two-step confidence: (1) moments-based design vs the
    # draws-based design on the rows held; (2) nm_conf_split() against
    # nm_confidence() (draws-based) on 3000 of those rows
    d <- designs$full
    lib <- colSums(counts[, d$animals])
    cd_m <- nm_conf_design(
      arr,
      seq_along(ctrl),
      lib,
      ids_d,
      pairs,
      mean = mom_mean[, ctrl],
      var = mom_var[, ctrl]
    )
    rv_a <- rowMeans(inferential_rv(list(
      mean = rowMeans(arr, dims = 2L),
      var = apply(arr, c(1, 2), stats::var)
    )))
    say(
      "moments vs draws on %d rows of U: max |InfRV diff| %.3g, assessed mismatches %d",
      length(ids_d),
      max(abs(cd_m$mean_inf_rv - rv_a)),
      sum(cd_m$assessed != (rowMeans(rowMeans(arr, dims = 2L)) >= 1))
    )
    in_b <- d$animals %in% d$splits[[1]]
    a1 <- nm_conf_split(cd_m, in_b)
    chk_ids <- ids_d[seq(1, length(ids_d), length.out = 3000)]
    sub_ids <- ids_d
    a2 <- nm_confidence(
      arr,
      seq_along(ctrl),
      in_b,
      lib,
      chk_ids,
      pairs[
        pairs$transcript_id_1 %in% sub_ids & pairs$transcript_id_2 %in% sub_ids,
      ]
    )
    cmp <- a1[match(chk_ids, a1$transcript_id), ]
    say(
      "two-step vs nm_confidence on 3000 transcripts (pairs within U): %d sign-agreement diffs > 1e-12, %d InfRV-flag mismatches",
      sum(abs(cmp$sign_agreement - a2$sign_agreement) > 1e-12),
      sum(cmp$uncertain != a2$uncertain)
    )
    rm(cd_m, rv_a)
    # Vectorised Simes against the package adjust_usage() on the first split
    g <- factor(ifelse(d$animals %in% d$splits[[1]], "b", "a"))
    u <- nm_usage(counts[, d$animals], gene_all, g)
    pk <- adjust_usage(tibble::as_tibble(u$transcripts[, c(
      "transcript_id",
      "usage_log2FC",
      "pvalue",
      "gene_id"
    )]))
    pg <- unique(pk[, c("gene_id", "gene_pvalue", "gene_qvalue")])
    mm <- match(pg$gene_id, u$genes$gene_id)
    say(
      "vectorised Simes vs adjust_usage(): %d genes, max |gene p diff| %.3g, max |gene q diff| %.3g",
      nrow(pg),
      max(abs(pg$gene_pvalue - u$genes$gene_pvalue[mm])),
      max(abs(pg$gene_qvalue - u$genes$gene_qvalue[mm]))
    )
    idx <- seq_len(pilot)
    tm <- list()
    for (blk in intersect(blocks, c("dte", "usage", "swish"))) {
      t0 <- Sys.time()
      r <- run_block(blk, "full", idx)
      el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
      secs <- do.call(rbind, lapply(r, `[[`, "secs"))
      tm[[blk]] <- data.frame(
        block = blk,
        splits = pilot,
        wall_min = el,
        per_split_min_parallel = el / pilot,
        t(colMeans(secs))
      )
      say(
        "pilot %s: %.2f min for %d splits (%d workers); mean secs per method: %s",
        blk,
        el,
        pilot,
        workers,
        paste(
          names(colMeans(secs)),
          sprintf("%.1f", colMeans(secs)),
          collapse = ", "
        )
      )
    }
    if ("m10" %in% blocks) {
      d <- designs$full
      g <- factor(ifelse(d$animals %in% d$splits[[1]], "b", "a"))
      for (bk in c("saturn", "dexseq")) {
        t0 <- Sys.time()
        x <- try(nm_usage(
          counts[, d$animals],
          gene_all,
          g,
          backend = bk,
          workers = workers
        ))
        say(
          "pilot M10 %s: %.2f min (%d workers)%s",
          bk,
          as.numeric(difftime(Sys.time(), t0, units = "mins")),
          workers,
          if (inherits(x, "try-error")) paste(" ERROR", x) else ""
        )
      }
    }
    say("peak RSS of the main process %.2f GB", ue_peak_rss_gb())
    quit(save = "no")
  }
  for (blk in intersect(blocks, c("dte", "usage", "swish"))) {
    run_block(blk, "full", seq_len(126L))
    for (a in seq_along(ctrl)) {
      dn <- paste0("jk", a)
      run_block(blk, dn, designs[[dn]]$sample28)
    }
  }
  if (m1_fulljk && "dte" %in% blocks) {
    for (a in seq_along(ctrl)) {
      run_block("m1jk", paste0("jk", a), seq_len(126L))
    }
  }
  # Exact jackknife: the 98 splits of each deleted animal that the 28-split
  # sample left out, into fits/<block>_rest/ (the analysis appends them)
  for (blk in intersect(jk_exact, c("dte", "swish", "usage"))) {
    for (a in seq_along(ctrl)) {
      dn <- paste0("jk", a)
      run_block(
        blk,
        dn,
        setdiff(seq_len(126L), designs[[dn]]$sample28),
        dir_block = paste0(blk, "_rest")
      )
    }
  }
  if ("m10" %in% blocks) {
    f <- file.path(outd, "fits", "m10", "full.rds")
    if (!file.exists(f)) {
      dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)
      set.seed(seed0 + 210L)
      j10 <- sort(sample.int(126L, 10L))
      d <- designs$full
      res <- lapply(j10, function(j) {
        g <- factor(ifelse(d$animals %in% d$splits[[j]], "b", "a"))
        out <- list()
        for (bk in c("saturn", "dexseq")) {
          t0 <- Sys.time()
          x <- nm_usage(
            counts[, d$animals],
            gene_all,
            g,
            backend = bk,
            workers = workers
          )
          gg <- x$genes
          names(gg) <- c("feature_id", "pvalue", "qvalue")
          out[[bk]] <- sparse(gg, bk, gene_U)
          out[[paste0(bk, "_min")]] <- as.numeric(difftime(
            Sys.time(),
            t0,
            units = "mins"
          ))
        }
        say(
          "M10 split %d: satuRn %.1f min, DEXSeq %.1f min",
          j,
          out$saturn_min,
          out$dexseq_min
        )
        out
      })
      names(res) <- as.character(j10)
      saveRDS(res, f)
    }
  }
  say("fit done; peak RSS of the main process %.2f GB", ue_peak_rss_gb())
}

# Analysis ------------------------------------------------------------------------------
if (stage %in% c("analyse", "all") && pilot == 0) {
  say("analysis start")
  # 0 (default): no Monte Carlo correction for 28-split units (uncorrected,
  # conservative intervals); exact units never need one
  B_inner <- as.integer(Sys.getenv("MOUSE_JK_BOOT", "0"))
  # Records: rec[[method]][[design]][[split]] = sparse() list
  rec <- list()
  add_block <- function(block, methods, dns) {
    for (dn in dns) {
      f <- file.path(outd, "fits", block, paste0(dn, ".rds"))
      if (!file.exists(f)) {
        next
      }
      r <- readRDS(f)
      # The remaining 98 splits of a deleted animal (exact jackknife)
      fr <- file.path(outd, "fits", paste0(block, "_rest"), paste0(dn, ".rds"))
      if (block != "m1jk" && dn != "full" && file.exists(fr)) {
        r <- c(r, readRDS(fr))
        r <- r[order(as.integer(names(r)))]
      }
      for (m in methods) {
        key <- if (block == "m1jk") "M1full" else m
        # jarl-ignore undesirable_operator: appends to a result list of the enclosing scope
        rec[[key]][[dn]] <<- lapply(r, `[[`, m)
      }
    }
  }
  dns <- names(designs)
  add_block("dte", c("M1", "M2", "M7", "M11"), dns)
  add_block("usage", c("M8", "M9", "M8_tx", "M9_tx"), dns)
  add_block("swish", "M6", dns)
  add_block("m1jk", "M1", setdiff(dns, "full"))
  if (!is.null(rec$M1full)) {
    rec$M1full$full <- rec$M1$full
  }
  # A method is exact when every deleted animal has all 126 of its splits
  is_exact <- function(m) {
    all(vapply(
      paste0("jk", seq_along(ctrl)),
      function(dn) length(rec[[m]][[dn]]) == 126L,
      TRUE
    ))
  }
  exact <- vapply(names(rec), is_exact, TRUE)
  say(
    "exact jackknife (126 splits per deleted animal): %s; 28 sampled splits, uncorrected: %s",
    paste(names(exact)[exact], collapse = ", "),
    paste(names(exact)[!exact], collapse = ", ")
  )
  methods_tx <- intersect(c("M1", "M2", "M6", "M7", "M1full"), names(rec))
  methods_gene <- intersect(c("M8", "M9", "M11"), names(rec))
  say("methods with fits: %s", paste(names(rec), collapse = ", "))

  # Strata of transcripts (U) and of genes (G_U)
  st_U$tertile <- factor(st_U$infrv_tertile, levels = c("low", "mid", "high"))
  st_U$tertile15 <- factor(
    st_U$infrv_tertile_d15,
    levels = c("low", "mid", "high")
  )
  st_U$decile <- factor(st_U$expr_decile)
  st_U$iso <- factor(st_U$isoform_class, levels = c("1", "2-4", "5-9", "10+"))
  lib_s <- colSums(counts_stre)
  gcpm <- rowsum(sweep(counts_stre, 2, lib_s / 1e6, "/"), gene_all)
  gene_tab <- data.frame(gene_id = gene_U)
  gene_tab$mean_cpm <- rowMeans(gcpm)[gene_tab$gene_id]
  gene_tab$max_infrv <- st_U$gene_max_infrv[match(
    gene_tab$gene_id,
    st_U$gene_id
  )]
  gene_tab$n_iso <- ann$n_isoforms_index[match(gene_tab$gene_id, ann$gene_id)]
  q3 <- stats::quantile(gene_tab$max_infrv, c(1 / 3, 2 / 3))
  gene_tab$tertile <- cut(
    gene_tab$max_infrv,
    c(-Inf, q3, Inf),
    labels = c("low", "mid", "high")
  )
  gene_tab$decile <- factor(cut(
    gene_tab$mean_cpm,
    c(-Inf, stats::quantile(gene_tab$mean_cpm, seq(0.1, 0.9, 0.1)), Inf),
    labels = FALSE
  ))
  gene_tab$iso <- factor(cut(
    gene_tab$n_iso,
    c(0, 1, 4, 9, Inf),
    labels = c("1", "2-4", "5-9", "10+")
  ))
  ue_write(gene_tab, file.path(outd, "gene_strata.parquet"), ue, seed0)
  universe_of <- function(m) if (m %in% methods_gene) gene_tab else st_U
  id_of <- function(m) if (m %in% methods_gene) "gene_id" else "transcript_id"

  # Tested set of each method and design (label-free filters: one set per
  # design, checked against the checksum of every split)
  tested_set <- list()
  for (m in names(rec)) {
    for (dn in names(rec[[m]])) {
      xs <- rec[[m]][[dn]]
      sums <- vapply(xs, function(x) x$tested_sum, numeric(1))
      has <- which(!vapply(xs, function(x) is.null(x$tested_idx), TRUE))
      first <- xs[[has[1]]]$tested_idx
      if (length(unique(sums)) != 1L || sum(as.numeric(first)) != sums[1]) {
        stop("tested set differs between splits of ", m, " ", dn, call. = FALSE)
      }
      tested_set[[m]][[dn]] <- first
    }
  }
  say("tested sets constant within every design and method (checksums)")

  # Per split counts and sparse significant tables
  per_split <- list()
  sig_rows <- list()
  for (m in names(rec)) {
    if (m %in% c("M8_tx", "M9_tx")) {
      next
    }
    uni <- universe_of(m)
    uid <- uni[[id_of(m)]]
    for (dn in names(rec[[m]])) {
      for (j in names(rec[[m]][[dn]])) {
        x <- rec[[m]][[dn]][[j]]
        sg <- x$sig
        inU <- sg$feature_id %in% uid
        tert <- uni$tertile[match(sg$feature_id, uid)]
        p01 <- !is.na(sg$pvalue) & sg$pvalue < 0.01
        tt <- uni$tertile[tested_set[[m]][[dn]]]
        per_split[[length(per_split) + 1]] <- data.frame(
          method = m,
          design = dn,
          split = as.integer(j),
          n_tested = x$n_tested,
          n_tested_U = length(tt),
          n_tested_low = sum(tt %in% "low"),
          n_tested_mid = sum(tt %in% "mid"),
          n_tested_high = sum(tt %in% "high"),
          n_q05 = sum(sg$qvalue < 0.05, na.rm = TRUE),
          n_q05_U = sum(sg$qvalue < 0.05 & inU, na.rm = TRUE),
          n_q10_U = sum(sg$qvalue < 0.10 & inU, na.rm = TRUE),
          n_p01_U = sum(p01 & inU),
          n_p001_U = sum(sg$pvalue < 0.001 & inU, na.rm = TRUE),
          n_p01_low = sum(p01 & tert %in% "low"),
          n_p01_mid = sum(p01 & tert %in% "mid"),
          n_p01_high = sum(p01 & tert %in% "high")
        )
        if (m != "M1full" && nrow(sg)) {
          sig_rows[[length(sig_rows) + 1]] <- data.frame(
            method = m,
            design = dn,
            split = as.integer(j),
            feature_id = sg$feature_id,
            pvalue = sg$pvalue,
            qvalue = sg$qvalue
          )
        }
      }
    }
  }
  per_split <- do.call(rbind, per_split)
  per_split <- merge(
    per_split,
    split_tab[, c("design", "split", "log2_depth_diff", "maprate_diff")],
    by = c("design", "split")
  )
  ue_write(per_split, file.path(outd, "per_split_counts.parquet"), ue, seed0)
  ue_write(
    do.call(rbind, sig_rows),
    file.path(outd, "per_split_sig.parquet"),
    ue,
    seed0
  )
  rm(sig_rows)

  # Keys: a split set is a data frame of (design, split)
  keys_full <- data.frame(design = "full", split = seq_len(126L))
  keys_units <- function(m, subset = "all") {
    lapply(seq_along(ctrl), function(a) {
      dn <- paste0("jk", a)
      k <- data.frame(design = dn, split = as.integer(names(rec[[m]][[dn]])))
      sub_keys(k, subset)
    })
  }
  # Batch subsets of a split set (see batch1 above)
  bal_key <- paste(split_tab$design, split_tab$split)
  sub_keys <- function(k, subset) {
    if (subset == "all") {
      return(k)
    }
    b <- split_tab$batch_balanced[match(paste(k$design, k$split), bal_key)]
    k[if (subset == "balanced") b else !b, , drop = FALSE]
  }
  ps_key <- paste(per_split$method, per_split$design, per_split$split)
  getps <- function(m, keys, col) {
    per_split[[col]][match(paste(m, keys$design, keys$split), ps_key)]
  }
  # Per-transcript number of splits with p < 0.01 over a split set
  sig_p01 <- function(m, keys) {
    uni <- universe_of(m)
    uid <- uni[[id_of(m)]]
    y <- integer(length(uid))
    for (i in seq_len(nrow(keys))) {
      sg <- rec[[m]][[keys$design[i]]][[as.character(keys$split[i])]]$sig
      f <- sg$feature_id[!is.na(sg$pvalue) & sg$pvalue < 0.01]
      w <- match(f, uid)
      w <- w[!is.na(w)]
      y[w] <- y[w] + 1L
    }
    y
  }
  # Number of splits of a split set in which each feature was tested
  tested_count <- function(m, keys) {
    uni <- universe_of(m)
    n <- integer(nrow(uni))
    for (dn in unique(keys$design)) {
      w <- tested_set[[m]][[dn]]
      n[w] <- n[w] + sum(keys$design == dn)
    }
    n
  }
  # p < 0.01 rate within an InfRV tertile among the tested features (the
  # WP-A definition) and per feature of the universe (untested count as p >= 0.01)
  rate_tert <- function(m, keys, lvl) {
    sum(getps(m, keys, paste0("n_p01_", lvl))) /
      sum(getps(m, keys, paste0("n_tested_", lvl)))
  }
  rate_tert_U <- function(m, keys, lvl) {
    uni <- universe_of(m)
    sum(getps(m, keys, paste0("n_p01_", lvl))) /
      (nrow(keys) * sum(uni$tertile %in% lvl))
  }
  adj_rr <- function(m, keys, tert_col = "tertile") {
    uni <- universe_of(m)
    d <- data.frame(
      y = sig_p01(m, keys),
      n = tested_count(m, keys),
      t = uni[[tert_col]],
      dec = uni$decile,
      iso = uni$iso
    )
    d <- d[d$n > 0, ]
    # Not estimable without events in both the high and the low tertile
    if (sum(d$y[d$t %in% "high"]) == 0 || sum(d$y[d$t %in% "low"]) == 0) {
      return(NA_real_)
    }
    fit <- suppressWarnings(stats::glm(
      y ~ t + dec + iso + offset(log(n)),
      family = stats::quasipoisson(),
      data = d
    ))
    unname(exp(stats::coef(fit)["thigh"]))
  }
  mh_rr <- function(m, keys) {
    uni <- universe_of(m)
    y <- sig_p01(m, keys)
    n <- tested_count(m, keys)
    k <- interaction(uni$decile, uni$iso, drop = TRUE)
    hi <- uni$tertile %in% "high"
    lo <- uni$tertile %in% "low"
    a <- tapply(y * hi, k, sum)
    b <- tapply(y * lo, k, sum)
    N1 <- tapply(n * hi, k, sum)
    N0 <- tapply(n * lo, k, sum)
    tot <- N1 + N0
    ok <- tot > 0
    sum((a * N0 / tot)[ok]) / sum((b * N1 / tot)[ok])
  }
  jk <- function(
    m,
    stat,
    log = FALSE,
    subset = "all",
    units = keys_units(m, subset),
    pool = 126,
    B = B_inner
  ) {
    # Exact units hold every split of their subset: pool = their size
    if (exact[[m]]) {
      pool <- max(vapply(units, nrow, 1L))
    }
    nm_jackknife(
      sub_keys(keys_full, subset),
      units,
      function(k) stat(m, k),
      pool = pool,
      B = B,
      log = log,
      seed = seed0 + 400L
    )
  }
  results <- list()
  put <- function(m, what, r, note = "") {
    # jarl-ignore undesirable_operator: appends to a result list of the enclosing scope
    results[[length(results) + 1]] <<- data.frame(
      method = m,
      statistic = what,
      r,
      note = note
    )
  }
  # Per-method results are cached (the jackknife of the model-based
  # statistics takes minutes per method); delete analysis_cache/ to redo
  cache_dir <- file.path(outd, "analysis_cache_v2")
  dir.create(cache_dir, showWarnings = FALSE)
  for (m in c(methods_tx, methods_gene)) {
    if (is.null(rec[[m]]$full) || length(rec[[m]]) < 11) {
      next
    }
    cf <- file.path(cache_dir, paste0(m, ".rds"))
    if (file.exists(cf)) {
      results <- c(results, readRDS(cf))
      next
    }
    n_before <- length(results)
    level <- if (m %in% methods_gene) "genes" else "transcripts"
    put(
      m,
      "false_calls_q05_U_per_split",
      jk(m, function(m, k) mean(getps(m, k, "n_q05_U"))),
      paste("q < 0.05, features in the universe;", level)
    )
    put(
      m,
      "false_calls_q05_all_per_split",
      jk(m, function(m, k) mean(getps(m, k, "n_q05"))),
      paste("q < 0.05, all tested", level)
    )
    put(
      m,
      "false_calls_q10_U_per_split",
      jk(m, function(m, k) mean(getps(m, k, "n_q10_U")))
    )
    put(
      m,
      "splits_with_any_call_q05_U",
      jk(m, function(m, k) mean(getps(m, k, "n_q05_U") > 0)),
      "family-wise: share of splits with at least one call"
    )
    put(
      m,
      "p01_rate_tested_U",
      jk(
        m,
        function(m, k) {
          sum(getps(m, k, "n_p01_U")) / sum(getps(m, k, "n_tested_U"))
        },
        log = TRUE
      ),
      "p < 0.01 rate among the tested features of the universe; 0.01 under a calibrated null"
    )
    put(
      m,
      "p01_rate_per_U",
      jk(
        m,
        function(m, k) mean(getps(m, k, "n_p01_U")) / nrow(universe_of(m)),
        log = TRUE
      ),
      "per feature of the universe; untested count as p >= 0.01"
    )
    put(
      m,
      "n_tested_U_per_split",
      jk(m, function(m, k) mean(getps(m, k, "n_tested_U")))
    )
    put(
      m,
      "n_tested_per_split",
      jk(m, function(m, k) mean(getps(m, k, "n_tested")))
    )
    for (lv in c("low", "mid", "high")) {
      put(
        m,
        paste0("p01_inflation_", lv),
        jk(m, function(m, k) rate_tert(m, k, lv) / 0.01, log = TRUE),
        "p < 0.01 rate / 0.01 among the tested features of the InfRV tertile (strata samples)"
      )
      put(
        m,
        paste0("p01_inflation_perU_", lv),
        jk(m, function(m, k) rate_tert_U(m, k, lv) / 0.01, log = TRUE),
        "per feature of the tertile; untested count as p >= 0.01"
      )
    }
    put(
      m,
      "p01_ratio_high_low",
      jk(
        m,
        function(m, k) rate_tert(m, k, "high") / rate_tert(m, k, "low"),
        log = TRUE
      )
    )
    if (m != "M1full") {
      put(
        m,
        "adjusted_rate_ratio_high_low",
        jk(m, adj_rr, log = TRUE),
        "quasi-Poisson: splits with p < 0.01 ~ InfRV tertile + expression decile + isoform class"
      )
      put(
        m,
        "mh_rate_ratio_high_low",
        jk(m, mh_rr, log = TRUE),
        "Mantel-Haenszel rate ratio, strata decile x isoform class"
      )
      if (m %in% methods_tx) {
        put(
          m,
          "adjusted_rate_ratio_high_low_d15",
          jk(m, function(m, k) adj_rr(m, k, "tertile15"), log = TRUE),
          "robustness: tertiles from draws 1-15"
        )
      }
    }
    if (m == "M1full") {
      put(
        m,
        "adjusted_rate_ratio_high_low",
        jk(m, adj_rr, log = TRUE),
        "quasi-Poisson: splits with p < 0.01 ~ InfRV tertile + expression decile + isoform class"
      )
      put(
        m,
        "mh_rate_ratio_high_low",
        jk(m, mh_rr, log = TRUE),
        "Mantel-Haenszel rate ratio, strata decile x isoform class"
      )
    }
    # Library batch: statistics over the batch-balanced and the batch-aligned
    # splits only (each deleted animal: its splits of the same kind)
    for (sb in c("balanced", "aligned")) {
      if (m == "M1full") {
        next
      }
      nt <- sprintf("batch-%s splits only (fragment-length batch)", sb)
      put(
        m,
        paste0("false_calls_q05_U_per_split_batch_", sb),
        jk(m, function(m, k) mean(getps(m, k, "n_q05_U")), subset = sb),
        nt
      )
      put(
        m,
        paste0("splits_with_any_call_q05_U_batch_", sb),
        jk(m, function(m, k) mean(getps(m, k, "n_q05_U") > 0), subset = sb),
        nt
      )
      put(
        m,
        paste0("p01_rate_tested_U_batch_", sb),
        jk(
          m,
          function(m, k) {
            sum(getps(m, k, "n_p01_U")) / sum(getps(m, k, "n_tested_U"))
          },
          log = TRUE,
          subset = sb
        ),
        nt
      )
      put(
        m,
        paste0("p01_ratio_high_low_batch_", sb),
        jk(
          m,
          function(m, k) rate_tert(m, k, "high") / rate_tert(m, k, "low"),
          log = TRUE,
          subset = sb
        ),
        nt
      )
      put(
        m,
        paste0("adjusted_rate_ratio_high_low_batch_", sb),
        jk(m, adj_rr, log = TRUE, subset = sb),
        nt
      )
    }
    saveRDS(results[seq.int(n_before + 1L, length(results))], cf)
    say("analysed %s", m)
  }
  # Differences between methods on the same splits (log rate ratio difference)
  pairs_cmp <- list(c("M1", "M2"), c("M1", "M6"), c("M1", "M7"), c("M8", "M9"))
  for (pc in pairs_cmp) {
    if (!all(pc %in% names(rec))) {
      next
    }
    pool_pc <- if (all(exact[pc])) 126 else 1e9
    for (sb in c("all", "balanced")) {
      units <- keys_units(pc[2], sb)
      sfx <- if (sb == "all") "" else "_batch_balanced"
      pool_sb <- if (all(exact[pc])) max(vapply(units, nrow, 1L)) else pool_pc
      put(
        paste(pc, collapse = "-"),
        paste0("log_adjusted_rr_difference", sfx),
        nm_jackknife(
          sub_keys(keys_full, sb),
          units,
          function(k) log(adj_rr(pc[1], k)) - log(adj_rr(pc[2], k)),
          pool = pool_sb,
          B = B_inner,
          seed = seed0 + 401L
        ),
        "log RR(first) - log RR(second); > 0 means the second reduces the excess"
      )
      put(
        paste(pc, collapse = "-"),
        paste0("false_calls_q05_U_difference", sfx),
        nm_jackknife(
          sub_keys(keys_full, sb),
          units,
          function(k) {
            mean(getps(pc[1], k, "n_q05_U")) - mean(getps(pc[2], k, "n_q05_U"))
          },
          pool = pool_sb,
          B = B_inner,
          seed = seed0 + 402L
        )
      )
    }
  }
  # CMH test (independence across splits assumed: anti-conservative; a check only)
  for (m in intersect(c("M1", "M2", "M6", "M7", "M11"), names(rec))) {
    uni <- universe_of(m)
    y <- sig_p01(m, keys_full)
    k <- interaction(uni$decile, uni$iso, drop = TRUE)
    sel <- uni$tertile %in% c("low", "high")
    ev <- tapply(
      y[sel],
      list(factor(uni$tertile[sel], levels = c("high", "low")), k[sel]),
      sum
    )
    ev[is.na(ev)] <- 0
    arr3 <- array(0, c(2, 2, ncol(ev)))
    arr3[, 1, ] <- ev
    nt <- tapply(
      tested_count(m, keys_full)[sel],
      list(factor(uni$tertile[sel], levels = c("high", "low")), k[sel]),
      sum
    )
    nt[is.na(nt)] <- 0
    arr3[, 2, ] <- nt - ev
    keepk <- apply(arr3, 3, function(z) {
      all(rowSums(z) > 0) && all(colSums(z) > 0)
    })
    ct <- stats::mantelhaen.test(arr3[,, keepk])
    put(
      m,
      "cmh_common_or_high_low",
      data.frame(
        estimate = unname(ct$estimate),
        se = NA,
        se_uncorrected = NA,
        ci_lo = ct$conf.int[1],
        ci_hi = ct$conf.int[2],
        ci_lo_uncorrected = NA,
        ci_hi_uncorrected = NA,
        n_units = NA,
        mc_share = NA,
        mc_corrected = NA
      ),
      sprintf(
        "CMH over transcript-splits, p = %.3g; treats splits as independent, so the CI is too narrow",
        ct$p.value
      )
    )
  }
  # Depth imbalance and calls
  for (m in intersect(c("M1", "M2", "M6"), names(rec))) {
    x <- per_split[per_split$method == m & per_split$design == "full", ]
    ct <- stats::cor.test(
      abs(x$log2_depth_diff),
      x$n_q05_U,
      method = "spearman",
      exact = FALSE
    )
    put(
      m,
      "spearman_calls_vs_abs_depth_diff",
      data.frame(
        estimate = unname(ct$estimate),
        se = NA,
        se_uncorrected = NA,
        ci_lo = NA,
        ci_hi = NA,
        ci_lo_uncorrected = NA,
        ci_hi_uncorrected = NA,
        n_units = NA,
        mc_share = NA,
        mc_corrected = NA
      ),
      sprintf("126 full splits, p = %.3g (splits not independent)", ct$p.value)
    )
  }
  # M10 (10 splits) against M9 on the same splits
  f10 <- file.path(outd, "fits", "m10", "full.rds")
  if (file.exists(f10)) {
    r10 <- readRDS(f10)
    gid <- gene_tab$gene_id
    m10 <- do.call(
      rbind,
      lapply(names(r10), function(j) {
        x <- r10[[j]]
        s9 <- rec$M9$full[[j]]$sig
        s8 <- rec$M8$full[[j]]$sig
        data.frame(
          split = as.integer(j),
          saturn_q05_U = sum(
            x$saturn$sig$qvalue < 0.05 & x$saturn$sig$feature_id %in% gid,
            na.rm = TRUE
          ),
          dexseq_q05_U = sum(
            x$dexseq$sig$qvalue < 0.05 & x$dexseq$sig$feature_id %in% gid,
            na.rm = TRUE
          ),
          M9_q05_U = sum(
            s9$qvalue < 0.05 & s9$feature_id %in% gid,
            na.rm = TRUE
          ),
          M8_q05_U = sum(
            s8$qvalue < 0.05 & s8$feature_id %in% gid,
            na.rm = TRUE
          ),
          saturn_min = x$saturn_min,
          dexseq_min = x$dexseq_min
        )
      })
    )
    ue_write(m10, file.path(outd, "m10_splits.parquet"), ue, seed0 + 210L)
  }
  res <- do.call(rbind, results)
  res$n_splits_full <- 126L
  ex_of <- function(m) {
    ms <- strsplit(m, "-")[[1]]
    all(ms %in% names(exact)) && all(exact[ms])
  }
  res$jackknife <- ifelse(
    vapply(res$method, ex_of, TRUE),
    "exact: all splits per deleted animal",
    ifelse(
      res$mc_corrected %in% TRUE,
      "28 sampled splits per deleted animal, Monte Carlo corrected",
      "28 sampled splits per deleted animal, uncorrected (conservative)"
    )
  )
  ue_write(res, file.path(outd, "results_mouse_null.parquet"), ue, seed0 + 400L)
  utils::write.table(
    res,
    file.path(outd, "results_mouse_null.tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  ue_write_session(outd, ue)
  say("analysis done")
}
