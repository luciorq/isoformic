# WP-A step 3: transcript usage, isoform switches and productive-fraction
# (NMD) shifts on the placenta null splits.
#
# Purpose: count false usage calls (every call is false: control versus
# control) with the package's usage backends on matrices, scale = FALSE
# against scale = TRUE (counts divided by the overdispersion of the split's
# own samples, as test_usage(scale = TRUE)), by gene-level ambiguity
# stratum; and, on the 7v7 leave-one-unit-out splits, false isoform switches
# (isoform_switches() rule on M1 calls against M2 calls, with the M11 gene
# calls as `deg`) and false productive-fraction shifts
# (test_productive_fraction() test on raw against divided counts).
#   diffsplice   all 7v7_louo (210) and 5v5 (100) splits, scale FALSE and TRUE
#   satuRn       30 7v7_louo splits (the first two per left-out unit)
#   DEXSeq       10 7v7_louo splits, only if one split takes less than the
#                pre-registered stop rule (prereg.yaml stop_rules, 20 min)
# Matrix code: the package internals filter_usage(), usage_diffsplice(),
# usage_saturn(), usage_dexseq_fit(), simes_pvalue() (vectorised here, the
# same formula), two_group_proportion_test() and productivity classes from
# the S0 store (transcript_productivity() of the snapshot).
#
# Inputs: <out>/store (S0), <out>/null_placenta/cache/dte_cache.rds and
# dte/main/tx, dte/main/gene (21_null_dte.R), splits (20_null_splits.R).
# Env vars:
#   ISOFORMIC_UE_OUT       output root (ue_common.R)
#   ISOFORMIC_UE_WORKERS   forked workers for diffsplice (default 4)
#   ISOFORMIC_UE_PART      diffsplice | saturn | dexseq | switches (default all
#                          four, in that order)
# Resources (measured, peak RSS of the largest process): diffsplice 44 min
# with 2 workers, 2.4 GB; satuRn 20 min, 2.6 GB; DEXSeq 35 min, 3.1 GB; switches and NMD
# 13 min, 2.4 GB. Run inside a memory cap of about 8 GB (see README.md).
# Run from the repository root:
#   ISOFORMIC_UE_WORKERS=4 pixi run --frozen Rscript data-raw/uncertainty_evidence/22_null_dtu.R
# Outputs: <out>/null_placenta/dtu/<backend>/{genes,summary}/<split>_<scale>.parquet,
# dtu/switches/<split>.parquet, dtu/nmd/<split>.parquet, dtu/timing.parquet.

source("data-raw/uncertainty_evidence/null_common.R")
ue <- ue_setup("22_null_dtu")
p <- wpa_paths()
dir.create(p$dtu, recursive = TRUE, showWarnings = FALSE)
seed <- ue$seed_base + 1000L * ue$prereg$seeds$package_index$A + 300L
set.seed(seed)
workers <- as.integer(Sys.getenv("ISOFORMIC_UE_WORKERS", "4"))
parts <- strsplit(
  Sys.getenv("ISOFORMIC_UE_PART", "diffsplice,saturn,dexseq,switches"),
  ",",
  fixed = TRUE
)[[1L]]

C <- readRDS(file.path(p$cache, "dte_cache.rds"))
strata <- wpa_read(file.path(p$store, "strata_tx.parquet"))
gstrata <- wpa_read(file.path(p$store, "strata_gene.parquet"))
prod <- wpa_read(file.path(p$store, "productivity.parquet"))
splits <- wpa_read(file.path(p$splits, "splits.parquet"))
split_samples <- wpa_read(file.path(p$splits, "split_samples.parquet"))
prim_units <- strata$transcript_id[strata$on_primary_chrom %in% TRUE]
gene_of_all <- C$gene_all

split_setup <- function(sid) {
  sm <- split_samples[
    split_samples$split_id == sid & split_samples$member_set == "primary",
  ]
  s <- c(sm$sample_id[sm$arm == "a"], sm$sample_id[sm$arm == "b"])
  groups <- factor(
    rep(c("a", "b"), c(sum(sm$arm == "a"), sum(sm$arm == "b"))),
    levels = c("a", "b")
  )
  od_tab <- overdispersion_from_moments(list(
    mean = C$mom_mean[, s],
    var = C$mom_var[, s],
    n = 30L
  ))
  od <- stats::setNames(od_tab$overdispersion, od_tab$transcript_id)
  list(s = s, groups = groups, od = od)
}

# Vectorised adjust_usage(): BH over transcripts, Simes per gene, BH over genes
adjust_fast <- function(res) {
  res$qvalue <- stats::p.adjust(res$pvalue, "BH")
  ok <- !is.na(res$pvalue)
  d <- res[ok, c("gene_id", "pvalue")]
  d <- d[order(d$gene_id, d$pvalue), ]
  n <- stats::ave(d$pvalue, d$gene_id, FUN = length)
  r <- stats::ave(d$pvalue, d$gene_id, FUN = seq_along)
  sp <- tapply(n * d$pvalue / r, d$gene_id, min)
  genes <- data.frame(gene_id = unique(res$gene_id))
  genes$gene_pvalue <- as.numeric(sp[genes$gene_id])
  genes$gene_qvalue <- stats::p.adjust(genes$gene_pvalue, "BH")
  merge(res, genes, by = "gene_id", all.x = TRUE)
}

usage_one <- function(sid, backend, scale) {
  t0 <- Sys.time()
  st <- split_setup(sid)
  counts <- C$counts[prim_units, st$s, drop = FALSE]
  if (isTRUE(scale)) {
    counts <- counts / st$od[prim_units]
  }
  gene_of <- gene_of_all[prim_units]
  kept <- filter_usage(counts, gene_of, st$groups, 10, 0.1)
  counts <- counts[kept, , drop = FALSE]
  gene_of <- gene_of[kept]
  contrasts <- resolve_contrasts(NULL, levels(st$groups))
  res <- switch(
    backend,
    diffsplice = usage_diffsplice(counts, gene_of, st$groups, contrasts),
    satuRn = usage_saturn(counts, gene_of, st$groups, contrasts, 1L),
    DEXSeq = usage_dexseq(counts, gene_of, st$groups, contrasts, 1L)
  )
  res <- as.data.frame(res)
  res$gene_id <- gene_of[match(res$transcript_id, rownames(counts))]
  res <- adjust_fast(res)
  genes <- res[
    !duplicated(res$gene_id),
    c("gene_id", "gene_pvalue", "gene_qvalue")
  ]
  genes$n_tx <- as.integer(table(res$gene_id)[genes$gene_id])
  genes$n_tx_q05 <- as.integer(tapply(
    res$qvalue < 0.05,
    res$gene_id,
    sum,
    na.rm = TRUE
  )[genes$gene_id])
  genes$split_id <- sid
  genes$backend <- backend
  genes$scaled <- scale
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  # Summary by gene stratum
  gt <- gstrata$max_infrv_tertile[match(genes$gene_id, gstrata$gene_id)]
  gn <- as.character(gstrata$any_non_identifiable[match(
    genes$gene_id,
    gstrata$gene_id
  )])
  gt[is.na(gt)] <- "not_in_U"
  gn[is.na(gn)] <- "not_in_U"
  sv <- list(
    all = rep("all", nrow(genes)),
    gene_max_infrv_tertile = gt,
    gene_any_non_identifiable = gn
  )
  rows <- list()
  for (v in names(sv)) {
    for (l in unique(sv[[v]])) {
      w <- sv[[v]] == l
      rows[[length(rows) + 1L]] <- data.frame(
        split_id = sid,
        backend = backend,
        scaled = scale,
        stratum = v,
        level = l,
        n_genes_tested = sum(w),
        n_genes_q05 = sum(genes$gene_qvalue[w] < 0.05, na.rm = TRUE),
        n_genes_q10 = sum(genes$gene_qvalue[w] < 0.10, na.rm = TRUE),
        n_genes_p01 = sum(genes$gene_pvalue[w] < 0.01, na.rm = TRUE),
        n_tx_tested = sum(res$gene_id %in% genes$gene_id[w]),
        n_tx_q05 = sum(
          res$qvalue[res$gene_id %in% genes$gene_id[w]] < 0.05,
          na.rm = TRUE
        ),
        n_tx_p01 = sum(
          res$pvalue[res$gene_id %in% genes$gene_id[w]] < 0.01,
          na.rm = TRUE
        )
      )
    }
  }
  summ <- do.call(rbind, rows)
  summ$seconds <- secs
  summ$peak_rss_gb <- ue_peak_rss_gb()
  tag <- paste0(sid, "_", if (scale) "scaled" else "raw")
  ue_write(
    genes,
    file.path(p$dtu, backend, "genes", paste0(tag, ".parquet")),
    ue,
    seed
  )
  ue_write(
    summ,
    file.path(p$dtu, backend, "summary", paste0(tag, ".parquet")),
    ue,
    seed
  )
  data.frame(
    split_id = sid,
    backend = backend,
    scaled = scale,
    seconds = secs,
    peak_rss_gb = ue_peak_rss_gb()
  )
}

run_backend <- function(ids, backend, n_workers) {
  jobs <- expand.grid(
    sid = ids,
    scale = c(FALSE, TRUE),
    stringsAsFactors = FALSE
  )
  tag <- paste0(jobs$sid, "_", ifelse(jobs$scale, "scaled", "raw"))
  done <- file.exists(file.path(
    p$dtu,
    backend,
    "summary",
    paste0(tag, ".parquet")
  ))
  jobs <- jobs[!done, ]
  wpa_log(backend, ": ", nrow(jobs), " jobs (", sum(done), " done)")
  out <- parallel::mclapply(
    seq_len(nrow(jobs)),
    function(k) {
      r <- try(usage_one(jobs$sid[[k]], backend, jobs$scale[[k]]))
      if (inherits(r, "try-error")) {
        cat("FAILED ", backend, jobs$sid[[k]], r, "\n")
      }
      r
    },
    mc.cores = n_workers,
    mc.preschedule = FALSE
  )
  ok <- !vapply(out, inherits, logical(1L), "try-error")
  if (any(ok)) {
    tm <- do.call(rbind, out[ok])
    f <- file.path(p$dtu, paste0("timing_", backend, ".parquet"))
    if (file.exists(f)) {
      tm <- rbind(wpa_read(f), tm)
    }
    ue_write(tm, f, ue, seed)
  }
  invisible(out)
}

louo <- splits[splits$design == "7v7_louo", ]
louo <- louo[order(louo$rank_in_unit, louo$left_out), ]
if ("diffsplice" %in% parts) {
  ids <- c(louo$split_id, splits$split_id[splits$design == "5v5"])
  run_backend(ids, "diffsplice", workers)
}
if ("saturn" %in% parts) {
  run_backend(louo$split_id[1:30], "satuRn", workers)
}
if ("dexseq" %in% parts) {
  # Stop rule: one split first, the rest only when it is under the limit
  limit <- ue$prereg$stop_rules$dexseq_minutes_per_split * 60
  first <- louo$split_id[[1L]]
  run_backend(first, "DEXSeq", 2L)
  tm <- wpa_read(file.path(p$dtu, "timing_DEXSeq.parquet"))
  worst <- max(tm$seconds[tm$split_id == first])
  wpa_log(
    "DEXSeq first split: ",
    round(worst / 60, 1),
    " min (limit ",
    limit / 60,
    ")"
  )
  if (worst <= limit) {
    run_backend(louo$split_id[1:10], "DEXSeq", workers)
  } else {
    ue_write(
      data.frame(
        backend = "DEXSeq",
        first_split_minutes = worst / 60,
        limit_minutes = limit / 60,
        decision = "stopped"
      ),
      file.path(p$dtu, "dexseq_stop_rule.parquet"),
      ue,
      seed
    )
  }
}

# Switches and productive-fraction shifts on the 7v7_louo splits -----------------
if ("switches" %in% parts) {
  prod_ids <- prod$transcript_id
  is_prod <- stats::setNames(
    prod$productivity %in% c("productive", "altered_productive"),
    prod$transcript_id
  )
  coding_genes <- unique(gene_of_all[prod$transcript_id[
    !prod$productivity %in% "noncoding"
  ]])
  sw_one <- function(sid) {
    tx <- wpa_read(file.path(p$dte, "main", "tx", paste0(sid, ".parquet")))
    ge <- wpa_read(file.path(p$dte, "main", "gene", paste0(sid, ".parquet")))
    st <- split_setup(sid)
    out <- list()
    for (arm in c("naive", "aware")) {
      m <- if (arm == "naive") "M1" else "M2"
      g <- if (arm == "naive") "M11raw" else "M11div"
      pv <- tx[[paste0(m, "_p")]]
      ok <- !is.na(pv)
      q <- rep(NA_real_, length(pv))
      q[ok] <- stats::p.adjust(pv[ok], "BH")
      lfc <- tx[[paste0(m, "_lfc")]]
      sig <- ok & q <= 0.05 & abs(lfc) >= 1
      gid <- gene_of_all[tx$transcript_id]
      up <- tapply(sig & lfc > 0, gid, sum)
      dn <- tapply(sig & lfc < 0, gid, sum)
      nsig <- tapply(sig, gid, sum)
      gp <- ge[[paste0(g, "_p")]]
      gok <- !is.na(gp)
      gq <- rep(NA_real_, length(gp))
      gq[gok] <- stats::p.adjust(gp[gok], "BH")
      deg <- ge$gene_id[gok & gq <= 0.05 & abs(ge[[paste0(g, "_lfc")]]) >= 1]
      genes <- names(nsig)[nsig > 0]
      type <- ifelse(
        up[genes] > 0 & dn[genes] > 0,
        "opposite_direction",
        ifelse(genes %in% deg, "consistent", "det_not_deg")
      )
      gt <- gstrata$max_infrv_tertile[match(genes, gstrata$gene_id)]
      gt[is.na(gt)] <- "not_in_U"
      out[[arm]] <- data.frame(
        split_id = sid,
        arm = arm,
        n_genes_with_sig_isoform = length(genes),
        n_opposite_direction = sum(type == "opposite_direction"),
        n_det_not_deg = sum(type == "det_not_deg"),
        n_consistent = sum(type == "consistent"),
        n_opposite_top_tertile = sum(
          type == "opposite_direction" & gt == "high"
        ),
        n_opposite_low_tertile = sum(
          type == "opposite_direction" & gt == "low"
        ),
        n_deg = length(deg)
      )
    }
    sw <- do.call(rbind, out)
    # Productive fraction (test_productive_fraction() defaults: productive +
    # altered_productive over coding genes, no count filter)
    nm <- list()
    for (arm in c("raw", "divided")) {
      cnt <- C$counts[prod_ids, st$s, drop = FALSE]
      if (arm == "divided") {
        cnt <- cnt / st$od[prod_ids]
      }
      gids <- gene_of_all[prod_ids]
      total <- rowsum(cnt, gids)
      productive <- rowsum(cnt * is_prod[prod_ids], gids)
      keep <- rownames(total) %in% coding_genes
      test <- two_group_proportion_test(
        productive[keep, , drop = FALSE],
        total[keep, , drop = FALSE],
        st$groups
      )
      gt <- gstrata$max_infrv_tertile[match(
        rownames(total)[keep],
        gstrata$gene_id
      )]
      gt[is.na(gt)] <- "not_in_U"
      nm[[arm]] <- data.frame(
        split_id = sid,
        arm = arm,
        n_genes_tested = sum(!is.na(test$pvalue)),
        n_padj05 = sum(test$padj < 0.05, na.rm = TRUE),
        n_padj10 = sum(test$padj < 0.10, na.rm = TRUE),
        n_p01 = sum(test$pvalue < 0.01, na.rm = TRUE),
        n_padj05_top_tertile = sum(
          test$padj < 0.05 & gt == "high",
          na.rm = TRUE
        ),
        n_tested_top_tertile = sum(!is.na(test$pvalue) & gt == "high")
      )
    }
    ue_write(
      sw,
      file.path(p$dtu, "switches", paste0(sid, ".parquet")),
      ue,
      seed
    )
    ue_write(
      do.call(rbind, nm),
      file.path(p$dtu, "nmd", paste0(sid, ".parquet")),
      ue,
      seed
    )
    sid
  }
  ids <- louo$split_id[file.exists(file.path(
    p$dte,
    "main",
    "tx",
    paste0(louo$split_id, ".parquet")
  ))]
  ids <- ids[!file.exists(file.path(p$dtu, "nmd", paste0(ids, ".parquet")))]
  wpa_log("switches / NMD: ", length(ids), " splits")
  r <- parallel::mclapply(
    ids,
    function(s) try(sw_one(s)),
    mc.cores = workers,
    mc.preschedule = FALSE
  )
  bad <- vapply(r, inherits, logical(1L), "try-error")
  if (any(bad)) print(r[bad][1:min(3, sum(bad))])
}
ue_write_session(p$dtu, ue)
wpa_log("done")
