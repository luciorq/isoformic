# WP-A step 4: gene-thinning arm on the placenta null splits (power and
# between-gene FDR).
#
# Purpose: plant known changes in real control-versus-control splits and
# measure empirical FDR and TPR of M1, M2, M6, M7, M9 and M10. In each split,
# 5% of the filtered multi-isoform genes (sampled within each tertile of the
# gene's maximum PE InfRV) plus as many single-isoform genes get a log2FC
# drawn uniformly from +/-[0.58, 2]; every transcript of such a gene is
# thinned binomially in one arm (the arm that goes down) with f = 2^-|log2FC|.
# Draws of the thinned rows: d' = f*m + sqrt(f)*(d - m) (m = mean of the
# sample's draws; keeps the ratio of draw variance to mean, so the
# overdispersion is unchanged), sensitivity variant d' = f*d (variance scaled
# by f^2). The moments of every folded transcript are transformed the same
# way, so M2's prior is recomputed on the thinned data.
# Label: power and between-gene FDR only. Whole-gene thinning keeps the
# isoform ratios, so this arm tests neither misattribution between isoforms
# nor transcript usage.
# Truth: a tested transcript is DE when its gene was thinned.
#
# Inputs: cache and run_split() of 21_null_dte.R (sourced with
# options(wpa.library_only = TRUE)), splits (20_null_splits.R), S0 strata.
# Env vars:
#   ISOFORMIC_UE_OUT       output root (ue_common.R)
#   ISOFORMIC_UE_WORKERS   forked workers (default 4)
#   ISOFORMIC_UE_THIN_DRAWS  affine (default) or scaled (d' = f*d)
#   ISOFORMIC_UE_THIN_N    splits per size (default 30)
# Resources (measured with 1 worker): affine 23 min, scaled 32 min, peak
# 6.0 GB each. Run inside a memory cap of about 8 GB (see
# README.md).
# Run from the repository root:
#   ISOFORMIC_UE_WORKERS=4 pixi run --frozen Rscript data-raw/uncertainty_evidence/23_thin_dte.R
# Outputs: <out>/null_placenta/thin/<draw rule>/{tx,truth,metrics}/<split>.parquet

options(wpa.library_only = TRUE)
Sys.setenv(
  ISOFORMIC_UE_VARIANT = "main",
  ISOFORMIC_UE_METHODS = "M1,M2,M6,M7,M9"
)
source("data-raw/uncertainty_evidence/21_null_dte.R")
ue$step <- "23_thin_dte"
rule <- Sys.getenv("ISOFORMIC_UE_THIN_DRAWS", "affine")
n_per <- as.integer(Sys.getenv("ISOFORMIC_UE_THIN_N", "30"))
tdir <- file.path(p$thin, rule)
dir.create(tdir, recursive = TRUE, showWarnings = FALSE)
gst <- wpa_read(file.path(p$store, "strata_gene.parquet"))

louo <- splits[splits$design == "7v7_louo", ]
louo <- louo[order(louo$rank_in_unit, louo$left_out), ]
ids <- c(
  louo$split_id[seq_len(n_per)],
  splits$split_id[splits$design == "5v5"][seq_len(n_per)]
)

# M6/M7 must run on the thinning splits whatever methods_for() says for 7v7_louo/5v5
thin_one <- function(sid) {
  k <- match(sid, splits$split_id)
  # Same seed for both draw rules: the same genes, fold-changes and thinned
  # counts, so the two rules differ only in the draws
  seed <- seed_pkg + 5000L + k
  set.seed(seed)
  sm <- split_samples[
    split_samples$split_id == sid & split_samples$member_set == "primary",
  ]
  s <- c(sm$sample_id[sm$arm == "a"], sm$sample_id[sm$arm == "b"])
  in_b <- c(rep(FALSE, sum(sm$arm == "a")), rep(TRUE, sum(sm$arm == "b")))
  # Filtered genes (design-free filter on the raw counts of the split)
  lib <- colSums(C$counts[, s])
  keep1 <- C$U[design_free_keep(C$counts[C$U, s], lib)]
  genes_f <- unique(C$gene_all[keep1])
  n_iso <- strata$n_isoforms_gene[match(genes_f, strata$gene_id)]
  multi <- genes_f[n_iso >= 2L]
  single <- genes_f[n_iso == 1L]
  tert <- gst$max_infrv_tertile[match(multi, gst$gene_id)]
  pick_multi <- unlist(lapply(split(multi, tert), function(g) {
    sample(g, round(0.05 * length(g)))
  }))
  pick_single <- sample(single, min(length(single), length(pick_multi)))
  picked <- c(pick_multi, pick_single)
  lfc <- stats::runif(length(picked), 0.58, 2) *
    sample(c(-1, 1), length(picked), replace = TRUE)
  names(lfc) <- picked
  # Copy of the cache restricted to the split's samples, thinned. Vectorised
  # over the cells of the thinned genes: factor f = 2^-|log2FC| in the arm that
  # goes down (b up, lfc > 0: thin arm a; b down: thin arm b), 1 elsewhere.
  CC <- C
  cnt <- C$counts[, s]
  mmean <- C$mom_mean[, s]
  mvar <- C$mom_var[, s]
  ar <- C$arr[, s, , drop = FALSE]
  factor_rows <- function(genes) {
    lf <- lfc[genes]
    f <- 2^-abs(lf)
    thin_b <- lf < 0
    fm <- matrix(1, length(genes), length(s))
    fm[, !in_b] <- ifelse(thin_b, 1, f)
    fm[, in_b] <- ifelse(thin_b, f, 1)
    fm
  }
  r_all <- which(C$gene_all[rownames(cnt)] %in% picked)
  F_all <- factor_rows(C$gene_all[rownames(cnt)[r_all]])
  y <- cnt[r_all, , drop = FALSE]
  cells <- which(F_all < 1)
  y[cells] <- stats::rbinom(length(cells), round(y[cells]), F_all[cells])
  cnt[r_all, ] <- y
  mmean[r_all, ] <- F_all * mmean[r_all, , drop = FALSE]
  mvar[r_all, ] <- (if (rule == "affine") F_all else F_all^2) *
    mvar[r_all, , drop = FALSE]
  r_U <- which(C$gene_all[C$U] %in% picked)
  if (length(r_U)) {
    F_U <- factor_rows(C$gene_all[C$U[r_U]])
    d <- ar[r_U, , , drop = FALSE]
    dm <- rowMeans(d, dims = 2L)
    for (k in seq_len(dim(d)[3L])) {
      d[,, k] <- if (rule == "affine") {
        F_U * dm + sqrt(F_U) * (d[,, k] - dm)
      } else {
        F_U * d[,, k]
      }
    }
    ar[r_U, , ] <- d
  }
  CC$counts <- cnt
  CC$mom_mean <- mmean
  CC$mom_var <- mvar
  CC$arr <- ar
  CC$efflen <- C$efflen[, s]
  CC$pc30 <- C$pc30[, s]
  rm(cnt, mmean, mvar, ar, y, d)
  # Gene moments are not used here (M11 is not run)
  r <- run_split(sid, CC = CC)
  tab <- r$tab
  tab$de_true <- C$gene_all[tab$transcript_id] %in% picked
  tab$true_log2fc <- unname(lfc[C$gene_all[tab$transcript_id]])
  tab$tertile <- strata$infrv_tertile[match(
    tab$transcript_id,
    strata$transcript_id
  )]
  tab$thin_rule <- rule
  ue_write(tab, file.path(tdir, "tx", paste0(sid, ".parquet")), ue, seed)
  ue_write(
    data.frame(
      split_id = sid,
      gene_id = picked,
      log2fc = unname(lfc),
      multi_isoform = picked %in% pick_multi
    ),
    file.path(tdir, "truth", paste0(sid, ".parquet")),
    ue,
    seed
  )
  # Metrics
  calls <- method_calls(tab)
  k_m2 <- sum(calls$M2$q < 0.05, na.rm = TRUE)
  out <- list()
  for (m in names(calls)) {
    pv <- calls[[m]]$p
    q <- calls[[m]]$q
    tested <- calls[[m]]$tested %||% !is.na(pv)
    for (sub in c("all", "top_tertile")) {
      w <- tested & (sub == "all" | tab$tertile %in% "high")
      for (thr in c(0.05, 0.10)) {
        call <- w & q < thr
        call[is.na(call)] <- FALSE
        out[[length(out) + 1L]] <- data.frame(
          split_id = sid,
          design = splits$design[k],
          method = m,
          subset = sub,
          rule = paste0("q<", thr),
          n_tested = sum(w),
          n_true = sum(w & tab$de_true),
          n_calls = sum(call),
          n_false = sum(call & !tab$de_true),
          n_true_calls = sum(call & tab$de_true)
        )
      }
      # FDR at matched k: the k smallest p-values, k = M2 calls at q < 0.05
      kk <- min(k_m2, sum(w))
      o <- order(pv[w])
      top <- which(w)[o[seq_len(kk)]]
      out[[length(out) + 1L]] <- data.frame(
        split_id = sid,
        design = splits$design[k],
        method = m,
        subset = sub,
        rule = "top_k_M2",
        n_tested = sum(w),
        n_true = sum(w & tab$de_true),
        n_calls = kk,
        n_false = sum(!tab$de_true[top]),
        n_true_calls = sum(tab$de_true[top])
      )
    }
  }
  met <- do.call(rbind, out)
  met$n_genes_thinned <- length(picked)
  met$thin_rule <- rule
  ue_write(met, file.path(tdir, "metrics", paste0(sid, ".parquet")), ue, seed)
  sid
}

# run_split() picks M6/M7 by design; make sure the thinning splits get them
methods_for <- function(design, split_id) c("M1", "M2", "M6", "M7", "M9")
todo <- ids[!file.exists(file.path(tdir, "metrics", paste0(ids, ".parquet")))]
wpa_log("thinning (", rule, "): ", length(todo), " splits")
res <- parallel::mclapply(
  todo,
  function(sid) {
    out <- try(thin_one(sid))
    if (inherits(out, "try-error")) {
      cat("FAILED ", sid, ": ", out, "\n")
    }
    out
  },
  mc.cores = workers,
  mc.preschedule = FALSE
)
wpa_log("done; failed: ", sum(vapply(res, inherits, logical(1L), "try-error")))
ue_write_session(tdir, ue)
