# 44_rhalf_replicate.R: replication of the flag result between technical read
# halves of the ESRP knockdown runs (WP-C2).
#
# Purpose: the read-half arm of the replication analysis. Each of the 12
# knockdown libraries (GSE129540 SRR8870780-85, ESRP2 siRNA, 3 control + 3
# knockdown; GSE75489 SRR2966412-17, ESRP1/2 shRNA, 3 + 3) was split into
# two random, disjoint halves of its read pairs (42_split_reads.sh) and each
# half quantified with Salmon and 30 Gibbs draws (43_rhalf_quant.sh). The two
# halves of one library differ only by read sampling and the EM, the
# quantity InfRV and the Gibbs draws estimate; biological and library
# variation is shared. This arm therefore measures sampling and EM noise
# only: it is a technical consistency check, not proof of quantification
# error, and it says nothing about biological replication.
#
# Per series and in both directions (A to B, B to A):
#   * Discovery on half D (3 vs 3): edgeR QL DTE on transcript point counts
#     (M1, primary) and on counts divided by the RTA overdispersion of half
#     D (M2, Baldoni et al. 2024). Calls: FDR < 0.10, with the WP-B fallback
#     to the top 500 by p when there are fewer than 100 calls.
#   * Flags from the draws of half D only (wpb_confidence(), the
#     add_uncertainty() + de_confidence() rules, as WP-B): primary flag =
#     class != high; secondary: sign agreement < 0.95, InfRV >= 0.5
#     (uncertain), top InfRV tertile of the expressed transcripts of half D,
#     ambiguous sibling.
#   * Validation on half V: the same edgeR QL test on its point counts. A
#     call replicates when the half-V logFC has the same sign and one-sided
#     p < 0.05 (prereg `primary`). M2 calls are validated by the same raw
#     test on half V, so the yardstick is identical for M1 and M2.
#   * Statistics pooled over both directions, per series: replication of
#     flagged and unflagged calls, deficit, deficit within joint deciles of
#     count and |z|, adjusted OR (logistic: replicated ~ flag + log2 count +
#     |z| + log2 isoforms + transcript type), delta-AUC, count-matched flag
#     C4 (per direction as many calls as the primary flag, the lowest mean
#     counts), placebo flag (200 permutations within count x |z| cells).
#     Calls that half V's own filterByExpr() drops have no validation
#     p-value and count as not replicated; the rates, deficits and counts
#     among calls tested in half V (*_testable, n_untestable_v) are reported
#     beside them. summary_c4 also gives the overlap of the unflagged sets
#     of the primary flag and C4.
#     CIs: bootstrap of genes (all calls of a resampled gene, in both
#     directions, kept together), 1000 resamples, percentile intervals.
#   * Permutation baseline: all 10 distinct 3 vs 3 partitions of the 6
#     libraries (the true one among them), each in both directions, with the
#     same labels in both halves. Because the halves share every library,
#     any between-library difference replicates; this baseline shows how
#     much of the replication is shared-library signal.
#   * Delta-AUC transfer: the logistic models of non-replication (count +
#     |z|, with and without the primary flag) fitted on one series and
#     evaluated on the other, both ways; CI from a gene bootstrap of the
#     evaluation series with the models fixed.
#   * Sensitivity: the existing consecutive halves (first and second 6 M
#     read pairs of the 12 M file) of GSE129540 in human/quants_half.
#   * Known ESRP targets (data-raw/build_knockdown_example.R) are listed
#     descriptively only; they never select transcripts or strata.
#
# Haplotype and patch copies are folded per draw onto their primary
# transcript with the S0 fold map (store/fold_map.parquet, prereg amendment
# A1); only folded units on the primary chromosomes are tested.
#
# Inputs:
#   $ISOFORMIC_UE_OUT/rhalf/quants/<GSE>_<SRR>_<A|B>/   (43_rhalf_quant.sh)
#   $ISOFORMIC_V50_DIR/human/quants_half/<SRR>_h1|h2 (consecutive halves)
#   $ISOFORMIC_UE_OUT/store/fold_map.parquet, store/tx2gene.parquet (S0)
#   prereg.yaml and the package snapshot (ue_common.R), 38_audit_common.R
# Env vars (besides ISOFORMIC_UE_OUT, ISOFORMIC_UE_PKG, ISOFORMIC_UE_THREADS):
#   RH_STAGE     pairs | summarise | all (default all)
#   RH_WORKERS   forked workers for the permutation fits and the bootstrap (2)
#   RH_BOOT      bootstrap resamples (1000)
#   ISOFORMIC_V50_DIR  work directory of the Salmon runs (required for the
#                consecutive halves of the pairs stage)
# Outputs: $ISOFORMIC_UE_OUT/rhalf/replication/
#   pairs/<halving>_<series>_<method>_<dir>.parquet  call tables
#   perm/<halving>_<series>_p<k>_<dir>.parquet         permutation call tables
#   summary_flags.parquet, summary_c4.parquet, summary_placebo.parquet,
#   summary_perm.parquet, summary_transfer.parquet, known_targets.parquet,
#   summary_info.parquet, 44_rhalf_replicate.log, session_info.parquet
#
# Seeds: base + 1000 * 3 (package C) + 600 (bootstrap), + 601 (placebo),
#   + 602 (transfer).
#
# Run from the repository root. Measured: pairs stage 4 to 8 min per
# series and halving, peak 4.4 GB with RH_WORKERS=1 (3 forked workers
# exceeded an 8 GB cap); summarise about 15 min, 0.7 GB. Run each stage
# inside a memory cap of 8 GB and 6 GB (see README.md):
#   RH_STAGE=pairs RH_WORKERS=1 \
#     pixi run --frozen Rscript data-raw/uncertainty_evidence/44_rhalf_replicate.R
#   RH_STAGE=summarise RH_WORKERS=3 \
#     pixi run --frozen Rscript data-raw/uncertainty_evidence/44_rhalf_replicate.R

source("data-raw/uncertainty_evidence/ue_common.R")
ue <- ue_setup("44_rhalf_replicate")
source("data-raw/uncertainty_evidence/38_audit_common.R")
seed0 <- ue$seed_base + 1000L * ue$prereg$seeds$package_index$C
stage <- Sys.getenv("RH_STAGE", "all")
workers <- as.integer(Sys.getenv("RH_WORKERS", "2"))
n_boot <- as.integer(Sys.getenv("RH_BOOT", "1000"))
root <- file.path(ue_out_root(), "rhalf")
out <- file.path(root, "replication")
pair_dir <- file.path(out, "pairs")
perm_dir <- file.path(out, "perm")
dir.create(pair_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(perm_dir, recursive = TRUE, showWarnings = FALSE)
say <- wpb_logger(file.path(out, "44_rhalf_replicate.log"))
say(
  "start stage=%s workers=%d boot=%d sha=%s",
  stage,
  workers,
  n_boot,
  ue$git_sha
)

series <- list(
  GSE129540 = list(
    control = c("SRR8870780", "SRR8870781", "SRR8870782"),
    knockdown = c("SRR8870783", "SRR8870784", "SRR8870785")
  ),
  GSE75489 = list(
    control = c("SRR2966415", "SRR2966416", "SRR2966417"),
    knockdown = c("SRR2966412", "SRR2966413", "SRR2966414")
  )
)
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
# Quant directories of a halving (random: 42/43; consecutive: quants_half)
quant_dir <- function(halving, gse, run, half) {
  if (halving == "random") {
    file.path(root, "quants", sprintf("%s_%s_%s", gse, run, half))
  } else {
    file.path(
      ue_v50_dir(),
      "human",
      "quants_half",
      sprintf("%s_%s", run, c(A = "h1", B = "h2")[[half]])
    )
  }
}
halvings <- list(random = names(series), consecutive = "GSE129540")

# Annotation and fold map (S0) ---------------------------------------------------------
store <- file.path(ue_out_root(), "store")
fm <- wpb_read(file.path(store, "fold_map.parquet"))
fold_of <- stats::setNames(fm$folded_id, fm$transcript_id)
t2g <- wpb_read(file.path(store, "tx2gene.parquet"))
t2g <- t2g[t2g$on_primary_chrom %in% TRUE, ]
type_class_of <- function(tt) {
  ifelse(
    tt %in%
      c(
        "protein_coding",
        "retained_intron",
        "nonsense_mediated_decay",
        "lncRNA",
        "protein_coding_CDS_not_defined"
      ),
    tt,
    "other"
  )
}

# Load one series of one halving: folded point counts and draws -------------------------
load_series <- function(halving, gse) {
  runs <- c(series[[gse]]$control, series[[gse]]$knockdown)
  dirs <- expand.grid(run = runs, half = c("A", "B"), stringsAsFactors = FALSE)
  dirs$dir <- mapply(quant_dir, halving, gse, dirs$run, dirs$half)
  stopifnot(file.exists(file.path(dirs$dir, "quant.sf")))
  meta <- do.call(
    rbind,
    lapply(seq_len(nrow(dirs)), function(i) {
      m <- jsonlite::read_json(file.path(
        dirs$dir[i],
        "aux_info",
        "meta_info.json"
      ))
      data.frame(
        halving = halving,
        series = gse,
        run = dirs$run[i],
        half = dirs$half[i],
        samp_type = m$samp_type,
        n_draws = m$num_bootstraps,
        num_processed = m$num_processed,
        num_mapped = m$num_mapped,
        percent_mapped = m$percent_mapped,
        salmon_version = m$salmon_version,
        start_time = m$start_time,
        end_time = m$end_time
      )
    })
  )
  stopifnot(meta$samp_type == "gibbs", meta$n_draws == 30L)
  pc <- lapply(dirs$dir, function(d) {
    q <- data.table::fread(
      file.path(d, "quant.sf"),
      select = c("Name", "NumReads")
    )
    f <- fold_of[q$Name]
    f[is.na(f)] <- q$Name[is.na(f)]
    rowsum(q$NumReads, f, reorder = TRUE)
  })
  ids <- rownames(pc[[1]])
  stopifnot(vapply(pc, function(x) identical(rownames(x), ids), TRUE))
  cnt <- do.call(cbind, pc)
  colnames(cnt) <- paste(dirs$run, dirs$half, sep = "_")
  keep <- rowSums(cnt) > 0 & ids %in% t2g$transcript_id
  kid <- ids[keep]
  out <- list(
    meta = meta,
    runs = runs,
    groups = factor(
      rep(c("control", "knockdown"), each = 3),
      c("control", "knockdown")
    )
  )
  for (h in c("A", "B")) {
    cols <- paste(runs, h, sep = "_")
    out$counts[[h]] <- cnt[kid, cols, drop = FALSE]
    colnames(out$counts[[h]]) <- runs
    arr <- array(
      0,
      c(length(kid), length(runs), 30L),
      dimnames = list(kid, runs, NULL)
    )
    mu <- v <- matrix(0, length(ids), length(runs), dimnames = list(ids, runs))
    for (j in seq_along(runs)) {
      b <- read_bootstrap_dir(dirs$dir[dirs$run == runs[j] & dirs$half == h])
      f <- fold_of[rownames(b$counts)]
      f[is.na(f)] <- rownames(b$counts)[is.na(f)]
      x <- rowsum(b$counts, f, reorder = TRUE)
      stopifnot(identical(rownames(x), ids), ncol(x) == 30L)
      mu[, j] <- rowMeans(x)
      v[, j] <- matrixStats::rowVars(x)
      arr[, j, ] <- x[kid, ]
      rm(b, x)
    }
    out$arr[[h]] <- arr
    od <- ue_overdispersion(mu, v, 30L)
    out$od[[h]] <- stats::setNames(od$overdispersion, od$transcript_id)
    out$od_prior[[h]] <- attr(od, "prior")
    rm(arr, mu, v)
  }
  out$ids <- kid
  out$gene_of <- stats::setNames(
    t2g$gene_id[match(kid, t2g$transcript_id)],
    kid
  )
  out$n_iso <- stats::setNames(
    t2g$n_isoforms_gene[match(kid, t2g$transcript_id)],
    kid
  )
  out$type <- stats::setNames(
    type_class_of(t2g$transcript_type[match(kid, t2g$transcript_id)]),
    kid
  )
  out$gene_name <- stats::setNames(
    t2g$gene_name[match(kid, t2g$transcript_id)],
    kid
  )
  invisible(gc(verbose = FALSE))
  out
}

calls_rule <- function(fit, q = 0.10, min_calls = 100L, top = 500L) {
  ok <- !is.na(fit$pvalue)
  hits <- fit$feature_id[ok & fit$qvalue < q]
  if (length(hits) >= min_calls) {
    return(list(ids = hits, fallback = FALSE))
  }
  o <- order(fit$pvalue)
  list(ids = fit$feature_id[o[seq_len(min(top, sum(ok)))]], fallback = TRUE)
}

# One directed analysis: discovery on half D, validation on half V ------------------------
analyse <- function(S, D, V, labels, methods = c("M1", "M2")) {
  mD <- S$counts[[D]]
  mV <- S$counts[[V]]
  fV <- wpb_edger(mV, labels)
  in_b <- labels == "knockdown"
  lib <- colSums(mD)
  res <- list()
  for (m in methods) {
    fD <- if (m == "M1") {
      wpb_edger(mD, labels)
    } else {
      wpb_edger(wpb_divide(mD, S$od[[D]]), labels)
    }
    cr <- calls_rule(fD)
    ids <- cr$ids
    cf <- wpb_confidence(
      S$arr[[D]],
      s = seq_along(S$runs),
      in_b = in_b,
      lib = lib,
      ids = ids,
      gene_of = unname(S$gene_of)
    )
    rv <- cf$mean_inf_rv_all
    terc <- stats::quantile(
      rv[cf$mean_count_all >= 1],
      c(1 / 3, 2 / 3),
      names = FALSE
    )
    fd <- fD[match(ids, fD$feature_id), ]
    fv <- fV[match(ids, fV$feature_id), ]
    tab <- data.frame(
      transcript_id = ids,
      gene_id = unname(S$gene_of[ids]),
      gene_name = unname(S$gene_name[ids]),
      method = m,
      fallback = cr$fallback,
      n_calls_q10 = sum(fD$qvalue < 0.10, na.rm = TRUE),
      d_logFC = fd$logFC,
      d_p = fd$pvalue,
      d_q = fd$qvalue,
      d_absz = abs(fd$z),
      d_mean_count = rowMeans(mD[ids, , drop = FALSE]),
      n_isoforms = unname(S$n_iso[ids]),
      type_class = unname(S$type[ids]),
      mean_inf_rv = cf$table$mean_inf_rv,
      sign_agreement = cf$table$sign_agreement,
      uncertain = cf$table$uncertain,
      ambiguous = cf$table$ambiguous,
      confidence = cf$table$confidence,
      infrv_top_tertile = unname(rv[ids] >= terc[2]),
      v_logFC = fv$logFC,
      v_p = fv$pvalue
    )
    tab$v_p1 <- wpb_one_sided(tab$v_logFC, tab$v_p, sign(tab$d_logFC))
    tab$replicated <- !is.na(tab$v_p1) & tab$v_p1 < 0.05
    res[[m]] <- tab
  }
  res
}

# Stage pairs ------------------------------------------------------------------------------
if (stage %in% c("pairs", "all")) {
  for (hv in names(halvings)) {
    for (gse in halvings[[hv]]) {
      done <- file.exists(file.path(
        perm_dir,
        sprintf("%s_%s_p10_BA.parquet", hv, gse)
      ))
      if (done) {
        next
      }
      runs_s <- c(series[[gse]]$control, series[[gse]]$knockdown)
      need <- unlist(lapply(c("A", "B"), function(h) {
        vapply(runs_s, function(r) quant_dir(hv, gse, r, h), "")
      }))
      if (
        !all(file.exists(file.path(
          need,
          "aux_info",
          "bootstrap",
          "bootstraps.gz"
        )))
      ) {
        say(
          "%s %s: quantification not complete yet, skipped (rerun later)",
          hv,
          gse
        )
        next
      }
      t0 <- Sys.time()
      S <- load_series(hv, gse)
      say(
        "%s %s loaded in %.1f min: %d folded units on primary chromosomes with reads; RTA prior A %.3f, B %.3f; peak RSS %.2f GB",
        hv,
        gse,
        as.numeric(difftime(Sys.time(), t0, units = "mins")),
        length(S$ids),
        S$od_prior$A,
        S$od_prior$B,
        ue_peak_rss_gb()
      )
      ue_write(
        S$meta,
        file.path(out, sprintf("meta_%s_%s.parquet", hv, gse)),
        ue
      )
      for (dd in list(c("A", "B"), c("B", "A"))) {
        dn <- paste0(dd, collapse = "")
        r <- analyse(S, dd[1], dd[2], S$groups)
        for (m in names(r)) {
          x <- r[[m]]
          x$halving <- hv
          x$series <- gse
          x$direction <- dn
          ue_write(
            x,
            file.path(pair_dir, sprintf("%s_%s_%s_%s.parquet", hv, gse, m, dn)),
            ue
          )
          say(
            "%s %s %s %s: %d calls (fallback %s), replicated %.3f, flagged %.1f%%",
            hv,
            gse,
            m,
            dn,
            nrow(x),
            x$fallback[1],
            mean(x$replicated),
            100 * mean(x$confidence != "high", na.rm = TRUE)
          )
        }
      }
      # Permutation baseline: the 10 distinct 3 vs 3 partitions
      cmb <- utils::combn(6L, 3L)
      cmb <- cmb[, cmb[1, ] == 1L, drop = FALSE]
      stopifnot(ncol(cmb) == 10L)
      jobs <- expand.grid(
        k = seq_len(10L),
        dir = c("AB", "BA"),
        stringsAsFactors = FALSE
      )
      invisible(parallel::mclapply(
        seq_len(nrow(jobs)),
        function(i) {
          k <- jobs$k[i]
          dn <- jobs$dir[i]
          lab <- factor(
            ifelse(seq_len(6L) %in% cmb[, k], "control", "knockdown"),
            c("control", "knockdown")
          )
          n_swapped <- sum(lab != S$groups)
          r <- tryCatch(
            analyse(
              S,
              substr(dn, 1, 1),
              substr(dn, 2, 2),
              lab,
              methods = "M1"
            )$M1,
            error = function(e) e
          )
          if (inherits(r, "error")) {
            say(
              "perm %s %s %d %s ERROR %s",
              hv,
              gse,
              k,
              dn,
              conditionMessage(r)
            )
            return(NULL)
          }
          r$halving <- hv
          r$series <- gse
          r$direction <- dn
          r$partition <- k
          r$control_set <- paste(S$runs[cmb[, k]], collapse = ";")
          r$true_partition <- n_swapped == 0L
          r$n_swapped <- n_swapped
          ue_write(
            r,
            file.path(
              perm_dir,
              sprintf("%s_%s_p%02d_%s.parquet", hv, gse, k, dn)
            ),
            ue
          )
          NULL
        },
        mc.cores = workers,
        mc.preschedule = FALSE
      ))
      say(
        "%s %s: permutations done; %.1f min in all",
        hv,
        gse,
        as.numeric(difftime(Sys.time(), t0, units = "mins"))
      )
      rm(S)
      invisible(gc(verbose = FALSE))
    }
  }
}

# Stage summarise -----------------------------------------------------------------------------
if (stage %in% c("summarise", "all")) {
  rd_all <- function(dir, pattern) {
    f <- list.files(dir, pattern = pattern, full.names = TRUE)
    do.call(rbind, lapply(f, wpb_read))
  }
  P_all <- rd_all(pair_dir, "\\.parquet$")
  PM_all <- rd_all(perm_dir, "\\.parquet$")
  P_all <- P_all[!is.na(P_all$confidence), ]
  prep <- function(P) {
    P$flag_primary <- P$confidence != "high"
    P$flag_unstable <- P$sign_agreement < 0.95
    P$flag_uncertain <- P$uncertain
    P$flag_infrv_tertile <- P$infrv_top_tertile
    P$flag_ambiguous <- P$ambiguous
    # Tested in half V: half V applies its own filterByExpr(), so a call it
    # filters out has no validation p-value and counts as not replicated
    P$testable_v <- !is.na(P$v_p)
    P$log2_count <- log2(P$d_mean_count + 1)
    P$log2_iso <- log2(pmax(P$n_isoforms, 1))
    P$type_class <- factor(P$type_class)
    P$cell <- interaction(
      dplyr::ntile(P$log2_count, 10),
      dplyr::ntile(P$d_absz, 10),
      drop = TRUE
    )
    # C4: per direction as many calls as the primary flag, the lowest counts
    P$flag_c4 <- FALSE
    for (ix in split(seq_len(nrow(P)), P$direction)) {
      k <- sum(P$flag_primary[ix])
      if (k > 0) P$flag_c4[ix[order(P$d_mean_count[ix])[seq_len(k)]]] <- TRUE
    }
    P
  }
  cov_formula <- "log2_count + d_absz + log2_iso + type_class"
  fit_or <- function(d, flag) {
    d <- d[!is.na(d[[flag]]), ]
    if (length(unique(d[[flag]])) < 2L) {
      return(NA_real_)
    }
    f <- stats::as.formula(paste("replicated ~", flag, "+", cov_formula))
    if (length(unique(d$type_class)) < 2L) {
      f <- stats::as.formula(paste(
        "replicated ~",
        flag,
        "+ log2_count + d_absz + log2_iso"
      ))
    }
    co <- stats::coef(suppressWarnings(stats::glm(
      f,
      data = d,
      family = stats::binomial()
    )))
    unname(exp(co[paste0(flag, "TRUE")]))
  }
  deficit <- function(d, flag) {
    f <- d[[flag]]
    mean(d$replicated[f %in% FALSE]) - mean(d$replicated[f %in% TRUE])
  }
  matched_diff <- function(d, flag) {
    f <- d[[flag]]
    num <- 0
    w <- 0
    for (ix in split(seq_len(nrow(d)), d$cell)) {
      a <- ix[f[ix] %in% TRUE]
      b <- ix[f[ix] %in% FALSE]
      if (length(a) && length(b)) {
        num <- num + length(a) * (mean(d$replicated[b]) - mean(d$replicated[a]))
        w <- w + length(a)
      }
    }
    num / w
  }
  auc_gain <- function(d, flag) {
    d <- d[!is.na(d[[flag]]), ]
    y <- !d$replicated
    m0 <- suppressWarnings(stats::glm(
      y ~ log2_count + d_absz,
      data = d,
      family = stats::binomial()
    ))
    m1 <- suppressWarnings(stats::glm(
      stats::as.formula(paste("y ~ log2_count + d_absz +", flag)),
      data = d,
      family = stats::binomial()
    ))
    wpb_auc(stats::fitted(m1), y) - wpb_auc(stats::fitted(m0), y)
  }
  # Gene bootstrap: all rows of a resampled gene (both directions) together
  boot_index <- function(P, B, seed) {
    set.seed(seed)
    by_gene <- split(seq_len(nrow(P)), P$gene_id)
    lapply(seq_len(B), function(b) {
      unlist(
        by_gene[sample.int(length(by_gene), replace = TRUE)],
        use.names = FALSE
      )
    })
  }
  par_boot <- function(P, idx, fun) {
    unlist(parallel::mclapply(
      idx,
      function(ix) tryCatch(fun(P[ix, ]), error = function(e) NA_real_),
      mc.cores = workers
    ))
  }
  flags <- c(
    "flag_primary",
    "flag_unstable",
    "flag_uncertain",
    "flag_infrv_tertile",
    "flag_ambiguous",
    "flag_c4"
  )
  rows <- list()
  c4_rows <- list()
  pl_rows <- list()
  sets <- unique(P_all[, c("halving", "series", "method")])
  rownames(sets) <- NULL
  for (i in seq_len(nrow(sets))) {
    st <- sets[i, ]
    P <- prep(P_all[
      P_all$halving == st$halving &
        P_all$series == st$series &
        P_all$method == st$method,
    ])
    idx <- boot_index(P, n_boot, seed0 + 600L)
    for (fl in flags) {
      f <- P[[fl]]
      bt_or <- par_boot(P, idx, function(x) fit_or(x, fl))
      bt_def <- par_boot(P, idx, function(x) deficit(x, fl))
      bt_md <- par_boot(P, idx, function(x) matched_diff(x, fl))
      bt_auc <- par_boot(P, idx, function(x) auc_gain(x, fl))
      rf <- par_boot(P, idx, function(x) mean(x$replicated[x[[fl]] %in% TRUE]))
      ru <- par_boot(P, idx, function(x) mean(x$replicated[x[[fl]] %in% FALSE]))
      rft <- par_boot(P, idx, function(x) {
        mean(x$replicated[x[[fl]] %in% TRUE & x$testable_v])
      })
      rut <- par_boot(P, idx, function(x) {
        mean(x$replicated[x[[fl]] %in% FALSE & x$testable_v])
      })
      dt <- par_boot(P, idx, function(x) deficit(x[x$testable_v, ], fl))
      Pt <- P[P$testable_v, ]
      rows[[length(rows) + 1L]] <- data.frame(
        st,
        flag = fl,
        n_calls = nrow(P),
        n_genes = length(unique(P$gene_id)),
        n_flagged = sum(f, na.rm = TRUE),
        pct_flagged = 100 * mean(f, na.rm = TRUE),
        replication_overall = mean(P$replicated),
        rate_flagged = mean(P$replicated[f %in% TRUE]),
        rate_flagged_ci = wpb_ci(rf),
        rate_unflagged = mean(P$replicated[f %in% FALSE]),
        rate_unflagged_ci = wpb_ci(ru),
        deficit = deficit(P, fl),
        deficit_ci = wpb_ci(bt_def),
        matched_deficit = matched_diff(P, fl),
        matched_deficit_ci = wpb_ci(bt_md),
        adjusted_or = fit_or(P, fl),
        adjusted_or_ci = wpb_ci(bt_or),
        delta_auc = auc_gain(P, fl),
        delta_auc_ci = wpb_ci(bt_auc),
        n_untestable_v = sum(!P$testable_v),
        n_untestable_v_flagged = sum(!P$testable_v & f %in% TRUE),
        rate_flagged_testable = mean(Pt$replicated[Pt[[fl]] %in% TRUE]),
        rate_flagged_testable_ci = wpb_ci(rft),
        rate_unflagged_testable = mean(Pt$replicated[Pt[[fl]] %in% FALSE]),
        rate_unflagged_testable_ci = wpb_ci(rut),
        deficit_testable = deficit(Pt, fl),
        deficit_testable_ci = wpb_ci(dt),
        per_direction_deficit = paste(
          vapply(
            split(P, P$direction),
            function(x) sprintf("%s %.3f", x$direction[1], deficit(x, fl)),
            ""
          ),
          collapse = "; "
        ),
        n_boot = n_boot
      )
    }
    dd <- function(x) deficit(x, "flag_primary") - deficit(x, "flag_c4")
    c4_rows[[i]] <- data.frame(
      st,
      comparison = "deficit(primary) - deficit(count-matched C4)",
      estimate = dd(P),
      ci = wpb_ci(par_boot(P, idx, dd)),
      deficit_primary = deficit(P, "flag_primary"),
      deficit_c4 = deficit(P, "flag_c4"),
      # Overlap of the two flags: a difference near 0 can come from a ceiling
      # (almost every unflagged call replicates under both flags), not from
      # identical sets
      n_unflagged_primary = sum(!P$flag_primary),
      n_unflagged_c4 = sum(!P$flag_c4),
      n_unflagged_both = sum(!P$flag_primary & !P$flag_c4),
      n_unflagged_primary_replicated = sum(P$replicated[!P$flag_primary]),
      n_unflagged_c4_replicated = sum(P$replicated[!P$flag_c4]),
      share_calls_same_flag = mean(P$flag_primary == P$flag_c4)
    )
    set.seed(seed0 + 601L)
    by_cell <- split(seq_len(nrow(P)), P$cell)
    placebo <- t(vapply(
      seq_len(200L),
      function(b) {
        f <- P$flag_primary
        for (ix in by_cell) {
          if (length(ix) > 1L) f[ix] <- f[ix][sample.int(length(ix))]
        }
        d <- P
        d$flag_placebo <- f
        c(deficit(d, "flag_placebo"), fit_or(d, "flag_placebo"))
      },
      numeric(2)
    ))
    obs <- c(deficit(P, "flag_primary"), fit_or(P, "flag_primary"))
    pl_rows[[i]] <- data.frame(
      st,
      statistic = c("deficit", "adjusted_or"),
      observed = obs,
      placebo_mean = colMeans(placebo, na.rm = TRUE),
      placebo_2.5 = apply(placebo, 2, wpb_q, 0.025),
      placebo_97.5 = apply(placebo, 2, wpb_q, 0.975),
      empirical_p_one_sided = c(
        (1 + sum(placebo[, 1] >= obs[1], na.rm = TRUE)) / 201,
        (1 + sum(placebo[, 2] <= obs[2], na.rm = TRUE)) / 201
      )
    )
    say(
      "%s %s %s: %d calls, primary flag %.1f%%, replication flagged %.3f vs unflagged %.3f, OR %.3f",
      st$halving,
      st$series,
      st$method,
      nrow(P),
      100 * mean(P$flag_primary),
      mean(P$replicated[P$flag_primary]),
      mean(P$replicated[!P$flag_primary]),
      fit_or(P, "flag_primary")
    )
  }
  flag_tab <- do.call(rbind, rows)
  ue_write(flag_tab, file.path(out, "summary_flags.parquet"), ue, seed0 + 600L)
  utils::write.table(
    flag_tab,
    file.path(out, "summary_flags.tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  ue_write(
    do.call(rbind, c4_rows),
    file.path(out, "summary_c4.parquet"),
    ue,
    seed0 + 600L
  )
  ue_write(
    do.call(rbind, pl_rows),
    file.path(out, "summary_placebo.parquet"),
    ue,
    seed0 + 601L
  )

  # Permutation baseline ------------------------------------------------------------------
  PM_all <- PM_all[!is.na(PM_all$confidence), ]
  PM_all$flag_primary <- PM_all$confidence != "high"
  pm <- do.call(
    rbind,
    lapply(
      split(
        PM_all,
        list(PM_all$halving, PM_all$series, PM_all$partition, PM_all$direction),
        drop = TRUE
      ),
      function(x) {
        data.frame(
          halving = x$halving[1],
          series = x$series[1],
          partition = x$partition[1],
          direction = x$direction[1],
          control_set = x$control_set[1],
          true_partition = x$true_partition[1],
          n_swapped = x$n_swapped[1],
          n_calls = nrow(x),
          fallback = x$fallback[1],
          n_calls_q10 = x$n_calls_q10[1],
          replication = mean(x$replicated),
          pct_flagged = 100 * mean(x$flag_primary),
          rate_flagged = mean(x$replicated[x$flag_primary]),
          rate_unflagged = mean(x$replicated[!x$flag_primary]),
          deficit = mean(x$replicated[!x$flag_primary]) -
            mean(x$replicated[x$flag_primary]),
          n_untestable_v = sum(is.na(x$v_p)),
          replication_testable = mean(x$replicated[!is.na(x$v_p)])
        )
      }
    )
  )
  ue_write(pm, file.path(out, "summary_perm.parquet"), ue)
  utils::write.table(
    pm,
    file.path(out, "summary_perm.tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  # Delta-AUC transfer between series (random halves, M1) -----------------------------------
  tr_rows <- list()
  pr <- lapply(names(series), function(g) {
    x <- P_all[
      P_all$halving == "random" & P_all$series == g & P_all$method == "M1",
    ]
    if (nrow(x)) prep(x) else x
  })
  names(pr) <- names(series)
  if (all(vapply(pr, nrow, 1L) > 0)) {
    for (k in 1:2) {
      tr_s <- names(series)[k]
      te_s <- names(series)[3 - k]
      tr <- pr[[tr_s]]
      te <- pr[[te_s]]
      tr$y <- !tr$replicated
      m0 <- suppressWarnings(stats::glm(
        y ~ log2_count + d_absz,
        data = tr,
        family = stats::binomial()
      ))
      m1 <- suppressWarnings(stats::glm(
        y ~ log2_count + d_absz + flag_primary,
        data = tr,
        family = stats::binomial()
      ))
      dauc <- function(d) {
        y <- !d$replicated
        wpb_auc(stats::predict(m1, d), y) - wpb_auc(stats::predict(m0, d), y)
      }
      idx <- boot_index(te, n_boot, seed0 + 602L)
      bt <- par_boot(te, idx, dauc)
      tr_rows[[k]] <- data.frame(
        train = tr_s,
        test = te_s,
        delta_auc_transfer = dauc(te),
        ci = wpb_ci(bt),
        auc_base_test = wpb_auc(stats::predict(m0, te), !te$replicated),
        auc_flag_test = wpb_auc(stats::predict(m1, te), !te$replicated),
        delta_auc_within_test = auc_gain(te, "flag_primary"),
        flag_coef_train = unname(stats::coef(m1)["flag_primaryTRUE"])
      )
    }
    tr_tab <- do.call(rbind, tr_rows)
    ue_write(
      tr_tab,
      file.path(out, "summary_transfer.parquet"),
      ue,
      seed0 + 602L
    )
  }

  # Known ESRP targets (descriptive only) ----------------------------------------------------
  kt <- P_all[
    P_all$gene_name %in% known_targets,
    c(
      "halving",
      "series",
      "method",
      "direction",
      "gene_name",
      "transcript_id",
      "d_logFC",
      "d_q",
      "d_mean_count",
      "confidence",
      "mean_inf_rv",
      "sign_agreement",
      "ambiguous",
      "v_logFC",
      "v_p",
      "v_p1",
      "replicated"
    )
  ]
  ue_write(kt, file.path(out, "known_targets.parquet"), ue)
  kts <- do.call(
    rbind,
    lapply(
      split(kt, list(kt$halving, kt$series, kt$method), drop = TRUE),
      function(x) {
        data.frame(
          halving = x$halving[1],
          series = x$series[1],
          method = x$method[1],
          n_calls = nrow(x),
          n_genes = length(unique(x$gene_name)),
          pct_flagged = 100 * mean(x$confidence != "high"),
          replication = mean(x$replicated),
          n_untestable_v = sum(is.na(x$v_p))
        )
      }
    )
  )
  ue_write(kts, file.path(out, "known_targets_summary.parquet"), ue)

  # Info: Salmon runs and split check ---------------------------------------------------------
  meta <- do.call(
    rbind,
    lapply(list.files(out, "^meta_.*\\.parquet$", full.names = TRUE), wpb_read)
  )
  qt <- utils::read.delim(file.path(root, "quant_times.tsv"))
  sc <- utils::read.delim(file.path(root, "split_check.tsv"))
  ue_write(meta, file.path(out, "summary_meta.parquet"), ue)
  ue_write(qt, file.path(out, "quant_times.parquet"), ue)
  ue_write(sc, file.path(out, "split_check.parquet"), ue)
  print(flag_tab[
    flag_tab$flag %in% c("flag_primary", "flag_c4"),
    c(
      "halving",
      "series",
      "method",
      "flag",
      "n_calls",
      "pct_flagged",
      "rate_flagged",
      "rate_unflagged",
      "deficit",
      "deficit_ci",
      "adjusted_or",
      "adjusted_or_ci"
    )
  ])
  say("summaries written")
}
ue_write_session(out, ue)
say("done; peak RSS %.2f GB", ue_peak_rss_gb())
