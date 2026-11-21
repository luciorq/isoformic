# Shared helpers of the WP-A scripts (genome-wide placenta null); sourced,
# not run.
#
# Purpose: paths of WP-A, the loader of the S0 store (10_placenta_prepare.R)
# into matrices, and the matrix-level implementations of the DTE methods that
# the per-split loops run (no IsoformicExperiment subsetting inside loops).
# Sourced by 20_null_splits.R, 21_null_dte.R, 22_null_dtu.R, 23_thin_dte.R
# and 29_null_summarise.R after ue_common.R (S0), which loads the frozen
# package snapshot and provides ue_write() with the provenance columns.
#
# Method numbering is the one of the WP-A work package (it differs from the
# `methods` list of prereg.yaml, see the mapping in 21_null_dte.R):
#   M1 edgeR QL raw, M2 edgeR QL divided, M3/M4 limma-voom raw/divided,
#   M5 DESeq2 Wald, M6 swish on the draws, M7 swish with every draw set to
#   the point counts, M8 Wilcoxon on TMM log CPM, M9 de_confidence over M1,
#   M10 M1 after removing as many lowest-expressed transcripts as M9 flags,
#   M11 gene level edgeR QL raw/divided.
#
# Env vars: ISOFORMIC_UE_OUT (ue_common.R); ISOFORMIC_V50_DIR is not read.

source("data-raw/uncertainty_evidence/ue_common.R")

wpa_paths <- function() {
  root <- ue_out_root()
  out <- file.path(root, "null_placenta")
  list(
    root = root,
    store = file.path(root, "store"),
    s0 = file.path(root, "S0"),
    out = out,
    cache = file.path(out, "cache"),
    splits = file.path(out, "splits"),
    dte = file.path(out, "dte"),
    dtu = file.path(out, "dtu"),
    thin = file.path(out, "thin"),
    summary = file.path(out, "summary"),
    logs = file.path(out, "logs")
  )
}

wpa_log <- function(...) {
  cat(format(Sys.time(), "%H:%M:%S"), " ", ..., "\n", sep = "")
}

wpa_read <- function(path, columns = NULL) {
  df <- as.data.frame(read_parquet(path, columns = columns))
  df[,
    setdiff(colnames(df), c("git_sha", "snapshot_sha256", "seed", "script")),
    drop = FALSE
  ]
}

#' DuckDB query on Parquet files of the store, returned as a data frame
wpa_query <- function(sql) {
  con <- DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  DBI::dbExecute(con, sprintf("SET threads = %d", ue_threads()))
  DBI::dbGetQuery(con, sql)
}

#' Long table (id, sample, value) to a matrix
long_to_matrix <- function(id, sample, value, ids = NULL, samples = NULL) {
  ids <- ids %||% sort(unique(id))
  samples <- samples %||% sort(unique(sample))
  m <- matrix(0, length(ids), length(samples), dimnames = list(ids, samples))
  i <- match(id, ids)
  j <- match(sample, samples)
  ok <- !is.na(i) & !is.na(j)
  m[cbind(i[ok], j[ok])] <- value[ok]
  m
}

# Statistics shared by the methods ----------------------------------------------

#' Row-wise Wilcoxon rank-sum p-values (two-sided) as stats::wilcox.test():
#' exact without ties (both groups < 50), normal approximation with
#' continuity and tie correction otherwise
row_wilcox <- function(x, in_b) {
  na <- sum(!in_b)
  nb <- sum(in_b)
  r <- matrixStats::rowRanks(x, ties.method = "average")
  w <- rowSums(r[, in_b, drop = FALSE]) - nb * (nb + 1) / 2
  rs <- t(apply(r, 1L, sort))
  ties <- if (ncol(rs) > 1L) {
    rowSums(rs[, -1L, drop = FALSE] == rs[, -ncol(rs), drop = FALSE]) > 0L
  } else {
    rep(FALSE, nrow(x))
  }
  p <- numeric(nrow(x))
  ex <- !ties
  if (any(ex)) {
    q <- w[ex]
    # wilcox.test: p = 2 * min(P(W <= q), P(W >= q)), capped at 1
    lo <- stats::pwilcox(q, nb, na)
    hi <- stats::pwilcox(q - 1, nb, na, lower.tail = FALSE)
    p[ex] <- pmin(1, 2 * pmin(lo, hi))
  }
  if (any(ties)) {
    n <- na + nb
    rr <- r[ties, , drop = FALSE]
    tie_term <- apply(rr, 1L, function(v) {
      t <- table(v)
      sum(t^3 - t)
    })
    sigma <- sqrt((na * nb / 12) * ((n + 1) - tie_term / (n * (n - 1))))
    z <- w[ties] - na * nb / 2
    z <- (z - sign(z) * 0.5) / sigma
    pp <- 2 * pmin(stats::pnorm(z), stats::pnorm(z, lower.tail = FALSE))
    pp[!is.finite(pp)] <- 1
    p[ties] <- pmin(pp, 1)
  }
  p
}

two_group_design <- function(in_b) {
  stats::model.matrix(~ factor(in_b, levels = c(FALSE, TRUE)))
}

#' edgeR v4 quasi-likelihood test (Baldoni et al. 2024 route: DGEList,
#' normLibSizes, glmQLFit, glmQLFTest) on the given rows
edger_ql <- function(counts, in_b) {
  y <- edgeR::normLibSizes(edgeR::DGEList(counts))
  fit <- edgeR::glmQLFit(y, two_group_design(in_b))
  tt <- edgeR::glmQLFTest(fit, coef = 2L)$table
  data.frame(id = rownames(counts), log2FC = tt$logFC, pvalue = tt$PValue)
}

limma_voom <- function(counts, in_b) {
  y <- edgeR::normLibSizes(edgeR::DGEList(counts))
  design <- two_group_design(in_b)
  v <- limma::voom(y, design)
  fit <- limma::eBayes(limma::lmFit(v, design))
  data.frame(
    id = rownames(counts),
    log2FC = fit$coefficients[, 2L],
    pvalue = fit$p.value[, 2L]
  )
}

wilcox_logcpm <- function(counts, in_b) {
  y <- edgeR::normLibSizes(edgeR::DGEList(counts))
  lc <- edgeR::cpm(y, log = TRUE)
  data.frame(
    id = rownames(counts),
    log2FC = rowMeans(lc[, in_b, drop = FALSE]) -
      rowMeans(lc[, !in_b, drop = FALSE]),
    pvalue = row_wilcox(lc, in_b)
  )
}

deseq2_wald <- function(counts, in_b) {
  cd <- data.frame(cond = factor(ifelse(in_b, "b", "a")))
  dds <- DESeq2::DESeqDataSetFromMatrix(round(counts), cd, ~cond)
  dds <- suppressMessages(DESeq2::DESeq(dds, quiet = TRUE))
  res <- DESeq2::results(dds, independentFiltering = FALSE)
  data.frame(
    id = rownames(counts),
    log2FC = res$log2FoldChange,
    pvalue = res$pvalue
  )
}

#' fishpond swish (scaleInfReps, labelKeep defaults, swish nperms) on a
#' transcript x sample x draw array; returns swish's own q-value as well
swish_test <- function(counts, arr, lengths, in_b, seed, nperms = 100L) {
  reps <- lapply(seq_len(dim(arr)[3L]), function(k) arr[,, k])
  names(reps) <- paste0("infRep", seq_along(reps))
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = c(list(counts = counts, length = lengths), reps),
    colData = S4Vectors::DataFrame(cond = factor(ifelse(in_b, "b", "a")))
  )
  se <- suppressMessages(fishpond::scaleInfReps(se, quiet = TRUE))
  se <- fishpond::labelKeep(se)
  se <- se[S4Vectors::mcols(se)$keep, ]
  set.seed(seed)
  se <- suppressMessages(fishpond::swish(
    se,
    x = "cond",
    nperms = nperms,
    quiet = TRUE
  ))
  mc <- S4Vectors::mcols(se)
  data.frame(
    id = rownames(se),
    log2FC = mc$log2FC,
    pvalue = mc$pvalue,
    swish_qvalue = mc$qvalue
  )
}

# Delete-one-unit jackknife ------------------------------------------------------

#' Jackknife over units of the mean of a per-split statistic. `value`: the
#' per-split values; `units_in`: list of the units each split uses; `units`:
#' all units. The leave-u-out estimate is the mean over the splits that do
#' not use unit u (for the 7v7 leave-one-unit-out design: the 14 splits that
#' leave u out). Monte Carlo component: each leave-u-out mean carries the
#' split-sampling variance s2_u / k_u; its expected contribution to the raw
#' jackknife variance, (n - 1)^2 / n * mean(s2_u / k_u), is subtracted.
jackknife_units <- function(value, units_in, units, stat = mean) {
  n <- length(units)
  theta <- stat(value)
  loo <- numeric(n)
  s2k <- numeric(n)
  for (j in seq_len(n)) {
    use <- !vapply(units_in, function(u) units[j] %in% u, logical(1L))
    v <- value[use]
    loo[j] <- stat(v)
    s2k[j] <- if (length(v) > 1L) stats::var(v) / length(v) else NA_real_
  }
  var_raw <- (n - 1) / n * sum((loo - mean(loo))^2)
  mc <- (n - 1)^2 / n * mean(s2k, na.rm = TRUE)
  var_cor <- max(var_raw - mc, 0)
  z <- stats::qnorm(0.975)
  data.frame(
    estimate = theta,
    se_raw = sqrt(var_raw),
    lo_raw = theta - z * sqrt(var_raw),
    hi_raw = theta + z * sqrt(var_raw),
    mc_var = mc,
    se_corrected = sqrt(var_cor),
    lo_corrected = theta - z * sqrt(var_cor),
    hi_corrected = theta + z * sqrt(var_cor),
    n_units = n,
    n_splits = length(value)
  )
}

boot_splits <- function(value, B = 2000L, seed = 1L, stat = mean) {
  set.seed(seed)
  b <- replicate(B, stat(sample(value, replace = TRUE)))
  q <- stats::quantile(b, c(0.025, 0.975), names = FALSE, na.rm = TRUE)
  data.frame(boot_lo = q[[1L]], boot_hi = q[[2L]])
}

wilson <- function(x, n) {
  if (n == 0) {
    return(c(NA_real_, NA_real_))
  }
  z <- stats::qnorm(0.975)
  p <- x / n
  den <- 1 + z^2 / n
  c0 <- (p + z^2 / (2 * n)) / den
  h <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / den
  c(c0 - h, c0 + h)
}
