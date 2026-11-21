# WP-B, post hoc checks (amendment A11): not pre-registered, added on review.
#
# Purpose: re-estimate the intervals of the held-out replication statistics
# under resampling schemes that respect the recurrence of genes across the
# directed pairs, and describe the pair-to-pair variation behind the pooled
# rates. No rule and no pre-registered number changes; R_B stays judged on
# the pre-registered two-level bootstrap of 31b_replication_summary.R.
#
# Why: 31b resamples the 30 half-splits, then the genes of each directed
# pair independently, so one gene gets unrelated bootstrap weights in each
# of the 60 pairs. All pairs are built from the same 23 samples, and the
# 167,318 calls come from about 20,000 transcripts in about 8,800 genes, so
# a gene's outcomes are correlated across pairs. Two schemes here:
#   gene:  genes resampled with replacement once, over all pairs (every
#          call of a gene, in every pair, moves together);
#   split_x_gene: the split weights of a resample of the 30 half-splits
#          multiplied by gene weights as above (crossed).
# Statistics (flags fixed on the full data, as in 31b): deficit and
# adjusted OR (the covariates of 31b: log2 count, |z|, log2 isoforms, type
# class) of flag_unstable and flag_primary, the adjusted OR of
# flag_unstable for sign concordance, and deficit(primary) - deficit(C4)
# (R_B's second criterion). Resampling uses frequency weights, which give
# the same estimates as stacking the resampled rows.
#
# Also, point estimates only:
#   - the adjusted ORs with directed-pair fixed effects (the 31b covariates
#     plus factor(pair_id));
#   - per directed pair: calls, replication rate, and the rates of flagged
#     and unflagged calls for flag_unstable; their ranges, the number of
#     pairs in which flagged calls replicate less, the equal-weight mean of
#     the pair rates, and the correlation of calls per pair with the pair
#     replication rate;
#   - recurrence: distinct transcripts and genes among all and flagged calls;
#   - read halves: calls at q < 0.10 and whether the top-500 fallback of
#     44_rhalf_replicate.R applied, for each series and direction (M1).
#
# Inputs ($ISOFORMIC_UE_OUT/): placenta_audit/replication/pooled_dte.parquet
#   (31b_replication_summary.R), rhalf/replication/pairs/*_M1_*.parquet
#   (44_rhalf_replicate.R).
# Outputs: $ISOFORMIC_UE_OUT/placenta_audit/post_hoc/
#   cluster_bootstrap.parquet  key, estimate, ci_lo, ci_hi, ci_type, n
#   cluster_bootstrap_draws.parquet  one row per resample and scheme
# Env vars: ISOFORMIC_UE_OUT, ISOFORMIC_UE_PKG (ue_common.R),
#   ISOFORMIC_UE_THREADS (workers, default 4), ISOFORMIC_UE_BOOT (resamples
#   per scheme, default 1,000).
# Resources (estimate): about 15 min with 4 workers, < 4 GB. Keep BLAS
# single-threaded (OPENBLAS_NUM_THREADS=1): the workers are forked, and a
# multithreaded BLAS in each makes the fits many times slower.
# Run from the repository root, inside a memory cap of about 6 GB (see
# README.md):
#   OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
#     pixi run --frozen Rscript data-raw/uncertainty_evidence/34_cluster_bootstrap.R

source("data-raw/uncertainty_evidence/ue_common.R")
ue <- ue_setup("34_cluster_bootstrap")
source("data-raw/uncertainty_evidence/38_audit_common.R")
seed <- wpb_seed(ue, 800L)
set.seed(seed)
out_dir <- wpb_dir("post_hoc")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
workers <- ue_threads()
n_boot <- as.integer(Sys.getenv("ISOFORMIC_UE_BOOT", "1000"))

rows <- list()
add <- function(
  key,
  estimate,
  ci_lo = NA_real_,
  ci_hi = NA_real_,
  ci_type = NA_character_,
  n = NA_real_
) {
  # jarl-ignore undesirable_operator: appends to a result list of the enclosing scope
  rows[[length(rows) + 1L]] <<- data.frame(
    key = key,
    estimate = as.numeric(estimate),
    ci_lo = as.numeric(ci_lo),
    ci_hi = as.numeric(ci_hi),
    ci_type = as.character(ci_type),
    n = as.numeric(n)
  )
}

P <- wpb_read(wpb_dir("replication", "pooled_dte.parquet"))
P$type_class <- factor(P$type_class)
y_rep <- as.numeric(P$replicated)
y_sc <- as.numeric(P$sign_concordant)

# Design matrices of the 31b model, built once per flag
cov_rhs <- "log2_count + d_absz + log2_iso + type_class"
X_of <- function(flag) {
  stats::model.matrix(stats::as.formula(paste("~", flag, "+", cov_rhs)), P)
}
X_unstable <- X_of("flag_unstable")
X_primary <- X_of("flag_primary")
wor <- function(X, y, w, flag) {
  fit <- suppressWarnings(stats::glm.fit(
    X,
    y,
    weights = w,
    family = stats::binomial()
  ))
  exp(fit$coefficients[[paste0(flag, "TRUE")]])
}
wdef <- function(f, w) {
  stats::weighted.mean(y_rep[!f], w[!f]) - stats::weighted.mean(y_rep[f], w[f])
}
stats_w <- function(w) {
  c(
    unstable_deficit = wdef(P$flag_unstable, w),
    unstable_or = wor(X_unstable, y_rep, w, "flag_unstable"),
    unstable_or_signconc = wor(X_unstable, y_sc, w, "flag_unstable"),
    primary_deficit = wdef(P$flag_primary, w),
    primary_or = wor(X_primary, y_rep, w, "flag_primary"),
    primary_minus_c4 = wdef(P$flag_primary, w) - wdef(P$flag_c4, w)
  )
}
est <- stats_w(rep(1, nrow(P)))

genes <- unique(P$gene_id)
gi <- match(P$gene_id, genes)
splits <- sort(unique(P$split))
si <- match(P$split, splits)
gene_w <- function() {
  tabulate(sample.int(length(genes), replace = TRUE), length(genes))
}
split_w <- function() {
  tabulate(sample.int(length(splits), replace = TRUE), length(splits))
}
schemes <- list(
  gene = function() gene_w()[gi],
  split_x_gene = function() split_w()[si] * gene_w()[gi]
)
draws <- list()
for (s in names(schemes)) {
  set.seed(seed + match(s, names(schemes)))
  # the weights are drawn in the parent process, so the result does not
  # depend on the number of workers
  W <- lapply(seq_len(n_boot), function(b) schemes[[s]]())
  t0 <- Sys.time()
  bt <- do.call(
    rbind,
    parallel::mclapply(
      W,
      function(w) {
        tryCatch(stats_w(w), error = function(e) rep(NA_real_, length(est)))
      },
      mc.cores = workers
    )
  )
  colnames(bt) <- names(est)
  draws[[s]] <- data.frame(scheme = s, b = seq_len(nrow(bt)), bt)
  for (k in names(est)) {
    q <- stats::quantile(bt[, k], c(0.025, 0.975), na.rm = TRUE, names = FALSE)
    add(
      sprintf("cluster.%s.%s", s, k),
      est[[k]],
      q[1],
      q[2],
      sprintf(
        "percentile bootstrap, %s resampling, %d resamples",
        s,
        sum(!is.na(bt[, k]))
      ),
      nrow(P)
    )
  }
  message(sprintf(
    "34: %s done (%d resamples, %.1f min)",
    s,
    nrow(bt),
    as.numeric(difftime(Sys.time(), t0, units = "mins"))
  ))
}

# Pair fixed effects ------------------------------------------------------------------
fe_or <- function(y, flag) {
  fit <- suppressWarnings(stats::glm(
    stats::as.formula(paste(y, "~", flag, "+", cov_rhs, "+ factor(pair_id)")),
    data = P,
    family = stats::binomial()
  ))
  exp(stats::coef(fit)[[paste0(flag, "TRUE")]])
}
for (flag in c("flag_unstable", "flag_primary")) {
  for (y in c("replicated", "sign_concordant")) {
    add(
      sprintf("pair_fe.%s.%s.adjusted_or", flag, y),
      fe_or(y, flag),
      ci_type = "point estimate; covariates of 31b plus directed-pair fixed effects",
      n = nrow(P)
    )
  }
}

# Per directed pair ----------------------------------------------------------------------
pp <- do.call(
  rbind,
  lapply(split(P, P$pair_id), function(d) {
    f <- d$flag_unstable
    data.frame(
      pair_id = d$pair_id[1],
      n_calls = nrow(d),
      rate = mean(d$replicated),
      rate_flagged = mean(d$replicated[f]),
      rate_unflagged = mean(d$replicated[!f])
    )
  })
)
n_pairs <- nrow(pp)
pt <- "per directed pair, point (range over pairs, not a CI)"
add("pairs.n", n_pairs, n = n_pairs)
add("pairs.calls_min", min(pp$n_calls), ci_type = pt, n = n_pairs)
add("pairs.calls_max", max(pp$n_calls), ci_type = pt, n = n_pairs)
add("pairs.calls_median", stats::median(pp$n_calls), ci_type = pt, n = n_pairs)
add("pairs.rate_min", min(pp$rate), ci_type = pt, n = n_pairs)
add("pairs.rate_max", max(pp$rate), ci_type = pt, n = n_pairs)
add(
  "pairs.cor_calls_rate",
  stats::cor(pp$n_calls, pp$rate),
  ci_type = "Pearson correlation over directed pairs",
  n = n_pairs
)
add(
  "pairs.unstable.rate_flagged_min",
  min(pp$rate_flagged),
  ci_type = pt,
  n = n_pairs
)
add(
  "pairs.unstable.rate_flagged_max",
  max(pp$rate_flagged),
  ci_type = pt,
  n = n_pairs
)
add(
  "pairs.unstable.rate_unflagged_min",
  min(pp$rate_unflagged),
  ci_type = pt,
  n = n_pairs
)
add(
  "pairs.unstable.rate_unflagged_max",
  max(pp$rate_unflagged),
  ci_type = pt,
  n = n_pairs
)
add(
  "pairs.unstable.rate_flagged_mean",
  mean(pp$rate_flagged),
  ci_type = "equal-weight mean over directed pairs",
  n = n_pairs
)
add(
  "pairs.unstable.rate_unflagged_mean",
  mean(pp$rate_unflagged),
  ci_type = "equal-weight mean over directed pairs",
  n = n_pairs
)
add(
  "pairs.unstable.n_flagged_lower",
  sum(pp$rate_flagged < pp$rate_unflagged),
  ci_type = "count of directed pairs",
  n = n_pairs
)
d_pp <- pp$rate_unflagged - pp$rate_flagged
add("pairs.unstable.deficit_min", min(d_pp), ci_type = pt, n = n_pairs)
add("pairs.unstable.deficit_max", max(d_pp), ci_type = pt, n = n_pairs)

# Recurrence -----------------------------------------------------------------------------
fu <- P$flag_unstable
add("recur.calls", nrow(P))
add("recur.transcripts", length(unique(P$transcript_id)))
add("recur.genes", length(unique(P$gene_id)))
add("recur.unstable.calls", sum(fu))
add("recur.unstable.transcripts", length(unique(P$transcript_id[fu])))
add("recur.unstable.genes", length(unique(P$gene_id[fu])))

# Read halves: calls at q < 0.10 and the fallback ----------------------------------------
rh <- list.files(
  file.path(ue_out_root(), "rhalf", "replication", "pairs"),
  pattern = "_M1_(AB|BA)\\.parquet$",
  full.names = TRUE
)
for (f in rh) {
  x <- wpb_read(f)
  tag <- sub("\\.parquet$", "", basename(f))
  add(paste0("rhalf.", tag, ".n_calls"), nrow(x))
  add(paste0("rhalf.", tag, ".n_calls_q10"), x$n_calls_q10[1])
  add(paste0("rhalf.", tag, ".fallback"), as.numeric(isTRUE(x$fallback[1])))
}

ue_write(
  do.call(rbind, rows),
  file.path(out_dir, "cluster_bootstrap.parquet"),
  ue,
  seed
)
ue_write(
  do.call(rbind, draws),
  file.path(out_dir, "cluster_bootstrap_draws.parquet"),
  ue,
  seed
)
ue_write_session(file.path(out_dir, "cluster_bootstrap_session"), ue)
message("34_cluster_bootstrap: done")
