# WP-B, post hoc checks (amendment A9): not pre-registered, descriptive only.
#
# Purpose: answer two review questions from outputs that already exist,
# without changing any rule or any pre-registered number.
#   1. Point estimate against the per-draw refits. For every F1 call (edgeR
#      QL on point counts, q < 0.05, full PE vs control), is the point
#      log2FC inside the 2.5-97.5% range of the 30 per-draw refits (the
#      range that 32_story.R reports for FLT1), and if not, by how much?
#      The distance is compared with FLT1-214's. The range is the spread of
#      the Gibbs draws refitted one at a time; it is not a confidence
#      interval for the EM point estimate, so the check also reports how
#      often the point Wald CI overlaps the range and how often the point
#      estimate is further from 0 than the whole range. By de_confidence()
#      class.
#   2. Power adjustment of the replication flags. 31b adjusts the flag ORs
#      linearly on the logit (log2 count, |z|, log2 isoforms, type class).
#      Here the same ORs with natural splines (5 df) on log2 count and |z|,
#      and with the spline tensor interaction, for flag_unstable and
#      flag_primary and both outcomes (replicated, sign concordant); the
#      sign-concordance deficit matched within count x |z| decile cells;
#      and sign concordance and replication among unflagged calls by count
#      and |z| decile. Point estimates only (no CI).
#   3. flag_unstable against flags of the same size in each directed pair
#      made of the lowest-count calls, the lowest-|z| calls, or the calls with
#      the lowest replication predicted from count and |z| (amendment A21):
#      deficit, share of the failures caught and the pairs won. Point
#      estimates only.
#
# Inputs ($ISOFORMIC_UE_OUT/placenta_audit/): prevalence/naive_f1_raw,
#   prevalence/perdraw_f1, prevalence/package_confidence (30_audit_prevalence.R),
#   replication/pooled_dte (31b_replication_summary.R), story/flt1_isoforms
#   (32_story.R).
# Outputs: $ISOFORMIC_UE_OUT/placenta_audit/post_hoc/
#   point_vs_draws.parquet     key, estimate, n (long)
#   power_sensitivity.parquet  key, estimate, n (long)
#   point_vs_draws_calls.parquet  one row per F1 call
#   unstable_same_size_controls.parquet  key, estimate, n (long; A21)
# Env vars: ISOFORMIC_UE_OUT, ISOFORMIC_UE_PKG (ue_common.R).
# Resources (estimate): about 2 min, < 3 GB.
# Run from the repository root, inside a memory cap of about 6 GB (see
# README.md):
#   pixi run --frozen Rscript data-raw/uncertainty_evidence/33_post_hoc_checks.R

source("data-raw/uncertainty_evidence/ue_common.R")
ue <- ue_setup("33_post_hoc_checks")
source("data-raw/uncertainty_evidence/38_audit_common.R")
seed <- wpb_seed(ue, 700L)
set.seed(seed)
out_dir <- wpb_dir("post_hoc")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
rd <- function(...) wpb_read(wpb_dir(...))
rows <- list()
add <- function(key, estimate, n = NA_real_) {
  # jarl-ignore undesirable_operator: appends to a result list of the enclosing scope
  rows[[length(rows) + 1L]] <<- data.frame(
    key = key,
    estimate = as.numeric(estimate),
    n = as.numeric(n)
  )
}

# 1. Point estimate against the per-draw refits ------------------------------------------
f1 <- rd("prevalence", "naive_f1_raw.parquet")
calls <- f1[!is.na(f1$qvalue) & f1$qvalue < 0.05, ]
pd <- rd("prevalence", "perdraw_f1.parquet")
pd <- pd[pd$feature_id %in% calls$feature_id, ]
n_draw <- max(pd$draw)
L <- matrix(NA_real_, nrow(calls), n_draw)
i <- match(pd$feature_id, calls$feature_id)
L[cbind(i, pd$draw)] <- pd$logFC
# the same quantiles as perdraw_stats() in 32_story.R
lo <- apply(L, 1L, wpb_q, 0.025)
hi <- apply(L, 1L, wpb_q, 0.975)
se <- abs(calls$logFC) / pmax(abs(calls$z), 1e-8)
conf <- rd("prevalence", "package_confidence.parquet")
conf <- conf[conf$feature_type == "transcript", ]
d <- data.frame(
  feature_id = calls$feature_id,
  point = calls$logFC,
  draw_lo = lo,
  draw_hi = hi,
  n_draws = rowSums(!is.na(L)),
  wald_lo = calls$logFC - 1.96 * se,
  wald_hi = calls$logFC + 1.96 * se,
  confidence = conf$confidence[match(calls$feature_id, conf$feature_id)]
)
d$outside_by <- pmax(d$draw_lo - d$point, d$point - d$draw_hi, 0)
d$outside <- d$outside_by > 0
d$wald_overlaps_range <- d$wald_lo <= d$draw_hi & d$draw_lo <= d$wald_hi
# point further from 0 than the whole range (on the side of its sign)
d$point_beyond_far_end <- d$outside &
  ((d$point > 0 & d$point > d$draw_hi) | (d$point < 0 & d$point < d$draw_lo))
fl <- wpb_read(wpb_dir("story", "flt1_isoforms.parquet"))
f214 <- fl[fl$transcript_name == "FLT1-214", ]
stopifnot(nrow(f214) == 1L)
ref <- d$outside_by[d$feature_id == f214$transcript_id]
stopifnot(length(ref) == 1L, ref > 0)
n <- nrow(d)
add("f1.n_calls", n, n)
add("f1.n_calls_without_draws", sum(d$n_draws == 0), n)
add("f1.outside_range.n", sum(d$outside), n)
add("f1.outside_range.share", mean(d$outside), n)
for (t in c(0.5, 1)) {
  add(sprintf("f1.outside_by_gt_%s.n", t), sum(d$outside_by > t), n)
  add(sprintf("f1.outside_by_gt_%s.share", t), mean(d$outside_by > t), n)
}
add("f1.flt1_214.outside_by", ref, 1)
add("f1.outside_by_ge_flt1_214.n", sum(d$outside_by >= ref), n)
add("f1.outside_by_ge_flt1_214.share", mean(d$outside_by >= ref), n)
add(
  "f1.outside.wald_overlaps_range.n",
  sum(d$outside & d$wald_overlaps_range),
  sum(d$outside)
)
add(
  "f1.outside.point_beyond_far_end.n",
  sum(d$point_beyond_far_end),
  sum(d$outside)
)
for (cl in sort(unique(stats::na.omit(d$confidence)))) {
  s <- d$confidence %in% cl
  add(sprintf("f1.class_%s.n", cl), sum(s), sum(s))
  add(
    sprintf("f1.class_%s.outside_range.share", cl),
    mean(d$outside[s]),
    sum(s)
  )
  add(
    sprintf("f1.class_%s.outside_by_ge_flt1_214.share", cl),
    mean(d$outside_by[s] >= ref),
    sum(s)
  )
}
# the named isoforms of the story (FLT1, LEP), from the same computation
named <- fl[
  fl$transcript_id %in% d$feature_id,
  c("transcript_id", "transcript_name")
]
txu <- wpb_read(file.path(ue_out_root(), "store", "strata_tx.parquet"))
lep_ids <- txu$transcript_id[
  txu$gene_name %in% "LEP" & txu$transcript_id %in% d$feature_id
]
named <- rbind(
  named,
  data.frame(
    transcript_id = lep_ids,
    transcript_name = txu$transcript_name[match(lep_ids, txu$transcript_id)]
  )
)
for (k in seq_len(nrow(named))) {
  r <- d[d$feature_id == named$transcript_id[k], ]
  nm <- named$transcript_name[k]
  for (cc in c(
    "point",
    "draw_lo",
    "draw_hi",
    "outside_by",
    "wald_lo",
    "wald_hi"
  )) {
    add(sprintf("f1.named.%s.%s", nm, cc), r[[cc]], 1)
  }
}
d$transcript_name <- txu$transcript_name[match(d$feature_id, txu$transcript_id)]
ue_write(d, file.path(out_dir, "point_vs_draws_calls.parquet"), ue, seed)
ue_write(
  do.call(rbind, rows),
  file.path(out_dir, "point_vs_draws.parquet"),
  ue,
  seed
)

# 2. Power adjustment of the replication flags --------------------------------------------
rows <- list()
P <- wpb_read(wpb_dir("replication", "pooled_dte.parquet"))
P$type_class <- factor(P$type_class)
P$cell <- interaction(
  dplyr::ntile(P$log2_count, 10),
  dplyr::ntile(P$d_absz, 10),
  drop = TRUE
)
lin <- "log2_count + d_absz + log2_iso + type_class"
spl <- "splines::ns(log2_count, 5) + splines::ns(d_absz, 5) + log2_iso + type_class"
int <- "splines::ns(log2_count, 5) * splines::ns(d_absz, 5) + log2_iso + type_class"
or_of <- function(y, flag, rhs) {
  dd <- P[!is.na(P[[flag]]), ]
  fit <- suppressWarnings(stats::glm(
    stats::as.formula(paste(y, "~", flag, "+", rhs)),
    data = dd,
    family = stats::binomial()
  ))
  exp(stats::coef(fit)[[paste0(flag, "TRUE")]])
}
matched <- function(y, flag) {
  f <- P[[flag]]
  num <- 0
  w <- 0
  for (ix in split(seq_len(nrow(P)), P$cell)) {
    a <- ix[f[ix] %in% TRUE]
    b <- ix[f[ix] %in% FALSE]
    if (length(a) && length(b)) {
      num <- num + length(a) * (mean(P[[y]][b]) - mean(P[[y]][a]))
      w <- w + length(a)
    }
  }
  num / w
}
for (flag in c("flag_unstable", "flag_primary")) {
  for (y in c("replicated", "sign_concordant")) {
    for (m in c("linear", "splines", "spline_interaction")) {
      rhs <- switch(m, linear = lin, splines = spl, spline_interaction = int)
      add(
        sprintf("power.%s.%s.or_%s", flag, y, m),
        or_of(y, flag, rhs),
        nrow(P)
      )
    }
    f <- P[[flag]]
    add(
      sprintf("power.%s.%s.raw_deficit", flag, y),
      mean(P[[y]][f %in% FALSE]) - mean(P[[y]][f %in% TRUE]),
      nrow(P)
    )
    add(
      sprintf("power.%s.%s.matched_deficit_10x10", flag, y),
      matched(y, flag),
      nrow(P)
    )
  }
}
# gradients among unflagged calls (flag_unstable FALSE)
u <- P[P$flag_unstable %in% FALSE, ]
for (v in c("log2_count", "d_absz")) {
  dec <- dplyr::ntile(u[[v]], 10)
  for (y in c("sign_concordant", "replicated")) {
    r <- tapply(u[[y]], dec, mean)
    for (k in names(r)) {
      add(
        sprintf("power.unflagged.%s.by_%s_decile.%s", y, v, k),
        r[[k]],
        sum(dec == as.integer(k))
      )
    }
  }
}
add(
  "power.median_log2_count.flag_unstable",
  stats::median(P$log2_count[P$flag_unstable %in% TRUE]),
  sum(P$flag_unstable %in% TRUE)
)
add(
  "power.median_log2_count.not_flag_unstable",
  stats::median(P$log2_count[P$flag_unstable %in% FALSE]),
  sum(P$flag_unstable %in% FALSE)
)
ue_write(
  do.call(rbind, rows),
  file.path(out_dir, "power_sensitivity.parquet"),
  ue,
  seed
)

# 3. flag_unstable against controls of the same size (amendment A21) ------------
# In each directed pair, flag as many calls as flag_unstable does, choosing the
# lowest mean count, the lowest |z|, or the lowest replication predicted from
# count and |z| (a logistic fit on all calls). Point estimates only (no CI).
rows <- list()
fit_cz <- suppressWarnings(stats::glm(
  replicated ~ log2_count + d_absz,
  data = P,
  family = stats::binomial()
))
pred_cz <- stats::fitted(fit_cz)
unstable <- P$flag_unstable %in% TRUE
by_pair <- split(seq_len(nrow(P)), P$pair_id)
same_size <- function(score) {
  out <- logical(nrow(P))
  for (ix in by_pair) {
    k <- sum(unstable[ix])
    if (k > 0L) {
      out[ix[order(score[ix])[seq_len(k)]]] <- TRUE
    }
  }
  out
}
deficit_of <- function(f, ix = seq_len(nrow(P))) {
  mean(P$replicated[ix][!f[ix]]) - mean(P$replicated[ix][f[ix]])
}
controls <- list(
  lowest_count = same_size(P$d_mean_count),
  lowest_absz = same_size(P$d_absz),
  lowest_predicted_count_absz = same_size(pred_cz)
)
add("same_size.flag_unstable.deficit", deficit_of(unstable), nrow(P))
add(
  "same_size.flag_unstable.sensitivity_nonreplication",
  mean(unstable[!P$replicated]),
  sum(!P$replicated)
)
for (nm in names(controls)) {
  f <- controls[[nm]]
  add(sprintf("same_size.%s.deficit", nm), deficit_of(f), nrow(P))
  add(
    sprintf("same_size.%s.sensitivity_nonreplication", nm),
    mean(f[!P$replicated]),
    sum(!P$replicated)
  )
  wins <- vapply(
    by_pair,
    function(ix) deficit_of(unstable, ix) > deficit_of(f, ix),
    logical(1L)
  )
  add(
    sprintf("same_size.%s.pairs_unstable_larger_deficit", nm),
    sum(wins, na.rm = TRUE),
    length(by_pair)
  )
}
ue_write(
  do.call(rbind, rows),
  file.path(out_dir, "unstable_same_size_controls.parquet"),
  ue,
  seed
)

ue_write_session(out_dir, ue)
message("33_post_hoc_checks: done")
