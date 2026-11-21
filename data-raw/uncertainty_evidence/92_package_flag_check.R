# WP-B, check of the package code (amendment A24): not pre-registered, added
# after the results, on review.
#
# Purpose: show that the `unstable` flag of the current de_confidence() is
# the flag that the evaluation scored. The evaluation ran on the frozen
# snapshot of prereg.yaml with its own implementation of the per-draw
# contrast (wpb_confidence(), 38_audit_common.R); the defaults and the help
# text of the package quote its results as properties of `unstable`. For
# each of the 60 directed half-split pairs of 31_halfsplit_replication.R,
# the discovery half is rebuilt as an IsoformicExperiment from the S0 store,
# the evaluation's edgeR calls are set with `de_tx<-` (every call counts as
# significant), and add_uncertainty(group = "condition") + de_confidence()
# of the package code recompute the sign agreement. The package values are
# compared with the evaluation's (sign_agreement < 0.95 among the assessed
# calls, mean InfRV, `uncertain`, `ambiguous`), and the package flag with
# `flag_unstable` of the pooled replication table. Nothing is re-estimated:
# when the flags agree on every call, every replication number of 31b and
# 34 describes the package flag as well.
#
# Stages (ISOFORMIC_UE_STAGE: pairs | summary | all, default all):
#   pairs    one table per directed pair (skips pairs already written)
#   summary  the long table of agreement counts
#
# Inputs ($ISOFORMIC_UE_OUT/): store/ (10_placenta_prepare.R),
#   placenta_audit/replication/pairs/<pair>_info.parquet and
#   <pair>_dte.parquet (31_halfsplit_replication.R),
#   placenta_audit/replication/pooled_dte.parquet (31b_replication_summary.R).
# Outputs: $ISOFORMIC_UE_OUT/placenta_audit/post_hoc/
#   package_check/pairs/<pair>.parquet  one row per call: evaluation and
#                                       package values side by side
#   package_flag_check.parquet          key, estimate, n (long), read by
#                                       90_article_tables.R as B.pkgcheck.*
# Env vars: ISOFORMIC_UE_OUT; ISOFORMIC_UE_PKG (ue_common.R; default here
#   ".", the repository code, which is the code being checked: run from the
#   repository root); ISOFORMIC_UE_ANNOT (annotation directory of GENCODE 50,
#   default the one in the isoformic cache); ISOFORMIC_UE_STAGE;
#   ISOFORMIC_UE_THREADS.
# Resources (measured 2026-11-20): pairs about 8 min (60 pairs, 2 to 10 s
#   each, after about 1 min to load the store), peak RSS 3.2 GB; summary
#   seconds. Run from the repository root, inside a memory cap of about 6 GB
#   (see README.md):
#   pixi run --frozen Rscript data-raw/uncertainty_evidence/92_package_flag_check.R

if (!nzchar(Sys.getenv("ISOFORMIC_UE_PKG"))) {
  Sys.setenv(ISOFORMIC_UE_PKG = ".")
}
source("data-raw/uncertainty_evidence/ue_common.R")
ue <- ue_setup("92_package_flag_check")
source("data-raw/uncertainty_evidence/38_audit_common.R")
stage <- Sys.getenv("ISOFORMIC_UE_STAGE", "all")
annot_dir <- Sys.getenv(
  "ISOFORMIC_UE_ANNOT",
  path.expand(
    "~/.cache/R/isoformic/gencode_v50_chr_patch_hapl_scaff_annotation_gff3"
  )
)
out_dir <- wpb_dir("post_hoc")
pair_out <- file.path(out_dir, "package_check", "pairs")
dir.create(pair_out, recursive = TRUE, showWarnings = FALSE)
eval_pair_dir <- wpb_dir("replication", "pairs")
say <- function(...) {
  message(format(Sys.time(), "%H:%M:%S "), sprintf(...))
}
# The evaluation's threshold of the unstable direction (prereg.yaml)
min_agreement <- 0.95

# Pairs: the package flag on every discovery half -------------------------------------
if (stage %in% c("pairs", "all")) {
  st <- wpb_load_store(draws = TRUE, moments = FALSE, say = say)
  runs <- st$runs
  arr <- st$arr
  st$arr <- NULL
  st$strata <- NULL
  iso <- IsoformicExperiment(
    experiment_name = "placenta_package_check",
    annot_path = annot_dir,
    assays = list(counts = st$counts),
    col_data = data.frame(
      sample_id = runs,
      condition = as.character(st$groups),
      library_size = unname(st$lib_size[runs])
    )
  )
  meta <- iso@metadata
  meta$inf_reps <- lapply(seq_len(dim(arr)[3L]), function(k) arr[,, k])
  iso@metadata <- meta
  rm(arr, meta)
  invisible(gc())
  say("object built, peak RSS %.2f GB", ue_peak_rss_gb())

  info_files <- list.files(
    eval_pair_dir,
    pattern = "_info\\.parquet$",
    full.names = TRUE
  )
  pair_ids <- sub("_info\\.parquet$", "", basename(info_files))
  for (pid in pair_ids) {
    dest <- file.path(pair_out, paste0(pid, ".parquet"))
    if (file.exists(dest)) {
      next
    }
    t0 <- Sys.time()
    info <- wpb_read(file.path(eval_pair_dir, paste0(pid, "_info.parquet")))
    discovery <- strsplit(info$D, ",", fixed = TRUE)[[1]]
    ev <- wpb_read(file.path(eval_pair_dir, paste0(pid, "_dte.parquet")))
    x <- iso[, discovery]
    de_tx(x) <- data.frame(
      transcript_id = ev$transcript_id,
      log2FC = ev$d_logFC,
      pvalue = ev$d_p,
      qvalue = ev$d_q
    )
    # every discovery call of the evaluation counts as significant
    de_cutoffs(x) <- list(log2fc = 0, pvalue = 1)
    x <- add_uncertainty(
      x,
      group = "condition",
      levels = c("control", "preeclampsia"),
      equivalence = FALSE,
      run_qc = FALSE,
      keep_draws = FALSE
    )
    conf <- as.data.frame(de_confidence(x))
    conf <- conf[conf$feature_type == "transcript", ]
    m <- conf[match(ev$transcript_id, conf$feature_id), ]
    res <- data.frame(
      pair_id = pid,
      transcript_id = ev$transcript_id,
      in_pkg = !is.na(m$feature_id),
      sign_agreement_ev = ev$sign_agreement,
      confidence_ev = ev$confidence,
      uncertain_ev = ev$uncertain,
      ambiguous_ev = ev$ambiguous,
      mean_inf_rv_ev = ev$mean_inf_rv,
      sign_agreement_pkg = m$sign_agreement,
      unstable_pkg = m$unstable,
      unstable_de_pkg = m$unstable_de,
      n_det_sig_pkg = m$n_det_sig,
      mean_inf_rv_pkg = m$mean_inf_rv,
      ambiguous_pkg = m$ambiguous
    )
    ue_write(res, dest, ue)
    say(
      "%s: %d calls, %d missing, %.1f s, peak RSS %.2f GB",
      pid,
      nrow(res),
      sum(!res$in_pkg),
      as.numeric(difftime(Sys.time(), t0, units = "secs")),
      ue_peak_rss_gb()
    )
  }
  rm(iso)
  invisible(gc())
}

# Summary: agreement counts ------------------------------------------------------------
if (stage %in% c("summary", "all")) {
  files <- list.files(pair_out, pattern = "\\.parquet$", full.names = TRUE)
  pk <- do.call(rbind, lapply(files, wpb_read))
  rows <- list()
  add <- function(key, estimate, n = NA_real_) {
    # jarl-ignore undesirable_operator: appends to a result list of the enclosing scope
    rows[[length(rows) + 1L]] <<- data.frame(
      key = key,
      estimate = as.numeric(estimate),
      n = as.numeric(n)
    )
  }
  # assessed by the evaluation: the calls with a confidence class
  assessed <- !is.na(pk$confidence_ev)
  n_calls <- nrow(pk)
  add("pairs", length(unique(pk$pair_id)))
  add("n_calls", n_calls)
  add("calls_missing", sum(!pk$in_pkg), n_calls)
  add(
    "assessed_mismatch",
    sum(assessed != !is.na(pk$sign_agreement_pkg)),
    n_calls
  )
  add("n_assessed", sum(assessed), n_calls)
  add(
    "unstable_disagree",
    sum(
      (pk$unstable_pkg != (pk$sign_agreement_ev < min_agreement))[assessed],
      na.rm = TRUE
    ),
    sum(assessed)
  )
  add(
    "max_abs_diff_sign_agreement",
    max(abs(pk$sign_agreement_pkg - pk$sign_agreement_ev)[assessed]),
    sum(assessed)
  )
  add(
    "max_abs_diff_mean_inf_rv",
    max(abs(pk$mean_inf_rv_pkg - pk$mean_inf_rv_ev), na.rm = TRUE),
    n_calls
  )
  add(
    "uncertain_disagree",
    sum((pk$mean_inf_rv_pkg >= 0.5) != pk$uncertain_ev, na.rm = TRUE),
    n_calls
  )
  add(
    "ambiguous_disagree",
    sum((pk$ambiguous_pkg != pk$ambiguous_ev)[assessed], na.rm = TRUE),
    sum(assessed)
  )
  # the flag scored by 31b and 34
  pooled <- wpb_read(wpb_dir("replication", "pooled_dte.parquet"))
  i <- match(
    paste(pooled$pair_id, pooled$transcript_id),
    paste(pk$pair_id, pk$transcript_id)
  )
  add("pooled_rows", nrow(pooled))
  add("pooled_unmatched", sum(is.na(i)), nrow(pooled))
  add(
    "pooled_flag_disagree",
    sum(pooled$flag_unstable != pk$unstable_pkg[i], na.rm = TRUE),
    nrow(pooled)
  )
  res <- do.call(rbind, rows)
  ue_write(res, file.path(out_dir, "package_flag_check.parquet"), ue)
  print(res, right = FALSE)
  say("summary done, peak RSS %.2f GB", ue_peak_rss_gb())
}
