# Does the quantification uncertainty predict what does not reproduce?
#
# Read-split experiment on GSE129540. Each library is quantified twice, on its
# first and second 6M read pairs (data-raw/knockdown/split_halves.sh), which
# gives technical replicates:
#   * two null contrasts, half 1 against half 2 of the same libraries (three
#     control and three knockdown samples): every difference is noise;
#   * the real contrast (control against knockdown), called in half 1 and
#     checked in half 2.
# Compared: swish, edgeR on the counts, edgeR on the counts divided by the
# overdispersion of the replicates (Baldoni et al. 2024), and the descriptive
# flags of de_confidence() (unstable direction, uncertain quantification,
# ambiguous with a sibling isoform).
#
# Usage:
#   ISOFORMIC_KNOCKDOWN_DIR=/path/to/workdir Rscript data-raw/knockdown/validate_uncertainty.R
#
# Writes the small result tables to inst/extdata/knockdown/read_split/.
#
# Memory: without the cached read_split_contrasts.rds the four contrasts are
# recomputed, and each holds several copies of 6 samples x 655k transcripts x
# 30 draws (about 0.94 GB per copy): the peak measured on GENCODE 50 is about
# 7.6 GiB of resident memory, so give that path a memory cap of at least
# 10 GB. Reading the cached .rds needs much less. A cache written before
# de_confidence() reported `unstable` (it has a `confidence` column instead)
# is refused: delete it to recompute.
suppressPackageStartupMessages({
  library(isoformic)
  library(dplyr)
  library(SummarizedExperiment)
  library(edgeR)
})

work_dir <- Sys.getenv("ISOFORMIC_KNOCKDOWN_DIR", unset = NA)
if (is.na(work_dir) || !dir.exists(work_dir)) {
  stop(
    "Set ISOFORMIC_KNOCKDOWN_DIR to the directory with quants_half/.",
    call. = FALSE
  )
}
quants_half <- file.path(work_dir, "quants_half")
rds <- file.path(work_dir, "read_split_contrasts.rds")
ctrl <- c("SRR8870780", "SRR8870781", "SRR8870782")
kd <- c("SRR8870783", "SRR8870784", "SRR8870785")

edger_test <- function(counts, group) {
  d <- DGEList(round(counts), group = group)
  d <- d[filterByExpr(d, group = group), , keep.lib.sizes = FALSE]
  d <- normLibSizes(d)
  design <- model.matrix(~group)
  d <- estimateDisp(d, design)
  tt <- topTags(glmQLFTest(glmQLFit(d, design), coef = 2), n = Inf)$table
  tibble::tibble(
    transcript_id = rownames(tt),
    logFC = tt$logFC,
    p = tt$PValue,
    q = tt$FDR
  )
}

contrast <- function(label, dirs, group, annot) {
  message("== ", label)
  names(dirs) <- paste0(group, "_", seq_along(dirs))
  samples <- data.frame(
    files = file.path(dirs, "quant.sf"),
    names = names(dirs),
    condition = factor(group, levels = unique(group))
  )
  se <- tximeta::tximeta(
    samples,
    type = "salmon",
    txOut = TRUE,
    skipMeta = TRUE,
    useHub = FALSE,
    ignoreAfterBar = TRUE,
    skipSeqinfo = TRUE,
    dropInfReps = FALSE
  )
  y <- fishpond::scaleInfReps(se)
  y <- fishpond::labelKeep(y, minCount = 10, minN = 3)
  # swish computes its p-values from permutations
  set.seed(20261029)
  y <- fishpond::swish(y, x = "condition")
  sw <- as.data.frame(mcols(y))
  sw <- tibble::tibble(
    transcript_id = rownames(sw),
    sw_lfc = sw$log2FC,
    sw_q = sw$qvalue
  )
  naive <- dplyr::rename(
    edger_test(SummarizedExperiment::assay(se, "counts"), samples$condition),
    ed_lfc = logFC,
    ed_p = p,
    ed_q = q
  )
  # edgeR 4.99 takes `sample.dirs` and returns a DGEList of divided counts by
  # default; earlier versions take `paths` and return undivided counts
  cs <- if ("sample.dirs" %in% names(formals(catchSalmon))) {
    catchSalmon(sample.dirs = unname(dirs), DGEList = FALSE, divide = FALSE)
  } else {
    catchSalmon(paths = unname(dirs))
  }
  scaled <- cs$counts / cs$annotation$Overdispersion
  rownames(scaled) <- sub("\\|.*", "", rownames(cs$counts))
  div <- dplyr::rename(
    edger_test(scaled, samples$condition),
    dv_lfc = logFC,
    dv_p = p,
    dv_q = q
  )
  keep <- rownames(se)[rowMeans(SummarizedExperiment::assay(se, "counts")) >= 1]
  tab <- function(a) {
    tibble::as_tibble(
      as.data.frame(SummarizedExperiment::assay(se, a)),
      rownames = "transcript_id"
    ) |>
      dplyr::filter(.data$transcript_id %in% keep)
  }
  iso <- IsoformicExperiment(
    experiment_name = label,
    annot_path = annot,
    assays = list(abundance = tab("abundance"), counts = tab("counts")),
    col_data = data.frame(
      sample_id = names(dirs),
      condition = as.character(group)
    )
  )
  iso <- add_uncertainty(
    iso,
    paths = dirs,
    group = "condition",
    levels = unique(as.character(group)),
    min_count = 5,
    run_qc = FALSE,
    equivalence = FALSE
  )
  conf <- de_confidence(iso) |>
    dplyr::transmute(
      transcript_id = feature_id,
      lfc_median,
      sign_agreement,
      unstable,
      mean_inf_rv,
      # NA for the transcripts below min_count, as `unstable`
      uncertain = dplyr::if_else(
        is.na(sign_agreement),
        NA,
        grepl("uncertain quantification", flags, fixed = TRUE)
      ),
      ambiguous,
      mean_count = mean_count
    )
  conf |>
    dplyr::left_join(sw, by = "transcript_id") |>
    dplyr::left_join(naive, by = "transcript_id") |>
    dplyr::left_join(div, by = "transcript_id")
}

if (file.exists(rds)) {
  res <- readRDS(rds)
  if (!all(c("unstable", "uncertain", "ambiguous") %in% colnames(res[[1L]]))) {
    stop(
      rds,
      " was written by an earlier version of this script; delete it to recompute.",
      call. = FALSE
    )
  }
} else {
  annot <- prepare_isoformic_annotation(download_reference(
    version = "50",
    file_type = "gff"
  ))
  halves <- function(runs, h) file.path(quants_half, paste0(runs, "_", h))
  halves_group <- rep(c("h1", "h2"), each = 3)
  res <- list(
    null_control = contrast(
      "null_control",
      c(halves(ctrl, "h1"), halves(ctrl, "h2")),
      halves_group,
      annot
    ),
    null_knockdown = contrast(
      "null_knockdown",
      c(halves(kd, "h1"), halves(kd, "h2")),
      halves_group,
      annot
    ),
    real_half1 = contrast(
      "real_half1",
      c(halves(ctrl, "h1"), halves(kd, "h1")),
      rep(c("control", "knockdown"), each = 3),
      annot
    ),
    real_half2 = contrast(
      "real_half2",
      c(halves(ctrl, "h2"), halves(kd, "h2")),
      rep(c("control", "knockdown"), each = 3),
      annot
    )
  )
  saveRDS(res, rds)
}

# 1. False positives of the null contrasts, by method ---------------------------
null_calls <- dplyr::bind_rows(lapply(
  c("null_control", "null_knockdown"),
  function(nm) {
    d <- res[[nm]]
    tibble::tibble(
      contrast = nm,
      method = c("swish", "edgeR", "edgeR, uncertainty divided out"),
      tested = c(sum(!is.na(d$sw_q)), sum(!is.na(d$ed_p)), sum(!is.na(d$dv_p))),
      q_below_0.05 = c(
        sum(d$sw_q < 0.05, na.rm = TRUE),
        sum(d$ed_q < 0.05, na.rm = TRUE),
        sum(d$dv_q < 0.05, na.rm = TRUE)
      ),
      p_below_0.001 = c(
        NA_integer_,
        sum(d$ed_p < 0.001, na.rm = TRUE),
        sum(d$dv_p < 0.001, na.rm = TRUE)
      )
    ) |>
      dplyr::mutate(expected_p_below_0.001 = round(.data$tested * 0.001, 1))
  }
))

# 2. Large changes that are stable over the draws in the null contrasts --------
null_stable <- dplyr::bind_rows(lapply(
  c("null_control", "null_knockdown"),
  function(nm) {
    d <- res[[nm]]
    s <- d[
      abs(d$lfc_median) > 1 &
        d$sign_agreement >= 0.95 &
        !is.na(d$sign_agreement),
    ]
    tibble::tibble(
      contrast = nm,
      stable_large_changes = nrow(s),
      uncertain = sum(s$uncertain),
      ambiguous = sum(s$ambiguous %in% TRUE),
      uncertain_or_ambiguous = sum(s$uncertain | s$ambiguous %in% TRUE),
      neither = sum(!s$uncertain & !(s$ambiguous %in% TRUE))
    )
  }
))

# 3. Technical replication of the real contrast ---------------------------------
h1 <- res$real_half1
h2 <- res$real_half2
both <- dplyr::inner_join(
  dplyr::select(
    h1,
    transcript_id,
    lfc1 = lfc_median,
    unstable1 = unstable,
    uncertain1 = uncertain,
    ambiguous1 = ambiguous,
    rv1 = mean_inf_rv,
    sa1 = sign_agreement,
    ed1 = ed_q
  ),
  dplyr::select(
    h2,
    transcript_id,
    lfc2 = lfc_median,
    ed2 = ed_q,
    n2 = mean_count
  ),
  by = "transcript_id"
)
# Each flag of half 1 on its own, among the calls of edgeR in half 1
calls1 <- dplyr::filter(both, .data$ed1 < 0.05)
replication_flag <- dplyr::bind_rows(
  lapply(
    c(
      `unstable direction` = "unstable1",
      `uncertain quantification` = "uncertain1",
      `ambiguous with a sibling` = "ambiguous1"
    ),
    function(column) {
      calls1 |>
        dplyr::filter(!is.na(.data[[column]])) |>
        dplyr::group_by(flagged = .data[[column]]) |>
        dplyr::summarise(
          calls = dplyr::n(),
          same_direction = round(mean(sign(.data$lfc1) == sign(.data$lfc2)), 3),
          significant_again = round(mean(.data$ed2 < 0.05, na.rm = TRUE), 3),
          .groups = "drop"
        )
    }
  ),
  .id = "flag"
)
replication_rv <- both |>
  dplyr::filter(abs(.data$lfc1) > 1, .data$sa1 >= 0.95, .data$n2 >= 5) |>
  dplyr::group_by(
    inferential_rv = cut(
      .data$rv1,
      c(-Inf, 0.5, 2, 10, Inf),
      labels = c("below 0.5", "0.5 to 2", "2 to 10", "above 10")
    )
  ) |>
  dplyr::summarise(
    transcripts = dplyr::n(),
    same_direction = round(mean(sign(.data$lfc1) == sign(.data$lfc2)), 3),
    .groups = "drop"
  ) |>
  dplyr::mutate(inferential_rv = as.character(.data$inferential_rv))

out <- fs::path("inst", "extdata", "knockdown", "read_split")
fs::dir_create(out)
write_parquet(null_calls, fs::path(out, "null_calls.parquet"))
write_parquet(null_stable, fs::path(out, "null_stable.parquet"))
write_parquet(replication_flag, fs::path(out, "replication_flag.parquet"))
write_parquet(replication_rv, fs::path(out, "replication_rv.parquet"))
print(null_calls)
print(null_stable)
print(replication_flag)
print(replication_rv)
