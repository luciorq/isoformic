# Compare isoform switch calls with IsoformSwitchAnalyzeR
#
# Runs IsoformSwitchAnalyzeR (DEXSeq test of isoform usage) on the full
# GSE129540 quantifications (control against ESRP2 knockdown) and compares its
# switching isoforms with test_usage() and with the stability and ambiguity of
# the same isoforms from add_uncertainty(). It reports agreement and
# disagreement; there is no ground truth for the real contrast, so nothing here
# says which tool is right.
#
# Usage:
#   ISOFORMIC_KNOCKDOWN_DIR=/path/to/workdir Rscript data-raw/knockdown/compare_isar.R
#
# Writes inst/extdata/knockdown/read_split/isar_comparison.parquet.
suppressPackageStartupMessages({
  library(isoformic)
  library(dplyr)
  library(IsoformSwitchAnalyzeR)
})

work_dir <- Sys.getenv("ISOFORMIC_KNOCKDOWN_DIR", unset = NA)
if (is.na(work_dir) || !dir.exists(work_dir)) {
  stop(
    "Set ISOFORMIC_KNOCKDOWN_DIR to the directory with quants_sub/.",
    call. = FALSE
  )
}
runs <- c(
  "SRR8870780",
  "SRR8870781",
  "SRR8870782",
  "SRR8870783",
  "SRR8870784",
  "SRR8870785"
)
sample_ids <- c(paste0("control_", 1:3), paste0("knockdown_", 1:3))
dirs <- stats::setNames(
  file.path(work_dir, "quants_sub", paste0("GSE129540_", runs)),
  sample_ids
)
design <- data.frame(
  sampleID = sample_ids,
  condition = rep(c("control", "knockdown"), each = 3)
)
rds <- file.path(work_dir, "isar_result.rds")

annot <- prepare_isoformic_annotation(download_reference(
  version = "50",
  file_type = "gff"
))
# IsoformSwitchAnalyzeR reads GTF, not the GFF3 used by isoformic; the file
# is the GENCODE v50 annotation of the transcripts that Salmon quantified
gtf <- file.path(
  work_dir,
  "ref",
  "gencode.v50.chr_patch_hapl_scaff.annotation.gtf.gz"
)
if (!file.exists(gtf)) {
  # The GTF is 123 MB: the default timeout of 60 s stops it partway
  options(timeout = max(3600, getOption("timeout")))
  download.file(
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/gencode.v50.chr_patch_hapl_scaff.annotation.gtf.gz",
    gtf,
    mode = "wb"
  )
}

if (file.exists(rds)) {
  features <- readRDS(rds)
} else {
  expression <- importIsoformExpression(
    sampleVector = stats::setNames(file.path(dirs, "quant.sf"), names(dirs)),
    ignoreAfterBar = TRUE,
    showProgress = FALSE
  )
  sw <- importRdata(
    isoformCountMatrix = expression$counts,
    isoformRepExpression = expression$abundance,
    designMatrix = design,
    isoformExonAnnoation = gtf,
    showProgress = FALSE
  )
  sw <- preFilter(
    sw,
    geneExpressionCutoff = 1,
    isoformExpressionCutoff = 0,
    IFcutoff = 0.01
  )
  sw <- isoformSwitchTestDEXSeq(sw, reduceToSwitchingGenes = FALSE)
  features <- sw$isoformFeatures
  saveRDS(features, rds)
}

# Uncertainty of the same isoforms, and the usage test of isoformic -------------
iso <- read_salmon(
  dirs,
  annot_path = annot,
  col_data = data.frame(sample_id = sample_ids, condition = design$condition),
  inferential = TRUE
)
iso <- add_uncertainty(
  iso,
  group = "condition",
  levels = c("control", "knockdown")
)
frac <- isoform_fraction(iso)
d_if <- rowMeans(frac[, design$condition == "knockdown", drop = FALSE]) -
  rowMeans(frac[, design$condition == "control", drop = FALSE])
d_if <- tibble::tibble(
  transcript_id = names(d_if),
  isoformic_dIF = unname(d_if)
)
usage <- function(scale) {
  test_usage(
    iso,
    group = "condition",
    levels = c("control", "knockdown"),
    method = "DEXSeq",
    scale = scale,
    workers = 8L
  ) |>
    dplyr::select("transcript_id", "qvalue") |>
    dplyr::left_join(d_if, by = "transcript_id")
}
unscaled <- usage(FALSE)
scaled <- usage(TRUE)

isar <- features |>
  dplyr::transmute(
    transcript_id = sub("\\|.*", "", isoform_id),
    isar_dIF = dIF,
    isar_q = isoform_switch_q_value
  ) |>
  dplyr::filter(!is.na(.data$isar_q))
call_sets <- list(
  `IsoformSwitchAnalyzeR` = isar$transcript_id[
    isar$isar_q < 0.05 & abs(isar$isar_dIF) >= 0.1
  ],
  `isoformic DEXSeq, raw counts` = unscaled$transcript_id[
    unscaled$qvalue < 0.05 & abs(unscaled$isoformic_dIF) >= 0.1
  ],
  `isoformic DEXSeq, uncertainty divided out` = scaled$transcript_id[
    scaled$qvalue < 0.05 & abs(scaled$isoformic_dIF) >= 0.1
  ]
)
tested_sets <- list(
  `IsoformSwitchAnalyzeR` = isar$transcript_id,
  `isoformic DEXSeq, raw counts` = unscaled$transcript_id,
  `isoformic DEXSeq, uncertainty divided out` = scaled$transcript_id
)
# `unstable` is the direction stability of the isoform's expression
# fold-change (de_confidence()), not of its usage change: for usage calls it
# mostly marks isoforms whose expression does not change.
conf <- de_confidence(iso)
unc <- tx_uncertainty(iso)
props <- function(ids) {
  u <- unc[match(ids, unc$transcript_id), ]
  cf <- conf[match(ids, conf$feature_id), ]
  tibble::tibble(
    uncertain = sum(u$uncertain %in% TRUE),
    ambiguous = sum(!is.na(u$ambiguity_group)),
    unstable = sum(cf$unstable %in% TRUE)
  )
}
summary_table <- dplyr::bind_rows(lapply(names(call_sets), function(nm) {
  dplyr::bind_cols(
    tibble::tibble(
      tool = nm,
      tested = length(tested_sets[[nm]]),
      calls = length(call_sets[[nm]]),
      shared_with_isar = length(intersect(call_sets[[nm]], call_sets[[1L]])),
      shared_with_scaled = length(intersect(call_sets[[nm]], call_sets[[3L]]))
    ),
    props(call_sets[[nm]])
  )
}))
out <- fs::path("inst", "extdata", "knockdown", "read_split")
fs::dir_create(out)
write_parquet(summary_table, fs::path(out, "isar_comparison.parquet"))
print(as.data.frame(summary_table))
versions <- tibble::tibble(
  package = c("IsoformSwitchAnalyzeR", "isoformic"),
  version = c(
    as.character(utils::packageVersion("IsoformSwitchAnalyzeR")),
    as.character(utils::packageVersion("isoformic"))
  )
)
print(versions)
