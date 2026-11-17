# Compare the groups of Terminus with the groups built from the correlation of
# the replicates, on the six GSE129540 samples
#
# Needs data-raw/knockdown/quantify_knockdown.sh (quants_sub/) and the
# Terminus run on Salmon output with --dumpEqWeights (terminus group and
# terminus collapse on all six samples, output folder terminus_all/out). The
# comparison reports how much the two groupings overlap and how the confidence
# classes change when the groups of Terminus are used; neither grouping is the
# truth.
#
# Usage:
#   ISOFORMIC_KNOCKDOWN_DIR=/path/to/workdir Rscript data-raw/knockdown/compare_terminus.R
#
# Writes inst/extdata/knockdown/read_split/terminus_comparison.parquet.
suppressPackageStartupMessages({
  library(isoformic)
  library(dplyr)
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
annot <- prepare_isoformic_annotation(download_reference(
  version = "50",
  file_type = "gff"
))
levels <- c("control", "knockdown")
cd <- data.frame(sample_id = sample_ids, condition = rep(levels, each = 3))

terminus <- read_terminus_groups(file.path(
  work_dir,
  "terminus_all",
  "out",
  runs[[1L]]
))
iso <- read_salmon(dirs, annot_path = annot, col_data = cd)
by_correlation <- add_uncertainty(
  iso,
  paths = dirs,
  group = "condition",
  levels = levels
)
by_terminus <- add_uncertainty(
  iso,
  paths = dirs,
  group = "condition",
  levels = levels,
  groups = terminus
)

unc_c <- tx_uncertainty(by_correlation)
unc_t <- tx_uncertainty(by_terminus)
assessed <- unc_c$transcript_id[unc_c$assessed %in% TRUE]
in_c <- unc_c$transcript_id[
  !is.na(unc_c$ambiguity_group) & unc_c$assessed %in% TRUE
]
in_t <- unc_t$transcript_id[
  !is.na(unc_t$ambiguity_group) & unc_t$assessed %in% TRUE
]

# Pairs of transcripts in the same correlation group that Terminus also groups
group_of_c <- stats::setNames(unc_c$ambiguity_group, unc_c$transcript_id)
group_of_t <- stats::setNames(unc_t$ambiguity_group, unc_t$transcript_id)
pairs_c <- dplyr::bind_rows(lapply(
  split(in_c, group_of_c[in_c]),
  function(ids) {
    if (length(ids) < 2L) {
      return(NULL)
    }
    as.data.frame(t(utils::combn(ids, 2L)), stringsAsFactors = FALSE)
  }
))
same_t <- group_of_t[pairs_c$V1] == group_of_t[pairs_c$V2]

classes <- function(obj) {
  conf <- de_confidence(obj)
  conf <- conf[abs(conf$lfc_median) > 1 & !is.na(conf$confidence), ]
  table(factor(conf$confidence, c("high", "moderate", "low")))
}
cc <- classes(by_correlation)
ct <- classes(by_terminus)

comparison <- tibble::tibble(
  quantity = c(
    "assessed transcripts",
    "transcripts in a correlation group",
    "transcripts in a Terminus group",
    "correlation-grouped transcripts also grouped by Terminus",
    "Terminus-grouped transcripts also in a correlation group",
    "correlation pairs also in the same Terminus group",
    "median InfRV, correlation-grouped",
    "median InfRV, Terminus-grouped",
    "median InfRV, in neither",
    "large changes (|log2FC| > 1): high, correlation groups",
    "large changes (|log2FC| > 1): moderate, correlation groups",
    "large changes (|log2FC| > 1): low, correlation groups",
    "large changes (|log2FC| > 1): high, Terminus groups",
    "large changes (|log2FC| > 1): moderate, Terminus groups",
    "large changes (|log2FC| > 1): low, Terminus groups"
  ),
  value = c(
    length(assessed),
    length(in_c),
    length(in_t),
    length(intersect(in_c, in_t)),
    length(intersect(in_t, in_c)),
    sum(same_t, na.rm = TRUE),
    stats::median(
      unc_c$mean_inf_rv[unc_c$transcript_id %in% in_c],
      na.rm = TRUE
    ),
    stats::median(
      unc_t$mean_inf_rv[unc_t$transcript_id %in% in_t],
      na.rm = TRUE
    ),
    stats::median(
      unc_c$mean_inf_rv[
        unc_c$transcript_id %in% setdiff(assessed, union(in_c, in_t))
      ],
      na.rm = TRUE
    ),
    cc[["high"]],
    cc[["moderate"]],
    cc[["low"]],
    ct[["high"]],
    ct[["moderate"]],
    ct[["low"]]
  )
)
comparison$value <- round(comparison$value, 3)
out <- fs::path("inst", "extdata", "knockdown", "read_split")
fs::dir_create(out)
write_parquet(comparison, fs::path(out, "terminus_comparison.parquet"))
print(as.data.frame(comparison))
cat("pairs in correlation groups:", nrow(pairs_c), "\n")
