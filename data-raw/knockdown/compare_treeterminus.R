# Compare three ways of grouping isoforms the reads cannot tell apart, on the
# six GSE129540 samples: the correlation of the replicates (the default of
# add_uncertainty()), the consensus groups of Terminus, and the TreeTerminus
# trees cut where the inferential variance falls below the cutoff.
#
# Needs the Salmon runs with --dumpEqWeights (quants_terminus/), the Terminus
# consensus (terminus_all/out, see run_terminus()) and the TreeTerminus trees
# of all samples (tt_mean/cluster_nwk.txt from `treeterminus group -d
# <folder of the six runs> -o tt_mean --mean_inf true`).
#
# Usage:
#   ISOFORMIC_KNOCKDOWN_DIR=/path/to/workdir Rscript data-raw/knockdown/compare_treeterminus.R
#
# Writes inst/extdata/knockdown/read_split/grouping_comparison.parquet.
suppressPackageStartupMessages({
  library(isoformic)
  library(dplyr)
})

work_dir <- Sys.getenv("ISOFORMIC_KNOCKDOWN_DIR", unset = NA)
if (is.na(work_dir) || !dir.exists(work_dir)) {
  stop("Set ISOFORMIC_KNOCKDOWN_DIR to the work directory.", call. = FALSE)
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
  file.path(work_dir, "quants_terminus", runs),
  sample_ids
)
annot <- prepare_isoformic_annotation(download_reference(
  version = "50",
  file_type = "gff"
))
levels <- c("control", "knockdown")
cd <- data.frame(sample_id = sample_ids, condition = rep(levels, each = 3))

iso <- read_salmon(dirs, annot_path = annot, col_data = cd)
base <- add_uncertainty(iso, paths = dirs, group = "condition", levels = levels)

terminus <- read_terminus_groups(file.path(
  work_dir,
  "terminus_all",
  "out",
  runs[[1L]]
))
t0 <- Sys.time()
trees <- read_treeterminus(
  file.path(work_dir, "tt_mean"),
  salmon_dir = dirs[[1L]]
)
cut <- cut_transcript_trees(trees, base)
seconds <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
message(
  "read and cut the trees in ",
  round(seconds),
  " s; ",
  length(unique(cut$group_id)),
  " groups"
)

objects <- list(
  correlation = base,
  Terminus = add_uncertainty(
    base,
    group = "condition",
    levels = levels,
    groups = terminus
  ),
  TreeTerminus = add_uncertainty(
    base,
    group = "condition",
    levels = levels,
    groups = cut
  )
)

summarise_grouping <- function(obj, name) {
  unc <- tx_uncertainty(obj)
  unc <- unc[unc$assessed %in% TRUE, ]
  is_grouped <- !is.na(unc$ambiguity_group)
  n_groups <- length(unique(unc$ambiguity_group[is_grouped]))
  n_left <- sum(unc$uncertain %in% TRUE & !is_grouped)
  collapsed <- de_confidence(
    obj,
    group = "condition",
    levels = levels,
    collapse_groups = TRUE
  )
  group_rows <- collapsed[collapsed$feature_type == "group", ]
  plain <- de_confidence(obj)
  large <- plain[abs(plain$lfc_median) > 1 & !is.na(plain$confidence), ]
  tibble::tibble(
    grouping = name,
    assessed = nrow(unc),
    grouped = sum(is_grouped),
    groups = n_groups,
    uncertain_left_alone = n_left,
    features_after_collapse = nrow(collapsed),
    median_group_inf_rv = stats::median(group_rows$mean_inf_rv, na.rm = TRUE),
    groups_below_cutoff = mean(group_rows$mean_inf_rv < 0.5, na.rm = TRUE),
    large_high = sum(large$confidence == "high"),
    large_moderate = sum(large$confidence == "moderate"),
    large_low = sum(large$confidence == "low")
  )
}
comparison <- dplyr::bind_rows(lapply(names(objects), function(nm) {
  summarise_grouping(objects[[nm]], nm)
}))
comparison$median_group_inf_rv <- round(comparison$median_group_inf_rv, 3)
comparison$groups_below_cutoff <- round(comparison$groups_below_cutoff, 3)
out <- fs::path("inst", "extdata", "knockdown", "read_split")
fs::dir_create(out)
write_parquet(comparison, fs::path(out, "grouping_comparison.parquet"))
print(as.data.frame(comparison))
