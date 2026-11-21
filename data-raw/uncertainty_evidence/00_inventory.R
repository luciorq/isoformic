# Inventory of the local quantifications that carry inferential replicates.
#
# Purpose: list, for every real dataset quantified in the work directories, what the
# uncertainty evidence experiments can use: organism, annotation (Salmon
# index), samples per group, read depth, library type, quantifier and
# version, number and kind of inferential replicates, and the run time of the
# quantifier. Everything is read from the quantifier's own output
# (aux_info/meta_info.json, cmd_info.json, lib_format_counts.json for Salmon;
# <run>.meta_info.json, <run>.ambig_info.tsv and <run>.infreps.pq for
# Oarfish); group labels come from the sample sheets in data-raw/ and from
# the GEO titles hard-coded below (GSE129540, GSE139962, GSE75489).
#
# Inputs (directories; missing subdirectories are reported, not fatal):
#   ISOFORMIC_V50_DIR        work directory of the GENCODE 50 Salmon and
#                            Oarfish runs (required; layout in README.md)
#   ISOFORMIC_KNOCKDOWN_DIR  work directory of data-raw/knockdown/ (required;
#                            its quants/ holds the GSE139962 runs)
#   ISOFORMIC_UE_OUT         output root (required)
#
# Resources: base R and jsonlite only, a few seconds, well under 1 GB.
# Run from the repository root:
#   pixi run --frozen Rscript data-raw/uncertainty_evidence/00_inventory.R
#
# Writes inventory_samples.tsv (one row per quantified sample) and
# inventory_datasets.tsv (one row per dataset) to ISOFORMIC_UE_OUT.

need_env <- function(name) {
  value <- Sys.getenv(name)
  if (!nzchar(value)) {
    stop(
      "Set the environment variable ",
      name,
      "; see data-raw/uncertainty_evidence/README.md",
      call. = FALSE
    )
  }
  value
}
v50 <- need_env("ISOFORMIC_V50_DIR")
kd <- need_env("ISOFORMIC_KNOCKDOWN_DIR")
out <- need_env("ISOFORMIC_UE_OUT")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

# Group labels -----------------------------------------------------------------
placenta <- utils::read.delim("data-raw/placenta/placenta_samples.tsv")
mouse <- utils::read.delim("data-raw/case_studies/mouse_bnst_samples.tsv")
prostate <- utils::read.delim("data-raw/case_studies/prostate_ont_samples.tsv")
geo <- c(
  # GSE129540 (ESRP2 knockdown, GEO sample titles)
  SRR8870780 = "control",
  SRR8870781 = "control",
  SRR8870782 = "control",
  SRR8870783 = "knockdown",
  SRR8870784 = "knockdown",
  SRR8870785 = "knockdown",
  # GSE75489 (ESRP1/2 knockdown, H358; GSM1957010-12 KD, 13-15 CTRL)
  SRR2966412 = "knockdown",
  SRR2966413 = "knockdown",
  SRR2966414 = "knockdown",
  SRR2966415 = "control",
  SRR2966416 = "control",
  SRR2966417 = "control"
)
labels <- c(
  stats::setNames(placenta$condition, placenta$run),
  stats::setNames(paste(mouse$sex, mouse$condition), mouse$run),
  stats::setNames(prostate$condition, prostate$run),
  geo
)

# Readers ----------------------------------------------------------------------
parse_time <- function(x) {
  if (is.null(x) || !nzchar(x)) {
    return(NA)
  }
  as.POSIXct(x, format = "%a %b %d %H:%M:%S %Y", tz = "UTC")
}

salmon_row <- function(dir, dataset) {
  meta <- jsonlite::read_json(file.path(dir, "aux_info", "meta_info.json"))
  cmd <- jsonlite::read_json(file.path(dir, "cmd_info.json"))
  lib_file <- file.path(dir, "lib_format_counts.json")
  lib <- if (file.exists(lib_file)) {
    jsonlite::read_json(lib_file)$expected_format
  } else {
    NA_character_
  }
  start <- parse_time(meta$start_time)
  end <- parse_time(meta$end_time)
  run <- sub("^GSE[0-9]+_", "", basename(dir))
  run <- sub("_h[12]$", "", run)
  data.frame(
    dataset = dataset,
    sample = basename(dir),
    group = unname(labels[run]) %||% NA_character_,
    quantifier = paste("salmon", meta$salmon_version),
    index = basename(cmd$index),
    library_type = lib,
    fragments = meta$num_processed,
    mapped = meta$num_mapped,
    percent_mapped = round(meta$percent_mapped, 1),
    replicate_kind = meta$samp_type,
    replicates = meta$num_bootstraps,
    replicate_file = file.exists(
      file.path(dir, "aux_info", "bootstrap", "bootstraps.gz")
    ),
    eq_classes = file.exists(file.path(dir, "aux_info", "eq_classes.txt.gz")) ||
      file.exists(file.path(dir, "aux_info", "eq_classes.txt")),
    threads = cmd$threads %||% NA,
    minutes = round(as.numeric(difftime(end, start, units = "mins")), 1)
  )
}

oarfish_row <- function(prefix, dataset) {
  meta <- jsonlite::read_json(paste0(prefix, ".meta_info.json"))
  ambig <- utils::read.delim(paste0(prefix, ".ambig_info.tsv"))
  start <- as.POSIXct(
    meta$start_time,
    format = "%a %b %d %H:%M:%S %Y",
    tz = "UTC"
  )
  end <- as.POSIXct(meta$end_time, format = "%a %b %d %H:%M:%S %Y", tz = "UTC")
  run <- basename(prefix)
  data.frame(
    dataset = dataset,
    sample = run,
    group = unname(labels[run]) %||% NA_character_,
    quantifier = paste("oarfish", meta$oarfish_version %||% NA),
    index = "GENCODE transcripts (built-in minimap2)",
    library_type = "ONT long read",
    # ambig_info.tsv has one row per transcript; the read total is the sum
    # of unique and ambiguous assignments, an upper bound on aligned reads
    # because an ambiguous read is counted for every transcript it hits
    fragments = NA_real_,
    mapped = sum(ambig$unique_reads),
    percent_mapped = NA_real_,
    replicate_kind = "bootstrap",
    replicates = meta$num_bootstraps,
    replicate_file = file.exists(paste0(prefix, ".infreps.pq")),
    eq_classes = FALSE,
    threads = meta$threads %||% NA,
    minutes = round(as.numeric(difftime(end, start, units = "mins")), 1)
  )
}

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

salmon_set <- function(root, dataset) {
  if (!dir.exists(root)) {
    message("missing: ", root)
    return(NULL)
  }
  dirs <- list.dirs(root, recursive = FALSE)
  dirs <- dirs[file.exists(file.path(dirs, "aux_info", "meta_info.json"))]
  do.call(rbind, lapply(dirs, salmon_row, dataset = dataset))
}

oarfish_set <- function(root, dataset) {
  files <- list.files(root, pattern = "\\.meta_info\\.json$", full.names = TRUE)
  prefixes <- sub("\\.meta_info\\.json$", "", files)
  do.call(rbind, lapply(prefixes, oarfish_row, dataset = dataset))
}

# Datasets ---------------------------------------------------------------------
sets <- list(
  salmon_set(file.path(v50, "human/placenta/quants"), "placenta_GSE148241"),
  salmon_set(
    file.path(v50, "cases/mouse_bnst/quants"),
    "mouse_bnst_PRJNA984605"
  ),
  salmon_set(file.path(v50, "human/quants_sub"), "esrp_kd_12M_v50"),
  salmon_set(file.path(v50, "human/quants_half"), "esrp2_kd_halves_v50"),
  salmon_set(file.path(v50, "human/g100"), "esrp2_kd_g100_v50"),
  salmon_set(file.path(v50, "human/quants_terminus"), "esrp2_kd_eqweights_v50"),
  salmon_set(file.path(kd, "quants"), "pc3_esrp1_kd_GSE139962_v49"),
  oarfish_set(
    file.path(v50, "cases/prostate_ont/quants"),
    "prostate_ont_GSE259334"
  )
)
samples <- do.call(rbind, sets)
samples$dataset[grepl("^GSE75489", samples$sample)] <- sub(
  "esrp_kd",
  "esrp12_kd_h358_GSE75489",
  samples$dataset[grepl("^GSE75489", samples$sample)],
  fixed = TRUE
)
samples$dataset[grepl("^GSE129540", samples$sample)] <- sub(
  "esrp_kd",
  "esrp2_kd_GSE129540",
  samples$dataset[grepl("^GSE129540", samples$sample)],
  fixed = TRUE
)
pc3 <- samples$dataset == "pc3_esrp1_kd_GSE139962_v49"
samples$group[pc3] <- ifelse(
  grepl("^siNC", samples$sample[pc3]),
  "control",
  "knockdown"
)

datasets <- do.call(
  rbind,
  lapply(split(samples, samples$dataset), function(d) {
    tab <- table(d$group, useNA = "ifany")
    data.frame(
      dataset = d$dataset[1],
      n_samples = nrow(d),
      groups = paste(names(tab), tab, sep = ": ", collapse = "; "),
      quantifier = paste(unique(d$quantifier), collapse = ","),
      index = paste(unique(d$index), collapse = ","),
      library_type = paste(unique(d$library_type), collapse = ","),
      fragments_median = stats::median(d$fragments),
      fragments_min = min(d$fragments),
      fragments_max = max(d$fragments),
      percent_mapped_range = paste(range(d$percent_mapped), collapse = "-"),
      replicates = paste(
        unique(paste(d$replicates, d$replicate_kind)),
        collapse = ","
      ),
      all_replicate_files = all(d$replicate_file),
      resampling_df = sum(d$replicates - 1),
      eq_classes = all(d$eq_classes),
      minutes_median = stats::median(d$minutes),
      threads = paste(unique(d$threads), collapse = ",")
    )
  })
)

utils::write.table(
  samples,
  file.path(out, "inventory_samples.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
utils::write.table(
  datasets,
  file.path(out, "inventory_datasets.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
print(datasets, right = FALSE)
