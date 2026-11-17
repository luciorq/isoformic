# Benchmark of test_usage() on the read-split data of GSE129540
#
# Needs the output of data-raw/knockdown/split_halves.sh. For every backend,
# with and without the division of the counts by the overdispersion of the
# replicates, it counts the false positives in the two null contrasts (half 1
# against half 2 of the same libraries) and checks how many calls of the real
# contrast (control against knockdown) in half 1 are found again in half 2.
#
# Usage:
#   ISOFORMIC_KNOCKDOWN_DIR=/path/to/workdir Rscript data-raw/knockdown/benchmark_usage.R
#
# Writes inst/extdata/knockdown/read_split/usage_null.parquet and
# usage_replication.parquet.
suppressPackageStartupMessages({
  library(isoformic)
  library(dplyr)
})

work_dir <- Sys.getenv("ISOFORMIC_KNOCKDOWN_DIR", unset = NA)
if (is.na(work_dir) || !dir.exists(work_dir)) {
  stop(
    "Set ISOFORMIC_KNOCKDOWN_DIR to the directory with quants_half/.",
    call. = FALSE
  )
}
workers <- as.integer(Sys.getenv("ISOFORMIC_WORKERS", unset = "8"))
quants_half <- file.path(work_dir, "quants_half")
ctrl <- c("SRR8870780", "SRR8870781", "SRR8870782")
kd <- c("SRR8870783", "SRR8870784", "SRR8870785")
annot <- prepare_isoformic_annotation(download_reference(
  version = "50",
  file_type = "gff"
))
halves <- function(runs, h) file.path(quants_half, paste0(runs, "_", h))

build <- function(dirs, group) {
  names(dirs) <- paste0(group, "_", seq_along(dirs))
  read_salmon(
    dirs,
    annot_path = annot,
    col_data = data.frame(sample_id = names(dirs), condition = group),
    inferential = TRUE,
    experiment_name = "usage"
  )
}

contrasts <- list(
  null_control = list(
    dirs = c(halves(ctrl, "h1"), halves(ctrl, "h2")),
    group = rep(c("h1", "h2"), each = 3)
  ),
  null_knockdown = list(
    dirs = c(halves(kd, "h1"), halves(kd, "h2")),
    group = rep(c("h1", "h2"), each = 3)
  ),
  real_half1 = list(
    dirs = c(halves(ctrl, "h1"), halves(kd, "h1")),
    group = rep(c("control", "knockdown"), each = 3)
  ),
  real_half2 = list(
    dirs = c(halves(ctrl, "h2"), halves(kd, "h2")),
    group = rep(c("control", "knockdown"), each = 3)
  )
)

rds <- file.path(work_dir, "usage_results.rds")
if (file.exists(rds)) {
  results <- readRDS(rds)
} else {
  results <- list()
  for (nm in names(contrasts)) {
    cc <- contrasts[[nm]]
    iso <- build(cc$dirs, cc$group)
    for (method in c("diffsplice", "satuRn", "DEXSeq")) {
      for (scale in c(FALSE, TRUE)) {
        key <- paste(nm, method, scale, sep = "|")
        message(key)
        t0 <- Sys.time()
        res <- suppressWarnings(test_usage(
          iso,
          group = "condition",
          method = method,
          scale = scale,
          workers = workers
        ))
        res$seconds <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
        results[[key]] <- res
      }
    }
    rm(iso)
    gc()
  }
  saveRDS(results, rds)
}

key_of <- function(x) do.call(rbind, strsplit(names(x), "|", fixed = TRUE))
keys <- key_of(results)
info <- tibble::tibble(
  contrast = keys[, 1],
  method = keys[, 2],
  scaled = keys[, 3] == "TRUE",
  key = names(results)
)

# 1. False positives in the null contrasts --------------------------------------
usage_null <- dplyr::bind_rows(lapply(seq_len(nrow(info)), function(i) {
  if (!startsWith(info$contrast[i], "null")) {
    return(NULL)
  }
  r <- results[[info$key[i]]]
  gene_rows <- r[!duplicated(r$gene_id), ]
  tibble::tibble(
    contrast = info$contrast[i],
    method = info$method[i],
    scaled = info$scaled[i],
    transcripts = nrow(r),
    q_below_0.05 = sum(r$qvalue < 0.05, na.rm = TRUE),
    p_below_0.001 = sum(r$pvalue < 0.001, na.rm = TRUE),
    expected_p_below_0.001 = round(nrow(r) * 0.001, 1),
    genes = nrow(gene_rows),
    genes_q_below_0.05 = sum(gene_rows$gene_qvalue < 0.05, na.rm = TRUE),
    seconds = round(r$seconds[1L])
  )
}))

# 2. Replication of the real contrast across the read halves --------------------
usage_replication <- dplyr::bind_rows(lapply(
  unique(paste(info$method, info$scaled)),
  function(m) {
    sel <- info[paste(info$method, info$scaled) == m, ]
    r1 <- results[[sel$key[sel$contrast == "real_half1"]]]
    r2 <- results[[sel$key[sel$contrast == "real_half2"]]]
    called <- r1[r1$qvalue < 0.05 & !is.na(r1$qvalue), ]
    j <- dplyr::inner_join(
      called[, c("transcript_id", "usage_log2FC")],
      r2[, c("transcript_id", "usage_log2FC", "qvalue")],
      by = "transcript_id",
      suffix = c("_1", "_2")
    )
    tibble::tibble(
      method = sel$method[1L],
      scaled = sel$scaled[1L],
      calls_half1 = nrow(called),
      tested_in_half2 = nrow(j),
      same_direction = round(
        mean(sign(j$usage_log2FC_1) == sign(j$usage_log2FC_2)),
        3
      ),
      significant_again = round(mean(j$qvalue < 0.05, na.rm = TRUE), 3)
    )
  }
))

out <- fs::path("inst", "extdata", "knockdown", "read_split")
fs::dir_create(out)
write_parquet(usage_null, fs::path(out, "usage_null.parquet"))
write_parquet(usage_replication, fs::path(out, "usage_replication.parquet"))
print(as.data.frame(usage_null))
print(as.data.frame(usage_replication))
