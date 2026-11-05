# Minimized ESRP knockdown examples (GSE75489 and GSE129540)
#
# Builds inst/extdata/knockdown/<series>/ and inst/extdata/knockdown_annotation/
# from full Salmon quantifications with 30 Gibbs draws. The Salmon output is
# large (about 100 MB per sample with the draws and the equivalence classes)
# and is not part of the package: it is rebuilt from the public reads with
# data-raw/knockdown/quantify_knockdown.sh and can be deposited on Zenodo.
#
# The swish test and the library sizes use the complete quantifications; the
# bundled tables keep only the transcripts of a gene panel (known ESRP targets
# plus the genes with the most significant transcripts), so q-values are those
# of the full analysis.
#
# Usage:
#   ISOFORMIC_KNOCKDOWN_DIR=/path/to/workdir Rscript data-raw/build_knockdown_example.R
#
# The work directory holds quants_sub/<series>_<run>/ (one Salmon output
# directory per run) as written by quantify_knockdown.sh.
suppressPackageStartupMessages({
  library(isoformic)
  library(dplyr)
  library(SummarizedExperiment)
})

work_dir <- Sys.getenv("ISOFORMIC_KNOCKDOWN_DIR", unset = NA)
if (is.na(work_dir) || !dir.exists(work_dir)) {
  stop(
    "Set ISOFORMIC_KNOCKDOWN_DIR to the directory with quants_sub/.",
    call. = FALSE
  )
}
extdata <- fs::path("inst", "extdata")
out_root <- fs::path(extdata, "knockdown")
n_top_genes <- 12L

series <- list(
  GSE75489 = list(
    description = paste(
      "ESRP1/2 shRNA knockdown in H358 lung epithelial cells",
      "(Yang et al. 2016, Mol Cell Biol)"
    ),
    runs = c(
      "SRR2966415",
      "SRR2966416",
      "SRR2966417",
      "SRR2966412",
      "SRR2966413",
      "SRR2966414"
    )
  ),
  GSE129540 = list(
    description = paste(
      "ESRP2 siRNA knockdown in prostate cancer cells",
      "(Munkley et al. 2019, eLife)"
    ),
    runs = c(
      "SRR8870780",
      "SRR8870781",
      "SRR8870782",
      "SRR8870783",
      "SRR8870784",
      "SRR8870785"
    )
  )
)

# Known ESRP targets (Warzecha 2009, 2010; Shapiro 2011; Yang 2016)
known_targets <- c(
  "ESRP1",
  "ESRP2",
  "CD44",
  "CTNND1",
  "SCRIB",
  "ENAH",
  "ARHGEF11",
  "MAP3K7",
  "FGFR2",
  "SLK",
  "EXOC7",
  "NUMB",
  "ITGA6",
  "FLNB",
  "KIF13A",
  "PLOD2",
  "CLSTN1",
  "TCF7L2",
  "MBNL1"
)

annot_dir <- prepare_isoformic_annotation(
  download_reference(version = "49", file_type = "gff")
)
tx_all <- read_parquet(fs::path(annot_dir, "annot_data_transcripts.parquet"))
genes_all <- read_parquet(fs::path(annot_dir, "annot_data_genes.parquet"))
tx_map <- tx_all |>
  dplyr::left_join(genes_all[, c("gene_id", "gene_name")], by = "gene_id")

# 1. Full analysis per series: swish, library sizes, run QC, equivalence -----
full <- list()
for (s in names(series)) {
  runs <- series[[s]]$runs
  sample_id <- c(paste0("control_", 1:3), paste0("knockdown_", 1:3))
  dirs <- stats::setNames(
    fs::path(work_dir, "quants_sub", paste0(s, "_", runs)),
    sample_id
  )
  sample_df <- data.frame(
    files = fs::path(dirs, "quant.sf"),
    names = sample_id,
    condition = factor(
      rep(c("control", "knockdown"), each = 3),
      levels = c("control", "knockdown")
    )
  )
  se <- tximeta::tximeta(
    sample_df,
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
  det <- as.data.frame(mcols(y))
  det$transcript_id <- rownames(y)
  full[[s]] <- list(
    se = se,
    dirs = dirs,
    det = det[, c("transcript_id", "log2FC", "pvalue", "qvalue")],
    library_size = colSums(assay(se, "counts"))
  )
  message(
    s,
    ": ",
    sum(det$qvalue < 0.05, na.rm = TRUE),
    " significant transcripts"
  )
}

# 2. Gene panel: known targets + genes with the most significant transcripts -
top_genes <- function(det, counts) {
  det |>
    dplyr::left_join(
      tx_map[, c("transcript_id", "gene_name")],
      by = "transcript_id"
    ) |>
    dplyr::filter(!is.na(qvalue), qvalue < 0.05, !is.na(gene_name)) |>
    dplyr::filter(
      transcript_id %in% rownames(counts)[rowMeans(counts) >= 20]
    ) |>
    dplyr::count(gene_name, sort = TRUE) |>
    dplyr::slice_head(n = n_top_genes) |>
    dplyr::pull(gene_name)
}
panel <- known_targets
for (s in names(full)) {
  panel <- union(panel, top_genes(full[[s]]$det, assay(full[[s]]$se, "counts")))
}
panel_genes <- genes_all$gene_id[genes_all$gene_name %in% panel]
panel_tx <- tx_all$transcript_id[tx_all$gene_id %in% panel_genes]
message(length(panel), " genes, ", length(panel_tx), " annotated transcripts")

# 3. Annotation subset (shared by both series) ---------------------------------
annot_out <- fs::path(extdata, "knockdown_annotation")
if (fs::dir_exists(annot_out)) {
  fs::dir_delete(annot_out)
}
fs::dir_create(annot_out)
write_parquet(
  dplyr::filter(genes_all, gene_id %in% panel_genes),
  fs::path(annot_out, "annot_data_genes.parquet")
)
write_parquet(
  dplyr::filter(tx_all, transcript_id %in% panel_tx),
  fs::path(annot_out, "annot_data_transcripts.parquet")
)
exons_all <- read_parquet(fs::path(annot_dir, "annot_data_exons.parquet"))
write_parquet(
  dplyr::filter(exons_all, transcript_id %in% panel_tx),
  fs::path(annot_out, "annot_data_exons.parquet")
)
meta_all <- read_parquet(fs::path(annot_dir, "annot_metadata.parquet"))
meta_all$value[meta_all$key == "annot_name"] <- "gencode_v49_knockdown_subset"
write_parquet(
  meta_all[!meta_all$key %in% c("source_path", "created_at"), ],
  fs::path(annot_out, "annot_metadata.parquet")
)
seq_all <- read_parquet(fs::path(annot_dir, "annot_seq_regions.parquet"))
write_parquet(
  dplyr::filter(
    seq_all,
    seqid %in% unique(tx_all$seqid[tx_all$transcript_id %in% panel_tx])
  ),
  fs::path(annot_out, "annot_seq_regions.parquet")
)

# 4. Per-series tables ---------------------------------------------------------
for (s in names(full)) {
  f <- full[[s]]
  out <- fs::path(out_root, s)
  if (fs::dir_exists(out)) {
    fs::dir_delete(out)
  }
  fs::dir_create(out)
  counts <- assay(f$se, "counts")
  use <- intersect(panel_tx, rownames(counts)[rowMeans(counts) >= 1])
  as_table <- function(m) {
    tibble::as_tibble(
      as.data.frame(round(m[use, ], 3)),
      rownames = "transcript_id"
    )
  }
  write_parquet(
    as_table(assay(f$se, "abundance")),
    fs::path(out, "abundance.parquet")
  )
  write_parquet(as_table(counts), fs::path(out, "counts.parquet"))
  write_parquet(
    dplyr::filter(f$det, transcript_id %in% use),
    fs::path(out, "det.parquet")
  )

  qc <- salmon_run_qc(f$dirs)
  col_data <- tibble::tibble(
    sample_id = names(f$dirs),
    condition = rep(c("control", "knockdown"), each = 3),
    library_size = unname(round(f$library_size)),
    run = series[[s]]$runs
  ) |>
    dplyr::left_join(qc, by = "sample_id")
  write_parquet(col_data, fs::path(out, "col_data.parquet"))

  draws <- read_inferential_replicates(f$dirs, transcripts = use, n_draws = 30)
  draws <- lapply(draws, function(m) round(m, 2))
  write_parquet(
    inferential_replicates_table(draws),
    fs::path(out, "draws.parquet")
  )

  write_parquet(
    equivalence_support(f$dirs, transcripts = use),
    fs::path(out, "equivalence.parquet")
  )
  writeLines(
    c(
      paste0("series: ", s),
      paste0("description: ", series[[s]]$description),
      "reads: first 12 million read pairs of each run",
      "reference: GENCODE v49 transcripts with the GRCh38 genome as decoy",
      "quantification: Salmon 2.8.0, 30 Gibbs draws, --dumpEq --seqBias --gcBias",
      paste0("runs: ", paste(series[[s]]$runs, collapse = ", ")),
      paste0("transcripts kept: ", length(use))
    ),
    fs::path(out, "README.txt")
  )
}
message(
  "Wrote ",
  fs::path(out_root),
  ": ",
  format(
    structure(
      sum(fs::file_size(fs::dir_ls(
        c(out_root, annot_out),
        recurse = TRUE,
        type = "file"
      ))),
      class = "object_size"
    ),
    units = "MB"
  )
)
