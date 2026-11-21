# Case study experiments for the articles
#
# Builds vignettes/articles/data/<study>/, one self-contained
# IsoformicExperiment directory per case study (write_isoformic() with the
# annotation copied inside), from the quantifications made by
# quantify_mouse_bnst.sh and quantify_prostate_ont.sh:
#
# * mouse_bnst: male adBNST, chronic stress vs control (PRJNA984605,
#   Degroat et al. 2024), Salmon against GENCODE vM39.
# * prostate_ont: 22Rv1 cells, lorecivivint vs DMSO, Oxford Nanopore cDNA
#   (GSE259334, Van Goubergen et al. 2024), Oarfish against GENCODE v50.
#
# The differential expression is computed on the complete quantifications:
# swish (fishpond) on the 30 inferential replicates for the transcripts and
# DESeq2 on the summed counts for the genes. The bundles keep the complete
# result tables (every tested transcript and gene, so the enrichment and the
# isoform switch summaries are those of the full analysis), but the assays
# and the exon structures only of a gene panel: the genes discussed in the
# articles plus those with the most significant transcripts. The articles are
# excluded from the package build, so the bundles are not installed.
#
# Usage:
#   ISOFORMIC_CASE_STUDIES_DIR=/path/to/workdir \
#     Rscript data-raw/case_studies/build_case_studies.R [mouse_bnst] [prostate_ont]
#
# The gene-level DESeq2 and apeglm results (dea/deg.parquet) are not
# bit-reproducible: they change at about 1e-5 with the number of BLAS/OpenMP
# threads (OPENBLAS_NUM_THREADS, OMP_NUM_THREADS, read when R starts) and the
# CPU. The calls and every number the articles state are unaffected, so a
# rebuild that changes only those digits need not be committed.
suppressPackageStartupMessages(library(isoformic))

work_dir <- Sys.getenv("ISOFORMIC_CASE_STUDIES_DIR", unset = NA)
if (is.na(work_dir) || !dir.exists(work_dir)) {
  stop(
    "Set ISOFORMIC_CASE_STUDIES_DIR to the work directory of the ",
    "quantify_*.sh scripts.",
    call. = FALSE
  )
}
studies <- commandArgs(trailingOnly = TRUE)
if (!length(studies)) {
  studies <- c("mouse_bnst", "prostate_ont")
}
out_root <- fs::path("vignettes", "articles", "data")
script_dir <- fs::path("data-raw", "case_studies")
n_top_genes <- 40L
n_draws <- 30L

study_config <- list(
  mouse_bnst = list(
    title = "Chronic stress in the male mouse adBNST (PRJNA984605)",
    reference = list(version = "M39", organism = "mouse"),
    reads = "Illumina paired end, all reads",
    quantification = paste(
      "Salmon 2.11.0 against GENCODE vM39 transcripts with the GRCm39",
      "genome as decoy, 30 Gibbs draws, --posBias --seqBias --gcBias"
    ),
    levels = c("control", "stressed"),
    # Genes discussed in the article (the first three come from the
    # original analysis of the dataset by the isoformic authors)
    genes = c("Trip12", "Meg3", "Map4", "Fkbp5", "Sgk1", "Npy", "Crh"),
    min_count = 10,
    min_n = 5,
    # The samples fall in two batches that dominate the first principal
    # component of the gene-level counts (run accessions SRR24948925-35 and
    # SRR24948940-55; the batch is not annotated in GEO). The batches are
    # balanced between conditions, so both tests adjust for them.
    covariate = "batch",
    # MSigDB mouse M2 canonical pathways from Reactome (CC BY 4.0)
    genesets = paste0(
      "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2025.1.Mm/",
      "m2.cp.reactome.v2025.1.Mm.symbols.gmt"
    )
  ),
  prostate_ont = list(
    title = "Lorecivivint in 22Rv1 prostate cancer cells, nanopore cDNA (GSE259334)",
    reference = list(version = "50", organism = "human"),
    reads = "Oxford Nanopore PCR-cDNA (SQK-PCS109), PromethION, all reads",
    quantification = paste(
      "Oarfish 0.11.0 (--seq-tech ont-cdna) against GENCODE v50",
      "transcripts, 30 bootstrap replicates"
    ),
    levels = c("vehicle", "lorecivivint"),
    genes = c("AR", "CLK1", "CLK2", "CLK3", "CLK4", "SRSF9", "SRSF3", "SCN8A"),
    min_count = 10,
    min_n = 3,
    # MSigDB human hallmark gene sets (CC BY 4.0)
    genesets = paste0(
      "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2025.1.Hs/",
      "h.all.v2025.1.Hs.symbols.gmt"
    )
  )
)

# Paths of the per-sample quantifications and the sample table
study_samples <- function(study) {
  if (identical(study, "mouse_bnst")) {
    sheet <- utils::read.delim(fs::path(script_dir, "mouse_bnst_samples.tsv"))
    sheet <- sheet[sheet$sex == "male", ]
    sheet$sample_id <- paste0(sheet$condition, "_", sheet$replicate)
    run_number <- as.integer(sub("^SRR", "", sheet$run))
    sheet$batch <- ifelse(run_number <= 24948935L, "A", "B")
    paths <- fs::path(work_dir, "mouse_bnst", "quants", sheet$run)
  } else {
    sheet <- utils::read.delim(fs::path(script_dir, "prostate_ont_samples.tsv"))
    sheet$sample_id <- paste0(sheet$condition, "_", sheet$replicate)
    paths <- fs::path(work_dir, "prostate_ont", "quants", sheet$run)
  }
  levels <- study_config[[study]]$levels
  sheet$condition <- factor(sheet$condition, levels = levels)
  sheet <- sheet[order(sheet$condition, sheet$replicate), ]
  paths <- paths[match(sheet$run, basename(paths))]
  list(
    paths = stats::setNames(as.character(paths), sheet$sample_id),
    col_data = tibble::as_tibble(sheet[, intersect(
      c("sample_id", "condition", "batch", "replicate", "run"),
      colnames(sheet)
    )])
  )
}

# Transcript-level swish on the inferential replicates
run_swish <- function(iso, paths, cfg) {
  se <- as_summarized_experiment(iso, level = "transcript")
  draws <- read_inferential_replicates(paths, n_draws = n_draws)
  keep <- intersect(rownames(se), rownames(draws[[1L]]))
  se <- se[keep, ]
  draws <- lapply(draws, function(m) m[keep, colnames(se), drop = FALSE])
  names(draws) <- paste0("infRep", seq_along(draws))
  SummarizedExperiment::assays(se) <- c(
    as.list(SummarizedExperiment::assays(se)),
    draws
  )
  se$condition <- factor(se$condition, levels = cfg$levels)
  if (!is.null(cfg$covariate)) {
    se[[cfg$covariate]] <- factor(se[[cfg$covariate]])
  }
  y <- fishpond::scaleInfReps(se, quiet = TRUE)
  y <- fishpond::labelKeep(y, minCount = cfg$min_count, minN = cfg$min_n)
  y <- y[SummarizedExperiment::mcols(y)$keep, ]
  # swish computes its p-values from permutations
  set.seed(20261109)
  # With a covariate, swish stratifies the permutations by its levels
  y <- fishpond::swish(y, x = "condition", cov = cfg$covariate, quiet = TRUE)
  SummarizedExperiment::mcols(y) |>
    as.data.frame() |>
    dplyr::mutate(transcript_id = rownames(y)) |>
    dplyr::select("transcript_id", "log2FC", "pvalue", "qvalue") |>
    tibble::as_tibble()
}

# Gene-level DESeq2 on the summed counts
run_deseq2 <- function(iso, cfg) {
  se <- as_summarized_experiment(iso, level = "gene")
  # Only `counts` is used: Oarfish lengths are not effective lengths
  SummarizedExperiment::assays(se) <- list(
    counts = round(SummarizedExperiment::assay(se, "counts"))
  )
  se$condition <- factor(se$condition, levels = cfg$levels)
  design <- stats::reformulate(c(cfg$covariate, "condition"))
  if (!is.null(cfg$covariate)) {
    se[[cfg$covariate]] <- factor(se[[cfg$covariate]])
  }
  dds <- DESeq2::DESeqDataSet(se, design = design)
  dds <- dds[rowSums(DESeq2::counts(dds) >= cfg$min_count) >= cfg$min_n, ]
  dds <- DESeq2::DESeq(dds, quiet = TRUE)
  coef <- paste0("condition_", cfg$levels[2], "_vs_", cfg$levels[1])
  res <- DESeq2::lfcShrink(dds, coef = coef, type = "apeglm", quiet = TRUE)
  deg <- as.data.frame(res) |>
    tibble::rownames_to_column("gene_id") |>
    tibble::as_tibble()
  # Principal components of the variance-stabilised counts (2000 most
  # variable genes), saved for the batch check of the article
  vst <- SummarizedExperiment::assay(DESeq2::vst(dds, blind = TRUE))
  top <- order(-apply(vst, 1L, stats::var))[seq_len(min(2000L, nrow(vst)))]
  pca <- stats::prcomp(t(vst[top, ]))
  var_explained <- pca$sdev^2 / sum(pca$sdev^2)
  pcs <- tibble::as_tibble(pca$x[, 1:4], rownames = "sample_id")
  attr(deg, "pca") <- list(pcs = pcs, var_explained = var_explained[1:4])
  deg
}

# Annotation directory restricted to the tested genes; exons for the panel
subset_annotation <- function(annot_dir, gene_ids, panel_gene_ids, out) {
  if (fs::dir_exists(out)) {
    fs::dir_delete(out)
  }
  fs::dir_create(out)
  genes <- read_parquet(fs::path(annot_dir, "annot_data_genes.parquet"))
  tx <- read_parquet(fs::path(annot_dir, "annot_data_transcripts.parquet"))
  exons <- read_parquet(fs::path(annot_dir, "annot_data_exons.parquet"))
  genes <- dplyr::filter(genes, .data$gene_id %in% gene_ids)
  tx <- dplyr::filter(tx, .data$gene_id %in% gene_ids)
  panel_tx <- tx$transcript_id[tx$gene_id %in% panel_gene_ids]
  write_parquet(genes, fs::path(out, "annot_data_genes.parquet"))
  write_parquet(tx, fs::path(out, "annot_data_transcripts.parquet"))
  write_parquet(
    dplyr::filter(exons, .data$transcript_id %in% panel_tx),
    fs::path(out, "annot_data_exons.parquet")
  )
  meta <- read_parquet(fs::path(annot_dir, "annot_metadata.parquet"))
  meta <- meta[!meta$key %in% c("source_path", "created_at"), ]
  meta$value[meta$key == "annot_name"] <- paste0(
    meta$value[meta$key == "annot_name"],
    "_case_study_subset"
  )
  write_parquet(meta, fs::path(out, "annot_metadata.parquet"))
  seqs <- read_parquet(fs::path(annot_dir, "annot_seq_regions.parquet"))
  write_parquet(
    dplyr::filter(seqs, .data$seqid %in% unique(tx$seqid)),
    fs::path(out, "annot_seq_regions.parquet")
  )
  invisible(out)
}

# Gene sets restricted to those with 15 to 500 tested genes, written as GMT
write_genesets <- function(url, tested_genes, genes_tbl, out) {
  gmt <- fs::path(tempdir(), basename(url))
  utils::download.file(url, gmt, quiet = TRUE, mode = "wb")
  sets <- read_gmt(gmt)
  tested <- unique(genes_tbl$gene_name[genes_tbl$gene_id %in% tested_genes])
  n <- vapply(sets, function(g) sum(g %in% tested), integer(1))
  sets <- sets[n >= 15L & n <= 500L]
  writeLines(
    vapply(
      names(sets),
      function(nm) paste(c(nm, basename(url), sets[[nm]]), collapse = "\t"),
      character(1)
    ),
    out
  )
  message(length(sets), " gene sets written to ", out)
}

build_study <- function(study) {
  cfg <- study_config[[study]]
  message("== ", study)
  samples <- study_samples(study)
  annot_dir <- prepare_isoformic_annotation(
    download_reference(
      version = cfg$reference$version,
      organism = cfg$reference$organism,
      file_type = "gff"
    )
  )
  reader <- if (identical(study, "mouse_bnst")) read_salmon else read_oarfish
  iso <- reader(
    samples$paths,
    annot_path = annot_dir,
    col_data = samples$col_data,
    experiment_name = study
  )
  det <- run_swish(iso, samples$paths, cfg)
  deg <- run_deseq2(iso, cfg)
  pca <- attr(deg, "pca")
  attr(deg, "pca") <- NULL
  message(
    nrow(det),
    " transcripts tested, ",
    sum(det$qvalue < 0.05, na.rm = TRUE),
    " with q < 0.05; ",
    nrow(deg),
    " genes tested, ",
    sum(deg$padj < 0.05, na.rm = TRUE),
    " with padj < 0.05"
  )

  # Gene panel: genes of the article + the genes with the most significant
  # transcripts (smallest q-value, then largest fold change), preferring
  # those whose gene-level result is not significant (isoform switches)
  tx_map <- annot_data_transcripts(iso)[, c("transcript_id", "gene_id")]
  genes_tbl <- annot_data_genes(iso)[, c("gene_id", "gene_name")]
  sig_deg <- deg$gene_id[!is.na(deg$padj) & deg$padj < 0.05]
  top <- det |>
    dplyr::filter(!is.na(.data$qvalue), .data$qvalue < 0.05) |>
    dplyr::left_join(tx_map, by = "transcript_id") |>
    dplyr::summarise(
      min_q = min(.data$qvalue),
      max_lfc = max(abs(.data$log2FC)),
      .by = "gene_id"
    ) |>
    dplyr::mutate(switch = !.data$gene_id %in% sig_deg) |>
    dplyr::arrange(
      dplyr::desc(.data$switch),
      .data$min_q,
      dplyr::desc(.data$max_lfc)
    ) |>
    dplyr::slice_head(n = n_top_genes) |>
    dplyr::pull("gene_id")
  panel <- union(genes_tbl$gene_id[genes_tbl$gene_name %in% cfg$genes], top)
  tested_genes <- union(
    deg$gene_id,
    tx_map$gene_id[tx_map$transcript_id %in% det$transcript_id]
  )
  message(
    length(panel),
    " panel genes, ",
    length(tested_genes),
    " tested genes"
  )

  out <- fs::path(out_root, study)
  annot_tmp <- fs::path(tempdir(), paste0(study, "_annotation"))
  subset_annotation(annot_dir, tested_genes, panel, annot_tmp)
  panel_tx <- tx_map$transcript_id[tx_map$gene_id %in% panel]
  keep_tx <- intersect(rownames(tx_assay(iso, "counts")), panel_tx)
  assays_small <- lapply(
    isoformic::tx_assays(iso)[c("counts", "abundance")],
    function(m) round(m[keep_tx, , drop = FALSE], 2)
  )
  small <- IsoformicExperiment(
    experiment_name = study,
    annot_path = annot_tmp,
    assays = assays_small,
    col_data = samples$col_data,
    metadata = list(
      title = cfg$title,
      reads = cfg$reads,
      quantification = cfg$quantification,
      transcript_de = paste0(
        "swish (fishpond), 30 inferential replicates",
        if (!is.null(cfg$covariate)) paste0(", stratified by ", cfg$covariate)
      ),
      gene_de = paste0(
        "DESeq2, apeglm-shrunken log2 fold changes",
        if (!is.null(cfg$covariate)) {
          paste0(", design ~ ", cfg$covariate, " + condition")
        }
      ),
      panel = "assays and exons of the gene panel only"
    )
  )
  de_tx(small) <- det
  de_gene(small) <- deg
  if (fs::dir_exists(out)) {
    fs::dir_delete(out)
  }
  write_isoformic(small, out, copy_annotation = TRUE)
  write_parquet(
    tibble::add_column(
      pca$pcs,
      !!!stats::setNames(
        as.list(round(pca$var_explained, 4)),
        paste0("var_PC", 1:4)
      )
    ),
    fs::path(out, "sample_pca.parquet")
  )
  if (!is.null(cfg$genesets)) {
    write_genesets(
      cfg$genesets,
      tested_genes,
      genes_tbl,
      fs::path(out, "genesets.gmt")
    )
  }
  size <- sum(fs::file_size(fs::dir_ls(out, recurse = TRUE, type = "file")))
  message(
    "Wrote ",
    out,
    " (",
    format(structure(size, class = "object_size"), units = "MB"),
    ")"
  )
}

for (study in studies) {
  build_study(study)
}
