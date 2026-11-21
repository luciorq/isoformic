# Regenerates the bundled example data under inst/extdata/, and the placenta
# tables the articles read from vignettes/articles/data/preeclampsia/, from
# the Salmon quantifications of the preeclampsia placenta dataset
# (PRJNA622166) against GENCODE release 50, made by
# data-raw/placenta/quantify_placenta.sh, and from the GENCODE 50 annotation.
# Run from the package root:
#   ISOFORMIC_PLACENTA_DIR=/path/to/workdir Rscript data-raw/prepare_extdata.R
#
# Writes to vignettes/articles/data/preeclampsia/ (only the articles use
# them; they are too large for the package):
#   pe_abundance.parquet, pe_det.parquet, pe_deg.parquet, pe_col_data.parquet
#     TPM of the tested transcripts in the 25 samples, swish results for
#     every tested transcript and gene, and the sample table
# and to inst/extdata/:
#   example_experiment/  69 genes (showcase genes + top 60 genes by q-value),
#     annotation with transcript sequences copied inside
#   example.gff3.gz, example_transcripts.fa.gz  three genes, for the annotation
#     and sequence examples
#   mini_quants/, mini_draws.parquet, mini_annotation/  six samples
#     quantified against the transcripts of 15 genes, their 20 Gibbs draws as
#     a long table, and the matching annotation subset
#   reactome_subset.gmt  the 40 Reactome pathways with most example genes of
#     the complete MSigDB collection, which the articles read from
#     vignettes/articles/data/preeclampsia/
suppressPackageStartupMessages({
  library(dplyr)
  library(SummarizedExperiment)
})
devtools::load_all(quiet = TRUE)

release <- "50"
work_dir <- Sys.getenv("ISOFORMIC_PLACENTA_DIR", unset = NA)
if (is.na(work_dir) || !dir.exists(fs::path(work_dir, "placenta", "quants"))) {
  stop(
    "Set ISOFORMIC_PLACENTA_DIR to the work directory of ",
    "data-raw/placenta/quantify_placenta.sh.",
    call. = FALSE
  )
}
extdata <- fs::path("inst", "extdata")
article_data <- fs::path("vignettes", "articles", "data", "preeclampsia")
fs::dir_create(article_data)
sheet <- utils::read.delim(fs::path(
  "data-raw",
  "placenta",
  "placenta_samples.tsv"
))

# 1. GENCODE annotation (cached in the user cache) ------------------------------
gff <- download_reference(version = release, file_type = "gff")
fasta <- download_reference(version = release, file_type = "fasta")
annot <- prepare_isoformic_annotation(gff)
genes_all <- read_parquet(fs::path(annot, "annot_data_genes.parquet"))
tx_all <- read_parquet(fs::path(annot, "annot_data_transcripts.parquet"))
exons_all <- read_parquet(fs::path(annot, "annot_data_exons.parquet"))
meta_all <- read_parquet(fs::path(annot, "annot_metadata.parquet"))
seq_all <- read_parquet(fs::path(annot, "annot_seq_regions.parquet"))
# Genes on the reference chromosomes, without the PAR copies on chrY: gene
# names also match the copies on patches and alternate haplotypes
primary_genes <- genes_all$gene_id[
  grepl("^chr([0-9]+|X|Y|M)$", genes_all$seqid) &
    !grepl("_PAR_Y$", genes_all$gene_id)
]

# Writes an annotation subset to `dir`
write_annotation_subset <- function(dir, tx_ids, name) {
  if (fs::dir_exists(dir)) {
    fs::dir_delete(dir)
  }
  fs::dir_create(dir)
  tx <- filter(tx_all, transcript_id %in% tx_ids)
  write_parquet(
    filter(genes_all, gene_id %in% tx$gene_id),
    fs::path(dir, "annot_data_genes.parquet")
  )
  write_parquet(tx, fs::path(dir, "annot_data_transcripts.parquet"))
  write_parquet(
    filter(exons_all, transcript_id %in% tx_ids),
    fs::path(dir, "annot_data_exons.parquet")
  )
  meta <- meta_all
  meta$value[meta$key == "annot_name"] <- name
  write_parquet(
    meta[!meta$key %in% c("source_path", "created_at"), ],
    fs::path(dir, "annot_metadata.parquet")
  )
  write_parquet(
    filter(seq_all, seqid %in% unique(tx$seqid)),
    fs::path(dir, "annot_seq_regions.parquet")
  )
  invisible(dir)
}

# 2. Complete placenta quantifications and swish -------------------------------
col_data <- tibble::tibble(
  sample_id = sheet$run,
  condition = factor(sheet$condition, levels = c("control", "preeclampsia"))
)
quant_dirs <- fs::path(work_dir, "placenta", "quants", sheet$run)
se <- tximeta::tximeta(
  data.frame(
    files = fs::path(quant_dirs, "quant.sf"),
    names = sheet$run,
    condition = col_data$condition
  ),
  type = "salmon",
  txOut = TRUE,
  skipMeta = TRUE,
  useHub = FALSE,
  ignoreAfterBar = TRUE,
  skipSeqinfo = TRUE,
  dropInfReps = FALSE
)
inf_names <- grep("^infRep", assayNames(se), value = TRUE)

# Transcripts: swish on the 30 Gibbs draws
y_tx <- fishpond::scaleInfReps(se, quiet = TRUE)
y_tx <- fishpond::labelKeep(y_tx, minCount = 10, minN = 3)
y_tx <- y_tx[mcols(y_tx)$keep, ]
set.seed(20261113) # swish computes its p-values from permutations
y_tx <- fishpond::swish(y_tx, x = "condition", quiet = TRUE)
det <- tibble::tibble(
  transcript_id = rownames(y_tx),
  log2FC = mcols(y_tx)$log2FC,
  pvalue = mcols(y_tx)$pvalue,
  qvalue = mcols(y_tx)$qvalue
)

# Genes: counts and draws summed over the transcripts of each gene, length
# weighted by abundance (as tximport does), then swish
tx2gene <- tx_all$gene_id[match(rownames(se), tx_all$transcript_id)]
known <- !is.na(tx2gene)
gene_sum <- function(m) rowsum(m[known, , drop = FALSE], tx2gene[known])
abund_g <- gene_sum(SummarizedExperiment::assay(se, "abundance"))
len_g <- gene_sum(
  SummarizedExperiment::assay(se, "abundance") *
    SummarizedExperiment::assay(se, "length")
) /
  abund_g
len_g[!is.finite(len_g)] <- NA
mean_len <- rowMeans(len_g, na.rm = TRUE)
for (j in seq_len(ncol(len_g))) {
  miss <- is.na(len_g[, j])
  len_g[miss, j] <- mean_len[miss]
}
# Genes without abundance in any sample have no length: one NA length makes
# every scaled draw NA in fishpond::scaleInfReps()
len_g[!is.finite(len_g)] <- stats::median(len_g, na.rm = TRUE)
gene_assays <- c(
  list(
    counts = gene_sum(SummarizedExperiment::assay(se, "counts")),
    abundance = abund_g,
    length = len_g
  ),
  lapply(stats::setNames(inf_names, inf_names), function(a) {
    gene_sum(SummarizedExperiment::assay(se, a))
  })
)
se_gene <- SummarizedExperiment(
  assays = gene_assays,
  colData = SummarizedExperiment::colData(se)
)
y_gene <- fishpond::scaleInfReps(se_gene, quiet = TRUE)
y_gene <- fishpond::labelKeep(y_gene, minCount = 10, minN = 3)
y_gene <- y_gene[mcols(y_gene)$keep, ]
set.seed(20261113)
y_gene <- fishpond::swish(y_gene, x = "condition", quiet = TRUE)
deg <- tibble::tibble(
  gene_id = rownames(y_gene),
  log10mean = mcols(y_gene)$log10mean,
  stat = mcols(y_gene)$stat,
  log2FC = mcols(y_gene)$log2FC,
  pvalue = mcols(y_gene)$pvalue,
  locfdr = mcols(y_gene)$locfdr,
  qvalue = mcols(y_gene)$qvalue
)

# Only the tested transcripts (the rows of pe_det): the whole annotation would
# make the table five times larger, mostly with transcripts without reads
abundance <- tibble::as_tibble(
  as.data.frame(SummarizedExperiment::assay(se, "abundance")),
  rownames = "transcript_id"
)
abundance <- abundance[abundance$transcript_id %in% det$transcript_id, ]
# TPM as 32-bit floats (seven significant digits; read back as double): half
# the size of doubles
write_parquet_float <- function(x, path) {
  con <- DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE))
  duckdb::duckdb_register(con, "x", x)
  num <- names(x)[vapply(x, is.double, logical(1))]
  cols <- ifelse(
    names(x) %in% num,
    paste0("CAST(\"", names(x), "\" AS FLOAT) AS \"", names(x), "\""),
    paste0("\"", names(x), "\"")
  )
  DBI::dbExecute(
    con,
    paste0(
      "COPY (SELECT ",
      paste(cols, collapse = ", "),
      " FROM x) TO '",
      path,
      "' (FORMAT parquet, COMPRESSION zstd)"
    )
  )
  invisible(path)
}
write_parquet_float(abundance, fs::path(article_data, "pe_abundance.parquet"))
write_parquet(det, fs::path(article_data, "pe_det.parquet"))
write_parquet(deg, fs::path(article_data, "pe_deg.parquet"))
write_parquet(
  dplyr::mutate(col_data, condition = as.character(.data$condition)),
  fs::path(article_data, "pe_col_data.parquet")
)

# 3. Example experiment: showcase genes + top differentially expressed genes ---
showcase <- c(
  "RBPJ",
  "EGFR",
  "PNCK",
  "FLT1",
  "ENG",
  "IL2RA",
  "XIST",
  "TSIX",
  "MTG2"
)
top_deg <- deg |>
  filter(gene_id %in% tx_all$gene_id) |>
  arrange(qvalue, pvalue) |>
  slice_head(n = 60) |>
  pull(gene_id)
gene_ids <- unique(c(
  genes_all$gene_id[
    genes_all$gene_name %in% showcase & genes_all$gene_id %in% primary_genes
  ],
  top_deg
))
tx_sel <- tx_all$transcript_id[tx_all$gene_id %in% gene_ids]
example_annot <- write_annotation_subset(
  fs::path(tempdir(), "example_annotation"),
  tx_sel,
  paste0("gencode_v", release, "_example_subset")
)
# Transcript sequences (for transcript_orfs(), coding_potential(), polya_sites())
extract_transcript_sequences(
  example_annot,
  transcript_fasta = fasta,
  overwrite = TRUE
)
in_example <- rownames(se) %in% tx_sel
iso <- IsoformicExperiment(
  experiment_name = "pe_placenta_example",
  annot_path = example_annot,
  assays = list(
    abundance = SummarizedExperiment::assay(se, "abundance")[in_example, ],
    counts = SummarizedExperiment::assay(se, "counts")[in_example, ]
  ),
  col_data = dplyr::mutate(col_data, condition = as.character(.data$condition)),
  metadata = list(
    source = paste(
      "PRJNA622166 (Landscape of dysregulated placental RNA editing",
      "associated with preeclampsia), 8 preeclampsia and 17 control placentas"
    ),
    annotation = paste0("GENCODE v", release, " subset"),
    quantification = paste(
      "Salmon 2.8.0 against GENCODE v50 transcripts with the GRCh38 genome",
      "as decoy, 30 Gibbs draws, --seqBias --gcBias; swish on the transcripts",
      "and on the genes (draws summed per gene)"
    )
  )
)
de_tx(iso) <- filter(det, transcript_id %in% tx_sel)
de_gene(iso) <- filter(deg, gene_id %in% gene_ids)
example_dir <- fs::path(extdata, "example_experiment")
if (fs::dir_exists(example_dir)) {
  fs::dir_delete(example_dir)
}
write_isoformic(iso, example_dir, copy_annotation = TRUE)

# 4. Tiny GFF3 for annotation examples ----------------------------------------
tiny_genes <- genes_all$gene_id[
  genes_all$gene_name %in%
    c("RBPJ", "PNCK", "FLT1") &
    genes_all$gene_id %in% primary_genes
]
lines <- readLines(gzfile(gff))
header <- lines[startsWith(lines, "#")]
header <- header[
  !startsWith(header, "##sequence-region") | grepl("chr4 |chrX |chr13 ", header)
]
pattern <- paste0("gene_id=(", paste(tiny_genes, collapse = "|"), ");")
body <- lines[grepl(pattern, lines)]
con <- gzfile(fs::path(extdata, "example.gff3.gz"), "w", compression = 9)
writeLines(c(header, body), con)
close(con)
# Matching transcript FASTA (sequences of the same transcripts)
tiny_tx <- tx_all$transcript_id[tx_all$gene_id %in% tiny_genes]
tiny_fasta <- fs::path(extdata, "example_transcripts.fa.gz")
con <- gzfile(tiny_fasta, "w")
stream_fasta(
  fasta,
  keep = function(name) {
    strsplit(name, "|", fixed = TRUE)[[1L]][[1L]] %in% tiny_tx
  },
  on_record = function(name, sequence) {
    writeLines(
      c(
        paste0(">", name),
        substring(
          sequence,
          seq(1, nchar(sequence), by = 60),
          pmin(seq(60, nchar(sequence) + 59, by = 60), nchar(sequence))
        )
      ),
      con
    )
  }
)
close(con)

# 5. Mini quantifications, their Gibbs draws and annotation --------------------
# The Salmon directories are copied without the draws (aux_info/bootstrap),
# which go to mini_draws.parquet, without the aux_info/*.gz files, and with
# the paths of the work directory
# replaced by /work.
mini_src <- fs::dir_ls(
  fs::path(work_dir, "placenta", "mini_quants"),
  type = "directory",
  regexp = "_quant$"
)
names(mini_src) <- sub("_quant$", "", fs::path_file(mini_src))
draws <- read_inferential_replicates(mini_src)
write_parquet(
  inferential_replicates_table(draws),
  fs::path(extdata, "mini_draws.parquet")
)
mini_dir <- fs::path(extdata, "mini_quants")
if (fs::dir_exists(mini_dir)) {
  fs::dir_delete(mini_dir)
}
fs::dir_create(mini_dir)
for (src in mini_src) {
  dest <- fs::path(mini_dir, fs::path_file(src))
  fs::dir_copy(src, dest)
  fs::dir_delete(fs::path(dest, "aux_info", "bootstrap"))
  # The bias and fragment length files (*.gz) are not read by the package
  # and are ignored by git: without them a local build and a build from a
  # clean clone ship the same files
  fs::file_delete(fs::dir_ls(fs::path(dest, "aux_info"), glob = "*.gz"))
  for (f in fs::dir_ls(dest, recurse = TRUE, type = "file")) {
    if (grepl("[.](json|log|txt)$", f)) {
      txt <- readLines(f, warn = FALSE)
      writeLines(gsub(work_dir, "/work", txt, fixed = TRUE), f)
    }
  }
}
quant_ids <- unique(sub(
  "[|].*$",
  "",
  read_tabular(
    fs::path(mini_dir, "SRR11498039_quant", "quant.sf")
  )$Name
))
write_annotation_subset(
  fs::path(extdata, "mini_annotation"),
  quant_ids,
  paste0("gencode_v", release, "_mini_subset")
)

# 6. Reactome gene set subset (pathways with example genes) -------------------
all_sets <- read_gmt(fs::path(
  article_data,
  "c2.cp.reactome.v2023.1.Hs.symbols.gmt.txt"
))
example_gene_names <- genes_all$gene_name[genes_all$gene_id %in% gene_ids]
hits <- vapply(all_sets, function(s) sum(s %in% example_gene_names), integer(1))
subset <- all_sets[order(-hits)][seq_len(40)]
writeLines(
  vapply(
    names(subset),
    function(n) paste(c(n, "na", subset[[n]]), collapse = "\t"),
    character(1)
  ),
  fs::path(extdata, "reactome_subset.gmt")
)
cat("done\n")
