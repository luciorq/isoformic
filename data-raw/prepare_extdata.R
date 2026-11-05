# Regenerates the files under inst/extdata/ from the original CSV tables of
# the preeclampsia placenta dataset (PRJNA622166, GENCODE v33) and the
# GENCODE v33 GFF3 annotation. Run from the package root:
#   Rscript data-raw/prepare_extdata.R
library(dplyr)
devtools::load_all(quiet = TRUE)

extdata <- fs::path("inst", "extdata")
raw_dir <- fs::path("data-raw", "pe_csv") # original CSV tables (not tracked)

# 1. GENCODE v33 annotation directory (cached in the user cache) ---------------
gff_v33 <- download_reference(version = "33", file_type = "gff")
fasta_v33 <- download_reference(version = "33", file_type = "fasta")
annot_v33 <- prepare_isoformic_annotation(gff_v33, overwrite = TRUE)
genes_v33 <- read_parquet(fs::path(annot_v33, "annot_data_genes.parquet"))
tx_v33 <- read_parquet(fs::path(annot_v33, "annot_data_transcripts.parquet"))
exons_v33 <- read_parquet(fs::path(annot_v33, "annot_data_exons.parquet"))

# 2. Full preeclampsia tables as Parquet ---------------------------------------
abundance <- read_tabular(fs::path(raw_dir, "PE_1_counts.csv"))
names(abundance)[1] <- "transcript_id"
det <- read_tabular(fs::path(raw_dir, "DETs_fixed.csv"))
deg <- read_tabular(fs::path(raw_dir, "DEGs_PE_fixed2.csv"))
sample_ids <- setdiff(colnames(abundance), "transcript_id")
col_data <- tibble::tibble(
  sample_id = sample_ids,
  condition = c(rep("preeclampsia", 8), rep("control", length(sample_ids) - 8))
)
write_parquet(abundance, fs::path(extdata, "pe_abundance.parquet"))
write_parquet(det, fs::path(extdata, "pe_det.parquet"))
write_parquet(deg, fs::path(extdata, "pe_deg.parquet"))
write_parquet(col_data, fs::path(extdata, "pe_col_data.parquet"))

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
  filter(gene_id %in% tx_v33$gene_id) |>
  arrange(qvalue) |>
  slice_head(n = 60) |>
  pull(gene_id)
gene_ids <- unique(c(
  genes_v33$gene_id[genes_v33$gene_name %in% showcase],
  top_deg
))
tx_sel <- tx_v33 |> filter(gene_id %in% gene_ids)
example_annot <- fs::path(tempdir(), "example_annotation")
fs::dir_create(example_annot)
genes_sel <- filter(genes_v33, gene_id %in% gene_ids)
exons_sel <- filter(exons_v33, transcript_id %in% tx_sel$transcript_id)
write_parquet(genes_sel, fs::path(example_annot, "annot_data_genes.parquet"))
write_parquet(tx_sel, fs::path(example_annot, "annot_data_transcripts.parquet"))
write_parquet(exons_sel, fs::path(example_annot, "annot_data_exons.parquet"))
meta <- read_parquet(fs::path(annot_v33, "annot_metadata.parquet"))
meta$value[meta$key == "annot_name"] <- "gencode_v33_example_subset"
meta <- meta[!meta$key %in% c("source_path", "created_at"), ]
write_parquet(meta, fs::path(example_annot, "annot_metadata.parquet"))
seq_regions <- read_parquet(fs::path(annot_v33, "annot_seq_regions.parquet"))
seq_regions_sel <- filter(seq_regions, seqid %in% unique(tx_sel$seqid))
write_parquet(
  seq_regions_sel,
  fs::path(example_annot, "annot_seq_regions.parquet")
)
# Transcript sequences (for transcript_orfs(), coding_potential(), polya_sites())
extract_transcript_sequences(
  example_annot,
  transcript_fasta = fasta_v33,
  overwrite = TRUE
)

abundance_sel <- abundance |> filter(transcript_id %in% tx_sel$transcript_id)
tpm <- as.matrix(abundance_sel[, sample_ids])
rownames(tpm) <- abundance_sel$transcript_id
counts <- round(tpm * 20) # pseudo-counts for a second assay (TPM x scaling)
iso <- IsoformicExperiment(
  experiment_name = "pe_placenta_example",
  annot_path = example_annot,
  assays = list(abundance = tpm, counts = counts),
  col_data = col_data,
  metadata = list(
    source = "PRJNA622166 (Landscape of dysregulated placental RNA editing associated with preeclampsia)",
    annotation = "GENCODE v33 subset",
    note = "counts are TPM x 20 rounded, for demonstration only"
  )
)
de_tx(iso) <- filter(det, transcript_id %in% tx_sel$transcript_id)
de_gene(iso) <- filter(deg, gene_id %in% gene_ids)
example_dir <- fs::path(extdata, "example_experiment")
if (fs::dir_exists(example_dir)) {
  fs::dir_delete(example_dir)
}
write_isoformic(iso, example_dir, copy_annotation = TRUE)

# 4. Tiny GFF3 for annotation examples ----------------------------------------
tiny_genes <- genes_v33$gene_id[
  genes_v33$gene_name %in% c("RBPJ", "PNCK", "FLT1")
]
lines <- readLines(gzfile(gff_v33))
header <- lines[startsWith(lines, "#")]
header <- header[
  !startsWith(header, "##sequence-region") | grepl("chr4 |chrX |chr13 ", header)
]
pattern <- paste0("gene_id=(", paste(tiny_genes, collapse = "|"), ");")
body <- lines[grepl(pattern, lines)]
writeLines(c(header, body), fs::path(extdata, "example.gff3"))
# Matching transcript FASTA (GENCODE v33 sequences of the same transcripts)
tiny_tx <- tx_v33$transcript_id[tx_v33$gene_id %in% tiny_genes]
tiny_fasta <- fs::path(extdata, "example_transcripts.fa.gz")
con <- gzfile(tiny_fasta, "w")
stream_fasta(
  fasta_v33,
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

# 5. GENCODE v49 annotation subset matching the bundled Salmon quantifications
quant_files <- fs::dir_ls(
  fs::path(extdata, "mini_quants"),
  recurse = TRUE,
  glob = "*quant.sf"
)
quant_ids <- unique(sub("\\|.*$", "", read_tabular(quant_files[1])$Name))
gff_v49 <- download_reference(version = "49", file_type = "gff")
annot_v49 <- prepare_isoformic_annotation(gff_v49, overwrite = TRUE)
tx49 <- read_parquet(fs::path(annot_v49, "annot_data_transcripts.parquet")) |>
  filter(transcript_id %in% quant_ids)
mini_annot <- fs::path(extdata, "mini_annotation")
if (fs::dir_exists(mini_annot)) {
  fs::dir_delete(mini_annot)
}
fs::dir_create(mini_annot)
genes49 <- read_parquet(fs::path(annot_v49, "annot_data_genes.parquet"))
genes49 <- filter(genes49, gene_id %in% tx49$gene_id)
exons49 <- read_parquet(fs::path(annot_v49, "annot_data_exons.parquet"))
exons49 <- filter(exons49, transcript_id %in% tx49$transcript_id)
write_parquet(genes49, fs::path(mini_annot, "annot_data_genes.parquet"))
write_parquet(tx49, fs::path(mini_annot, "annot_data_transcripts.parquet"))
write_parquet(exons49, fs::path(mini_annot, "annot_data_exons.parquet"))
meta49 <- read_parquet(fs::path(annot_v49, "annot_metadata.parquet"))
meta49$value[meta49$key == "annot_name"] <- "gencode_v49_mini_subset"
write_parquet(
  meta49[!meta49$key %in% c("source_path", "created_at"), ],
  fs::path(mini_annot, "annot_metadata.parquet")
)
seq49 <- read_parquet(fs::path(annot_v49, "annot_seq_regions.parquet"))
seq49 <- filter(seq49, seqid %in% unique(tx49$seqid))
write_parquet(seq49, fs::path(mini_annot, "annot_seq_regions.parquet"))

# 6. Reactome gene set subset (pathways with example genes) -------------------
all_sets <- read_gmt(fs::path(
  extdata,
  "c2.cp.reactome.v2023.1.Hs.symbols.gmt.txt"
))
example_gene_names <- genes_v33$gene_name[genes_v33$gene_id %in% gene_ids]
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
