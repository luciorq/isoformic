# AlphaGenome pilot benchmark, step 1: choose the genes and write the targets
#
# Question of the pilot: do sequence-based splicing predictions of AlphaGenome
# (splice-site probability and usage, junction counts, RNA-seq coverage) tell
# productive from unproductive isoforms beyond what the structure rules of
# transcript_productivity() already give?
#
# This step needs no AlphaGenome access. It writes, under
# <work_dir>/targets/:
#   genes.parquet        one row per gene: the 1 Mb window to predict (0-based,
#                        half-open, as AlphaGenome intervals) and the stratum
#   transcripts.parquet  transcripts of those genes with the truth labels
#                        (biotype, Karousis 2021 validated NMD pairs) and the
#                        calls of the rules (use_biotype = FALSE,
#                        use_evidence = FALSE)
#   exons.parquet        exons of those transcripts (1-based, inclusive)
#   introns.parquet      transcript x intron table (1-based, inclusive)
#   sites.parquet        unique splice sites per gene (first and last intron
#                        base, 1-based)
#   regions.parquet      unique exons and introns per gene, to average the
#                        predicted RNA-seq coverage over
#
# Usage (from the package root):
#   ISOFORMIC_ALPHAGENOME_DIR=/path/to/workdir Rscript data-raw/alphagenome/01_targets.R
# Optional environment variables:
#   ISOFORMIC_ALPHAGENOME_GENCODE  GENCODE release (default: the package default)
#   ISOFORMIC_ALPHAGENOME_N        genes per stratum (default 250)
#   ISOFORMIC_ALPHAGENOME_SEED     sampling seed (default 1)

devtools::load_all(".", quiet = TRUE)

work_dir <- Sys.getenv("ISOFORMIC_ALPHAGENOME_DIR", unset = NA)
if (is.na(work_dir) || !nzchar(work_dir)) {
  stop("Set ISOFORMIC_ALPHAGENOME_DIR to a work directory.", call. = FALSE)
}
release <- Sys.getenv("ISOFORMIC_ALPHAGENOME_GENCODE", unset = "")
release <- if (nzchar(release)) release else NULL
n_stratum <- as.integer(Sys.getenv("ISOFORMIC_ALPHAGENOME_N", unset = "250"))
seed <- as.integer(Sys.getenv("ISOFORMIC_ALPHAGENOME_SEED", unset = "1"))
out_dir <- fs::dir_create(fs::path(work_dir, "targets"))

window <- 1048576L # AlphaGenome SEQUENCE_LENGTH_1MB
max_gene_length <- 900000L
chromosomes <- paste0("chr", c(1:22, "X"))
unversioned <- function(x) sub("\\..*$", "", x)

# 1. Annotation with transcript sequences ------------------------------------
gff <- download_reference(version = release, file_type = "gff")
annot_dir <- prepare_isoformic_annotation(gff)
if (!has_sequences(annot_dir)) {
  genome <- download_reference(
    version = release,
    file_type = "genome_fasta",
    region = "primary"
  )
  extract_transcript_sequences(annot_dir, genome_fasta = genome)
}
tx <- read_parquet(fs::path(annot_dir, "annot_data_transcripts.parquet"))
genes <- read_parquet(fs::path(annot_dir, "annot_data_genes.parquet"))
exons <- read_parquet(fs::path(annot_dir, "annot_data_exons.parquet"))

# 2. Truth independent of the annotation: Karousis 2021 validated pairs ------
nmd <- read_frozen_data("nmd_exons")
tx_key <- unversioned(tx$transcript_id)
pairs <- nmd |>
  dplyr::transmute(
    nmd_exon = paste0(.data$seqid, ":", .data$start_pos, "-", .data$end_pos),
    nmd_tx = strsplit(.data$nmd_transcript_ids, "[,;] *"),
    non_nmd_tx = strsplit(.data$non_nmd_transcript_ids, "[,;] *")
  ) |>
  tidyr::unnest("nmd_tx") |>
  tidyr::unnest("non_nmd_tx") |>
  dplyr::mutate(
    nmd_tx = tx$transcript_id[match(unversioned(.data$nmd_tx), tx_key)],
    non_nmd_tx = tx$transcript_id[match(unversioned(.data$non_nmd_tx), tx_key)]
  ) |>
  dplyr::filter(!is.na(.data$nmd_tx), !is.na(.data$non_nmd_tx)) |>
  dplyr::distinct()
pairs$gene_id <- tx$gene_id[match(pairs$nmd_tx, tx$transcript_id)]
pairs <- pairs[
  pairs$gene_id == tx$gene_id[match(pairs$non_nmd_tx, tx$transcript_id)],
]

# 3. Gene strata --------------------------------------------------------------
per_gene <- tx |>
  dplyr::summarise(
    n_tx = dplyr::n(),
    n_coding = sum(.data$transcript_type == "protein_coding"),
    n_nmd = sum(.data$transcript_type == "nonsense_mediated_decay"),
    n_ri = sum(.data$transcript_type == "retained_intron"),
    has_canonical = any(.data$is_canonical %in% TRUE),
    .by = "gene_id"
  )
eligible <- genes |>
  dplyr::filter(
    .data$gene_type == "protein_coding",
    .data$seqid %in% chromosomes,
    is.na(.data$haplotype),
    .data$end_pos - .data$start_pos + 1L <= max_gene_length
  ) |>
  dplyr::inner_join(per_gene, by = "gene_id") |>
  dplyr::filter(.data$has_canonical, .data$n_coding >= 2L)

set.seed(seed)
take <- function(ids, n) sort(sample(ids, min(n, length(ids))))
karousis_genes <- take(intersect(eligible$gene_id, pairs$gene_id), n_stratum)
biotype_pool <- eligible$gene_id[eligible$n_nmd >= 1L & eligible$n_ri >= 1L]
biotype_genes <- take(setdiff(biotype_pool, karousis_genes), n_stratum)
sel <- eligible |>
  dplyr::filter(.data$gene_id %in% c(karousis_genes, biotype_genes)) |>
  dplyr::mutate(
    stratum = dplyr::if_else(
      .data$gene_id %in% karousis_genes,
      "karousis",
      "biotype"
    )
  )

# 4. Windows: 1 Mb centred on the gene, kept inside the chromosome -----------
sizes <- read_parquet(fs::path(annot_dir, "annot_seq_regions.parquet")) |>
  dplyr::transmute(.data$seqid, chrom_length = .data$end_pos)
stopifnot(all(chromosomes %in% sizes$seqid))
sel <- dplyr::left_join(sel, sizes, by = "seqid")
mid <- (sel$start_pos + sel$end_pos) %/% 2L
w_start <- pmax(0L, mid - window %/% 2L)
w_start <- pmin(w_start, sel$chrom_length - window)
gene_tbl <- sel |>
  dplyr::transmute(
    .data$gene_id,
    .data$gene_name,
    .data$seqid,
    .data$strand,
    .data$start_pos,
    .data$end_pos,
    window_start = w_start,
    window_end = w_start + window,
    .data$stratum,
    .data$n_tx
  )
stopifnot(
  all(gene_tbl$window_start < gene_tbl$start_pos),
  all(gene_tbl$window_end >= gene_tbl$end_pos)
)

# 5. Transcripts, labels and the calls of the rules ---------------------------
sel_tx <- tx[tx$gene_id %in% gene_tbl$gene_id, ]
iso <- IsoformicExperiment(
  experiment_name = "alphagenome_pilot",
  annot_path = annot_dir,
  assays = list(
    abundance = matrix(
      1,
      nrow = nrow(sel_tx),
      ncol = 1L,
      dimnames = list(sel_tx$transcript_id, "s1")
    )
  ),
  col_data = data.frame(sample_id = "s1")
)
rules <- transcript_productivity(
  iso,
  use_biotype = FALSE,
  use_evidence = FALSE
)
tx_tbl <- sel_tx |>
  dplyr::select(
    "transcript_id",
    "transcript_name",
    "transcript_type",
    "gene_id",
    "is_canonical",
    "tags",
    "strand",
    "n_exons",
    "tx_length"
  ) |>
  dplyr::left_join(
    dplyr::transmute(
      rules,
      .data$transcript_id,
      rule_productivity = as.character(.data$productivity),
      rule_reason = .data$reason,
      rule_source = .data$source
    ),
    by = "transcript_id"
  ) |>
  dplyr::mutate(
    karousis = dplyr::case_when(
      .data$transcript_id %in% pairs$nmd_tx ~ "nmd",
      .data$transcript_id %in% pairs$non_nmd_tx ~ "non_nmd",
      .default = NA_character_
    )
  )

# 6. Exons, introns, splice sites and regions ---------------------------------
ex_tbl <- exons |>
  dplyr::filter(.data$transcript_id %in% tx_tbl$transcript_id) |>
  dplyr::select("transcript_id", "seqid", "strand", "start_pos", "end_pos") |>
  dplyr::arrange(.data$transcript_id, .data$start_pos)
ex_tbl$gene_id <- tx_tbl$gene_id[match(
  ex_tbl$transcript_id,
  tx_tbl$transcript_id
)]
intron_tbl <- ex_tbl |>
  dplyr::mutate(
    intron_start = .data$end_pos + 1L,
    intron_end = dplyr::lead(.data$start_pos) - 1L,
    .by = "transcript_id"
  ) |>
  dplyr::filter(!is.na(.data$intron_end)) |>
  dplyr::select(
    "gene_id",
    "transcript_id",
    "seqid",
    "strand",
    "intron_start",
    "intron_end"
  )
stopifnot(all(intron_tbl$intron_end >= intron_tbl$intron_start))

# Donor: first intron base in the direction of transcription
plus <- intron_tbl$strand == "+"
site_tbl <- dplyr::bind_rows(
  tibble::tibble(
    gene_id = intron_tbl$gene_id,
    seqid = intron_tbl$seqid,
    strand = intron_tbl$strand,
    site = "donor",
    pos = ifelse(plus, intron_tbl$intron_start, intron_tbl$intron_end)
  ),
  tibble::tibble(
    gene_id = intron_tbl$gene_id,
    seqid = intron_tbl$seqid,
    strand = intron_tbl$strand,
    site = "acceptor",
    pos = ifelse(plus, intron_tbl$intron_end, intron_tbl$intron_start)
  )
) |>
  dplyr::distinct()
region_tbl <- dplyr::bind_rows(
  dplyr::transmute(
    ex_tbl,
    .data$gene_id,
    .data$seqid,
    .data$strand,
    region = "exon",
    .data$start_pos,
    .data$end_pos
  ),
  dplyr::transmute(
    intron_tbl,
    .data$gene_id,
    .data$seqid,
    .data$strand,
    region = "intron",
    start_pos = .data$intron_start,
    end_pos = .data$intron_end
  )
) |>
  dplyr::distinct()

write_parquet(gene_tbl, fs::path(out_dir, "genes.parquet"))
write_parquet(tx_tbl, fs::path(out_dir, "transcripts.parquet"))
write_parquet(ex_tbl, fs::path(out_dir, "exons.parquet"))
write_parquet(intron_tbl, fs::path(out_dir, "introns.parquet"))
write_parquet(site_tbl, fs::path(out_dir, "sites.parquet"))
write_parquet(region_tbl, fs::path(out_dir, "regions.parquet"))
write_parquet(pairs, fs::path(out_dir, "karousis_pairs.parquet"))

cli::cli_h2("Targets")
cli::cli_inform(c(
  i = "Annotation: {.path {annot_dir}}",
  i = "Genes: {nrow(gene_tbl)} ({sum(gene_tbl$stratum == 'karousis')} with validated NMD pairs)",
  i = "Transcripts: {nrow(tx_tbl)}; unique introns: {nrow(dplyr::distinct(intron_tbl, .data$gene_id, .data$intron_start, .data$intron_end))}",
  i = "Karousis labels: {sum(tx_tbl$karousis %in% 'nmd')} NMD, {sum(tx_tbl$karousis %in% 'non_nmd')} non-NMD"
))
print(table(tx_tbl$transcript_type, tx_tbl$rule_productivity, useNA = "ifany"))
