evidence_experiment <- function(envir = parent.frame()) {
  annot <- mini_annotation("mini_cds.gff3", envir = envir)
  extract_transcript_sequences(
    annot,
    genome_fasta = mini_genome_fasta(envir = envir)
  )
  mat <- matrix(
    c(10, 5, 2, 1, 3, 10, 4, 3, 2, 3, 2, 1, 8, 6, 3, 1, 2, 9, 7, 3),
    nrow = 5,
    dimnames = list(c("tA1", "tA2", "tA3", "tA4", "tB1"), paste0("s", 1:4))
  )
  IsoformicExperiment(
    experiment_name = "evidence",
    annot_path = annot,
    assays = list(abundance = mat),
    col_data = data.frame(
      sample_id = paste0("s", 1:4),
      condition = c("a", "a", "b", "b")
    )
  )
}

testthat::test_that("add_appris joins principal isoform labels and scores", {
  iso <- evidence_experiment()
  appris <- tibble::tibble(
    transcript_id = c("tA1", "tA3", "zzz"),
    gene_id = "gA",
    appris_annotation = c("PRINCIPAL:1", "MINOR", "PRINCIPAL:2"),
    appris_principal = c(TRUE, FALSE, TRUE),
    appris_score = c(10, 2, 5),
    firestar_residues = c(3, 0, 1),
    spade_domain_score = c(100, 20, 50)
  )
  testthat::expect_message(add_appris(iso, appris), "Added")
  rd <- row_data(iso)
  testthat::expect_identical(
    rd$appris_annotation[rd$transcript_id == "tA1"],
    "PRINCIPAL:1"
  )
  testthat::expect_true(is.na(rd$appris_annotation[rd$transcript_id == "tA2"]))
  p <- transcript_productivity(iso)
  testthat::expect_identical(
    p$reason[p$transcript_id == "tA1"],
    "APPRIS principal isoform"
  )
  testthat::expect_error(
    add_appris(iso, appris[3, ]),
    class = "isoformic_invalid_evidence"
  )
})

riboseq_table <- function() {
  tibble::tibble(
    orf_id = c("main", "uorf", "listed", "skip"),
    phase = c(2L, 1L, 1L, 2L),
    seqid = "chr1",
    start_pos = c(1100, 1010, 3000, 5000),
    end_pos = c(8300, 1021, 3200, 5200),
    strand = "+",
    n_blocks = c(4L, 1L, 1L, 1L),
    block_starts = c("1100;3000;5000;8000", "1010", "3000", "5000"),
    block_ends = c("1200;3200;5200;8300", "1021", "3200", "5200"),
    orf_type = c(NA, "uORF", "intORF", NA),
    gene_id = "gA",
    gene_name = "GENEA",
    transcript_biotype = "protein_coding",
    transcript_ids = c(NA, NA, "tA2.1", NA),
    replicated = c(NA, TRUE, TRUE, NA),
    studies = NA_character_,
    protein = NA_character_
  )
}

testthat::test_that("Ribo-seq ORFs are matched by exon-block compatibility and listed transcripts", {
  iso <- evidence_experiment()
  exons <- iso@annot_data_exons
  pairs <- compatible_orf_transcripts(riboseq_table(), exons)
  main <- pairs$transcript_id[pairs$orf_id == "main"]
  # tA2 skips exon 3, tA3 has an extra exon between the blocks, tA4 retains the intron
  testthat::expect_identical(main, "tA1")
  uorf <- sort(pairs$transcript_id[pairs$orf_id == "uorf"])
  testthat::expect_identical(uorf, c("tA1", "tA2", "tA3", "tA4"))
  skip <- sort(pairs$transcript_id[pairs$orf_id == "skip"])
  testthat::expect_identical(skip, c("tA1", "tA3", "tA4"))
  add_riboseq_orfs(iso, riboseq_table())
  rd <- row_data(iso)
  rd <- rd[match(c("tA1", "tA2", "tA3", "tA4", "tB1"), rd$transcript_id), ]
  testthat::expect_identical(
    rd$riboseq_main_orf,
    c(TRUE, FALSE, FALSE, FALSE, NA)
  )
  testthat::expect_identical(rd$riboseq_uorf, c(TRUE, TRUE, TRUE, TRUE, NA))
  testthat::expect_equal(rd$riboseq_n_orfs, c(4L, 2L, 3L, 3L, NA))
  testthat::expect_match(rd$riboseq_orf_types[2], "intORF")
  p <- transcript_productivity(iso)
  testthat::expect_identical(p$source[p$transcript_id == "tA1"], "evidence")
  testthat::expect_match(p$reason[p$transcript_id == "tA1"], "Ribo-seq")
  # Without sequences the annotated CDS defines the main ORF
  annot2 <- mini_annotation("mini_cds.gff3")
  add_riboseq_orfs(annot2, riboseq_table())
  tx2 <- read_parquet(fs::path(annot2, "annot_data_transcripts.parquet"))
  testthat::expect_true(tx2$riboseq_main_orf[tx2$transcript_id == "tA1"])
  testthat::expect_error(
    add_riboseq_orfs(annot2, riboseq_table()[0, ]),
    class = "isoformic_invalid_evidence"
  )
})

testthat::test_that("add_nmd_exons flags transcripts with validated NMD exons", {
  iso <- evidence_experiment()
  nmd <- tibble::tibble(
    seqid = c("1", "1"),
    start_pos = c(4000L, 123L),
    end_pos = c(4100L, 456L),
    strand = "+",
    nmd_transcript_ids = c("ENST1", "tA4.3"),
    non_nmd_transcript_ids = c("ENST2", "x"),
    gene_id = "gA",
    gene_name = "GENEA"
  )
  add_nmd_exons(iso, nmd)
  rd <- row_data(iso)
  testthat::expect_true(rd$nmd_validated_exon[rd$transcript_id == "tA3"])
  testthat::expect_true(is.na(rd$nmd_validated_exon[rd$transcript_id == "tA1"]))
  # tA4.3 is not an Ensembl identifier, so its version is not stripped: no match
  testthat::expect_true(is.na(rd$nmd_validated_exon[rd$transcript_id == "tA4"]))
  p <- transcript_productivity(iso)
  testthat::expect_identical(p$source[p$transcript_id == "tA3"], "evidence")
  testthat::expect_error(
    add_nmd_exons(iso, nmd[2, ]),
    class = "isoformic_invalid_evidence"
  )
})

testthat::test_that("domain_retention compares isoform domains with the canonical isoform", {
  iso <- evidence_experiment()
  domains <- tibble::tibble(
    transcript_id = c("tA1", "tA1", "tA3"),
    gene_id = "gA",
    seqid = "chr1",
    start_pos = c(3000L, 5000L, 3000L),
    end_pos = c(3200L, 5200L, 3200L),
    strand = "+",
    score = 10,
    pfam_acc = c("PF1", "PF2", "PF1"),
    pfam_name = c("DomA", "DomB", "DomA"),
    evalue = 1e-5,
    pep_start = c(40L, 110L, 40L),
    pep_end = c(100L, 170L, 100L)
  )
  dr <- domain_retention(iso, domains = domains)
  dr <- dr[match(c("tA1", "tA2", "tA3", "tA4", "tB1"), dr$transcript_id), ]
  testthat::expect_identical(
    dr$in_domain_table,
    c(TRUE, FALSE, TRUE, FALSE, FALSE)
  )
  testthat::expect_equal(dr$n_domains_canonical, c(2L, 2L, 2L, 2L, 0L))
  testthat::expect_equal(dr$n_retained, c(2L, 1L, 1L, 1L, 0L))
  testthat::expect_equal(dr$n_lost, c(0L, 1L, 1L, 0L, 0L))
  testthat::expect_equal(dr$n_partial, c(0L, 0L, 0L, 1L, 0L))
  testthat::expect_identical(dr$domains_lost, c(NA, "DomB", "DomB", NA, NA))
  testthat::expect_equal(dr$domain_retention_fraction, c(1, 0.5, 0.5, 0.5, NA))
  testthat::expect_identical(
    domain_retention(iso, gene_name = "GENEB", domains = domains)$transcript_id,
    "tB1"
  )
})

testthat::test_that("gtex_reference derives the healthy-tissue state per gene", {
  iso <- evidence_experiment()
  gtex <- tibble::tibble(
    transcript_id = c("tA1", "tA2", "tA3", "tB1", "tA1", "tA3"),
    gene_id = c("gA", "gA", "gA", "gB", "gA", "gA"),
    tissue = c(rep("Tissue X", 4), "Tissue Y", "Tissue Y"),
    median_tpm = c(8, 2, 10, 1, 1, 9),
    n_samples = 5L
  )
  ref <- gtex_reference(iso, tissue = "Tissue X", gtex = gtex)
  a <- ref[ref$gene_id == "gA", ]
  testthat::expect_equal(a$reference_gene_tpm, 20)
  testthat::expect_equal(a$reference_productive_fraction, 0.5)
  testthat::expect_identical(a$reference_dominant_transcript_id, "tA3")
  testthat::expect_identical(a$reference_dominant_productivity, "unproductive")
  testthat::expect_equal(a$n_reference_transcripts, 3L)
  b <- ref[ref$gene_id == "gB", ]
  testthat::expect_true(is.na(b$reference_productive_fraction))
  strict <- gtex_reference(
    iso,
    tissue = "Tissue X",
    gtex = gtex,
    include_altered = FALSE
  )
  testthat::expect_equal(
    strict$reference_productive_fraction[strict$gene_id == "gA"],
    0.4
  )
  testthat::expect_error(
    gtex_reference(iso, tissue = "Nope", gtex = gtex),
    class = "isoformic_unknown_tissue"
  )
  tissues <- gtex_tissues(gtex)
  testthat::expect_identical(tissues$tissue, c("Tissue X", "Tissue Y"))
  pf <- productive_fraction(iso, group = "condition", reference = ref)
  testthat::expect_equal(pf$reference_productive_fraction, 0.5)
})
