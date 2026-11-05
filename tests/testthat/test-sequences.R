orf_value <- function(orfs, id, column) {
  orfs[[column]][orfs$transcript_id == id]
}

testthat::test_that("sequences are extracted from a genome FASTA", {
  annot <- mini_annotation("mini_cds.gff3")
  fasta <- mini_genome_fasta()
  testthat::expect_false(has_sequences(annot))
  testthat::expect_error(
    transcript_sequences(annot),
    class = "isoformic_no_sequences"
  )
  out <- extract_transcript_sequences(annot, genome_fasta = fasta)
  testthat::expect_true(fs::file_exists(out))
  testthat::expect_true(has_sequences(annot))
  seqs <- transcript_sequences(annot)
  tx <- read_parquet(fs::path(annot, "annot_data_transcripts.parquet"))
  testthat::expect_setequal(seqs$transcript_id, tx$transcript_id)
  testthat::expect_equal(
    seqs$seq_length[match(tx$transcript_id, seqs$transcript_id)],
    tx$tx_length
  )
  genome <- mini_genome_seq()
  # Plus strand: exon 1 of tA1 is genomic 1000-1200
  seq_a1 <- seqs$sequence[seqs$transcript_id == "tA1"]
  testthat::expect_identical(
    substr(seq_a1, 1L, 201L),
    substr(genome, 1000L, 1200L)
  )
  testthat::expect_identical(substr(seq_a1, 101L, 103L), "ATG")
  # Minus strand: tB1 starts with the reverse complement of 20800-21000
  seq_b1 <- seqs$sequence[seqs$transcript_id == "tB1"]
  testthat::expect_identical(
    substr(seq_b1, 1L, 201L),
    reverse_complement(substr(genome, 20800L, 21000L))
  )
  meta <- kv_to_list(read_parquet(fs::path(annot, "annot_metadata.parquet")))
  testthat::expect_match(meta$sequences_source, "^genome:")
  # Existing sequences are kept unless overwrite = TRUE
  testthat::expect_message(
    extract_transcript_sequences(annot, genome_fasta = fasta),
    "already exists"
  )
  testthat::expect_error(
    extract_transcript_sequences(annot, overwrite = TRUE),
    class = "isoformic_sequence_source"
  )
})

testthat::test_that("genome records match seqids with or without chr prefix, gzip too", {
  annot <- mini_annotation("mini_cds.gff3")
  fasta <- mini_genome_fasta(seq_name = "1", gzip = TRUE)
  extract_transcript_sequences(annot, genome_fasta = fasta)
  seqs <- transcript_sequences(annot, transcript_id = "tA1")
  testthat::expect_identical(substr(seqs$sequence, 101L, 103L), "ATG")
  bad <- mini_genome_fasta(seq_name = "chrZ")
  testthat::expect_error(
    extract_transcript_sequences(annot, genome_fasta = bad, overwrite = TRUE),
    class = "isoformic_sequence_source"
  )
})

testthat::test_that("sequences are read from a transcript FASTA by identifier", {
  annot <- mini_annotation("mini_cds.gff3")
  path <- withr::local_tempfile(fileext = ".fa")
  writeLines(
    c(
      ">tA1.2|gene|extra",
      "ACGTACGTAC",
      "GT",
      ">tB1",
      "GGGG",
      ">other",
      "TTTT"
    ),
    path
  )
  testthat::expect_warning(
    extract_transcript_sequences(annot, transcript_fasta = path),
    "No sequence found for 3"
  )
  seqs <- transcript_sequences(annot)
  testthat::expect_identical(seqs$transcript_id, c("tA1", "tB1"))
  testthat::expect_identical(seqs$sequence[[1L]], "ACGTACGTACGT")
  testthat::expect_identical(seqs$seq_length, c(12L, 4L))
})

testthat::test_that("sequence utilities translate and reverse complement", {
  testthat::expect_identical(reverse_complement("ATGCN"), "NGCAT")
  testthat::expect_identical(translate_dna("ATGGCATAA"), "MA*")
  testthat::expect_identical(translate_dna("ATGNNN"), "MX")
  testthat::expect_identical(translate_dna("AT"), "")
})

testthat::test_that("ORFs are projected from the canonical start codon", {
  annot <- mini_annotation("mini_cds.gff3")
  extract_transcript_sequences(annot, genome_fasta = mini_genome_fasta())
  orfs <- transcript_orfs(annot)
  testthat::expect_identical(
    orfs$orf_source,
    c("annotation", "annotation", "annotation", "projected", "longest")
  )
  # Canonical: complete ORF 1100..8300 = 804 nt, stop TAA, no NMD
  testthat::expect_equal(orf_value(orfs, "tA1", "orf_length"), 804)
  testthat::expect_identical(orf_value(orfs, "tA1", "stop_codon"), "TAA")
  testthat::expect_identical(orf_value(orfs, "tA1", "nmd_class"), "none")
  testthat::expect_equal(orf_value(orfs, "tA1", "orf_genomic_start"), 1100)
  testthat::expect_equal(orf_value(orfs, "tA1", "orf_genomic_end"), 8300)
  testthat::expect_equal(orf_value(orfs, "tA1", "n_uorf"), 1L)
  testthat::expect_equal(orf_value(orfs, "tA1", "utr5_length"), 100)
  testthat::expect_true(orf_value(orfs, "tA1", "annotation_cds_agrees"))
  # In-frame exon skipping keeps start and stop, frame preserved
  testthat::expect_equal(orf_value(orfs, "tA2", "orf_length"), 603)
  testthat::expect_true(orf_value(orfs, "tA2", "start_same_as_canonical"))
  testthat::expect_true(orf_value(orfs, "tA2", "stop_same_as_canonical"))
  testthat::expect_identical(
    orf_value(orfs, "tA2", "frame_vs_canonical"),
    "in_frame"
  )
  # Cassette exon with a premature stop: annotation corrected, PTC called
  testthat::expect_false(orf_value(orfs, "tA3", "annotation_cds_agrees"))
  testthat::expect_equal(orf_value(orfs, "tA3", "orf_genomic_end"), 4048)
  testthat::expect_identical(orf_value(orfs, "tA3", "nmd_class"), "ptc")
  testthat::expect_identical(
    orf_value(orfs, "tA3", "frame_vs_canonical"),
    "novel_sequence"
  )
  # Retained intron without annotated CDS: projected ORF ends in the intron
  testthat::expect_identical(orf_value(orfs, "tA4", "orf_source"), "projected")
  testthat::expect_equal(orf_value(orfs, "tA4", "orf_genomic_end"), 3303)
  testthat::expect_identical(orf_value(orfs, "tA4", "nmd_class"), "ptc")
  testthat::expect_true(orf_value(orfs, "tA4", "stop_to_last_junction") > 50)
  # lncRNA: longest ORF only, low coding probability
  testthat::expect_identical(orf_value(orfs, "tB1", "orf_source"), "longest")
  testthat::expect_lt(orf_value(orfs, "tB1", "coding_probability"), 0.364)
  testthat::expect_gt(orf_value(orfs, "tA1", "coding_probability"), 0.364)
  testthat::expect_identical(
    orfs$canonical_transcript_id,
    c("tA1", "tA1", "tA1", "tA1", "tB1")
  )
  # Filtering and the coding potential view
  testthat::expect_identical(
    transcript_orfs(annot, gene_name = "GENEB")$transcript_id,
    "tB1"
  )
  cp <- coding_potential(annot)
  testthat::expect_true(cp$is_coding[cp$transcript_id == "tA1"])
  testthat::expect_false(cp$is_coding[cp$transcript_id == "tB1"])
})

testthat::test_that("start-proximal stops and long 3' UTRs are flagged", {
  annot <- mini_annotation("mini_cds.gff3")
  extract_transcript_sequences(annot, genome_fasta = mini_genome_fasta())
  orfs <- transcript_orfs(annot, start_proximal = 400, long_utr3 = 500)
  testthat::expect_identical(
    orf_value(orfs, "tA3", "nmd_class"),
    "ptc_start_proximal"
  )
  testthat::expect_identical(orf_value(orfs, "tA1", "nmd_class"), "long_3utr")
})

testthat::test_that("ORF helpers handle edge cases", {
  stops <- find_stop_positions(c("ATGTAATAG", "", "AAAA"))
  testthat::expect_identical(stops[[1L]], c(4L, 7L))
  testthat::expect_identical(stops[[2L]], integer())
  testthat::expect_identical(
    first_inframe_stop(stops, c(1, NA, 1)),
    c(6, NA, NA)
  )
  longest <- longest_orfs(
    c("CCATGAAATAACC", "ATGAAA"),
    find_stop_positions(c("CCATGAAATAACC", "ATGAAA")),
    min_orf_length = 6L
  )
  testthat::expect_equal(longest$start, c(3, 1))
  testthat::expect_equal(longest$end, c(11, 6))
  testthat::expect_equal(longest$length, c(9, 6))
  u <- count_uorfs(
    c("ATGAAATAACCATGCCC", "ATGCCATGCCC"),
    find_stop_positions(c("ATGAAATAACCATGCCC", "ATGCCATGCCC")),
    c(12, 6)
  )
  testthat::expect_identical(u$n_uorf, c(1L, 0L))
  testthat::expect_identical(u$n_overlapping_uorf, c(0L, 1L))
})

testthat::test_that("transcript_orfs works on IsoformicExperiment objects and is cached", {
  iso <- mini_experiment(with_dea = FALSE)
  testthat::expect_error(transcript_orfs(iso), class = "isoformic_no_sequences")
  annot <- mini_annotation("mini_cds.gff3")
  extract_transcript_sequences(annot, genome_fasta = mini_genome_fasta())
  mat <- matrix(
    1,
    nrow = 5,
    ncol = 2,
    dimnames = list(paste0("t", c("A1", "A2", "A3", "A4", "B1")), c("s1", "s2"))
  )
  iso <- IsoformicExperiment(
    experiment_name = "cds",
    annot_path = annot,
    assays = list(abundance = mat),
    col_data = data.frame(sample_id = c("s1", "s2"), condition = c("a", "b"))
  )
  orfs <- transcript_orfs(iso)
  testthat::expect_equal(nrow(orfs), 5L)
  testthat::expect_identical(transcript_orfs(iso), orfs)
  testthat::expect_identical(
    transcript_orfs(iso, gene_name = "GENEA")$transcript_id,
    paste0("tA", 1:4)
  )
})

testthat::test_that("build_coding_model writes a usable model", {
  set.seed(7)
  random_seq <- function(n, alphabet = c("A", "C", "G", "T")) {
    paste(sample(alphabet, n, replace = TRUE), collapse = "")
  }
  cds <- vapply(
    1:20,
    function(i) paste0("ATG", random_seq(150, c("A", "C", "G")), "TAA"),
    character(1L)
  )
  mrna <- paste0(
    vapply(1:20, function(i) random_seq(30), character(1L)),
    cds,
    "AAAAAA"
  )
  non <- vapply(1:20, function(i) random_seq(200), character(1L))
  dir <- withr::local_tempdir()
  testthat::expect_message(
    res <- suppressWarnings(build_coding_model(
      cds,
      mrna,
      non,
      output_dir = dir
    )),
    "Coding model written"
  )
  testthat::expect_true(fs::file_exists(fs::path(
    dir,
    "hexamer_logodds.parquet"
  )))
  testthat::expect_true(fs::file_exists(fs::path(
    dir,
    "coding_model_coefficients.parquet"
  )))
  testthat::expect_equal(nrow(res$hexamers), 4096L)
  testthat::expect_gt(res$accuracy, 0.9)
  # Restore the shipped model for the other tests
  rm("hexamer_table", "coding_model", envir = annotation_cache)
})

testthat::test_that("stream_fasta keeps records longer than 1e6 characters intact across chunks", {
  # substring()'s default `last = 1000000L` once truncated sequences after a header
  set.seed(3)
  long_seq <- paste(
    sample(c("A", "C", "G", "T"), 2500000, replace = TRUE),
    collapse = ""
  )
  path <- withr::local_tempfile(fileext = ".fa")
  lines <- substring(
    long_seq,
    seq(1, nchar(long_seq), by = 60),
    pmin(seq(60, nchar(long_seq) + 59, by = 60), nchar(long_seq))
  )
  writeLines(c(">big", lines, ">small", "ACGT", ">big2", lines), path)
  got <- new.env()
  stream_fasta(
    path,
    keep = function(name) TRUE,
    on_record = function(name, sequence) assign(name, sequence, envir = got),
    chunk_bytes = 3e6
  )
  testthat::expect_identical(got$big, long_seq)
  testthat::expect_identical(got$big2, long_seq)
  testthat::expect_identical(got$small, "ACGT")
})
