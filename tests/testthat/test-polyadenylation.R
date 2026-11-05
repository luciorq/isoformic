apa_experiment <- function(envir = parent.frame(), with_sequences = TRUE) {
  annot <- mini_annotation("mini_apa.gff3", envir = envir)
  if (isTRUE(with_sequences)) {
    extract_transcript_sequences(
      annot,
      genome_fasta = mini_genome_fasta(envir = envir)
    )
  }
  ids <- c("tA1", "tA2", "tA3", "tA4", "tA5", "tA6", "tB1", "tB2")
  mat <- matrix(
    c(
      10,
      5,
      2,
      1,
      4,
      2,
      3,
      1,
      10,
      4,
      3,
      2,
      2,
      2,
      3,
      2,
      2,
      1,
      8,
      6,
      1,
      3,
      3,
      4,
      1,
      2,
      9,
      7,
      2,
      1,
      3,
      5
    ),
    nrow = 8,
    dimnames = list(ids, paste0("s", 1:4))
  )
  IsoformicExperiment(
    experiment_name = "apa",
    annot_path = annot,
    assays = list(abundance = mat, counts = mat * 10),
    col_data = data.frame(
      sample_id = paste0("s", 1:4),
      condition = c("a", "a", "b", "b")
    )
  )
}

testthat::test_that("polya_sites ranks 3' ends and classifies APA events", {
  iso <- apa_experiment()
  pa <- polya_sites(iso)
  pa <- pa[match(c("tA1", "tA5", "tA6", "tB1", "tB2"), pa$transcript_id), ]
  testthat::expect_equal(pa$tes_pos, c(9000, 8400, 6500, 20000, 20100))
  testthat::expect_equal(pa$n_tes_groups, c(3L, 3L, 3L, 2L, 2L))
  # Plus strand: proximal = smallest position; minus strand: proximal = largest
  testthat::expect_equal(pa$tes_rank, c(3L, 2L, 1L, 2L, 1L))
  testthat::expect_identical(pa$is_distal, c(TRUE, FALSE, FALSE, TRUE, FALSE))
  testthat::expect_identical(pa$is_proximal, c(FALSE, FALSE, TRUE, FALSE, TRUE))
  testthat::expect_identical(
    pa$apa_type,
    c("same_end", "tandem_utr", "intronic", "same_end", "tandem_utr")
  )
  testthat::expect_equal(pa$utr3_length[1:2], c(700, 100))
  testthat::expect_identical(pa$seqid, rep("chr1", 5L))
  testthat::expect_true(all(
    c("polya_signal", "polya_signal_distance") %in% colnames(pa)
  ))
  # The synthetic 3' UTR of tA1 carries no planted signal: NA or a weak hexamer
  testthat::expect_true(
    is.na(pa$polya_signal[1]) ||
      pa$polya_signal[1] %in% polya_signal_hexamers()
  )
  testthat::expect_identical(
    polya_sites(iso, gene_name = "GENEB")$transcript_id,
    c("tB1", "tB2")
  )
})

testthat::test_that("polya_sites works without sequences and joins an atlas", {
  iso <- apa_experiment(with_sequences = FALSE)
  pa <- polya_sites(iso)
  testthat::expect_true(all(is.na(pa$polya_signal)))
  testthat::expect_equal(pa$utr3_length[pa$transcript_id == "tA5"], 100)
  atlas <- tibble::tibble(
    site_id = c("s1", "s2", "s3"),
    seqid = c("1", "1", "1"),
    cleavage_pos = c(9010L, 8390L, 6000L),
    strand = c("+", "+", "+"),
    site_class = c("TE", "TE", "IN"),
    tpm = c(5, 1, 0.2)
  )
  pa <- polya_sites(iso, polya_atlas = atlas, window = 20)
  pa <- pa[match(c("tA1", "tA5", "tA6"), pa$transcript_id), ]
  testthat::expect_identical(pa$pas_site_id, c("s1", "s2", NA))
  testthat::expect_identical(pa$pas_supported, c(TRUE, TRUE, FALSE))
  testthat::expect_equal(pa$pas_distance, c(10L, 10L, NA))
  testthat::expect_identical(pa$pas_class, c("TE", "TE", NA))
  # Minimal atlas columns are completed
  pa2 <- polya_sites(
    iso,
    polya_atlas = atlas[, c("seqid", "cleavage_pos", "strand")]
  )
  testthat::expect_true(pa2$pas_supported[pa2$transcript_id == "tA1"])
  testthat::expect_error(
    polya_sites(iso, polya_atlas = data.frame(x = 1)),
    class = "isoformic_invalid_evidence"
  )
  # A GENCODE poly(A) GFF3 is parsed into the same shape
  gff <- withr::local_tempfile(fileext = ".gff3")
  writeLines(
    c(
      "##gff-version 3",
      "chr1\tHAVANA\tpolyA_site\t9005\t9006\t.\t+\t.\tID=1;gene_type=polya_site",
      "chr1\tHAVANA\tpolyA_signal\t8970\t8975\t.\t+\t.\tID=2;gene_type=polya_signal",
      "chr1\tHAVANA\tpolyA_site\t19995\t19996\t.\t-\t.\tID=3;gene_type=polya_site"
    ),
    gff
  )
  sites <- read_polya_gff(gff)
  testthat::expect_equal(nrow(sites), 3L)
  testthat::expect_equal(sites$cleavage_pos, c(9006L, 8975L, 19995L))
  pa3 <- polya_sites(iso, polya_atlas = gff)
  testthat::expect_identical(pa3$pas_site_id[pa3$transcript_id == "tA1"], "1")
  testthat::expect_identical(pa3$pas_site_id[pa3$transcript_id == "tB1"], "3")
})

testthat::test_that("find_polya_signal reports the strongest hexamer closest to the end", {
  found <- find_polya_signal(c(
    "CCCAATAAACCCCCATTAAACCC",
    "GGGGGG",
    "ATTAAACC",
    NA
  ))
  testthat::expect_identical(found$signal, c("AATAAA", NA, "ATTAAA", NA))
  testthat::expect_equal(found$distance, c(19L, NA, 7L, NA))
})

testthat::test_that("apa_usage quantifies distal usage and 3' UTR length per gene", {
  iso <- apa_experiment()
  usage <- apa_usage(iso, group = "condition")
  testthat::expect_setequal(usage$gene_id, c("gA", "gB"))
  a <- usage[usage$gene_id == "gA", ]
  # Distal group of gA: tA1..tA4 end at 9000; tA5 and tA6 are proximal
  s1_total <- 10 + 5 + 2 + 1 + 4 + 2
  s1_distal <- 10 + 5 + 2 + 1
  mat <- apa_usage(iso, as_matrix = TRUE)
  testthat::expect_equal(unname(mat["gA", "s1"]), s1_distal / s1_total)
  testthat::expect_equal(a$n_tes_groups, 3L)
  testthat::expect_match(a$apa_types, "intronic")
  testthat::expect_match(a$apa_types, "tandem_utr")
  testthat::expect_equal(
    a$delta_distal_fraction,
    a$distal_fraction_b - a$distal_fraction_a
  )
  testthat::expect_true(all(
    c("weighted_utr3_length_a", "delta_weighted_utr3_length") %in% colnames(a)
  ))
  # tA1 (700) and tA5 (100) carry the annotated 3' UTRs of the coding isoforms
  testthat::expect_true(
    a$weighted_utr3_length > 100 && a$weighted_utr3_length < 3000
  )
  b <- usage[usage$gene_id == "gB", ]
  testthat::expect_equal(b$distal_fraction, mean(c(3 / 4, 3 / 5, 3 / 7, 3 / 8)))
  testthat::expect_true(is.na(b$weighted_utr3_length))
  per_sample <- apa_usage(iso)
  testthat::expect_false("delta_distal_fraction" %in% colnames(per_sample))
  floored <- apa_usage(iso, min_expr = 3, as_matrix = TRUE)
  testthat::expect_equal(unname(floored["gA", "s1"]), 15 / 19)
  testthat::expect_identical(apa_usage(iso, gene_name = "GENEB")$gene_id, "gB")
})

testthat::test_that("test_apa tests the shift in distal usage between groups", {
  iso <- apa_experiment()
  res <- test_apa(iso, group = "condition")
  testthat::expect_setequal(res$gene_id, c("gA", "gB"))
  a <- res[res$gene_id == "gA", ]
  testthat::expect_equal(a$distal_fraction_a, (180 + 190) / (240 + 230))
  testthat::expect_true(all(
    c("log_odds_ratio", "pvalue", "padj", "dispersion") %in% colnames(res)
  ))
  testthat::expect_true(all(res$pvalue >= 0 & res$pvalue <= 1, na.rm = TRUE))
  testthat::expect_identical(
    test_apa(iso, group = "condition", gene_name = "GENEA")$gene_id,
    "gA"
  )
  testthat::expect_error(
    test_apa(iso, group = "sample_id"),
    class = "isoformic_group_levels"
  )
})

testthat::test_that("test_apa contrasts two of three groups", {
  iso <- apa_experiment()
  sample_table <- col_data(iso)
  sample_table$condition <- c("a", "a", "b", "other")
  col_data(iso) <- sample_table
  testthat::expect_error(
    test_apa(iso, group = "condition"),
    class = "isoformic_group_levels"
  )
  res <- test_apa(iso, group = "condition", levels = c("a", "b"))
  a <- res[res$gene_id == "gA", ]
  testthat::expect_equal(a$distal_fraction_a, (180 + 190) / (240 + 230))
})
