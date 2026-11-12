testthat::test_that("splice-site windows are cut in the direction of transcription", {
  annot <- mini_annotation("mini_cds.gff3")
  testthat::expect_false(has_splice_sites(annot))
  testthat::expect_error(
    splice_site_strength(annot),
    class = "isoformic_no_splice_sites"
  )
  extract_transcript_sequences(annot, genome_fasta = mini_genome_fasta())
  testthat::expect_true(has_splice_sites(annot))
  sites <- read_parquet(fs::path(annot, "annot_data_splice_sites.parquet"))
  genome <- mini_genome_seq()
  # Unique introns: 1201-2999, 3201-4999, 3201-3999, 4101-4999, 3201-7999,
  # 5201-7999 on the plus strand and 20301-20799 on the minus strand
  testthat::expect_equal(nrow(sites), 7L)
  testthat::expect_true(all(nchar(sites$donor_seq) == 9L))
  testthat::expect_true(all(nchar(sites$acceptor_seq) == 23L))
  first <- sites[sites$intron_start == 1201L, ]
  testthat::expect_identical(first$donor_seq, substr(genome, 1198L, 1206L))
  testthat::expect_identical(first$acceptor_seq, substr(genome, 2980L, 3002L))
  minus <- sites[sites$strand == "-", ]
  testthat::expect_equal(minus$intron_start, 20301L)
  testthat::expect_identical(
    minus$donor_seq,
    reverse_complement(substr(genome, 20794L, 20802L))
  )
  testthat::expect_identical(
    minus$acceptor_seq,
    reverse_complement(substr(genome, 20298L, 20320L))
  )
  meta <- kv_to_list(read_parquet(fs::path(annot, "annot_metadata.parquet")))
  testthat::expect_match(meta$splice_sites_source, "^genome:")
  testthat::expect_identical(meta$splice_sites_n, "7")
})

testthat::test_that("extract_splice_sites adds the windows without the sequences", {
  annot <- mini_annotation("mini_cds.gff3")
  out <- extract_splice_sites(annot, genome_fasta = mini_genome_fasta())
  testthat::expect_true(fs::file_exists(out))
  testthat::expect_false(has_sequences(annot))
  testthat::expect_message(
    extract_splice_sites(annot, genome_fasta = mini_genome_fasta()),
    "already exists"
  )
  # A window that runs off the sequence is not scored
  introns <- tibble::tibble(
    seqid = "chr1",
    strand = c("+", "+"),
    intron_start = c(2L, 30L),
    intron_end = c(25L, 60L)
  )
  chr <- paste(rep("ACGT", 20L), collapse = "")
  windows <- splice_site_windows(introns, chr)
  testthat::expect_true(is.na(windows$donor_seq[1L]))
  testthat::expect_identical(windows$donor_seq[2L], substr(chr, 27L, 35L))
  testthat::expect_identical(windows$acceptor_seq[2L], substr(chr, 41L, 63L))
  # A sequence region without introns
  none <- splice_site_windows(introns[0L, ], chr)
  testthat::expect_equal(nrow(none), 0L)
  testthat::expect_named(none, c(names(introns), "donor_seq", "acceptor_seq"))
})

testthat::test_that("position weight matrices score the consensus highest", {
  seqs <- c(rep("AAGGTAAGT", 8L), "CAGGTGAGT", "TTGGCACGA")
  background <- rep(0.25, 4L)
  pwm <- position_weight_matrix(seqs, background)
  testthat::expect_identical(dim(pwm), c(4L, 9L))
  testthat::expect_equal(pwm[["G", 4L]], log2((10 + 0.5) / 12 / 0.25))
  testthat::expect_equal(pwm[["C", 4L]], log2(0.5 / 12 / 0.25))
  scores <- score_pwm(c(seqs[c(1L, 9L, 10L)], "AAGGTAAGN", "AAG", NA), pwm)
  testthat::expect_gt(scores[1L], scores[2L])
  testthat::expect_gt(scores[2L], scores[3L])
  testthat::expect_true(all(is.na(scores[4:6])))
  testthat::expect_equal(
    scores[1L],
    sum(pwm[cbind(match(strsplit(seqs[1L], "")[[1L]], rownames(pwm)), 1:9)])
  )
})

testthat::test_that("splice_site_strength reports the weakest site of a transcript", {
  annot <- mini_annotation("mini_cds.gff3")
  extract_splice_sites(annot, genome_fasta = mini_genome_fasta())
  by_intron <- splice_site_strength(annot, by = "intron")
  # 3 + 2 + 4 + 2 + 1 introns
  testthat::expect_equal(nrow(by_intron), 12L)
  testthat::expect_identical(attr(by_intron, "training_set"), "all")
  testthat::expect_equal(attr(by_intron, "n_training"), 7L)
  testthat::expect_false(anyNA(by_intron$donor_score))
  by_tx <- splice_site_strength(annot)
  testthat::expect_equal(nrow(by_tx), 5L)
  one <- function(id, col) by_tx[[col]][by_tx$transcript_id == id]
  testthat::expect_equal(one("tA3", "n_introns"), 4L)
  testthat::expect_equal(one("tB1", "n_introns"), 1L)
  a3 <- by_intron[by_intron$transcript_id == "tA3", ]
  testthat::expect_equal(one("tA3", "weakest_donor"), min(a3$donor_score))
  testthat::expect_equal(one("tA3", "mean_acceptor"), mean(a3$acceptor_score))
  testthat::expect_equal(
    one("tA3", "n_noncanonical"),
    sum(!a3$canonical_dinucleotides)
  )
  # The same intron has the same score in every transcript that uses it
  shared <- by_intron[by_intron$intron_start == 1201L, ]
  testthat::expect_equal(nrow(shared), 4L)
  testthat::expect_length(unique(shared$donor_score), 1L)
  testthat::expect_identical(
    splice_site_strength(annot, gene_name = "GENEB")$transcript_id,
    "tB1"
  )
  testthat::expect_error(
    splice_site_strength(annot, by = "exon"),
    class = "rlang_error"
  )
})

testthat::test_that("transcripts without introns have no score", {
  exons <- tibble::tibble(
    transcript_id = c("t1", "t2", "t2"),
    seqid = "chr1",
    strand = "+",
    start_pos = c(1L, 1L, 11L),
    end_pos = c(50L, 10L, 20L)
  )
  # Adjacent exons leave no intron
  testthat::expect_equal(nrow(transcript_introns(exons)), 0L)
  testthat::expect_equal(nrow(unique_introns(exons[1L, ])), 0L)
})
