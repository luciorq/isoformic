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
  out <- suppressMessages(
    extract_transcript_sequences(annot, genome_fasta = fasta)
  )
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

testthat::test_that("transcript_sequences selects transcripts while reading", {
  annot <- mini_annotation("mini_cds.gff3")
  suppressMessages(extract_transcript_sequences(
    annot,
    genome_fasta = mini_genome_fasta()
  ))
  all_seqs <- transcript_sequences(annot)
  ids <- c(rev(all_seqs$transcript_id), all_seqs$transcript_id[1L], "nope")
  testthat::expect_identical(
    transcript_sequences(annot, transcript_id = ids),
    all_seqs[all_seqs$transcript_id %in% ids, ]
  )
  testthat::expect_identical(
    transcript_sequences(annot, transcript_id = "tB1"),
    all_seqs[all_seqs$transcript_id %in% "tB1", ]
  )
  testthat::expect_identical(
    transcript_sequences(annot, transcript_id = character()),
    all_seqs[0L, ]
  )
})

testthat::test_that("genome records match seqids with or without chr prefix, gzip too", {
  annot <- mini_annotation("mini_cds.gff3")
  fasta <- mini_genome_fasta(seq_name = "1", gzip = TRUE)
  suppressMessages(extract_transcript_sequences(annot, genome_fasta = fasta))
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
    suppressMessages(
      extract_transcript_sequences(annot, transcript_fasta = path)
    ),
    "No sequence found for 3"
  )
  seqs <- transcript_sequences(annot)
  testthat::expect_identical(seqs$transcript_id, c("tA1", "tB1"))
  testthat::expect_identical(seqs$sequence[[1L]], "ACGTACGTACGT")
  testthat::expect_identical(seqs$seq_length, c(12L, 4L))
  # Kept sequences from the same source: a message only
  testthat::expect_no_warning(
    suppressMessages(
      extract_transcript_sequences(annot, transcript_fasta = path)
    )
  )
  # From another source: a warning with the transcripts left without one
  n_tx <- nrow(read_parquet(fs::path(annot, "annot_data_transcripts.parquet")))
  testthat::expect_warning(
    suppressMessages(
      extract_transcript_sequences(annot, genome_fasta = mini_genome_fasta())
    ),
    class = "isoformic_sequence_source_mismatch"
  )
  testthat::expect_warning(
    suppressMessages(
      extract_transcript_sequences(annot, genome_fasta = mini_genome_fasta())
    ),
    paste0(n_tx - 2L, " of ", n_tx, " transcripts")
  )
  testthat::expect_identical(transcript_sequences(annot), seqs)
})

testthat::test_that("sequence utilities translate and reverse complement", {
  testthat::expect_identical(reverse_complement("ATGCN"), "NGCAT")
  testthat::expect_identical(translate_dna("ATGGCATAA"), "MA*")
  testthat::expect_identical(translate_dna("ATGNNN"), "MX")
  testthat::expect_identical(translate_dna("AT"), "")
})

testthat::test_that("ORFs are projected from the canonical start codon", {
  annot <- mini_annotation("mini_cds.gff3")
  suppressMessages(extract_transcript_sequences(
    annot,
    genome_fasta = mini_genome_fasta()
  ))
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

testthat::test_that("ORFs of selected genes of a directory match the full table", {
  annot <- mini_annotation("mini_cds.gff3")
  suppressMessages(extract_transcript_sequences(
    annot,
    genome_fasta = mini_genome_fasta()
  ))
  withr::local_options(isoformic.orf_cache = FALSE)
  full <- transcript_orfs(annot)
  gene_b <- full$gene_id[full$gene_name == "GENEB"][1L]
  for (genes in list("GENEA", c(gene_b, "GENEA"), c("GENEB", "NOPE"), "NOPE")) {
    testthat::expect_identical(
      transcript_orfs(annot, gene_name = genes),
      filter_genes(full, genes)
    )
  }
  # Computed for one gene first, the stored table gives the same rows
  withr::local_options(isoformic.orf_cache = TRUE)
  testthat::expect_identical(
    transcript_orfs(annot, gene_name = "GENEB"),
    filter_genes(full, "GENEB")
  )
  testthat::expect_identical(transcript_orfs(annot), full)
})

testthat::test_that("start-proximal stops and long 3' UTRs are flagged", {
  annot <- mini_annotation("mini_cds.gff3")
  suppressMessages(extract_transcript_sequences(
    annot,
    genome_fasta = mini_genome_fasta()
  ))
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

testthat::test_that("the ORF and coding potential helpers match their direct definitions", {
  withr::local_seed(11)
  random_seq <- function(n, alphabet = c("A", "C", "G", "T")) {
    paste(sample(alphabet, n, replace = TRUE), collapse = "")
  }
  seqs <- c(
    vapply(c(3, 7, 60, 301, 1000, 2500), random_seq, character(1L)),
    random_seq(400, c("A", "C", "G", "T", "N")),
    "",
    "AC",
    "acgtacgtacgtATGAAATAA",
    "ATGAAATAAATGCCCTGA"
  )
  stops <- find_stop_positions(seqs)

  # Fickett TESTCODE (Fickett 1982, the tables of CPAT), base by base
  fickett_reference <- function(s) {
    position_prob <- list(
      A = c(0.94, 0.68, 0.84, 0.93, 0.58, 0.68, 0.45, 0.34, 0.20, 0.22),
      C = c(0.80, 0.70, 0.70, 0.81, 0.66, 0.48, 0.51, 0.33, 0.30, 0.23),
      G = c(0.90, 0.88, 0.74, 0.64, 0.53, 0.48, 0.27, 0.16, 0.08, 0.08),
      T = c(0.97, 0.97, 0.91, 0.68, 0.69, 0.44, 0.54, 0.20, 0.09, 0.09)
    )
    position_weight <- c(A = 0.26, C = 0.18, G = 0.31, T = 0.33)
    position_para <- c(1.9, 1.8, 1.7, 1.6, 1.5, 1.4, 1.3, 1.2, 1.1, 0.0)
    content_prob <- list(
      A = c(0.28, 0.49, 0.44, 0.55, 0.62, 0.49, 0.54, 0.52, 0.36, 0.51),
      C = c(0.82, 0.64, 0.51, 0.64, 0.59, 0.59, 0.56, 0.49, 0.26, 0.36),
      G = c(0.40, 0.54, 0.47, 0.64, 0.64, 0.73, 0.41, 0.41, 0.33, 0.29),
      T = c(0.28, 0.24, 0.39, 0.40, 0.55, 0.75, 0.56, 0.69, 0.51, 0.58)
    )
    content_weight <- c(A = 0.11, C = 0.12, G = 0.15, T = 0.14)
    content_para <- c(0.33, 0.31, 0.29, 0.27, 0.25, 0.23, 0.21, 0.17, 0.17, 0)
    lookup <- function(value, para, prob) {
      i <- which(value >= para)[1L]
      if (is.na(i)) prob[[length(prob)]] else prob[[i]]
    }
    n <- nchar(s)
    if (n < 3L) {
      return(0)
    }
    ch <- strsplit(s, "", fixed = TRUE)[[1L]]
    frame <- (seq_along(ch) - 1L) %% 3L
    score <- 0
    for (b in c("A", "C", "G", "T")) {
      counts <- vapply(0:2, function(f) sum(ch[frame == f] == b), 0)
      asym <- max(counts) / (min(counts) + 1)
      comp <- sum(ch == b) / n
      score <- score +
        lookup(asym, position_para, position_prob[[b]]) *
          position_weight[[b]] +
        lookup(comp, content_para, content_prob[[b]]) * content_weight[[b]]
    }
    score
  }
  testthat::expect_identical(
    fickett_score(seqs),
    vapply(seqs, fickett_reference, numeric(1L), USE.NAMES = FALSE)
  )

  # Hexamer score: mean log-odds of the in-frame hexamers, unknown ones left
  # out
  table <- hexamer_table()
  hexamer_reference <- function(s) {
    n <- nchar(s)
    if (n < 6L) {
      return(0)
    }
    starts <- seq.int(1L, n - 5L, by = 3L)
    v <- table[substring(s, starts, starts + 5L)]
    v <- v[!is.na(v)]
    if (length(v)) mean(v) else 0
  }
  testthat::expect_identical(
    hexamer_score(seqs),
    vapply(seqs, hexamer_reference, numeric(1L), USE.NAMES = FALSE)
  )

  # A multibyte character counts as one position, as in nchar() and
  # substring(), so the frame of the bases after it does not shift
  odd <- c(
    "ACGT\u00c4ACGTACGTACGTATGGCC",
    "ATGAAAC\u00e9GGTTTCCCAAAGGGTAA",
    "\u4e2dATGCCCGGGAAATTTTAG"
  )
  testthat::expect_identical(
    hexamer_score(odd),
    vapply(odd, hexamer_reference, numeric(1L), USE.NAMES = FALSE)
  )
  testthat::expect_identical(
    fickett_score(odd),
    vapply(odd, fickett_reference, numeric(1L), USE.NAMES = FALSE)
  )

  # Longest ORF, frame by frame; ties go to the lower frame, then to the
  # upstream ORF
  longest_reference <- function(s, stops, min_orf_length) {
    a <- as.integer(gregexpr("ATG", s, fixed = TRUE)[[1L]])
    a <- a[a > 0L]
    best <- c(NA_real_, NA_real_, 0)
    for (f in 0:2) {
      for (x in a[a %% 3L == f]) {
        nxt <- stops[stops > x & (stops - x) %% 3L == 0L]
        end <- if (length(nxt)) {
          nxt[[1L]] + 2
        } else {
          x + floor((nchar(s) - x + 1) / 3) * 3 - 1
        }
        if (end - x + 1 > best[[3L]]) {
          best <- c(x, end, end - x + 1)
        }
      }
    }
    if (best[[3L]] < min_orf_length) {
      best <- rep(NA_real_, 3L)
    }
    best
  }
  for (min_len in c(6L, 90L)) {
    res <- longest_orfs(seqs, stops, min_orf_length = min_len)
    ref <- mapply(
      longest_reference,
      seqs,
      stops,
      MoreArgs = list(min_orf_length = min_len),
      USE.NAMES = FALSE
    )
    testthat::expect_identical(res$start, ref[1L, ])
    testthat::expect_identical(res$end, ref[2L, ])
    testthat::expect_identical(res$length, ref[3L, ])
  }
  tie <- longest_orfs(
    "ATGAAATAAATGCCCTGA",
    find_stop_positions("ATGAAATAAATGCCCTGA"),
    min_orf_length = 6L
  )
  testthat::expect_identical(c(tie$start, tie$end), c(1, 9))

  # uORFs: the first in-frame stop after each upstream ATG
  uorf_reference <- function(s, stops, p) {
    utr <- substring(s, 1L, p - 1)
    a <- as.integer(gregexpr("ATG", utr, fixed = TRUE)[[1L]])
    a <- a[a > 0L]
    stop_end <- vapply(
      a,
      function(x) {
        sx <- stops[stops > x & (stops - x) %% 3 == 0]
        if (length(sx)) sx[[1L]] + 2 else NA_real_
      },
      numeric(1L)
    )
    complete <- !is.na(stop_end) & stop_end < p & (stop_end - a + 1) >= 9L
    overlapping <- (p - a) %% 3 != 0 & (is.na(stop_end) | stop_end >= p)
    c(sum(complete), sum(overlapping))
  }
  long <- nchar(seqs) > 100L
  orf_start <- c(50, 120, 250, 333, 101)[seq_len(sum(long))]
  u <- count_uorfs(seqs[long], stops[long], orf_start)
  ref <- mapply(uorf_reference, seqs[long], stops[long], orf_start)
  testthat::expect_identical(u$n_uorf, unname(ref[1L, ]))
  testthat::expect_identical(u$n_overlapping_uorf, unname(ref[2L, ]))
})

testthat::test_that("transcript_orfs works on IsoformicExperiment objects and is cached", {
  iso <- mini_experiment(with_dea = FALSE)
  testthat::expect_error(transcript_orfs(iso), class = "isoformic_no_sequences")
  annot <- mini_annotation("mini_cds.gff3")
  suppressMessages(extract_transcript_sequences(
    annot,
    genome_fasta = mini_genome_fasta()
  ))
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
  orfs <- suppressWarnings(transcript_orfs(iso))
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

testthat::test_that("stream_fasta drops the carriage returns of CRLF files", {
  # writeLines() writes CRLF on Windows; the binary reader once kept the `\r`
  path <- withr::local_tempfile(fileext = ".fa")
  con <- file(path, open = "wb")
  writeLines(
    c(">chr1 synthetic", "ACGTACGTAC", "GTACGT", ">chr2", "TTTT"),
    con,
    sep = "\r\n"
  )
  close(con)
  got <- new.env()
  stream_fasta(
    path,
    keep = function(name) TRUE,
    on_record = function(name, sequence) assign(name, sequence, envir = got),
    chunk_bytes = 7
  )
  testthat::expect_setequal(ls(got), c("chr1", "chr2"))
  testthat::expect_identical(got$chr1, "ACGTACGTACGTACGT")
  testthat::expect_identical(got$chr2, "TTTT")
})

testthat::test_that("transcript_orfs keeps its rows in the annotation directory", {
  annot <- mini_annotation("mini_cds.gff3")
  suppressMessages(extract_transcript_sequences(
    annot,
    genome_fasta = mini_genome_fasta()
  ))
  file <- fs::path(annot, "annot_data_orfs.parquet")
  ids <- paste0("t", c("A1", "A2", "A3", "A4", "B1"))
  experiment <- function(ids) {
    IsoformicExperiment(
      experiment_name = "cds",
      annot_path = annot,
      assays = list(
        abundance = matrix(
          1,
          nrow = length(ids),
          ncol = 2,
          dimnames = list(ids, c("s1", "s2"))
        )
      ),
      col_data = data.frame(sample_id = c("s1", "s2"), condition = c("a", "b"))
    )
  }
  # The reference: nothing read, nothing written
  reference <- function(ids, ...) {
    withr::local_options(isoformic.orf_cache = FALSE)
    transcript_orfs(experiment(ids), ...)
  }
  expected <- reference(ids)
  expected_rules <- reference(ids, start_proximal = 400, long_utr3 = 500)
  testthat::expect_false(identical(
    expected_rules$nmd_class,
    expected$nmd_class
  ))
  testthat::expect_false(fs::file_exists(file))

  first <- transcript_orfs(experiment(ids))
  testthat::expect_identical(first, expected)
  testthat::expect_true(fs::file_exists(file))
  stored <- read_parquet(file)
  testthat::expect_setequal(stored$transcript_id, ids)

  # A new object reads the rows: with the computation taken away it still
  # answers, and the file is not written again
  written <- fs::file_info(file)$modification_time
  testthat::local_mocked_bindings(
    compute_transcript_orfs = function(...) {
      stop("computed again", call. = FALSE)
    }
  )
  testthat::expect_identical(transcript_orfs(experiment(ids)), expected)
  testthat::expect_identical(transcript_orfs(experiment(rev(ids))), {
    expected[match(rev(ids), expected$transcript_id), ]
  })
  testthat::expect_identical(transcript_orfs(annot), expected)
  testthat::expect_identical(fs::file_info(file)$modification_time, written)
  # The NMD class follows the arguments of the call, not those of the table
  testthat::expect_identical(
    transcript_orfs(experiment(ids), start_proximal = 400, long_utr3 = 500),
    expected_rules
  )
  testthat::expect_error(
    transcript_orfs(experiment(ids), min_orf_length = 30L),
    "computed again"
  )
})

testthat::test_that("the stored ORF rows are not used for another canonical isoform or another annotation", {
  annot <- mini_annotation("mini_cds.gff3")
  suppressMessages(extract_transcript_sequences(
    annot,
    genome_fasta = mini_genome_fasta()
  ))
  file <- fs::path(annot, "annot_data_orfs.parquet")
  transcripts <- read_parquet(fs::path(annot, "annot_data_transcripts.parquet"))
  transcripts <- transcripts[transcripts$gene_id %in% transcripts$gene_id[1L], ]
  canonical <- canonical_transcripts(transcripts)[[1L]]
  ids <- transcripts$transcript_id
  testthat::expect_gt(length(ids), 2L)
  experiment <- function(ids) {
    IsoformicExperiment(
      experiment_name = "cds",
      annot_path = annot,
      assays = list(
        abundance = matrix(
          1,
          nrow = length(ids),
          ncol = 2,
          dimnames = list(ids, c("s1", "s2"))
        )
      ),
      col_data = data.frame(sample_id = c("s1", "s2"), condition = c("a", "b"))
    )
  }
  reference <- function(ids, ...) {
    withr::local_options(isoformic.orf_cache = FALSE)
    transcript_orfs(experiment(ids), ...)
  }
  whole <- transcript_orfs(experiment(ids))
  expected_whole <- reference(ids)
  expected_others <- reference(setdiff(ids, canonical))
  # Without the canonical isoform of the annotation the gene has another one
  others <- setdiff(ids, canonical)
  partial <- transcript_orfs(experiment(others))
  testthat::expect_identical(partial, expected_others)
  testthat::expect_false(canonical %in% partial$canonical_transcript_id)
  # Its rows are not stored: the others still read those of the annotation
  testthat::expect_identical(
    unique(read_parquet(file)$canonical_transcript_id),
    canonical
  )
  written <- fs::file_info(file)$modification_time
  testthat::expect_identical(transcript_orfs(experiment(ids)), whole)
  testthat::expect_identical(fs::file_info(file)$modification_time, written)
  testthat::expect_identical(whole, expected_whole)

  # A table of other inputs, or one that cannot be read, is replaced
  testthat::expect_identical(
    transcript_orfs(experiment(ids), min_orf_length = 30L),
    reference(ids, min_orf_length = 30L)
  )
  testthat::expect_identical(transcript_orfs(experiment(ids)), whole)
  writeLines("not parquet", file)
  testthat::expect_identical(transcript_orfs(experiment(ids)), whole)
  testthat::expect_setequal(read_parquet(file)$transcript_id, ids)
  stale <- read_parquet(file)
  stale$orf_length <- -1
  stale$table_signature <- "another annotation"
  write_parquet(stale, file)
  testthat::expect_identical(transcript_orfs(experiment(ids)), whole)

  # Never into the installed package, and not with the option off
  testthat::expect_false(
    orf_table_writable(system.file("extdata", package = "isoformic"))
  )
  fs::file_delete(file)
  withr::with_options(list(isoformic.orf_cache = FALSE), {
    testthat::expect_identical(transcript_orfs(experiment(ids)), whole)
  })
  testthat::expect_false(fs::file_exists(file))
})

testthat::test_that("stored and computed ORF rows are mixed in the order of the call", {
  annot <- mini_annotation("mini_cds.gff3")
  suppressMessages(extract_transcript_sequences(
    annot,
    genome_fasta = mini_genome_fasta()
  ))
  file <- fs::path(annot, "annot_data_orfs.parquet")
  uncached <- function(...) {
    withr::local_options(isoformic.orf_cache = FALSE)
    transcript_orfs(...)
  }
  all_rows <- uncached(annot)
  gene_a <- all_rows$transcript_id[all_rows$gene_name == "GENEA"]
  gene_b <- setdiff(all_rows$transcript_id, gene_a)
  canonical_a <- unique(all_rows$canonical_transcript_id[
    all_rows$transcript_id %in% gene_a
  ])
  transcripts <- read_parquet(fs::path(annot, "annot_data_transcripts.parquet"))
  subset_orfs <- function(ids) {
    tx <- transcripts[match(ids, transcripts$transcript_id), ]
    tx$gene_name <- all_rows$gene_name[match(ids, all_rows$transcript_id)]
    stored_transcript_orfs(
      annot_dir = annot,
      transcripts = tx,
      exons = function() {
        read_parquet(fs::path(annot, "annot_data_exons.parquet"))
      },
      annotation = function() transcripts
    )
  }
  # Gene B first: its rows are stored, those of gene A are not
  subset_orfs(gene_b)
  testthat::expect_setequal(read_parquet(file)$transcript_id, gene_b)
  # Gene A is computed, gene B read, and the rows come back as asked
  seen <- new.env(parent = emptyenv())
  seen$computed <- character()
  compute <- compute_transcript_orfs
  testthat::local_mocked_bindings(
    compute_transcript_orfs = function(transcripts, ...) {
      seen$computed <- c(seen$computed, transcripts$transcript_id)
      compute(transcripts, ...)
    }
  )
  ids <- rev(all_rows$transcript_id)
  mixed <- subset_orfs(ids)
  testthat::expect_setequal(seen$computed, gene_a)
  testthat::expect_identical(
    as.data.frame(mixed),
    as.data.frame(all_rows[match(ids, all_rows$transcript_id), ])
  )
  # Writing gene A kept the rows of gene B
  testthat::expect_setequal(read_parquet(file)$transcript_id, ids)
  # A gene without its canonical isoform is computed although its rows are
  # stored, and the stored rows stay
  seen$computed <- character()
  others <- c(setdiff(gene_a, canonical_a), gene_b)
  subset_orfs(others)
  testthat::expect_setequal(seen$computed, setdiff(gene_a, canonical_a))
  testthat::expect_identical(
    unique(read_parquet(file)$canonical_transcript_id[
      read_parquet(file)$transcript_id %in% gene_a
    ]),
    canonical_a
  )
})

testthat::test_that("the stored ORF rows follow the annotation and are written only where allowed", {
  annot <- mini_annotation("mini_cds.gff3")
  suppressMessages(extract_transcript_sequences(
    annot,
    genome_fasta = mini_genome_fasta()
  ))
  file <- fs::path(annot, "annot_data_orfs.parquet")
  uncached <- function(...) {
    withr::local_options(isoformic.orf_cache = FALSE)
    transcript_orfs(...)
  }
  first <- transcript_orfs(annot)
  testthat::expect_identical(first, uncached(annot))

  # The names follow the genes table, which is not part of the signature
  genes_file <- fs::path(annot, "annot_data_genes.parquet")
  genes <- read_parquet(genes_file)
  genes$gene_name[genes$gene_name == "GENEA"] <- "RENAMED"
  write_parquet(genes, genes_file)
  renamed <- transcript_orfs(annot)
  testthat::expect_identical(renamed, uncached(annot))
  testthat::expect_identical(
    nrow(transcript_orfs(annot, gene_name = "RENAMED")),
    sum(first$gene_name == "GENEA")
  )

  # Other sequences, another signature
  before <- orf_table_signature(annot, 90L)
  sequences_file <- fs::path(annot, "annot_data_sequences.parquet")
  sequences <- read_parquet(sequences_file)
  sequences$sequence[[1L]] <- chartr("ACGT", "TGCA", sequences$sequence[[1L]])
  Sys.sleep(1.1)
  write_parquet(sequences, sequences_file)
  testthat::expect_false(identical(orf_table_signature(annot, 90L), before))
  testthat::expect_false(identical(orf_table_signature(annot, 30L), before))
  changed <- transcript_orfs(annot)
  testthat::expect_identical(changed, uncached(annot))
  testthat::expect_false(identical(changed, renamed))

  # Where writing is not allowed the rows are computed and nothing is written
  fs::file_delete(file)
  testthat::local_mocked_bindings(orf_table_writable = function(...) FALSE)
  testthat::expect_identical(transcript_orfs(annot), changed)
  testthat::expect_false(fs::file_exists(file))
})

testthat::test_that("the ORF table of the example goes to its private copy", {
  iso <- example_isoformic()
  orfs <- suppressWarnings(transcript_orfs(iso))
  testthat::expect_gt(nrow(orfs), 0L)
  testthat::expect_true(
    fs::file_exists(fs::path(iso@annot_path, "annot_data_orfs.parquet"))
  )
  testthat::expect_false(fs::file_exists(system.file(
    "extdata",
    "example_experiment",
    "annotation",
    "annot_data_orfs.parquet",
    package = "isoformic"
  )))
})
