# Full-data tests: real downloads and external tools. They run only with
# `ISOFORMIC_FULL=true` (see `skip_if_not_full()` in helper-fixtures.R).

testthat::test_that("GENCODE downloads and prepares into an annotation directory", {
  skip_if_not_full()
  testthat::skip_if_offline()
  out <- withr::local_tempdir()
  gff <- download_reference(
    version = "33",
    file_type = "gff",
    output_path = out
  )
  testthat::expect_true(fs::file_exists(gff))
  annot_dir <- suppressMessages(prepare_isoformic_annotation(
    gff,
    output_path = out
  ))
  testthat::expect_true(is_annotation_dir(annot_dir))
  tx <- read_parquet(fs::path(annot_dir, "annot_data_transcripts.parquet"))
  testthat::expect_gt(nrow(tx), 200000L)
  testthat::expect_true(all(
    c("protein_coding", "retained_intron", "nonsense_mediated_decay") %in%
      tx$transcript_type
  ))
})

testthat::test_that("Salmon indexes and quantifies the example transcripts", {
  skip_if_not_full()
  testthat::skip_if_offline()
  testthat::skip_if_not_installed("condathis")
  testthat::skip_on_os("windows")
  extdata <- fs::path_package("isoformic", "extdata")
  fasta <- fs::path(extdata, "example_transcripts.fa.gz")
  root <- withr::local_tempdir()

  # Simulate 200 single-end 76 bp reads from each transcript
  fa_lines <- readLines(gzfile(fasta))
  headers <- grepl("^>", fa_lines)
  seqs <- vapply(
    split(fa_lines[!headers], cumsum(headers)[!headers]),
    paste,
    character(1L),
    collapse = ""
  )
  seqs <- seqs[nchar(seqs) >= 100L]
  withr::local_seed(1L)
  reads <- unlist(lapply(seqs, function(s) {
    starts <- sample.int(nchar(s) - 75L, 200L, replace = TRUE)
    substring(s, starts, starts + 75L)
  }))
  fastq <- fs::path(root, "sample.fastq")
  writeLines(
    paste0(
      "@r",
      seq_along(reads),
      "\n",
      reads,
      "\n+\n",
      strrep("I", 76L)
    ),
    fastq
  )

  index <- fs::path(root, "salmon_index")
  suppressMessages(salmon_index(
    fasta_path = fasta,
    index_path = index,
    is_gencode = TRUE,
    num_threads = 2
  ))
  testthat::expect_true(fs::dir_exists(index))
  quant <- fs::path(root, "sample_quant")
  suppressMessages(salmon_quant(
    input_r1 = fastq,
    index_path = index,
    output_dir = quant,
    num_threads = 2,
    num_gibbs = 5
  ))
  testthat::expect_true(fs::file_exists(fs::path(quant, "quant.sf")))

  annot_dir <- suppressMessages(prepare_isoformic_annotation(
    fs::path(extdata, "example.gff3"),
    output_path = root
  ))
  iso <- suppressMessages(read_salmon(
    c(sample = quant),
    annot_path = annot_dir,
    inferential = TRUE
  ))
  testthat::expect_equal(col_names(iso), "sample")
  testthat::expect_gt(sum(assay(iso, "counts")), 0.9 * length(reads))
})
