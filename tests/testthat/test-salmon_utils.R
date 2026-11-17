testthat::test_that("salmon_index passes the user-supplied env_name through to condathis::run() (regression test)", {
  testthat::skip_if_not_installed("condathis")

  captured_args <- new.env(parent = emptyenv())
  captured_args$value <- NULL
  mock_run <- function(...) {
    captured_args$value <- list(...)
    invisible(NULL)
  }

  testthat::local_mocked_bindings(
    env_exists = function(...) TRUE,
    run = mock_run,
    .package = "condathis"
  )

  salmon_index(
    fasta_path = "transcripts.fa",
    index_path = "salmon_index",
    env_name = "my-custom-env"
  )

  testthat::expect_equal(captured_args$value$env_name, "my-custom-env")
})

testthat::test_that("salmon_quant passes the user-supplied env_name through to condathis::run() (regression test)", {
  testthat::skip_if_not_installed("condathis")

  captured_args <- new.env(parent = emptyenv())
  captured_args$value <- NULL
  mock_run <- function(...) {
    captured_args$value <- list(...)
    invisible(NULL)
  }

  index_dir <- withr::local_tempdir()

  testthat::local_mocked_bindings(
    env_exists = function(...) TRUE,
    run = mock_run,
    .package = "condathis"
  )

  salmon_quant(
    input_r1 = "reads_1.fq",
    index_path = index_dir,
    env_name = "another-custom-env"
  )

  testthat::expect_equal(captured_args$value$env_name, "another-custom-env")
})

testthat::test_that("salmon_index defaults to salmon-env when env_name is not supplied", {
  testthat::skip_if_not_installed("condathis")

  captured_args <- new.env(parent = emptyenv())
  captured_args$value <- NULL
  mock_run <- function(...) {
    captured_args$value <- list(...)
    invisible(NULL)
  }

  testthat::local_mocked_bindings(
    env_exists = function(...) TRUE,
    run = mock_run,
    .package = "condathis"
  )

  salmon_index(fasta_path = "transcripts.fa", index_path = "salmon_index")

  testthat::expect_equal(captured_args$value$env_name, "salmon-env")
})

testthat::test_that("salmon_quant builds the read arguments for single-end and paired-end input", {
  testthat::skip_if_not_installed("condathis")

  captured <- new.env()
  mock_run <- function(...) {
    assign("args", unlist(list(...)), envir = captured)
    invisible(NULL)
  }

  index_dir <- withr::local_tempdir()

  testthat::local_mocked_bindings(
    env_exists = function(...) TRUE,
    run = mock_run,
    .package = "condathis"
  )

  salmon_quant(input_r1 = "reads.fq", index_path = index_dir, num_gibbs = 20)
  args <- captured$args
  testthat::expect_identical(
    unname(args[which(args == "--unmatedReads") + 1L]),
    "reads.fq"
  )
  testthat::expect_false("--mates2" %in% args)
  testthat::expect_identical(
    unname(args[which(args == "--numGibbsSamples") + 1L]),
    "20"
  )
  # Flags that Salmon 2 ignores are no longer passed
  testthat::expect_false(any(
    c(
      "--no-version-check",
      "--disableChainingHeuristic",
      "--useVBOpt",
      "--validateMappings"
    ) %in%
      args
  ))

  salmon_quant(input_r1 = "r1.fq", input_r2 = "r2.fq", index_path = index_dir)
  args <- captured$args
  testthat::expect_identical(
    unname(args[which(args == "--mates1") + 1L]),
    "r1.fq"
  )
  testthat::expect_identical(
    unname(args[which(args == "--mates2") + 1L]),
    "r2.fq"
  )
  testthat::expect_false("--unmatedReads" %in% args)

  salmon_index(
    fasta_path = "transcripts.fa",
    index_path = "salmon_index",
    is_gencode = TRUE
  )
  args <- captured$args
  testthat::expect_true("--gencode" %in% args)
  testthat::expect_false(any(
    c("--type", "puff", "--no-version-check") %in% args
  ))
})

testthat::test_that("salmon_index builds a gentrome and decoys.txt from decoy_fasta", {
  testthat::skip_if_not_installed("condathis")

  tmp <- withr::local_tempdir()
  tx_fasta <- file.path(tmp, "transcripts.fa")
  writeLines(c(">tx1|gene1|", "ACGTACGT", ">tx2|gene2|", "GGGGCCCC"), tx_fasta)
  decoy_plain <- file.path(tmp, "genome.fa")
  writeLines(
    c(">chr1 1 dna:chromosome", "NNNNACGT", ">chrM", "TTTT"),
    decoy_plain
  )
  decoy_gz <- file.path(tmp, "genome.fa.gz")
  gz_con <- gzfile(decoy_gz, "wb")
  writeLines(readLines(decoy_plain), gz_con)
  close(gz_con)

  captured <- new.env()
  mock_run <- function(...) {
    args <- unlist(list(...))
    assign("args", args, envir = captured)
    if (!"--decoys" %in% args) {
      return(invisible(NULL))
    }
    gentrome <- args[[which(args == "--transcripts") + 1L]]
    decoys <- args[[which(args == "--decoys") + 1L]]
    assign("gentrome_name", basename(gentrome), envir = captured)
    assign("gentrome_is_gz", is_gzip_file(gentrome), envir = captured)
    assign("gentrome", readLines(gentrome), envir = captured)
    assign("decoys", readLines(decoys), envir = captured)
    invisible(NULL)
  }
  testthat::local_mocked_bindings(
    env_exists = function(...) TRUE,
    run = mock_run,
    .package = "condathis"
  )

  expected_gentrome <- c(
    ">tx1|gene1|",
    "ACGTACGT",
    ">tx2|gene2|",
    "GGGGCCCC",
    ">chr1 1 dna:chromosome",
    "NNNNACGT",
    ">chrM",
    "TTTT"
  )

  # Plain transcripts, gzipped genome: the gentrome is gzipped
  salmon_index(fasta_path = tx_fasta, decoy_fasta = decoy_gz)
  testthat::expect_true("--decoys" %in% captured$args)
  testthat::expect_false(tx_fasta %in% captured$args)
  testthat::expect_identical(captured$decoys, c("chr1", "chrM"))
  testthat::expect_identical(captured$gentrome, expected_gentrome)
  testthat::expect_identical(captured$gentrome_name, "gentrome.fa.gz")
  testthat::expect_true(captured$gentrome_is_gz)

  # Plain genome: the gentrome is plain
  salmon_index(fasta_path = tx_fasta, decoy_fasta = decoy_plain)
  testthat::expect_identical(captured$gentrome, expected_gentrome)
  testthat::expect_identical(captured$gentrome_name, "gentrome.fa")
  testthat::expect_false(captured$gentrome_is_gz)

  # The temporary files are removed once the index is built
  testthat::expect_false(file.exists(
    captured$args[[which(captured$args == "--decoys") + 1L]]
  ))

  # Without decoys nothing changes
  salmon_index(fasta_path = tx_fasta)
  testthat::expect_false("--decoys" %in% captured$args)
  testthat::expect_identical(
    unname(captured$args[[which(captured$args == "--transcripts") + 1L]]),
    tx_fasta
  )
})

testthat::test_that("salmon_index validates decoy_fasta before running Salmon", {
  testthat::skip_if_not_installed("condathis")

  tmp <- withr::local_tempdir()
  tx_fasta <- file.path(tmp, "transcripts.fa")
  writeLines(c(">tx1", "ACGT"), tx_fasta)
  no_header <- file.path(tmp, "no_header.fa")
  writeLines("ACGT", no_header)

  called <- new.env()
  called$run <- FALSE
  testthat::local_mocked_bindings(
    env_exists = function(...) TRUE,
    run = function(...) {
      called$run <- TRUE
      invisible(NULL)
    },
    .package = "condathis"
  )

  testthat::expect_error(
    salmon_index(fasta_path = tx_fasta, decoy_fasta = file.path(tmp, "x.fa")),
    class = "isoformic_file_not_found"
  )
  testthat::expect_error(
    salmon_index(fasta_path = file.path(tmp, "x.fa"), decoy_fasta = tx_fasta),
    class = "isoformic_file_not_found"
  )
  testthat::expect_error(
    salmon_index(fasta_path = tx_fasta, decoy_fasta = no_header),
    class = "isoformic_gentrome_error"
  )
  testthat::expect_false(called$run)
})
