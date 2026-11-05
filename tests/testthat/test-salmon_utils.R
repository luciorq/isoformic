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
