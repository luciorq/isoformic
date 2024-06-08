# tests/testthat/test-check_se_metadata_level.R
library(testthat)

# Create a dummy SummarizedExperiment object
se <- SummarizedExperiment::SummarizedExperiment(
  assays = list(counts = matrix(1:4, ncol = 2)),
  colData = S4Vectors::DataFrame(condition = c("A", "B")),
  metadata = list(level = "transcript")
)

test_that("check_se_metadata_level checks for 'level' in metadata", {
  se_no_level <- se
  se_no_level@metadata <- list()

  expect_error(check_se_metadata_level(se_no_level), class = "isoformic_se_metadata_without_level")
  expect_error(check_se_metadata_level(se), NA)
})

test_that("check_se_metadata_level_type checks if 'level' in metadata is a character", {
  se_wrong_type <- se
  se_wrong_type@metadata$level <- 123

  expect_error(check_se_metadata_level_type(se_wrong_type), class = "isoformic_se_metadata_level_not_char")
  expect_error(check_se_metadata_level_type(se), NA)
})


test_that("check_se_metadata_level_length checks if 'level' in metadata is of length 1", {
  se_wrong_length <- se
  se_wrong_length@metadata$level <- c("transcript", "gene")

  expect_error(check_se_metadata_level_length(se_wrong_length), class = "isoformic_se_level_length")
  expect_error(check_se_metadata_level_length(se), NA)
})
