# tests/testthat/test-se_experiment_level.R
library(testthat)

# Create a dummy SummarizedExperiment object
se <- SummarizedExperiment::SummarizedExperiment(
  assays = list(counts = matrix(1:4, ncol = 2)),
  colData = S4Vectors::DataFrame(condition = c("A", "B")),
  metadata = list(level = "txp")
)

test_that("se_experiment_level returns the correct level", {
  expect_equal(se_experiment_level(se), "txp")

  se@metadata$level <- "gene"
  expect_equal(se_experiment_level(se), "gene")
})

test_that("se_experiment_level handles invalid SummarizedExperiment objects", {
  expect_error(se_experiment_level(list()), class = "isoformic_obj_not_se")

  se_invalid <- se
  se_invalid@metadata <- list()
  expect_error(se_experiment_level(se_invalid), class = "isoformic_se_metadata_without_level")

  se_invalid@metadata$level <- c("txp", "gene")
  expect_error(se_experiment_level(se_invalid), class = "isoformic_se_level_length")

  se_invalid@metadata$level <- 123
  expect_error(se_experiment_level(se_invalid), class = "isoformic_se_metadata_level_not_char")
})

test_that("se_experiment_level handles invalid level values", {
  se_invalid <- se
  se_invalid@metadata$level <- "invalid_level"
  expect_error(se_experiment_level(se_invalid), class = "isoformic_se_level_not_txp_or_gene")
})
