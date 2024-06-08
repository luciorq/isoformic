# tests/testthat/test-check_se.R
library(testthat)

# Create a dummy SummarizedExperiment object
se <- SummarizedExperiment::SummarizedExperiment(
  assays = list(counts = matrix(1:4, ncol = 2)),
  colData = S4Vectors::DataFrame(condition = c("A", "B"))
)

test_that("is_se correctly identifies SummarizedExperiment objects", {
  expect_true(is_se(se))
  expect_false(is_se(list()))
  expect_false(is_se(data.frame()))
})

test_that("check_se correctly handles SummarizedExperiment objects", {
  expect_error(check_se(list()), class = "isoformic_obj_not_se")
  expect_error(check_se(data.frame()), class = "isoformic_obj_not_se")
  expect_error(check_se(se), NA)
})
