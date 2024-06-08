# tests/testthat/test-check_mae.R
library(testthat)

# Create a dummy MultiAssayExperiment object
mae <- MultiAssayExperiment::MultiAssayExperiment(
  experiments = list(
    transcript = SummarizedExperiment::SummarizedExperiment(
      assays = list(counts = matrix(1:8, ncol = 2, dimnames = list(
        c("tx1", "tx2", "tx3", "tx4"), c("txseq1", "txseq2")
      ))),
      colData = S4Vectors::DataFrame(condition = c("A", "B"))
      ),
    gene = SummarizedExperiment::SummarizedExperiment(
      assays = list(counts = matrix(5:8, ncol = 2, dimnames = list(
        c("gene1", "gene2"), c("geneseq1", "geneseq2")
      ))),
      colData = S4Vectors::DataFrame(condition = c("A", "B")))
  ),
  colData = S4Vectors::DataFrame(
    condition = c("A", "B"),
    row.names = c("sample1", "sample2"))
  ,
  sampleMap = S4Vectors::DataFrame(
    assay = as.factor(rep(c("transcript", "gene"), each = 2)),
    primary = c("sample1", "sample2", "sample1", "sample2"),
    colname = c("txseq1", "txseq2", "geneseq1", "geneseq2")
  )
)

test_that("is_mae correctly identifies MultiAssayExperiment objects", {
  expect_true(is_mae(mae))
  expect_false(is_mae(list()))
  expect_false(is_mae(data.frame()))
})

test_that("check_mae correctly handles MultiAssayExperiment objects", {
  expect_error(check_mae(list()), class = "isoformic_obj_not_mae")
  expect_error(check_mae(data.frame()), class = "isoformic_obj_not_mae")
  expect_error(check_mae(mae), NA)
})

test_that("check_mae_assay_rownames checks for matching rownames in assays", {
  # Correct case
  expect_error(check_mae_assay_rownames(mae), NA)

  # Incorrect case
  mae_wrong <- mae
  mae_wrong@sampleMap$primary <- c("sample1", "sample2", "sample3", "sample4")

  expect_error(check_mae_assay_rownames(mae_wrong), class = "isoformic_mae_experiment_rownames")
})
