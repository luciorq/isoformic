# tests/testthat/test-check_mae_isoformic.R
library(testthat)

# Create a dummy `tx_to_gene` slot
tx_to_gene_df <- S4Vectors::DataFrame(
  gene_id = c("gene1", "gene2", "gene1", "gene2"),
  tx_id = c("tx1", "tx2", "tx3", "tx4"),
  row.names = c("tx1", "tx2", "tx3", "tx4")
)

# Create a dummy `MultiAssayExperiment` object
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
  ),
  metadata = list(isoformic = list(tx_to_gene = tx_to_gene_df))
)

test_that("check_mae_isoformic_is_list checks if isoformic is a list in metadata", {
  expect_error(check_mae_isoformic_is_list(mae), NA)

  mae_wrong <- mae
  mae_wrong@metadata <- list(isoformic = "not_a_list")
  expect_error(check_mae_isoformic_is_list(mae_wrong), class = "isoformic_mae_metadata_is_list")
})

test_that("check_mae_isoformic_tx_to_gene checks if tx_to_gene is in isoformic metadata", {
  expect_error(check_mae_isoformic_tx_to_gene(mae), NA)

  mae_wrong <- mae
  mae_wrong@metadata <- list(isoformic = list())
  expect_error(check_mae_isoformic_tx_to_gene(mae_wrong), class = "isoformic_mae_metadata_no_tx_to_gene")
})

test_that("validate_isoformic_mae run all checks", {
  expect_error(check_mae_isoformic_is_list(mae), NA)

  mae_wrong <- mae
  mae_wrong@metadata <- list(isoformic = "not_a_list")
  expect_error(check_mae_isoformic_is_list(mae_wrong), class = "isoformic_mae_metadata_is_list")
})

test_that("check_mae_isoformic_tx_to_gene checks if tx_to_gene is in isoformic metadata", {
  expect_error(check_mae_isoformic_tx_to_gene(mae), NA)

  mae_wrong <- mae
  mae_wrong@metadata <- list(isoformic = list())
  expect_error(check_mae_isoformic_tx_to_gene(mae_wrong), class = "isoformic_mae_metadata_no_tx_to_gene")
})
