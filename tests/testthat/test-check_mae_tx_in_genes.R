# tests/testthat/test-check_mae_tx_in_genes.R
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
      colData = S4Vectors::DataFrame(condition = c("A", "B"))
    )
  ),
  colData = S4Vectors::DataFrame(
    condition = c("A", "B"),
    row.names = c("sample1", "sample2")
  ),
  sampleMap = S4Vectors::DataFrame(
    assay = as.factor(rep(c("transcript", "gene"), each = 2)),
    primary = c("sample1", "sample2", "sample1", "sample2"),
    colname = c("txseq1", "txseq2", "geneseq1", "geneseq2")
  ),
  metadata = list(isoformic = list(tx_to_gene = tx_to_gene_df))
)

test_that("check_mae_tx_in_genes checks if all transcripts have a gene equivalent", {
  # Correct case
  expect_error(check_mae_tx_in_genes(mae), NA)

  # Incorrect case
  # se_tx_wrong <- se_tx
  mae_wrong <- mae
  rownames(mae_wrong@ExperimentList[["transcript"]]) <- c("5", "7", "2", "tx5")
  expect_error(check_mae_tx_in_genes(mae_wrong), class = "isoformic_mae_missing_genes")
})
