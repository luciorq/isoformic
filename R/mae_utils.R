# ==============================================================================
# Checks for `MultiAssayExperiment` objects
# ==============================================================================

# Check if object is of class `MultiAssayExperiment`.
check_mae <- function(mae) {
  if (!is_mae(mae)) {
    cli::cli_abort(
      message = c(
        x = "The object is not of type {.code MultiAssayExperiment}."
      ),
      class = "isoformic_obj_not_mae"
    )
  }
}

# Check if object is of class `MultiAssayExperiment`.
is_mae <- function(mae) {
  if (isTRUE(inherits(x = mae, what = c("MultiAssayExperiment")))) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

# Check if the `rownames` of each Experiment in the `mae` obj are the same
check_mae_assay_rownames <- function(mae) {
  tx_sample_name <- MultiAssayExperiment::sampleMap(mae)[MultiAssayExperiment::sampleMap(mae)[["assay"]] %in% "transcript", ][["primary"]]
  gene_sample_name <- MultiAssayExperiment::sampleMap(mae)[MultiAssayExperiment::sampleMap(mae)[["assay"]] %in% "gene", ][["primary"]]
  if (!isTRUE(identical(tx_sample_name, gene_sample_name))) {
    cli::cli_abort(
      message = c(
        `x` = "{.code se_tx} and {.code se_gene} do {.strong not} contain the same {.code rownames(colData(se))}.",
        "{.var se_tx}: {.var {tx_sample_name}}",
        "{.var se_tx}: {.var {gene_sample_name}}"
      ),
      class = "isoformic_mae_experiment_rownames"
    )
  }
}

# Check if all transcripts have a gene equivalent
check_mae_tx_in_genes <- function(mae) {
  se_tx_rownames <- SummarizedExperiment::rownames(mae@ExperimentList[["transcript"]])
  se_gene_rownames <- SummarizedExperiment::rownames(mae@ExperimentList[["gene"]])
  tx_to_gene_df <- mae@metadata[["isoformic"]]$tx_to_gene
  tx_in_genes_vector <- tx_to_gene_df[tx_to_gene_df$gene_id %in% se_gene_rownames,]$tx_id
  if (!isTRUE(all(se_tx_rownames %in% tx_in_genes_vector))) {
    cli::cli_abort(
      message = c(
        `x` = "There are genes missing from {.code rownames(experiments(mae)$gene)}."
      ),
      class = "isoformic_mae_missing_genes"
    )
  }
}
