#' Create Isoformic Object from SummarizedExperiment
#'
#' @param se_tx `SummarizedExperiment` object containing **transcript-level**
#' differential expression analysis results.
#'
#' @param se_gene `SummarizedExperiment` object containing **gene-level**
#' differential expression analysis results.
#'
#' @export
create_isoformic_mae_from_se <- function(
    se_tx,
    se_gene,
    tx_to_gene = NULL) {
  # This condition is automatically satisfied if `se_gene` was generated
  # + with `tximport::summarizeToGene()`.
  if (isTRUE(
    identical(
      SummarizedExperiment::colData(se_tx),
      SummarizedExperiment::colData(se_gene)
    )
  )) {
    mae_coldata <- SummarizedExperiment::colData(se_tx)
  }
  if (!isTRUE(
    all(rownames(
      SummarizedExperiment::colData(se_tx)
    ) %in% rownames(SummarizedExperiment::colData(se_gene)))
  )) {
    cli::cli_abort(
      message = c(
        `x` = "{.code se_tx} and {.code se_gene} do {.strong not} contain the same {.code rownames(colData(se))}."
      ),
      class = "isoformic_se_experiment_rownames"
    )
  }
  mae_assays_list <- list(
    transcript = se_tx,
    gene = se_gene
  )
  mae <- MultiAssayExperiment::MultiAssayExperiment(
    experiments = mae_assays_list,
    colData = mae_coldata
  )
  check_mae_assay_rownames(mae)

  # Check for `tx_to_gene`
  row_data_tx_se <- SummarizedExperiment::rowData(se_tx)
  tx_id_vector <- rownames(row_data_tx_se)
  gene_id_vector <- NULL
  if (isTRUE(is.null(tx_to_gene))) {
    if (isTRUE("gene_id" %in% colnames(row_data_tx_se))) {
      gene_id_vector <- unlist(row_data_tx_se$gene_id)
    }
  }
  if (isTRUE(is.null(gene_id_vector))) {
    cli::cli_abort(
      message = c(
        `x` = "Variable {.var gene_id_vector} is not available."
      ),
      class = "isoformic_se_gene_id_vector_null"
    )
  }
  tx_to_gene_df <- data.frame(
    tx_id = tx_id_vector,
    gene_id = gene_id_vector
  )
  rownames(tx_to_gene_df) <- tx_to_gene_df$tx_id
  mae@metadata[["isoformic"]][["tx_to_gene"]] <- S4Vectors::DataFrame(tx_to_gene_df)
  check_mae_tx_in_genes(mae)
  return(mae)
}
