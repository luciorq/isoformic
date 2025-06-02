
#' Run Gene Set Enrichment Analysis for Different Transcript Types
#'
#' Performs gene set enrichment analysis (GSEA) on differential expression results for various transcript types,
#' using the `fgsea` package. The function iterates over specified transcript types, filters the data accordingly,
#' and runs GSEA for each type.
#'
#' @param det_df A `data.frame` or `tibble` containing transcript-level differential expression results,
#'   including `transcript_type`, `log2FC`, and `gene_name` columns.
#' @param genesets_list A list of gene sets to be used in the enrichment analysis.
#' @param pval_cutoff A numeric value specifying the p-value cutoff for the enrichment results. Default is `0.05`.
#' @param lfc_cutoff A numeric value specifying the log2 fold-change cutoff for filtering transcripts. Default is `1`.
#'
#' @return A `tibble` containing the enrichment analysis results for each transcript type, including pathway names,
#'   p-values, adjusted p-values, and the transcript type (experiment).
#'
#' @details
#' The function defines a list of transcript types and their corresponding labels.
#' It then filters the input differential expression data for each transcript type, ranks the genes by log2 fold-change,
#' and performs GSEA using the `fgsea` package.
#'
#' @examples
#' # Sample differential expression data
#' det_df <- data.frame(
#'   gene_name = c("GeneA", "GeneB", "GeneC", "GeneD"),
#'   transcript_type = c(
#'     "protein_coding", "retained_intron", "protein_coding_CDS_not_defined",
#'     "processed_transcript", "nonsense_mediated_decay"
#'   ),
#'   log2FC = c(1.5, -2.0, 0.8, -1.2)
#' )
#'
#' # Sample gene sets
#' genesets_list <- list(
#'   Pathway1 = c("GeneA", "GeneC"),
#'   Pathway2 = c("GeneB", "GeneD")
#' )
#'
#' # Run enrichment analysis
#' fgsea_results_df <- run_enrichment(
#'   det_df = det_df,
#'   genesets_list = genesets_list,
#'   pval_cutoff = 0.05,
#'   lfc_cutoff = 1
#' )
#'
#' # View the results
#' print(fgsea_results_df)
#'
#' @export
run_enrichment <- function(
    det_df,
    genesets_list,
    tx_to_gene,
    pval_cutoff = 0.05,
    lfc_cutoff = 1) {
  .data <- rlang::.data
  processed_or_cds <- ifelse(
    test = sum(tx_to_gene$transcript_type == "processed_transcript") > 1500,
    yes = "processed_transcript",
    no = "protein_coding_CDS_not_defined"
  )
  tx_types_list <- list(
    "protein_coding",
    c(
      "retained_intron",
      processed_or_cds,
      "nonsense_mediated_decay"
    ),
    "retained_intron",
    processed_or_cds,
    "nonsense_mediated_decay"
  )

  tx_type_names <- c(
    "protein_coding",
    "unproductive",
    "retained_intron",
    processed_or_cds,
    "nonsense_mediated_decay"
  )
  fgsea_results_df <- base::seq_along(tx_type_names) |>
    purrr::map(
      .f = function(x) {
        type_vec <- tx_types_list[[x]]
        type_name <- tx_type_names[x]
        res_fgsea <- det_df |>
          dplyr::filter(transcript_type %in% type_vec) |>
          dplyr::filter(!is.na(log2FC)) |>
          dplyr::arrange(log2FC)
        # TODO: @luciorq add option for log2FC * log10pvalue
        ranks <- res_fgsea$log2FC
        names(ranks) <- res_fgsea$transcript_name
        genesets_list <- tibble::enframe(genesets_list, name = "term", value = "gene_name") |>
          tidyr::unnest(cols = gene_name)
        genesets_list <- genesets_list |>
          dplyr::left_join(tx_to_gene |> dplyr::select(transcript_name, gene_name), by="gene_name")
        genesets_list <- split(genesets_list$transcript_name,
                               genesets_list$term)
        fgsea_results <- fgsea::fgseaMultilevel(
          pathways = genesets_list,
          stats = ranks,
          eps = 0
        ) |>
          #tibble::as_tibble() |>
          dplyr::filter(pval < pval_cutoff) |>
          dplyr::mutate(experiment = type_name)
        return(fgsea_results)
      }
    ) |> purrr::set_names(tx_type_names)

  fgsea_combined <- dplyr::bind_rows(fgsea_results_df)

  return(fgsea_combined)
}
