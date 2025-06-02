#' Annotate Transcripts with Differential Gene Expression Significance
#'
#' Adds a column to a transcript-level differential expression table indicating whether each transcript
#' originates from a gene that is significantly differentially expressed.
#'
#' @param DegsigVector A character vector containing the names of transcripts from significantly differentially expressed genes.
#' @param DET_table A `data.frame` or `tibble` containing transcript-level differential expression results,
#'   including a `transcript_name` column.
#'
#' @return A `tibble` with an additional column `DEG_sig` indicating whether the transcript is from a significantly
#'   differentially expressed gene (`"YES"` or `"NO"`).
#'
#' @examples
#' # Sample data
#' significant_transcripts <- c("transcript1", "transcript3")
#' DET_table <- data.frame(
#'   transcript_name = c("transcript1", "transcript2", "transcript3", "transcript4"),
#'   log2FC = c(2.5, -1.2, 0.8, -0.5),
#'   pvalue = c(0.01, 0.2, 0.03, 0.6)
#' )
#'
#' # Annotate transcripts with DEG significance
#' DET_table_annotated <- is_deg_sig(DegsigVector = significant_transcripts, DET_table = DET_table)
#'
#' # View the result
#' print(DET_table_annotated)
#'
#' @export
is_deg_sig <- function(DegsigVector, DET_table) {
  DETs_DEGs <- DET_table |>
    dplyr::filter(transcript_name %in% DegsigVector)
  DETs_DEGs$DEG_sig <- "YES"
  DETs_notDEGs <- DET_table |>
    dplyr::filter(!transcript_name %in% DegsigVector)
  DETs_notDEGs$DEG_sig <- "NO"
  DET_table_final <- dplyr::bind_rows(DETs_DEGs, DETs_notDEGs)
  return(DET_table_final)
}


#' Merge Gene and Transcript Level Differential Expression Tables
#'
#' Combines gene-level and transcript-level differential expression results into a single table,
#' annotates the combined data with significance labels based on specified cutoffs, and filters
#' transcripts based on their types.
#'
#' @param DEG_tab A `data.frame` or `tibble` containing gene-level differential expression results,
#'   including `gene_id`, `gene_name`, `log2FC`, and `pvalue` columns.
#' @param DET_final_tab A `data.frame` or `tibble` containing transcript-level differential expression results,
#'   including `transcript_id`, `transcript_name`, `transcript_type`, `log2FC`, and `pvalue` columns.
#' @param logfc_cut A numeric value specifying the absolute log2 fold-change cutoff for significance.
#' @param pval_cut A numeric value specifying the p-value cutoff for significance.
#'
#' @return A `tibble` combining gene and transcript differential expression results, with additional columns:
#'   - `id`: gene or transcript ID.
#'   - `name`: gene or transcript name.
#'   - `transcript_type`: type of transcript or `"gene"` for gene-level entries.
#'   - `abs_log2FC`: absolute value of log2 fold-change.
#'   - `significance`: `"sig"` if significant based on cutoffs, `"not_sig"` otherwise.
#'
#' @examples
#' # Sample gene-level data
#' DEG_tab <- data.frame(
#'   gene_id = c("gene1", "gene2"),
#'   gene_name = c("GeneA", "GeneB"),
#'   log2FC = c(1.5, -2.0),
#'   pvalue = c(0.01, 0.04)
#' )
#'
#' # Sample transcript-level data
#' DET_final_tab <- data.frame(
#'   transcript_id = c("tx1", "tx2", "tx3"),
#'   transcript_name = c("Transcript1", "Transcript2", "Transcript3"),
#'   transcript_type = c("protein_coding", "lncRNA", "processed_transcript"),
#'   log2FC = c(1.2, -1.8, 0.5),
#'   pvalue = c(0.02, 0.03, 0.2)
#' )
#'
#' # Merge and annotate differential expression results
#' DEGs_DETs_table <- join_DEG_DET(
#'   DEG_tab = DEG_tab,
#'   DET_final_tab = DET_final_tab,
#'   logfc_cut = 1,
#'   pval_cut = 0.05
#' )
#'
#' # View the result
#' print(DEGs_DETs_table)
#'
#' @export
join_DEG_DET <- function(DEG_tab, DET_final_tab, logfc_cut, pval_cut) {
  DEG_tab_mod <- DEG_tab |>
    dplyr::rename(id = gene_id)
  DEG_tab_mod <- DEG_tab_mod |>
    dplyr::rename(name = gene_name)
  DEG_tab_mod <- DEG_tab_mod |>
    dplyr::mutate(transcript_type = "gene")
  DEG_tab_mod <- DEG_tab_mod |>
    dplyr::mutate(gene_name = DEG_tab_mod$name)
  drop_columns <- c("DEG_sig")
  DET_final_tab <- DET_final_tab |>
    dplyr::filter(
      transcript_type %in% c(
        "retained_intron", "protein_coding_CDS_not_defined",
        "processed_transcript",
        "nonsense_mediated_decay",
        "lncRNA", "protein_coding",
        "pseudogene", "non_stop_decay",
        "processed_pseudogene",
        "transcribed_unprocessed_pseudogene",
        "transcribed_unitary_pseudogene",
        "unprocessed_pseudogene",
        "unitary_pseudogene"
      )
    )
  DET_final_tab_mod <- DET_final_tab |>
    dplyr::select(-dplyr::one_of(drop_columns))
  DET_final_tab_mod <- DET_final_tab_mod |>
    dplyr::rename(id = transcript_id)
  DET_final_tab_mod <- DET_final_tab_mod |>
    dplyr::rename(name = transcript_name)
  DEG_tab_mod <- DEG_tab_mod[
    colnames(DEG_tab_mod)[colnames(DEG_tab_mod)
    %in% colnames(DET_final_tab_mod)]
  ]
  DEGs_DETs_table <- dplyr::bind_rows(DEG_tab_mod, DET_final_tab_mod)
  DEGs_DETs_table$significance <- c()
  DEGs_DETs_table$abs_log2FC <- base::abs(DEGs_DETs_table$log2FC)
  DEGs_DETs_table$significance <- "not_sig"
  DEGs_DETs_table$significance[DEGs_DETs_table$abs_log2FC > logfc_cut &
    DEGs_DETs_table$pvalue < pval_cut] <- "sig"
  return(DEGs_DETs_table)
}




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
run_enrichment <- function(det_df,
                           genesets_list,
                           tx_to_gene,
                           pval_cutoff = 0.05,
                           lfc_cutoff = 1) {
  processed_or_cds <- ifelse(sum(tx_to_gene$transcript_type == "processed_transcript") > 1500,
                             "processed_transcript",
                             "protein_coding_CDS_not_defined")
  tx_types_list <- list(
    "protein_coding",
    c(
      "retained_intron", processed_or_cds,
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
