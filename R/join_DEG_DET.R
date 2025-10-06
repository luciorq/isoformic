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
join_DEG_DET <- function(
    DEG_tab,
    DET_final_tab,
    logfc_cut,
    pval_cut) {
  .data <- rlang::.data
  DEG_tab_mod <- DEG_tab |>
    dplyr::rename(id = "gene_id")
  DEG_tab_mod <- DEG_tab_mod |>
    dplyr::rename(name = "gene_name")
  DEG_tab_mod <- DEG_tab_mod |>
    dplyr::mutate(transcript_type = "gene")
  DEG_tab_mod <- DEG_tab_mod |>
    dplyr::mutate(gene_name = .data[["name"]])

  DET_final_tab <- DET_final_tab |>
    dplyr::filter(
      .data[["transcript_type"]] %in% c(
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

  drop_columns <- c("DEG_sig")
  if (any(colnames(DET_final_tab) %in% drop_columns)) {
    DET_final_tab_mod <- DET_final_tab |>
      dplyr::select(-dplyr::one_of(drop_columns))
  } else {
    DET_final_tab_mod <- DET_final_tab
  }

  DET_final_tab_mod <- DET_final_tab_mod |>
    dplyr::rename(id = "transcript_id")
  DET_final_tab_mod <- DET_final_tab_mod |>
    dplyr::rename(name = "transcript_name")
  DEG_tab_mod <- DEG_tab_mod[
    colnames(DEG_tab_mod)[colnames(DEG_tab_mod)
    %in% colnames(DET_final_tab_mod)]
  ]

  DEGs_DETs_table <- dplyr::bind_rows(DEG_tab_mod, DET_final_tab_mod)
  DEGs_DETs_table$significance <- c()
  DEGs_DETs_table$abs_log2FC <- base::abs(DEGs_DETs_table$log2FC)
  DEGs_DETs_table$significance <- "not_sig"
  DEGs_DETs_table$DEG_sig <- "NO"
  DEGs_DETs_table$significance[DEGs_DETs_table$abs_log2FC > logfc_cut &
    DEGs_DETs_table$pvalue < pval_cut] <- "sig"
  DEGs_DETs_table$DEG_sig[DEGs_DETs_table$abs_log2FC > logfc_cut &
    DEGs_DETs_table$pvalue < pval_cut] <- "YES"
  return(DEGs_DETs_table)
}
