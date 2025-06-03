#' Summarize Transcript-level Expression to Gene-level
#'
#' This function aggregates transcript-level expression data to gene-level by calculating the mean expression of transcripts belonging to the same gene.
#'
#' @param txi_transcript A `tibble` containing transcript-level expression abundances.
#' @param tx_to_gene A `data.frame` or `tibble` containing transcript-to-gene mapping information, including `transcript_id` and `gene_id` columns.
#'
#' @return A `tibble` containing gene-level expression abundances.
#'
#' @keywords internal
#'
#' @noRd
summarize_to_gene <- function(txi_transcript, tx_to_gene) {
  .data <- rlang::.data
  id_df <- tx_to_gene |>
    dplyr::select("transcript_id", "gene_id") |>
    dplyr::distinct()

  txi_gene <- id_df |>
    dplyr::left_join(txi_transcript, by = "transcript_id") |>
    dplyr::select(-c("transcript_id", "transcript_name")) |>
    tidyr::pivot_longer(
      -c("gene_id"),
      names_to = "samples",
      values_to = "expr"
    ) |>
    dplyr::group_by(gene_id, samples) |>
    dplyr::summarise(
      mean_expr = base::mean(.data$expr, na.rm = TRUE)
    ) |>
    dplyr::ungroup() |>
    dplyr::filter(!is.na(.data$mean_expr))
  txi_gene <- txi_gene |>
    tidyr::pivot_wider(
      names_from = "samples",
      values_from = "mean_expr"
    )
  return(txi_gene)
}
