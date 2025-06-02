#' Plot Log2 Fold-Change Results for Selected Genes
#'
#' Creates a bar plot of log2 fold-change values for transcripts of a selected gene,
#' differentiating transcript types and significance levels.
#'
#' @param DEG_DET_table A `data.frame` or `tibble` containing combined gene and transcript differential expression results,
#'   including `name`, `log2FC`, `transcript_type`, `significance`, and `gene_name` columns.
#' @param selected_gene A character string specifying the gene name to plot.
#' @param custom_colors An optional named vector of colors for different transcript types.
#'
#' @return A `ggplot2` object representing the bar plot.
#'
#' @details
#' The function filters the input table for the selected gene and creates a bar plot of log2 fold-change values.
#' If all transcripts are significant, it plots without adjusting alpha transparency; otherwise, it adjusts alpha
#' based on significance. The function uses predefined colors for transcript types, which can be overridden
#' by providing `custom_colors`.
#'
#' @examples
#' # Sample data
#' DEGs_DETs_table <- data.frame(
#'   name = c("Transcript1", "Transcript2", "GeneA"),
#'   log2FC = c(1.5, -2.0, 0.8),
#'   transcript_type = c("protein_coding", "lncRNA", "gene"),
#'   significance = c("sig", "not_sig", "sig"),
#'   gene_name = c("GeneA", "GeneA", "GeneA")
#' )
#'
#' # Plot log2 fold-change for the selected gene
#' plot_obj <- plot_log2FC(
#'   DEG_DET_table = DEGs_DETs_table,
#'   selected_gene = "GeneA"
#' )
#'
#' # Display the plot
#' print(plot_obj)
#'
#' @export
plot_log2FC <- function(
    DEG_DET_table,
    selected_gene,
    custom_colors = NULL) {
  .data <- rlang::.data
  if (
    all(
      DEG_DET_table$significance[DEG_DET_table$gene_name %in% selected_gene] == "sig"
    )
  ) {
    plot_obj <- ggplot2::ggplot(
      data = DEG_DET_table[DEG_DET_table$gene_name %in% selected_gene, ],
      mapping = ggplot2::aes(
        x = .data$name,
        y = .data$log2FC,
        fill = .data$transcript_type
      )
    ) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::scale_fill_manual(values = tx_type_palette()) +
      ggplot2::theme_bw()
  } else {
    plot_obj <- ggplot2::ggplot(
      data = DEG_DET_table[DEG_DET_table$gene_name %in% selected_gene, ],
      mapping = ggplot2::aes(
        x = .data$name,
        y = .data$log2FC,
        fill = .data$transcript_type,
        alpha = factor(
          .data$significance,
          levels = c("not_sig", "sig")
        )
      )
    ) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::scale_fill_manual(values = tx_type_palette()) +
      ggplot2::scale_alpha_manual(
        values = c(0.2, 1),
        labels = c("not_sig", "sig")
      ) +
      ggplot2::labs(alpha = "significance") +
      ggplot2::theme_bw()
  }
  return(plot_obj)
  # ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 0.5))
}
