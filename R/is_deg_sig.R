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
