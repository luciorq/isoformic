#' Convert to a tibble with transcript_id column
#'
#' @keywords internal
#' @noRd
convert_to_isoformic_tibble <- function(txi_transcript) {
  if (isFALSE(inherits(x = txi_transcript, what = c("tbl_df")))) {
    if (
      isTRUE(inherits(x = txi_transcript, what = c("matrix"))) ||
        isTRUE(inherits(x = txi_transcript, what = c("data.frame")))
    ) {
      txi_transcript <- as.data.frame(txi_transcript)
    }
    if (isFALSE("transcript_id" %in% colnames(txi_transcript))) {
      txi_transcript <- txi_transcript |>
        tibble::rownames_to_column(var = "transcript_id") |>
        tibble::as_tibble()
    }
  }
  return(txi_transcript)
}
