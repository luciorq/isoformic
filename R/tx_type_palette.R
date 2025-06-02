#' @internal
tx_type_palette <- function() {
  fixed_tx_biotypes <- c(
    "gene", "protein_coding", "retained_intron",
    "protein_coding_CDS_not_defined", "nonsense_mediated_decay",
    "lncRNA", "processed_pseudogene",
    "transcribed_unprocessed_pseudogene",
    "unprocessed_pseudogene", "non_stop_decay",
    "transcribed_unitary_pseudogene",
    "pseudogene", "unitary_pseudogene", "processed_transcript"
  )
  tx_type_color_names <- c(
    "#fb8072", "#a6d854", "#8da0cb", "#fc8d62",
    "#66c2a5", "#e78ac3", "#ffd92f", "#e5c494",
    "#d9d9d9", "#d9d9d9", "#d9d9d9", "#ffffb3",
    "#d9d9d9", "#d9d9d9"
  )
  names(tx_type_color_names) <- fixed_tx_biotypes
  return(tx_type_color_names)
}
