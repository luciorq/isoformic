#' Check if Isoformic `mae` have a list `isoformic` slot
check_mae_isoformic_is_list <- function(mae) {
  if (!isTRUE(inherits(mae@metadata[["isoformic"]], "list"))) {
    cli::cli_abort(
      message = c(
        `x` = "{.var MultiAssayExperiment} object do {.strong not} contain {.var isoformic} in the metadata slot."
      ),
      class = "isoformic_mae_metadata_is_list"
    )
  }
}

#' Check if Isoformic `mae` have a `tx_to_gene`
check_mae_isoformic_tx_to_gene <- function(mae) {
  if (!isTRUE(("tx_to_gene" %in% names(mae@metadata[["isoformic"]])))) {
    cli::cli_abort(
      message = c(
        `x` = "{.var MultiAssayExperiment} object do {.strong not} contain {.var tx_to_gene} in the {.var metadata(mae)$isoformic} slot."
      ),
      class = "isoformic_mae_metadata_no_tx_to_gene"
    )
  }
}

