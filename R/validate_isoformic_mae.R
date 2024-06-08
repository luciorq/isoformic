# Validate Isoformic Object
#' Validate if `{isoformic}` mae has been run for the object.
#' @param obj mae
#' @export
validate_isoformic_mae <- function(mae) {
  check_mae(mae)
  experiment_names_vec <- c("transcript", "gene")
  if (!isTRUE(all(experiment_names_vec %in% names(mae@ExperimentList)))) {
    cli::cli_abort(
      message = c(
        x = "{.field Names} of the {.code Experiments} slot of the {.code MultiAssayExperiment} object are {.strong not} contained in {.code {experiment_names_vec}}."
      ),
      class = "isoformic_mae_experiment_type"
    )
  }

  # Check if underlying `SummarizedExperiment`s have metadata level in
  # + c("txp", "gene")
  for (i in c("transcript", "gene")) {
    se_experiment_level(mae@ExperimentList[[i]])
  }
  check_mae_isoformic_is_list(mae)
  check_mae_isoformic_tx_to_gene(mae)
}
