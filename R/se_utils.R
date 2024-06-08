# ==============================================================================
# Checks for `SummarizedExperiment` objects
# ==============================================================================

# Check if object is of class `SummarizedExperiment`.
check_se <- function(se) {
  if (!is_se(se)) {
    cli::cli_abort(
      message = c(
        x = "The object is not of type {.code SummarizedExperiment}."
      ),
      class = "isoformic_obj_not_se"
    )
  }
}

# Check if object is of class `SummarizedExperiment`.
is_se <- function(se) {
  if (isTRUE(inherits(x = se, what = c("SummarizedExperiment")))) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

# Check if `se` has a `metadata` element named `"level"`.
check_se_metadata_level <- function(se) {
  if (!isTRUE("level" %in% names(se@metadata))) {
    cli::cli_abort(
      message = c(
        x = "The {.code metadata} slot of the {.code SummarizedExperiment} object do {.strong not} contain a {.field level} named field."
      ),
      class = "isoformic_se_metadata_without_level"
    )
  }
}

#' Check if content of metadata level slot is a string.
check_se_metadata_level_type <- function(se) {
  if (!isTRUE(inherits(se@metadata$level, what = "character"))) {
    cli::cli_abort(
      message = c(
        x = "{.field level} of the {.code metadata} slot of the {.code SummarizedExperiment} object is {.strong not} a {.type char}."
      ),
      class = "isoformic_se_metadata_level_not_char"
    )
  }
}
#' Check if content of metadata level slot is a string.
check_se_metadata_level_length <- function(se) {
  if (!isTRUE(length(se@metadata$level) == 1)) {
    cli::cli_abort(
      message = c(
        x = "{.field level} of the {.code metadata} slot of the {.code SummarizedExperiment} object is {.strong not} of length {.strong 1}."
      ),
      class = "isoformic_se_level_length"
    )
  }
}

#' Return `SummarizedExperiment` Experiment Level
#'
#' Check level of quantification for
#'
#' @param se object of class `SummarizedExperiment`.
#'
#' @param type String containing one of `c("txp", "gene")`.
se_experiment_level <- function(se) {
  check_se(se)
  check_se_metadata_level(se)
  check_se_metadata_level_length(se)
  check_se_metadata_level_type(se)

  level_str <- se@metadata$level
  level_vec <- c("txp", "gene")
  if (!isTRUE(level_str %in% level_vec)) {
    cli::cli_abort(
      message = c(
        x = "Experiment level should be one of {.var {level_vec}}"
      ),
      class = "isoformic_se_level_not_txp_or_gene"
    )
  }
  return(level_str)
}
