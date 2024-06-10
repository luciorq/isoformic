#' Run Swish Differential Expression
#'
#' Run swish method for a `SummarizedExperiment` with inferential replicates.
#'
#' @param se `SummarizedExperiment` object.
#' @param contrast_var column name from `colData(se)`.
#' @param ... Additional arguments passed to `fishpond::swish()`.
#'
#' @return `SummarizedExperiment` object.
#'
#' @export
run_swish_pairwise <- function(se, contrast_var = "condition", ...) {
  check_se(se)
  # TODO: @luciorq Create a check for scaled inf reps before
  # + running `fishpond::scaleInfReps()`.
  if (isTRUE(se@metadata[["infRepsScaled"]])) {
    se <- fishpond::scaleInfReps(se)
  }
  se <- fishpond::labelKeep(se, minCount = 10, minN = 3, x = contrast_var)
  se <- fishpond::swish(
    y = se,
    x = contrast_var,
    ...
  )

  se@metadata[["isoformic"]][["dea"]] <- S4Vectors::mcols(se)[, c("stat", "log2FC", "locfdr", "qvalue")]
  return(se)
}
