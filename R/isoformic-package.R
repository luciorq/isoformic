#' @keywords internal
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL

.onLoad <- function(...) {
  # S7::S4_register(IsoformicExperiment)
  S7::methods_register()
}

# enable usage of <S7_object>@name in package code
#' @rawNamespace if (getRversion() < "4.3.0") importFrom("S7", "@")
NULL
