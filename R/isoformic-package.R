#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom rlang .data
#' @importFrom rlang .env
#' @importFrom rlang :=
#' @importFrom rlang %||%
## usethis namespace: end
NULL

.onLoad <- function(...) {
  S7::methods_register()
}
