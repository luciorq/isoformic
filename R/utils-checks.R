#' Assert That a File Path Exists
#'
#' @param path Character string with the path to check.
#' @param arg Name of the argument being checked, used in the error message.
#'
#' @keywords internal
#' @noRd
assert_file_exists <- function(path, arg = rlang::caller_arg(path)) {
  if (!rlang::is_string(path) || !isTRUE(fs::file_exists(path))) {
    cli::cli_abort(
      message = c(
        x = "{.arg {arg}} must be the path to an existing file.",
        i = "Received: {.path {path}}."
      ),
      class = "isoformic_annot_file_dont_exist"
    )
  }
  return(invisible(TRUE))
}

#' Assert That a Directory Exists
#' @keywords internal
#' @noRd
assert_dir_exists <- function(path, arg = rlang::caller_arg(path)) {
  if (!rlang::is_string(path) || !isTRUE(fs::dir_exists(path))) {
    cli::cli_abort(
      message = c(
        x = "{.arg {arg}} must be the path to an existing directory.",
        i = "Received: {.path {path}}."
      ),
      class = "isoformic_dir_dont_exist"
    )
  }
  return(invisible(TRUE))
}

#' Check Whether Optional Dependencies Are Installed
#'
#' Iterates over packages and checks whether they are installed.
#' In a non-interactive session it returns `TRUE` or `FALSE`.
#' In an interactive session it prompts the user to install missing packages
#' and returns `TRUE` if the installation succeeded.
#'
#' @keywords internal
#' @noRd
check_installed <- function(pkgs) {
  missing_pkgs <- pkgs[!vapply(pkgs, rlang::is_installed, logical(1L))]
  if (identical(length(missing_pkgs), 0L)) {
    return(invisible(TRUE))
  }
  return_res <- FALSE
  if (rlang::is_interactive()) {
    return_res <- rlang::try_fetch(
      expr = {
        rlang::check_installed(
          pkg = missing_pkgs,
          call = rlang::caller_env(n = 2L)
        )
        TRUE
      },
      error = function(cnd) FALSE
    )
  }
  return(invisible(return_res))
}

#' Coerce to a Scalar String or Abort
#' @keywords internal
#' @noRd
check_string <- function(x, arg = rlang::caller_arg(x), allow_null = FALSE) {
  if (is.null(x) && isTRUE(allow_null)) {
    return(invisible(NULL))
  }
  if (!rlang::is_string(x)) {
    cli::cli_abort(
      message = c(x = "{.arg {arg}} must be a single character string."),
      class = "isoformic_not_string"
    )
  }
  return(invisible(x))
}
