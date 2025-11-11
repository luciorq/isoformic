#' Check of Transient Dependencies are installed
#'
#' Iterate over packages and check if they are installed.
#' If running on non interactive session it will return just TRUE or FALSE.
#' If running on interactive session it will prompt the user to install missing
#' packages.
#' And if installation is succesfull it will return TRUE, otherwise FALSE.
#'
#' @keywords internal
#' @noRd
check_installed <- function(pkgs) {
  missing_pkgs <- character(0L)
  for (pkg in pkgs) {
    if (!rlang::is_installed(pkg)) {
      missing_pkgs <- c(missing_pkgs, pkg)
    }
  }
  if (isTRUE(length(missing_pkgs) == 0L)) {
    return(invisible(TRUE))
  } else {
    return_res <- FALSE
  }

  if (rlang::is_interactive()) {
    return_res <- rlang::try_fetch(
      expr = {
        rlang::check_installed(
          pkg = missing_pkgs,
          call = rlang::caller_env(n = 2L),
        )
        TRUE
      },
      error = function(cnd) {
        FALSE
      }
    )
  }

  return(invisible(return_res))
}
