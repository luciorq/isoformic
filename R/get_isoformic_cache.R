#' Retrieve System-Dependent Cache Path for Isoformic
#'
#' Determines the appropriate user cache directory for the `isoformic` package
#' based on the operating system.
#' On macOS, it avoids using paths with spaces and follows the XDG base
#' directory specification.
#'
#' @details
#' This function uses [tools::R_user_dir()] to determine the user cache
#' directory, so `R_USER_CACHE_DIR` and `XDG_CACHE_HOME` move it (see
#' [isoformic_options]). The directory is not created.
#'
#' @param ... Additional path components to append to the cache directory.
#' @param ext An optional file extension (e.g., "rds", "csv") to append to
#' the final path.
#'
#' @returns
#' A path character string representing the path to the user cache
#' directory for the `isoformic` package.
#'
#' @examples
#' get_isoformic_cache()
#' get_isoformic_cache("frozen")
#'
#' @export
get_isoformic_cache <- function(..., ext = NULL) {
  rlang::check_dots_unnamed()
  if (isFALSE(rlang::is_null(ext))) {
    ext <- sub("^\\.", "", ext)
    ext <- tolower(ext)
  } else {
    ext <- ""
  }
  if (
    !nzchar(Sys.getenv(x = "XDG_CACHE_HOME")) &&
      isTRUE(startsWith(get_sys_arch(), "Darwin"))
  ) {
    withr::local_envvar(
      .new = list(
        `XDG_CACHE_HOME` = fs::path_home(".cache")
      )
    )
  }
  dir_path <- tools::R_user_dir(package = "isoformic", which = "cache")
  return(fs::path(dir_path, ..., ext = ext))
}

#' Retrieve Operating System and CPU Architecture
#' @keywords internal
#' @noRd
get_sys_arch <- function() {
  os <- base::Sys.info()["sysname"]
  cpu_arch <- base::Sys.info()["machine"]
  return(base::paste0(os, "-", cpu_arch))
}
