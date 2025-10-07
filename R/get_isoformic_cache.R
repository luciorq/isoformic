#' Retrieve System-Dependent Cache Path for Isoformic
#'
#' Determines the appropriate user cache directory for the `isoformic` package
#' based on the operating system.
#' On macOS, it avoids using paths with spaces and follows the XDG base
#' directory specification.
#'
#' @details
#' This function uses the `[tools::R_user_dir()]` function to determine the
#' user cache directory.
#'
#' @return
#' A path character string representing the path to the user cache
#' directory for the `isoformic` package.
#'
#' @export
get_isoformic_cache <- function(..., ext = NULL) {
  rlang::check_dots_unnamed()
  if (isFALSE(rlang::is_null(ext))) {
    ext <- stringr::str_remove(ext, pattern = "^\\.")
    ext <- stringr::str_to_lower(ext)
  } else {
    ext <- ""
  }
  if (
    isTRUE(Sys.getenv(x = "XDG_CACHE_HOME") == "") &&
      isTRUE(stringr::str_detect(get_sys_arch(), pattern = "^Darwin"))
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
