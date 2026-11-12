#!/usr/bin/env Rscript
# Windows only: install the R packages that the pixi environment cannot take
# from conda. bioconda has no Windows builds and conda-forge has no current
# Windows build of duckdb. Every package in `Imports` or `Suggests` of
# DESCRIPTION that the environment lacks, or has in a version older than the
# stated minimum, is installed as a binary package into the library of the
# environment, with the dependencies it still needs. The sources are:
#
#   - Bioconductor 3.22, the release for R 4.5 and the one bioconda builds
#     for Linux and macOS, so all platforms use the same Bioconductor
#     versions;
#   - https://luciorq.r-universe.dev, for condathis;
#   - CRAN, for duckdb and for dependencies that are not in the environment.
#
# `pixi.lock` does not record these packages. A Bioconductor release no
# longer changes once the next one is out, so the result is still stable.
# Set `ISOFORMIC_BIOC_REPOS` to a comma separated list of repositories to
# replace the Bioconductor ones, for example
# `https://bioc-release.r-universe.dev`.
#
# `pixi run r-packages` runs this script and the other Windows tasks depend
# on it. With nothing left to install it only reads the library.

# Keep in step with `r-base` in `pixi.toml`.
bioc_version <- "3.22"

# Package names and minimum versions from `Imports` and `Suggests`.
read_requirements <- function(path = "DESCRIPTION") {
  fields <- read.dcf(path, fields = c("Imports", "Suggests"))[1L, ]
  entries <- trimws(unlist(strsplit(fields[!is.na(fields)], ",", fixed = TRUE)))
  entries <- entries[nzchar(entries)]
  has_minimum <- grepl(">=", entries, fixed = TRUE)
  minimum <- rep(NA_character_, length(entries))
  minimum[has_minimum] <- sub(
    "^.*>=\\s*([[:digit:].-]+).*$",
    "\\1",
    entries[has_minimum]
  )
  data.frame(
    package = sub("\\s*\\(.*$", "", entries),
    minimum = minimum
  )
}

# Packages that are not installed in `lib`, or are older than the minimum.
find_unmet <- function(requirements, lib) {
  installed <- utils::installed.packages(lib.loc = lib)[, "Version"]
  unmet <- vapply(
    seq_len(nrow(requirements)),
    function(i) {
      version <- installed[requirements$package[i]]
      minimum <- requirements$minimum[i]
      is.na(version) ||
        (!is.na(minimum) && package_version(version) < minimum)
    },
    logical(1L)
  )
  requirements$package[unmet]
}

binary_repositories <- function() {
  bioc <- Sys.getenv("ISOFORMIC_BIOC_REPOS", unset = "")
  if (nzchar(bioc)) {
    bioc <- trimws(strsplit(bioc, ",", fixed = TRUE)[[1L]])
  } else {
    bioc <- sprintf(
      "https://bioconductor.org/packages/%s/%s",
      bioc_version,
      c("bioc", "data/annotation", "data/experiment")
    )
  }
  c(bioc, "https://luciorq.r-universe.dev", "https://cloud.r-project.org")
}

if (.Platform$OS.type != "windows") {
  stop(
    "This script is for Windows. On Linux and macOS `pixi install` takes ",
    "all packages from conda.",
    call. = FALSE
  )
}

lib <- .Library
requirements <- read_requirements()
to_install <- find_unmet(requirements, lib)

if (length(to_install) == 0L) {
  message("All packages in `Imports` and `Suggests` are installed.")
} else {
  message("Installing into ", lib, ": ", toString(to_install))
  options(timeout = max(600, getOption("timeout")))
  utils::install.packages(
    to_install,
    lib = lib,
    repos = binary_repositories(),
    type = "win.binary"
  )

  # `install.packages()` only warns when a package fails, and a binary
  # package built by CRAN or Bioconductor may not load in the R of conda.
  still_unmet <- find_unmet(requirements, lib)
  loads <- vapply(
    setdiff(to_install, still_unmet),
    requireNamespace,
    logical(1L),
    lib.loc = lib,
    quietly = TRUE
  )
  failed <- c(still_unmet, names(loads)[!loads])
  if (length(failed) > 0L) {
    stop(
      "Not installed or not loadable: ",
      toString(failed),
      call. = FALSE
    )
  }
  message("Installed and loaded: ", toString(to_install))
}
