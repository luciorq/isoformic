# Assemble the snapshot of processed files deposited on Zenodo.
#
# The bundle holds every frozen evidence dataset (Parquet), the coding
# potential model, an MD5SUMS file and a README with versions and sources.
# Run from the package root after `data-raw/build_frozen_data.R`:
#   Rscript data-raw/zenodo_deposit.R [output_dir]
# Default output: get_isoformic_cache("zenodo_deposit"). Then:
#   1. create a new version of the Zenodo deposit and upload the files of
#      the bundle directory (keep the file names; the package looks them up
#      by name: https://zenodo.org/records/<id>/files/<file>?download=1)
#   2. publish and copy the record id of the new version
#   3. set the id in frozen_zenodo_record() (R/frozen-data.R), cite the DOI
#      in NEWS.md and in inst/CITATION, bump the package version
# `fetch_frozen_data()` then falls back to the Zenodo files whenever the
# GitHub release asset is missing, and users can pin a record with
# options(isoformic.frozen_mirrors = "https://zenodo.org/records/<id>/files/").

devtools::load_all(".", quiet = TRUE)
args <- commandArgs(trailingOnly = TRUE)
output_dir <- if (length(args)) {
  args[[1L]]
} else {
  get_isoformic_cache("zenodo_deposit")
}
bundle_frozen_data(output_dir, build = "never")
cat("Deposit bundle:", output_dir, "\n")
print(readLines(fs::path(output_dir, "MD5SUMS")))
