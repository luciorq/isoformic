# Build (or rebuild) every frozen dataset from its sources and print the MD5
# checksums to record in R/frozen-data.R (frozen_registry()).
#
# Run from the package root after updating a source version in the registry:
#   Rscript data-raw/build_frozen_data.R
# Steps to publish (see data-raw/README.md):
#   1. update `version` and `sources` of the dataset in frozen_registry()
#   2. run this script; it downloads the sources into
#      get_isoformic_cache("sources", <name>, <version>) and writes the Parquet files
#      into get_isoformic_cache("frozen")
#   3. copy the printed MD5 values into frozen_registry()
#   4. upload the Parquet files as assets of the GitHub release tagged
#      `frozen-data` (gh release upload frozen-data <files> --clobber)
#   5. deposit the snapshot on Zenodo: `Rscript data-raw/zenodo_deposit.R`
#      writes the bundle, upload it and record the record id in
#      frozen_zenodo_record() (R/frozen-data.R)
#   6. run the test suite and record the change in NEWS.md
# Requirements: internet access, ~4 GB of disk for the GTEx source file and
# the `readxl` package for the NMD exon table.

devtools::load_all(".", quiet = TRUE)
args <- commandArgs(trailingOnly = TRUE)
datasets <- if (length(args)) args else frozen_datasets()$name
paths <- vapply(
  datasets,
  function(name) {
    cli::cli_h2(name)
    fetch_frozen_data(name, build = "always", verify = FALSE)
  },
  character(1L)
)
md5 <- tools::md5sum(paths)
cat("\nMD5 checksums to record in frozen_registry():\n")
for (name in names(paths)) {
  cat(sprintf(
    "  %-22s %s  (%s, %.1f MB)\n",
    name,
    unname(md5[[paths[[name]]]]),
    fs::path_file(paths[[name]]),
    fs::file_size(paths[[name]]) / 1e6
  ))
}
