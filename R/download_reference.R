#' Download Reference Files From GENCODE
#'
#' This function facilitates the downloading of reference files from the GENCODE database.
#' It supports downloading GTF, GFF, and transcriptome FASTA files for human and mouse genomes.
#' The function ensures that the correct version and file type are downloaded and handles directory creation
#' and file existence checks to avoid redundant downloads.
#'
#' @param version Character with the version string. For mouse references,
#'   the letter 'M' in the version string is mandatory.
#'
#' @param reference Character indicating the source of the reference file.
#'   One of "gencode" or "mane". Defaults to "gencode".
#'   **NOTE:** "mane" is not implemented yet.
#'
#' @param organism Character indicating the organism.
#'   For GENCODE, this can only `"human"` or `"mouse"`.
#'
#' @param file_type Character indicating the type of file to download.
#'   One of `"gtf"`, `"gff"`, or `"fasta"`. Defaults to `"gtf"`.
#'   **NOTE:** `"fasta"` refers to the transcriptome FASTA.
#'
#' @param output_path Character specifying the directory where the
#'   downloaded file will be saved. Defaults to `"data-raw"`.
#'
#' @param timeout_limit Numeric value specifying the time in seconds for the
#'   download limit. This argument takes precedence over
#'   `base::options("timeout")`. Defaults to 3600 seconds (1 Hour).
#'
#' @param method Character specifying the method used by
#'   `utils::download.file()`. Defaults to `"auto"`
#'
#' @return A character string with the full path to the downloaded file.
#'
#' @export
download_reference <- function(version = "46",
                               reference = "gencode",
                               organism = c("human", "mouse"),
                               file_type = c("gtf", "gff", "fasta"),
                               output_path = "data-raw",
                               timeout_limit = 3600,
                               method = "auto") {
  if (requireNamespace("curl", quietly = TRUE)) {
    if (!isTRUE(curl::has_internet())) {
      cli::cli_abort(
        message = c(
          "x" = "No internet connection available."
        ),
        class = "isoformic_download_reference_no_internet"
      )
    }
  }
  reference <- stringr::str_to_lower(reference)
  reference <- rlang::arg_match(reference)
  organism <- stringr::str_to_lower(organism)
  organism <- rlang::arg_match(organism)
  file_type <- stringr::str_to_lower(file_type)
  file_type <- rlang::arg_match(file_type)
  version <- as.character(version)
  base_url <- switch(organism,
    human = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_",
    mouse = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_"
  )
  file_string <- switch(file_type,
    gff = ".annotation.gff3.gz",
    gtf = ".annotation.gtf.gz",
    fasta = ".transcripts.fa.gz"
  )
  download_url <- paste0(
    base_url, version, "/gencode.v", version, file_string
  )
  full_output_path <- paste0(
    output_path, "/gencode.v", version, file_string
  )
  if (!isTRUE(fs::dir_exists(output_path))) {
    fs::dir_create(output_path)
  }
  if (isTRUE(fs::file_exists(full_output_path))) {
    cli::cli_inform(
      c(
        `i` = "{.path {full_output_path}} already exists."
      )
    )
    return(full_output_path)
  }
  withr::with_options(
    new = base::list(
      timeout = base::max(
        base::unlist(base::options("timeout")),
        timeout_limit
      )
    ),
    code = {
      rlang::try_fetch(
        expr = {
          download_res <- utils::download.file(
            url = download_url,
            destfile = full_output_path,
            method = method
          )
        },
        error = function(cnd) {
          cli::cli_abort(
            message = c(
              `x` = "Failed to download file from {.url {download_url}}.",
              `x` = "Error: {cnd$message}"
            ),
            class = "isoformic_download_reference_error"
          )
        }
      )
    }
  )
  if (isTRUE(download_res == 0)) {
    cli::cli_inform(
      c(
        `v` = "{.path {full_output_path}} successfully downloaded."
      )
    )
  }
  return(full_output_path)
}
