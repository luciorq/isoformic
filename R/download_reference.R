#' Download Reference Files from GENCODE
#'
#' Downloads reference annotation files from the GENCODE database for human or mouse genomes.
#' Supports downloading GTF, GFF, and transcriptome FASTA files. The function handles directory
#' creation and checks for existing files to avoid redundant downloads.
#'
#' @param version A character string specifying the GENCODE release version.
#'   For mouse references, include the letter 'M' in the version string (e.g., `"M32"`).
#'   Default is `"46"`.
#' @param reference A character string specifying the source of the reference file.
#'   Currently, only `"gencode"` is supported. Default is `"gencode"`.
#' @param organism A character string specifying the organism.
#'   Valid options are `"human"` or `"mouse"`.
#' @param file_type A character string specifying the type of file to download.
#'   Valid options are `"gtf"`, `"gff"`, or `"fasta"`. Defaults to `"gtf"`.
#'   **Note:** `"fasta"` refers to the transcriptome FASTA file.
#' @param output_path A character string specifying the directory where the
#'   downloaded file will be saved. Defaults to `"data-raw"`.
#' @param timeout_limit A numeric value specifying the maximum time in seconds for the
#'   download to complete. This argument takes precedence over `options("timeout")`.
#'   Defaults to `3600` seconds (1 hour).
#' @param method A character string specifying the method used by
#'   `utils::download.file()`. Defaults to `"auto"`.
#'
#' @return A character string with the full path to the downloaded file.
#'
#' @details
#' The function constructs the appropriate download URL based on the specified organism,
#' version, and file type, and downloads the file to the specified output path.
#' If the file already exists in the output directory, the function will not download it again
#' and will return the existing file path. The function requires an internet connection and
#' handles timeout settings to prevent download interruptions.
#'
#' @note Currently, only `"gencode"` reference files are supported.
#'   The `"mane"` reference is not implemented yet.
#'
#' @examples
#' # Download human GTF file for GENCODE release 43
#' gtf_file <- download_reference(
#'   version = "43",
#'   organism = "human",
#'   file_type = "gtf",
#'   output_path = "data-raw"
#' )
#'
#' # Download mouse GTF file for GENCODE release M32
#' gtf_file_mouse <- download_reference(
#'   version = "M32",
#'   organism = "mouse",
#'   file_type = "gtf",
#'   output_path = "data-raw"
#' )
#'
#' # Download human transcriptome FASTA file for GENCODE release 43
#' fasta_file <- download_reference(
#'   version = "43",
#'   organism = "human",
#'   file_type = "fasta",
#'   output_path = "data-raw"
#' )
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
  invisible(full_output_path)
}
