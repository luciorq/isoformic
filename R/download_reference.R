#' Download Reference Files from GENCODE
#'
#' Downloads reference annotation files from the GENCODE database
#' for human or mouse genomes.
#' Supports downloading GFF, GTF, transcriptome FASTA,
#' and genome FASTA files.
#' The function handles directory creation and checks for existing files
#' to avoid redundant downloads.
#'
#' @param version A character string specifying the GENCODE release version.
#' For mouse references, include the letter 'M' in the version string
#' (e.g., `"M38"`). Default is `"49"`.
#' @param reference A character string specifying the source of the reference
#' file.
#' Currently, only `"gencode"` is supported. Default is `"gencode"`.
#' @param organism A character string specifying the organism.
#' Valid options are `"human"` or `"mouse"`.
#' @param file_type A character string specifying the type of file to download.
#' Valid options are `"gff"`, `"gtf"`, `"fasta"` or `"genome_fasta"`.
#' Defaults to `"gff"`.
#' **Note:** `"fasta"` refers to the transcriptome FASTA file.
#' `"genome_fasta"` refers to the whole genome sequence FASTA file.
#' @param output_path A character string specifying the directory where the
#' downloaded file will be saved. Defaults to `":cache:"`.
#' Cache path is defined by the `[isoformic::get_isoformic_cache()]` function.
#' @param timeout_limit A numeric value specifying the maximum time in seconds
#' for the download to complete. This argument takes precedence over
#' `options("timeout")`. Defaults to `3600` seconds (1 hour).
#' @param method A character string specifying the method used by
#' `utils::download.file()`. Defaults to `"auto"`.
#'
#' @return A character string with the full path to the downloaded file.
#'
#' @details
#' The function constructs the appropriate download URL based on the
#' specified organism, version, and file type, and downloads the file to
#' the specified output path, being a user cache by default.
#' If the file already exists in the output directory, the function will not
#' download it again and will return the existing file path. The function
#' requires an internet connection and handles timeout settings to prevent
#' download interruptions.
#'
#' @note Currently, only `"gencode"` reference files are supported.
#'   The `"mane"` reference is not implemented yet.
#'
#' @examples
#' \dontrun{
#' # Download human GFF file for GENCODE release 49
#' gff_file <- download_reference(
#'   version = "49",
#'   organism = "human",
#'   file_type = "gff",
#'   output_path = ":cache:"
#' )
#'
#' # Download mouse GFF file for GENCODE release M38
#' gff_file_mouse <- download_reference(
#'   version = "M38",
#'   organism = "mouse",
#'   file_type = "gff",
#'   output_path = ":cache:"
#' )
#'
#' # Download human transcriptome FASTA file for GENCODE release 49
#' fasta_file <- download_reference(
#'   version = "49",
#'   organism = "human",
#'   file_type = "fasta",
#'   output_path = ":cache:"
#' )
#' }
#'
#' @export
download_reference <- function(
  version = "49",
  reference = "gencode",
  organism = c("human", "mouse"),
  file_type = c("gff", "gtf", "fasta", "genome_fasta"),
  output_path = ":cache:",
  timeout_limit = 3600,
  method = "auto"
) {
  if (requireNamespace("curl", quietly = TRUE)) {
    if (!isTRUE(curl::has_internet())) {
      cli::cli_abort(
        message = c(
          "x" = "No internet connection available."
        ),
        class = "isoformic_download_reference_no_internet_error"
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
  if (identical(output_path, ":cache:")) {
    output_path <- get_isoformic_cache()
  }

  base_url <- switch(
    organism,
    human = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_",
    mouse = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_"
  )

  # Example URLs for reference:
  # https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.chr_patch_hapl_scaff.annotation.gff3.gz
  # https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/gencode.v34.chr_patch_hapl_scaff.annotation.gff3.gz
  # For the Genome sequence FASTA file
  # "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/GRCh38.p14.genome.fa.gz"
  # "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M38/GRCm39.genome.fa.gz"
  if (identical(file_type, "genome_fasta")) {
    file_string <- switch(
      organism,
      human = {
        version_num <- as.numeric(version)
        if (version_num >= 44) {
          "GRCh38.p14.genome.fa.gz"
        } else if (version_num >= 32) {
          "GRCh38.p13.genome.fa.gz"
        } else if (version_num >= 28) {
          "GRCh38.p12.genome.fa.gz"
        } else if (version_num >= 26) {
          "GRCh38.p10.genome.fa.gz"
        } else if (version_num == 25) {
          "GRCh38.p7.genome.fa.gz"
        } else if (version_num == 24) {
          "GRCh38.p5.genome.fa.gz"
        } else if (version_num == 23) {
          "GRCh38.p3.genome.fa.gz"
        } else if (version_num == 22) {
          "GRCh38.p2.genome.fa.gz"
        } else if (version_num %in% c(20, 21)) {
          "GRCh38.genome.fa.gz"
        } else if (version_num == 19) {
          "GRCh37.p13.genome.fa.gz"
        } else if (version_num == 18) {
          "GRCh37.p12.genome.fa.gz"
        } else if (version_num == 17) {
          "GRCh37.p11.genome.fa.gz"
        } else {
          cli::cli_abort(
            message = c(
              "x" = "GENCODE version {version} not supported for human genome FASTA download.",
              "i" = "Supported versions are 17 and above."
            ),
            class = "isoformic_download_reference_unsupported_version_error"
          )
        }
      },
      mouse = {
        version_num <- stringr::str_remove(version, "^M")
        version_num <- as.numeric(version_num)
        if (version_num >= 26) {
          "GRCm39.genome.fa.gz"
        } else if (version_num >= 17) {
          "GRCm38.p6.genome.fa.gz"
        } else if (version_num >= 12) {
          "GRCm38.p5.genome.fa.gz"
        } else if (version_num >= 6) {
          "GRCm38.p4.genome.fa.gz"
        } else if (version_num >= 3) {
          "GRCm38.p3.genome.fa.gz"
        } else if (version_num == 2) {
          "GRCm38.p2.genome.fa.gz"
        } else if (version_num == 1) {
          "NCBIM37.genome.fa.gz"
        } else {
          cli::cli_abort(
            message = c(
              "x" = "GENCODE version {version} not supported for mouse genome FASTA download.",
              "i" = "Supported versions are 27 and above."
            ),
            class = "isoformic_download_reference_unsupported_version_error"
          )
        }
      }
    )
    download_url <- paste0(base_url, version, "/", file_string)
    full_output_path <- paste0(output_path, "/", file_string)
  } else {
    file_string <- switch(
      file_type,
      # CHR - Chr only annot
      # gff = ".annotation.gff3.gz",
      # gtf = ".annotation.gtf.gz",
      # ALL - Full annot
      gff = ".chr_patch_hapl_scaff.annotation.gff3.gz",
      gtf = ".chr_patch_hapl_scaff.annotation.gtf.gz",
      fasta = ".transcripts.fa.gz"
    )
    download_url <- paste0(
      base_url,
      version,
      "/gencode.v",
      version,
      file_string
    )
    full_output_path <- paste0(
      output_path,
      "/gencode.v",
      version,
      file_string
    )
  }

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
