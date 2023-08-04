#' Download Reference Files From GENCODE
#' @param reference Source of the reference file
#'   One of "gencode" or "mane". Defaults to "gencode".
#' @param file_type One of "gff", "fasta", or "gtf". Defaults to "gff".
#'
#' @export
download_reference <- function(version = "43",
                               reference = "gencode",
                               file_type = "gff",
                               output_path = "data-raw",
                               timeout_limit = 3600,
                               method = "auto") {
  base_url <- "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/"
  if (isTRUE(file_type == "gff")) {
    file_string <- ".chr_patch_hapl_scaff.annotation.gff3.gz"
  } else if (isTRUE(file_type == "gtf")) {
    file_string <- ".chr_patch_hapl_scaff.annotation.gtf.gz"
  } else if (isTRUE(file_type == "fasta")) {
    file_string <- ".transcripts.fa.gz"
  }
  gff_url <- paste0(
    base_url, "release_", version, "/gencode.v", version, file_string
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
      dl_res <- utils::download.file(
        url = gff_url,
        destfile = full_output_path,
        method = method
      )
    }
  )
  if (isTRUE(dl_res == 0)) {
    cli::cli_inform(
      c(
        `v` = "{.path {full_output_path}} successfully downloaded."
      )
    )
  }
  return(full_output_path)
}
