#' Extract Attribute from Column Vector
#'
#' Extracts a specific attribute from a column vector based on the file type.
#'
#' @keywords internal
#' @noRd
extract_attribute <- function(column_vector, attr_name, file_type) {
  if (isTRUE(file_type == "gtf")) {
    attr_vector <- column_vector |>
      stringr::str_extract(pattern = paste0(attr_name, ".*?;")) |>
      stringr::str_remove(paste0("^", attr_name, " \"")) |>
      stringr::str_remove("\";$") |>
      stringr::str_squish()
  } else if (isTRUE(file_type == "gff")) {
    attr_vector <- column_vector |>
      stringr::str_extract(pattern = paste0(attr_name, "=.*?;")) |>
      stringr::str_remove(paste0("^", attr_name, "=")) |>
      stringr::str_remove(";$") |>
      stringr::str_squish()
  }
  return(attr_vector)
}
