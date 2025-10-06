#' Parse Annotation Attributes
#'
#' Parses the attributes of a specific feature type from a column vector.
#'
#' @keywords internal
#' @noRd
parse_annot_attributes <- function(
    column_vector,
    feature_type = c("gene", "transcript", "exon"),
    file_type = c("gtf", "gff")) {
  feature_type <- rlang::arg_match(feature_type)
  file_type <- rlang::arg_match(file_type)

  if (isTRUE(feature_type == "gene")) {
    gene_id_vector <- extract_attribute(
      column_vector = column_vector,
      attr_name = "gene_id",
      file_type = file_type
    )
    gene_name_vector <- extract_attribute(
      column_vector = column_vector,
      attr_name = "gene_name",
      file_type = file_type
    )
    gene_type_vector <- extract_attribute(
      column_vector = column_vector,
      attr_name = "gene_type",
      file_type = file_type
    )
    attr_df <- tibble::tibble(
      gene_id = gene_id_vector,
      gene_name = gene_name_vector,
      gene_type = gene_type_vector
    )
  } else if (isTRUE(feature_type == "exon")) {
    exon_id_vector <- extract_attribute(
      column_vector = column_vector,
      attr_name = "exon_id",
      file_type = file_type
    )
    exon_number_vector <- extract_attribute(
      column_vector = column_vector,
      attr_name = "exon_number",
      file_type = file_type
    )
    transcript_id_vector <- extract_attribute(
      column_vector = column_vector,
      attr_name = "transcript_id",
      file_type = file_type
    )
    attr_df <- tibble::tibble(
      exon_id = exon_id_vector,
      exon_number = exon_number_vector,
      transcript_id = transcript_id_vector
    )
  } else if (isTRUE(feature_type == "transcript")) {
    transcript_id_vector <- extract_attribute(
      column_vector = column_vector,
      attr_name = "transcript_id",
      file_type = file_type
    )
    transcript_name_vector <- extract_attribute(
      column_vector = column_vector,
      attr_name = "transcript_name",
      file_type = file_type
    )
    transcript_type_vector <- extract_attribute(
      column_vector = column_vector,
      attr_name = "transcript_type",
      file_type = file_type
    )
    gene_id_vector <- extract_attribute(
      column_vector = column_vector,
      attr_name = "gene_id",
      file_type = file_type
    )
    attr_df <- tibble::tibble(
      transcript_id = transcript_id_vector,
      transcript_name = transcript_name_vector,
      transcript_type = transcript_type_vector,
      gene_id = gene_id_vector
    )
  }
  return(attr_df)
}
