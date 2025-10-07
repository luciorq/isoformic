# ContextData-class

#' ContextData Class
#'
#' The `ContextData` class holds information on the genomic context of
#' transcripts and their annotation in a `IsoformicExperiment` analysis.
#'
#' The preferred way to construct an object of this class is through the
#' [`create_context_data()`] function.
#'
#' @param gene_name Character vector of gene names to print.
#' @param gff_file Character string specifying the path to a GFF file
#'  containing the gene annotation.
#' @param txdb A `TxDb` object containing the transcript database.
#' @param annotation_table A data frame containing the transcript annotation.
#' @param annotation_name Character string specifying the name of the
#' annotation (e.g., "Ensembl_v104").
#' @param assembly_name Character string specifying the name of the
#' genome assembly (e.g., "GRCh38").
#' @param ideogram_assembly Character string specifying the name of the
#' genome assembly for ideogram plotting (e.g., "hg38").
#' @param organism Character string specifying the organism name
#' (e.g., "Homo sapiens").
#' @param orgdb_package Character string specifying the name of the
#' organism database package (e.g., "org.Hs.eg.db").
#' @param bsgenome_package Character string specifying the name of the
#' BSgenome package (e.g., "BSgenome.Hsapiens.UCSC.hg38").
#' @param tx_type_palette Character vector specifying the color palette
#' for transcript types.
#'
#' @keywords internal
#' @export
ContextData <- S7::new_class(
  name = "ContextData",
  package = "isoformic",
  properties = list(
    gene_name = S7::class_character,
    gff_file = S7::class_character,
    txdb = S7::class_any,
    annotation_table = S7::class_any,
    annotation_name = S7::class_character,
    assembly_name = S7::class_character,
    ideogram_assembly = S7::class_character,
    organism = S7::class_character,
    orgdb_package = S7::class_character,
    bsgenome_package = S7::class_character,
    tx_type_palette = S7::class_character
  )
)

#' Create ContextData Object
#'
#' This function instantiates a `ContextData` object containing information
#' about the genomic context of transcripts.
#' It requires a GFF file for gene annotation and constructs a `TxDb`
#' object from it. The function also prepares an annotation table and
#' updates the transcript names in the `TxDb` object.
#' The `ContextData` object can then be used in conjunction with
#' `IsoformicExperiment` for transcriptomic analyses.
#'
#' @param gff_file Character string specifying the path to a GFF file
#' containing the gene annotation.
#' @param organism Character string specifying the organism name
#' (e.g., "Homo sapiens").
#' @param orgdb_package Character string specifying the name of the
#' organism database package (e.g., "org.Hs.eg.db").
#' @param bsgenome_package Character string specifying the name of the
#' BSgenome package (e.g., "BSgenome.Hsapiens.UCSC.hg38").
#' @param tx_type_palette Named character vector specifying the color palette
#' for transcript types.
#'
#' @inheritParams rlang::args_dots_empty
#'
#' @rdname create_context_data
#'
#' @export
create_context_data <- function(
    gff_file,
    ...,
    organism,
    orgdb_package,
    bsgenome_package,
    tx_type_palette = NULL) {
  rlang::check_dots_empty()
  .data <- rlang::.data
  context_data <- ContextData()
  context_data@gff_file <- gff_file
  context_data@organism <- organism
  context_data@orgdb_package <- orgdb_package
  context_data@bsgenome_package <- bsgenome_package

  if (rlang::is_null(tx_type_palette)) {
    context_data@tx_type_palette <- tx_type_palette()
  } else {
    context_data@tx_type_palette <- tx_type_palette
  }

  context_data@txdb <- txdbmaker::makeTxDbFromGFF(
    file = context_data@gff_file,
    dataSource = "Custom Isoformic TxDb",
    organism = context_data@organism
  )

  context_data@annotation_table <- prepare_annotation(
    context_data@gff_file,
    file_type = "gff"
  )$transcript |>
    dplyr::select(
      "transcript_id",
      "gene_id",
      "transcript_type",
      "transcript_name"
    ) |>
    dplyr::rename(
      tx_name = "transcript_id",
      gene_name = "gene_id",
      tx_biotype = "transcript_type",
      tx_symbol = "transcript_name"
    )
  tx_replace_df <- DBI::dbGetQuery(
    context_data@txdb$conn, "SELECT * FROM transcript;"
  ) |>
    tibble::as_tibble() |>
    dplyr::left_join(
      context_data@annotation_table,
      by = "tx_name",
      relationship = "many-to-many"
    ) |>
    dplyr::distinct()

  tx_replace_df2 <- tx_replace_df |>
    dplyr::rename(
      `tx_ensemblid` = "tx_name"
    ) |>
    dplyr::mutate(
      `tx_name` = .data$tx_symbol
    )

  DBI::dbWriteTable(
    context_data@txdb$conn,
    "transcript",
    tx_replace_df2,
    overwrite = TRUE
  )
  return(context_data)
}
