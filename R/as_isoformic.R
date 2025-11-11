#' Convert a SummarizedExperiment to an IsoformicExperiment Object
#'
#' This function converts a `SummarizedExperiment` object to an `IsoformicExperiment` object.
#' It extracts the assays, row data, column data, and metadata from the input object
#' and uses them to create a new `IsoformicExperiment` object.
#'
#' @param se A `SummarizedExperiment` object to be converted.
#' @param annot_path Path to the annotation file. This can be a GFF file or the
#' path pre-built annotation database created with `[prepare_isoformic_annotation()]`.
#' @param annot_type Type of the annotation file provided.
#' Options are "gff" for GFF files and "annot_db" for pre-built annotation
#' databases.
#'
#' @export
as_isoformic <- function(se, annot_path, annot_type = c("gff", "annot_db")) {
  # rlang::check_installed("SummarizedExperiment")
  rlang::check_required(se)
  rlang::check_required(annot_path)
  annot_type <- rlang::arg_match(annot_type)

  if (S7::S7_inherits(x = se, class = isoformic::IsoformicExperiment)) {
    return(se)
  }

  if (isFALSE(inherits(x = se, what = "SummarizedExperiment"))) {
    cli::cli_abort(
      c(
        "The input object must be a SummarizedExperiment.",
        "x" = "An object of class {.cls {class(se)}} was provided."
      ),
      class = "isoformic_convert_se_error"
    )
  }

  if (isFALSE(fs::file_exists(annot_path))) {
    cli::cli_abort(
      c(
        "Annotation file not found.",
        "x" = "The file {.file {annot_path}} does not exist."
      ),
      class = "isoformic_convert_se_error"
    )
  }

  if (identical(annot_type, "gff")) {
    annot_db_path <- prepare_isoformic_annotation(
      input_path = annot_path,
      output_path = get_isoformic_cache(),
      file_type = "gff"
    )
    annot_metadata <- get_annot_metadata(gff_file = annot_path)
  } else if (identical(annot_type, "annot_db")) {
    annot_db_path <- annot_path
    annot_metadata <- NULL
  }

  experiment_name <- rlang::hash(se)
  data_path <- get_isoformic_cache(experiment_name)
  if (!fs::dir_exists(data_path)) {
    fs::dir_create(data_path, recurse = TRUE)
  }
  col_data <- se@colData |>
    base::as.data.frame() |>
    tibble::rownames_to_column("sample_id") |>
    dplyr::mutate(dplyr::across(dplyr::where(is.factor), as.character)) |>
    tibble::as_tibble()

  iso_obj <- IsoformicExperiment(
    experiment_name = experiment_name,
    data_path = data_path,
    annot_path = annot_db_path,
    col_data = col_data
  )

  is_gencode_tx_id <- grepl("^ENST", se@NAMES[1])
  is_gencode_pipe_names <- grepl("\\|", se@NAMES[1])
  is_gencode_pipe_assays <- grepl(
    "\\|",
    attr(se@assays@data@listData[[1]], "dimnames")[[1]][1]
  )

  if (
    isTRUE(is_gencode_tx_id) &&
      isTRUE((is_gencode_pipe_names || is_gencode_pipe_assays))
  ) {
    se@NAMES <- stringr::str_remove(se@NAMES, "\\|.*")
    for (i in names(se@assays@data@listData)) {
      attr(se@assays@data@listData[[i]], "dimnames")[[1]] <- se@NAMES
    }
  }

  # TODO: @luciorq - Verify if this approach still works for
  # + RangedSummarizedExperiment and DelayedArray based objects assays
  # if (inherits(se@assays@data@listData, "list")) {
  assay_list <- se@assays@data@listData
  # } else {
  # assay_list <- list(counts = se@assays@data@listData)
  # }

  iso_obj@assay <- assay_list
  iso_obj@annot_metadata <- annot_metadata

  return(iso_obj)
}
