# IsoformicExperiment-class

#' IsoformicExperiment Class
#'
#' The `IsoformicExperiment` class encapsulates the core data structure
#' for transcriptomic analyses in the `isoformic` package.
#' It holds the path to the dataset, sample metadata, and provides
#' access to transcript, gene, and exon annotations through properties.
#' The preferred way to construct an object of this class is through the
#' [`isoformic_experiment()`] function.
#'
#' @param path Character string specifying the path to the data directory.
#' @param annot_path Character string specifying the path to the annotation
#' database directory.
#' @param col_data A data frame containing sample metadata.
#' @param row_data_transcripts A property that retrieves transcript
#' annotation data.
#' @param row_data_genes A property that retrieves gene annotation data.
#' @param row_data_exons A property that retrieves exon annotation data.
#' @param row_data A property that aggregates transcript, gene,
#' and exon annotation data.
#'
#' @keywords internal
#' @export
IsoformicExperiment <- S7::new_class(
  "IsoformicExperiment",
  properties = list(
    experiment_name = S7::class_character,
    data_path = S7::class_character,
    annot_path = S7::class_character,
    col_data = S7::class_any,
    row_data_transcripts = S7::new_property(
      default = NULL,
      getter = function(self) {
        get_row_data_type(self, "transcripts")
      }
    ),
    row_data_genes = S7::new_property(
      default = NULL,
      getter = function(self) {
        get_row_data_type(self, "genes")
      }
    ),
    row_data_exons = S7::new_property(
      default = NULL,
      getter = function(self) {
        get_row_data_type(self, "exons")
      }
    ),
    row_data = S7::new_property(
      default = NULL,
      getter = function(self) {
        get_row_data(self)
      }
    ),
    annot_metadata = S7::class_any
  )
)

#' Read Sample Metadata
col_data <- S7::new_generic("col_data", "self")

S7::method(col_data, IsoformicExperiment) <- function(self) {
  return(self@col_data)
}

#' Read Transcript Annotation
row_data_transcripts <- S7::new_generic("row_data_transcripts", "self")

S7::method(row_data_transcripts, IsoformicExperiment) <- function(self) {
  return(self@row_data_transcripts)
}

#' Read Gene Annotation
row_data_genes <- S7::new_generic("row_data_genes", "self")

S7::method(row_data_genes, IsoformicExperiment) <- function(self) {
  return(self@row_data_genes)
}

#' Read Exon Annotation
row_data_exons <- S7::new_generic("row_data_exons", "self")

S7::method(row_data_exons, IsoformicExperiment) <- function(self) {
  return(self@row_data_exons)
}

#' Aggregate Whole Annotation Table
row_data <- S7::new_generic("row_data", "self")

S7::method(row_data, IsoformicExperiment) <- function(self, compute = FALSE) {
  return(get_row_data(self, compute))
}

S7::method(print, IsoformicExperiment) <- function(x, ...) {
  cat("<IsoformicExperiment>\n")
  cat(" Experiment Name: ", x@experiment_name, "\n", sep = "")
  cat(" Data Path: ", x@data_path, "\n", sep = "")
  cat(" Annotation Name: ", fs::path_file(x@annot_path), "\n", sep = "")
  cat(" Samples: ", nrow(x@col_data), "\n", sep = "")
  cat(" Transcripts: ", nrow(x@row_data_transcripts), "\n", sep = "")
  cat(" Genes: ", nrow(x@row_data_genes), "\n", sep = "")
  cat(" Exons: ", nrow(x@row_data_exons), "\n", sep = "")
  cat(" Annotation Metadata: ", length(x@annot_metadata), " slots\n", sep = "")
}

# =========================================================================

# Utils
#' Get Row Data Annotation by Type
get_row_data_type <- function(self, type) {
  if (!(type %in% c("transcripts", "genes", "exons"))) {
    cli::cli_abort(
      message = c(
        x = "{.var type} must be one of {.val transcripts}, {.val genes}, or {.val exons}."
      ),
      class = "isoformic_annot_invalid_type"
    )
  }
  if (!isTRUE(fs::dir_exists(self@annot_path))) {
    cli::cli_abort(
      message = c(
        x = "The annotation path {.path {self@annot_path}} does {.strong not} exist."
      ),
      class = "isoformic_annot_path_dont_exist"
    )
  }
  arrow::read_parquet(
    file = fs::path(
      self@annot_path, paste0("row_data_", type), ext = "parquet"
    ),
    as_data_frame = FALSE
  )
}

get_row_data <- function(self, compute = FALSE) {
  query_res <- self@row_data_transcripts |>
    dplyr::left_join(
      self@row_data_genes,
      by = "gene_id",
      suffix = c("_tx", "_gene")
    ) |>
    dplyr::left_join(
      self@row_data_exons,
      by = "transcript_id",
      suffix = c("", "_exon")
    ) |>
    dplyr::select(
      -c(dplyr::starts_with("id"))
    )
  if (isTRUE(compute)) {
    return(dplyr::collect(query_res))
  } else {
    return(query_res)
  }
}

# read_rds | load

# write_rds | save

#' Instantiate an IsoformicExperiment Object
#'
#' This function creates an instance of the IsoformicExperiment class,
#' which is used to manage and analyze transcriptomic data.
#'
#' It requires the path to a dataset that contains the necessary
#' data files, including sample metadata and annotation tables.
#'
#' @param db_path Character string specifying the path to the dataset.
#'
#' @rdname isoformic_experiment
#'
#' @export
isoformic_experiment <- function(db_path) {
  IsoformicExperiment(
    path = db_path
  )
}
