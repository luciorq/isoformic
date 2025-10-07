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
#' @param path Character string specifying the path to the dataset.
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
    path = S7::class_character,
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
    )
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

S7::method(row_data, IsoformicExperiment) <- function(self) {
  return(self@row_data)
}

# Utils
#' Get Row Data Annotation by Type
get_row_data_type <- function(self, type) {
  arrow::open_dataset(
    sources = fs::path(self@path, paste0("row_data_", type), ext = "parquet"),
    format = "parquet"
  )
}

get_row_data <- function(self) {
  self@row_data_transcripts |>
    dplyr::left_join(
      self@row_data_genes,
      by = "gene_id",
      suffix = c("_tx", "_gene")
    ) |>
    dplyr::left_join(
      self@row_data_exons,
      by = "transcript_id",
      suffix = c("", "_exon")
    )
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
