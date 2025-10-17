# IsoformicExperiment-class

#' IsoformicExperiment Class
#'
#' The `IsoformicExperiment` class encapsulates the core data structure
#' for transcriptomic analyses in the `isoformic` package.
#' It holds the path to the dataset, sample metadata, and provides
#' access to transcript, gene, and exon annotations through properties.
#' The preferred way to construct an object of this class is through the
#' [`IsoformicExperiment()`] constructor.
#'
#' @param experiment_name Character string specifying the name of the experiment.
#' This name is used for caching the assays experiment.
#' If a name is not provided a random identifier is used.
#' @param data_path Character string specifying the path to the data directory.
#' @param annot_path Character string specifying the path to the annotation
#' database directory.
#' @param assay A list of matrices or data frames containing assay data,
#' with transcript IDs as row names and sample IDs as column names.
#' Each element of the list represents a different assay (e.g., TPM, counts).
#' @param annot_metadata A list containing metadata about the annotation,
#' such as source, version, and date.
#' @param dea A list containing differential expression analysis results
#' for transcripts and genes.
#' @param gsea A list containing gene set enrichment analysis results.
#' @param tx_type_palette A named character vector specifying the color
#' palette for different transcript types.
#' @param col_data A data frame containing sample metadata.
#' First column must be `sample_id` matching the column names of the assays.
#'
#' @param annot_data_transcripts A property that retrieves transcript
#' annotation data.
#' @param annot_data_genes A property that retrieves gene annotation data.
#' @param annot_data_exons A property that retrieves exon annotation data.
#' @param annot_data A property that aggregates transcript, gene,
#' and exon annotation data.
#'
#' @param self An `IsoformicExperiment` object.
#'
#' @param ... Additional arguments passed to methods.
#'
#' @rdname IsoformicExperiment
#'
#' @export
IsoformicExperiment <- S7::new_class(
  "IsoformicExperiment",
  properties = list(
    experiment_name = S7::new_property(
      class = S7::class_character,
      default = NA_character_
    ),
    data_path = S7::new_property(
      class = S7::class_any,
      default = NULL
    ),
    annot_path = S7::new_property(
      class = S7::class_any,
      default = NULL
    ),
    assay = S7::new_property(
      class = S7::class_any,
      default = NULL
    ),
    col_data = S7::new_property(
      class = S7::class_any,
      default = NULL
    ),
    annot_data_transcripts = S7::new_property(
      class = S7::class_any,
      default = NULL,
      getter = function(self) {
        get_annot_data_type(self, "transcripts")
      }
    ),
    annot_data_genes = S7::new_property(
      class = S7::class_any,
      default = NULL,
      getter = function(self) {
        get_annot_data_type(self, "genes")
      }
    ),
    annot_data_exons = S7::new_property(
      class = S7::class_any,
      default = NULL,
      getter = function(self) {
        get_annot_data_type(self, "exons")
      }
    ),
    annot_data = S7::new_property(
      class = S7::class_any,
      default = NULL,
      getter = function(self) {
        get_annot_data(self)
      }
    ),
    annot_metadata = S7::class_any,
    dea = S7::class_any,
    gsea = S7::class_any,
    tx_type_palette = S7::new_property(
      class = S7::class_character,
      default = NULL,
      getter = function(self) {
        if (isTRUE(rlang::is_null(self@tx_type_palette))) {
          return(tx_type_palette())
        } else {
          return(self@tx_type_palette)
        }
      },
      setter = function(self, value) {
        if (!length(value)) {
          return(self)
        }

        message(names(value))
        if (
          !is.character(value) ||
            rlang::is_null(names(value))
        ) {
          cli::cli_abort(
            message = c(
              x = "{.var tx_type_palette} must be a named character vector."
            ),
            class = "isoformic_invalid_tx_type_palette"
          )
        }
        self@tx_type_palette <- value
        return(self)
      }
    )
  ),
  validator = function(self) {
    if (isTRUE(length(self@experiment_name) > 1)) {
      cli::cli_abort(
        message = c(
          x = "{.var experiment_name} must be a single character string."
        ),
        class = "isoformic_invalid_experiment_name"
      )
    }

    if (isFALSE(rlang::is_null(self@data_path))) {
      if (isFALSE(is.character(self@data_path)) || isTRUE(length(self@data_path) != 1)) {
        cli::cli_abort(
          message = c(
            x = "{.var data_path} must be a single character string."
          ),
          class = "isoformic_invalid_data_path"
        )
      }
      if (isFALSE(rlang::is_null(self@data_path))) {
        if (isFALSE(fs::dir_exists(self@data_path))) {
          cli::cli_abort(
            message = c(
              x = "The data path {.path {self@data_path}} does {.strong not} exist."
            ),
            class = "isoformic_data_path_dont_exist"
          )
        }
      }
    }

    if (isFALSE(rlang::is_null(self@annot_path))) {
      if (isFALSE(fs::dir_exists(self@annot_path))) {
        cli::cli_abort(
          message = c(
            x = "The annotation path {.path {self@annot_path}} does {.strong not} exist."
          ),
          class = "isoformic_annot_path_dont_exist"
        )
      }
    }


    if (isTRUE(length(self@assay) > 0)) {
      for (assay_name in names(self@assay)) {
        validate_assay_rownames(self, assay_name)
        validate_assay_colnames(self, assay_name)
      }
    }
    if (isFALSE(rlang::is_null(self@dea)) && isTRUE(length(self@dea) > 0)) {
      validate_dea(self)
    }
    return(NULL)
  }
)

# =============================================================================
# IsoformicExperiment Accessors - New S7 Generics with S7 Methods
# + Each one of those should correspond to one property in the class definition
# =============================================================================

#' Read Sample Metadata
#' @rdname IsoformicExperiment
col_data <- S7::new_generic("col_data", "self")

S7::method(col_data, IsoformicExperiment) <- function(self) {
  return(self@col_data)
}

#' Read Transcript Annotation
#' @rdname IsoformicExperiment
annot_data_transcripts <- S7::new_generic("annot_data_transcripts", "self")

S7::method(annot_data_transcripts, IsoformicExperiment) <- function(self) {
  return(self@annot_data_transcripts)
}

#' Read Gene Annotation
#' @rdname IsoformicExperiment
annot_data_genes <- S7::new_generic("annot_data_genes", "self")

S7::method(annot_data_genes, IsoformicExperiment) <- function(self) {
  return(self@annot_data_genes)
}

#' Read Exon Annotation
#' @rdname IsoformicExperiment
annot_data_exons <- S7::new_generic("annot_data_exons", "self")

S7::method(annot_data_exons, IsoformicExperiment) <- function(self) {
  return(self@annot_data_exons)
}

#' Aggregate Whole Annotation Table
#' @rdname IsoformicExperiment
annot_data <- S7::new_generic("annot_data", "self")

S7::method(annot_data, IsoformicExperiment) <- function(self, compute = FALSE) {
  return(get_annot_data(self, compute))
}

#' Return Row Names for Transcript Level Annotation
#' @rdname IsoformicExperiment
annot_row_names <- S7::new_generic("annot_row_names", "self")

S7::method(annot_row_names, IsoformicExperiment) <- function(self) {
  self@annot_data_transcripts |>
    dplyr::pull("transcript_id", as_vector = TRUE)
}

#' Return Column Names for Transcript Level Assays
#' @rdname IsoformicExperiment
col_names <- S7::new_generic("col_names", "self")

S7::method(col_names, IsoformicExperiment) <- function(self) {
  self@col_data |>
    dplyr::pull("sample_id")
}

#' Return Row Names for Transcript Level Assays
#' @rdname IsoformicExperiment
row_names <- S7::new_generic("row_names", "self")
S7::method(row_names, IsoformicExperiment) <- function(self) {
  if (length(self@assay) == 0L) {
    return(0L)
  }
  assay_name <- names(self@assay)[1]
  attr(self@assay[[assay_name]], "dimnames")[[1]]
}
# =========================================================================
# S7 Methods for S4 Generics - Any really needed?
# + If S4 generics are needed, add `S7::S4_register(IsoformicExperiment)`
# + to the `.onLoad` function in `R/isoformic-package.R` file.
# =========================================================================

# #' Return Row Names for Transcript Level Assays
# #' @rdname IsoformicExperiment
# S7::new_generic("rownames", "x")

# S7::method(rownames, IsoformicExperiment) <- function(x, do.NULL, prefix) {
# rownames.IsoformicExperiment <- function(x, do.NULL, prefix) {
#  rlang::check_required(x)
#  if (length(x@assay) == 0L) {
#    return(0L)
#  }
#  assay_name <- names(x@assay)[1]
#  return(attr(x@assay[[assay_name]], "dimnames")[[1]])
# }

# =============================================================================
# S7 Methods for common S3 Generics
# =============================================================================

# S3 Methods to implement:
# load - read_rds
# + Check if unique name already exists in cache/path, and warn / ask about
#  + overwriting.
# + potential arguments data_path, experiment_name, overwrite, etc.
# save - write_rds
# + Compress cache files and serialize properly.

# Print IsoformicExperiment Object Summary
S7::method(print, IsoformicExperiment) <- function(x, ...) {
  # Main Print
  cat("<IsoformicExperiment>\n")
  cat(" Experiment Name: ", x@experiment_name, "\n", sep = "")
  cat(" Data Path: ", x@data_path, "\n", sep = "")

  # Assay Block
  cat(" Assays: ", paste(names(x@assay), collapse = ", "), "\n", sep = "")
  cat(" Samples: ", nrow(x@col_data), "\n", sep = "")
  cat(" Transcripts in Assay: ", length(rownames(x)), "\n", sep = "")

  # Annotation Block
  cat(" Annotation Name: ", fs::path_file(x@annot_path), "\n", sep = "")
  cat(" Annotation Path: ", x@annot_path, "\n", sep = "")
  cat(
    " Transcripts in Annotation: ",
    nrow(x@annot_data_transcripts),
    "\n",
    sep = ""
  )
  cat(" Genes in Annotation: ", nrow(x@annot_data_genes), "\n", sep = "")
  cat(" Exons in Annotation: ", nrow(x@annot_data_exons), "\n", sep = "")
  cat(" Annotation Metadata: ", length(x@annot_metadata), " slots\n", sep = "")
}

S7::method(dimnames, IsoformicExperiment) <- function(x) {
  if (length(x@assay) == 0L) {
    return(list(NULL, NULL))
  }
  assay_name <- names(x@assay)[1]
  # message("Using assay: ", assay_name)
  return(dimnames(x@assay[[assay_name]]))
}

S7::method(dim, IsoformicExperiment) <- function(x) {
  if (length(x@assay) == 0L) {
    return(c(0L, 0L))
  }
  assay_name <- names(x@assay)[1]
  # message("Using assay: ", assay_name)
  return(dim(x@assay[[assay_name]]))
}

# =============================================================================
# New S7 Generics and S7 Methods - Convenience Functions to Access Common Data
# =============================================================================

#' Retrieve Transcript to Gene Mapping Table from Annotation
#' @rdname IsoformicExperiment
tx_to_gene <- S7::new_generic("tx_to_gene", "self")

S7::method(tx_to_gene, IsoformicExperiment) <- function(self) {
  self@annot_data_transcripts |>
    dplyr::select("transcript_id", "gene_id") |>
    dplyr::distinct() |>
    dplyr::collect()
}

#' Retrieve Transcript Annotation Table
#' @rdname IsoformicExperiment
tx_annot <- S7::new_generic("tx_annot", "self")

S7::method(tx_annot, IsoformicExperiment) <- function(self) {
  cols_to_remove <- c(
    "seqid",
    "start_pos",
    "end_pos",
    "strand"
  )
  self@annot_data_transcripts |>
    dplyr::select(-dplyr::any_of(cols_to_remove)) |>
    # dplyr::select("transcript_id", "gene_id") |>
    dplyr::distinct() |>
    dplyr::left_join(
      y = self@annot_data_genes |>
        dplyr::select(
          dplyr::any_of(c("gene_id", "gene_name", "gene_type"))
        ) |>
        dplyr::distinct(),
      by = "gene_id"
    ) |>
    dplyr::distinct() |>
    dplyr::collect()
}

#' Retrieve Differential Expression Results for Transcripts
#' @rdname IsoformicExperiment
de_tx <- S7::new_generic("de_tx", "self")

S7::method(de_tx, IsoformicExperiment) <- function(self, de_type = "det") {
  get_dea_results(self, de_type = de_type)
}

#' Retrieve Differential Expression Results for Genes
#' @rdname IsoformicExperiment
de_gene <- S7::new_generic("de_gene", "self")

S7::method(de_gene, IsoformicExperiment) <- function(self, de_type = "deg") {
  get_dea_results(self, de_type = de_type)
}


get_dea_results <- function(self, de_type = c("det", "deg")) {
  .data <- rlang::.data
  .env <- rlang::.env
  `:=` <- rlang::`:=`
  de_type <- rlang::arg_match(de_type)
  if (isTRUE(rlang::is_null(self@dea[[de_type]]))) {
    cli::cli_abort(
      message = c(
        x = "Differential expression results not found in {.cls IsoformicExperiment} object.",
        i = "Please run differential expression analysis and combine results before accessing."
      ),
      class = "isoformic_missing_dea_results"
    )
  }
  self@dea[["sig_cutoff"]] <- list(
    det = list(log2FC = 1, pvalue = 0.05, fdr = "qvalue"),
    deg = list(log2FC = 1, pvalue = 0.05, fdr = "qvalue")
  )

  log2fc_cutoff <- self@dea$sig_cutoff[[de_type]]$log2FC
  pvalue_cutoff <- self@dea$sig_cutoff[[de_type]]$pvalue

  pvalue_col_to_use <- self@dea$sig_cutoff[[de_type]]$fdr
  if (
    rlang::is_null(pvalue_col_to_use) ||
      identical(pvalue_col_to_use, "") ||
      identical(pvalue_col_to_use, "none")
  ) {
    pvalue_col_to_use <- "pvalue"
  }

  self@dea[[de_type]] |>
    dplyr::mutate(
      "{de_type}_sig" := dplyr::case_when(
        abs(.data$log2FC) >= .env$log2fc_cutoff &
          .data[[pvalue_col_to_use]] <= .env$pvalue_cutoff ~
          "yes",
        TRUE ~ "no"
      )
    ) |>
    dplyr::arrange(dplyr::desc(abs(.data$log2FC))) |>
    dplyr::collect()
}


# ==============================================================================
# Utility Getter and Setter Functions for IsoformicExperiment Class
# ==============================================================================

#' Retrieve Annotation Data by Feature Type
#' @keywords internal
#' @noRd
get_annot_data_type <- function(iso_obj, type) {
  if (!(type %in% c("transcripts", "genes", "exons"))) {
    cli::cli_abort(
      message = c(
        x = "{.var type} must be one of {.val transcripts}, {.val genes}, or {.val exons}."
      ),
      class = "isoformic_annot_invalid_type"
    )
  }
  if (!isTRUE(fs::dir_exists(iso_obj@annot_path))) {
    cli::cli_abort(
      message = c(
        x = "The annotation path {.path {iso_obj@annot_path}} does {.strong not} exist."
      ),
      class = "isoformic_annot_path_dont_exist"
    )
  }

  arrow::read_parquet(
    file = fs::path(
      iso_obj@annot_path,
      paste0("annot_data_", type),
      ext = "parquet"
    ),
    as_data_frame = FALSE
  )
}

#' Retrieve Comprehensive Annotation Data for Transcript Level Annotation
#' @keywords internal
#' @noRd
get_annot_data <- function(self, compute = FALSE) {
  query_res <- self@annot_data_transcripts |>
    dplyr::left_join(
      self@annot_data_genes,
      by = "gene_id",
      suffix = c("_tx", "_gene")
    ) |>
    dplyr::left_join(
      self@annot_data_exons,
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

# TODO: @luciorq Check list of scientist names from moby for random names.
#' Set a Random Experiment Name
#'
#' This function generates a random experiment name for the IsoformicExperiment
#' class if `@experiment_name` is not already set.
#'
#' @param data_path Character string specifying the data path where the
#'   experiment data is stored.
#' @param experiment_name Character string specifying the name of the experiment.
#'
#' @keywords internal
#' @noRd
set_random_experiment_name <- function(
    data_path = get_isoformic_cache(),
    experiment_name = NULL) {
  rlang::try_fetch(
    expr = {
      if (isFALSE(fs::dir_exists(data_path))) {
        fs::dir_create(data_path, recurse = TRUE)
      }
    },
    error = function(e) {
      cli::cli_abort(
        message = c(
          `x` = "Could not create or access the data path {.path {data_path}}.",
          `!` = "Please check your file permissions and disk space.",
          `i` = c(
            "If the problem persists, consider using a {.field @data_path} directory on a different filesystem."
          )
        ),
        class = "isoformic_data_path_access_error"
      )
    }
  )
  if (
    isFALSE(rlang::is_null(experiment_name)) &&
      isFALSE(identical(experiment_name, ""))
  ) {
    return(experiment_name)
  }

  repeat {
    random_name <- rlang::hash(base::as.character(base::Sys.time()))
    experiment_dir <- fs::path(data_path, random_name)
    if (!isTRUE(fs::dir_exists(experiment_dir))) {
      fs::dir_create(experiment_dir, recurse = TRUE)
      return(random_name)
    }
  }
}



# ==============================================================================
# Validation Functions for IsoformicExperiment Class
# ==============================================================================
validate_assay_rownames <- function(self, assay_name) {
  annot_tx_rownames <- self@annot_data_transcripts |>
    # dplyr::collect() |>
    dplyr::pull("transcript_id", as_vector = TRUE)
  assay_rownames <- attr(self@assay[[assay_name]], "dimnames")[[1]]

  if (!all(assay_rownames %in% annot_tx_rownames)) {
    cli::cli_abort(
      c(
        "Assay {.val {assay_name}} contains transcript IDs not present in the IsoformicExperiment annotation data.",
        "i" = "Please ensure that all transcript IDs in the assay are present in the used annotation."
      ),
      class = "isoformic_invalid_assay_rownames"
    )
  }

  return(invisible(TRUE))
}

validate_assay_colnames <- function(self, assay_name) {
  col_data_ids <- self@col_data$sample_id
  assay_colnames <- attr(self@assay[[assay_name]], "dimnames")[[2]]
  if (!all(assay_colnames %in% col_data_ids)) {
    cli::cli_abort(
      c(
        "Assay {.val {assay_name}} contains sample IDs not present in the IsoformicExperiment colData.",
        `i` = "Please ensure that all sample IDs in the assay are present in the colData."
      ),
      class = "isoformic_invalid_assay_colnames"
    )
  }
  return(invisible(TRUE))
}

validate_dea <- function(self) {
  # if (is.null(self@dea[["method"]])) {
  #   cli::cli_abort(
  #     c(
  #       `x` = "DEA method is not specified in the '@dea' slot."
  #     ),
  #     class = "isoformic_dea_missing_method"
  #   )
  # }
  # if (is.null(self@dea[["det"]]) && is.null(self@dea[["deg"]])) {
  #   cli::cli_abort(
  #     c(
  #       `x` = "DEA results are missing in the '@dea' slot."
  #     ),
  #     class = "isoformic_dea_missing_results"
  #   )
  # }

  if (!rlang::is_null(self@dea[["det"]])) {
    det <- self@dea[["det"]]
    # required_cols <- c("transcript_id", "log2FC", "pvalue", "padj")
    required_cols <- c("transcript_id", "log2FC", "pvalue", "qvalue")
    missing_cols <- setdiff(required_cols, colnames(det))
    if (length(missing_cols) > 0) {
      cli::cli_abort(
        c(
          `x` = paste(
            "The following required columns are missing in the DET results:",
            paste(missing_cols, collapse = ", ")
          )
        ),
        class = "isoformic_dea_invalid_det"
      )
    }
    annot_tx_ids <- self@annot_data_transcripts |>
      dplyr::collect() |>
      dplyr::pull("transcript_id")
    if (!all(det$transcript_id %in% annot_tx_ids)) {
      missing_tx <- det$transcript_id[!det$transcript_id %in% annot_tx_ids]
      cli::cli_abort(
        c(
          `x` = paste(
            "There are missing transcript_ids in the DET results that are not present in the annotation:",
            paste(utils::head(missing_tx, 10), collapse = ", "),
            ifelse(length(missing_tx) > 10, " ...", "")
          )
        ),
        class = "isoformic_dea_invalid_det"
      )
    }
  }
  if (!rlang::is_null(self@dea[["deg"]])) {
    deg <- self@dea[["deg"]]
    # required_cols <- c("gene_id", "log2FC", "pvalue", "padj")
    required_cols <- c("gene_id", "log2FC", "pvalue", "qvalue")
    missing_cols <- setdiff(required_cols, colnames(deg))
    if (length(missing_cols) > 0) {
      cli::cli_abort(
        c(
          `x` = paste(
            "The following required columns are missing in the DEG results:",
            paste(missing_cols, collapse = ", ")
          )
        ),
        class = "isoformic_dea_invalid_deg"
      )
    }
    annot_gene_ids <- self@annot_data_genes |>
      dplyr::collect() |>
      dplyr::pull("gene_id")
    if (!all(deg$gene_id %in% annot_gene_ids)) {
      missing_genes <- deg$gene_id[!deg$gene_id %in% annot_gene_ids]
      cli::cli_abort(
        c(
          `x` = paste(
            "There are missing gene_ids in the DEG results that are not present in the annotation:",
            paste(utils::head(missing_genes, 10), collapse = ", "),
            ifelse(length(missing_genes) > 10, " ...", "")
          )
        ),
        class = "isoformic_dea_invalid_deg"
      )
    }
  }
  return(invisible(TRUE))
}
