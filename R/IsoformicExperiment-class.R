# IsoformicExperiment class ---------------------------------------------------

#' IsoformicExperiment: Transcript, Gene and Exon Level Data in One Object
#'
#' @description
#' `IsoformicExperiment` is an S7 class modelled after
#' `SummarizedExperiment`, but organised around *transcriptomic diversity*:
#' the same samples are described at three levels that are kept consistent
#' by construction.
#'
#' * **Transcript level**: assays are numeric matrices with transcript IDs
#'   as row names and sample IDs as column names (for example `abundance`,
#'   `counts`, `length`). All assays share the same dimnames.
#' * **Gene level**: gene assays are derived from the transcript assays by
#'   summing the isoforms of each gene (abundance-weighted means for
#'   `length`), cached, and can be overridden with [gene_assays<-].
#' * **Exon level**: exon structures come from the annotation and are
#'   available through [annot_data_exons()], [exon_table()] and
#'   [exon_diversity()].
#'
#' Annotation is never copied into the object. The `annot_path` property
#' points to an *annotation directory* of Parquet files created by
#' [prepare_isoformic_annotation()] (a GFF3/GTF file is accepted and
#' converted on the fly). Tables are read with DuckDB on first access and
#' cached in memory keyed by file path, size and modification time, so a
#' re-created annotation is picked up automatically and several objects
#' sharing one annotation share one copy of the tables.
#'
#' Derived results (gene assays, row data, isoform fractions, ORFs,
#' transcript features and productivity classes) are cached per object in
#' the `.cache` environment, keyed by a hash of their inputs, and therefore
#' survive copies but never go stale. See [clear_cache()].
#'
#' Objects persist to disk as a directory of Parquet files with
#' [write_isoformic()] and are restored with [read_isoformic()].
#'
#' @section Properties:
#' \describe{
#'   \item{`experiment_name`}{Single string. Generated when not supplied.}
#'   \item{`annot_path`}{Path to the annotation directory, or `NULL` for an
#'     empty object.}
#'   \item{`data_path`}{Directory the object was written to or read from
#'     (see [write_isoformic()]), or `NULL`.}
#'   \item{`assays`}{Named list of transcript-level matrices.}
#'   \item{`gene_assays`}{Named list of gene-level matrices. Derived from
#'     `assays` when not set explicitly.}
#'   \item{`col_data`}{Tibble of sample metadata; its `sample_id` column
#'     matches the assay column names.}
#'   \item{`dea`}{List with differential expression results: `det`
#'     (transcripts), `deg` (genes) and `cutoffs`, and the uncertainty
#'     tables of [add_uncertainty()]. Use [de_tx<-] and [de_gene<-] to set
#'     the results.}
#'   \item{`gsea`}{Named list of enrichment result tables.}
#'   \item{`metadata`}{Free-form named list.}
#'   \item{`annot_metadata`}{Named list describing the annotation source,
#'     read from the annotation directory unless set explicitly.}
#'   \item{`tx_type_palette`}{Named character vector of colours per
#'     transcript type, used by all plots. Defaults to [tx_type_palette()].}
#'   \item{`annot_data_transcripts`, `annot_data_genes`,
#'     `annot_data_exons`}{Read-only, cached annotation tables.}
#' }
#'
#' @section Subsetting:
#' `x[i, j]` keeps the transcripts `i` and the samples `j`, given as IDs,
#' positions or logical vectors; a missing `i` or `j` keeps all of them, and
#' the subset must hold at least one transcript and one sample. Everything
#' stored by transcript or by sample follows: the assays, the sample table,
#' stored gene assays (the genes of the kept transcripts) and the
#' inferential replicates in `x@metadata$inf_reps`. The uncertainty tables
#' of [add_uncertainty()] are restricted to the kept transcripts (a
#' collapsed ambiguity group stays only with all its members); when samples
#' are dropped they are not recomputed, and a message says so. Differential
#' expression and enrichment results are kept as they are.
#'
#' The inferential replicates are not copied: the `metadata$inf_reps` of a
#' subset holds the draws of `x` and a `selection` attribute, a list with
#' the IDs of the kept `transcripts` and `samples`. The functions of the
#' package (and [inferential_replicates_table()]) read the draws through
#' that selection, and [write_isoformic()] writes only the selected rows and
#' columns, so the results are those of a copy of the kept draws, while the
#' memory of a subset does not grow with the draws. The trade-off: the draws
#' of `x` (every transcript and sample) stay in memory as long as a subset
#' made from it exists, even after `x` itself is removed, and they go with
#' every serialised copy of the subset: [saveRDS()], [serialize()] and the
#' objects sent to parallel workers that do not share memory (such as
#' `future` or `mirai` workers) carry every transcript and sample of the
#' draws of `x`, while [write_isoformic()] writes only the selected ones.
#' Code that indexes the matrices of `metadata$inf_reps` directly sees the
#' draws of `x`. A list operation on `metadata$inf_reps` (such as
#' [lapply()] or `[`) drops the `selection` attribute and returns the draws
#' of `x`: the functions that take the object still read only its kept
#' transcripts and samples, while code that gets only the list, such as
#' [inferential_replicates_table()] without its `x` argument, reads every
#' transcript and sample of `x`. When no kept
#' transcript has draws, the subset keeps only the replicate moments (with
#' zero-row draws), so that [scale_counts()] and
#' [inferential_overdispersion()] still work, but the functions that need
#' the draws treat it as an object without inferential replicates; without
#' moments, `metadata$inf_reps` is `NULL`. To replace the draws by plain matrices of
#' the kept draws (and free those of `x`), write the subset with
#' [write_isoformic()] and read it back; to drop them, set
#' `metadata$inf_reps <- NULL`. To compare two groups of a larger object,
#' pass `levels` to the function instead of subsetting: the samples of the
#' other groups are ignored. A subset with other transcripts
#' starts with an empty cache of derived tables (the transcript annotation
#' join is shared).
#'
#' @param experiment_name Single string naming the experiment; a random name
#'   is generated when `NULL`.
#' @param annot_path Annotation directory or GFF3/GTF file.
#' @param assays Named list of transcript-level matrices (or data frames
#'   with a `transcript_id` column).
#' @param col_data Data frame of sample metadata with a `sample_id` column
#'   (or sample IDs as row names / first column). `NULL` builds it from the
#'   column names of the assays.
#' @param gene_assays Optional named list of gene-level matrices.
#' @param dea Optional list of differential expression results.
#' @param gsea Optional named list of enrichment results.
#' @param metadata Optional named list of free-form metadata.
#' @param data_path Optional on-disk location of the experiment.
#' @param tx_type_palette Optional named colour vector per transcript type.
#'
#' @returns An `IsoformicExperiment` object.
#'
#' @seealso [example_isoformic()] for a ready-made object,
#'   [read_salmon()] and [as_isoformic()] to create one from quantifications,
#'   [IsoformicExperiment-accessors] for the accessor generics.
#'
#' @examples
#' iso <- example_isoformic()
#' iso
#' dim(iso)
#' assay_names(iso)
#' head(row_data(iso))
#' head(gene_assay(iso, "abundance"))
#'
#' # Subset like a SummarizedExperiment. Stored gene assays, inferential
#' # replicates and uncertainty tables follow the transcripts and samples kept
#' iso[1:10, c("SRR11498039", "SRR11498063")]
#'
#' # Build one from matrices
#' annot <- system.file("extdata", "example_experiment", "annotation", package = "isoformic")
#' tpm <- tx_assay(iso, "abundance")[1:5, 1:3]
#' small <- IsoformicExperiment(
#'   experiment_name = "tiny",
#'   annot_path = annot,
#'   assays = list(abundance = tpm),
#'   col_data = data.frame(sample_id = colnames(tpm), condition = c("a", "a", "b"))
#' )
#' small
#'
#' @export
IsoformicExperiment <- S7::new_class(
  "IsoformicExperiment",
  package = "isoformic",
  properties = list(
    experiment_name = S7::new_property(
      class = S7::class_character,
      default = NA_character_
    ),
    annot_path = S7::new_property(
      class = S7::new_union(NULL, S7::class_character),
      default = NULL,
      setter = function(self, value) {
        if (!is.null(value)) {
          value <- resolve_annot_path(value)
        }
        attr(self, "annot_path") <- value
        self
      }
    ),
    data_path = S7::new_property(
      class = S7::new_union(NULL, S7::class_character),
      default = NULL
    ),
    assays = S7::new_property(
      class = S7::class_list,
      default = list(),
      setter = function(self, value) {
        value <- standardize_assays(value)
        attr(self, "assays") <- value
        # Hashed once here: the derived-value caches key on it at every access
        attr(self, "assays_hash") <- rlang::hash(value)
        self
      }
    ),
    gene_assays = S7::new_property(
      class = S7::class_list,
      default = list(),
      getter = function(self) {
        stored <- prop_raw(self, "gene_assays")
        if (length(stored) > 0L) {
          return(stored)
        }
        derive_gene_assays(self)
      },
      setter = function(self, value) {
        attr(self, "gene_assays") <- standardize_assays(value, level = "gene")
        self
      }
    ),
    col_data = S7::new_property(
      class = S7::new_union(NULL, S7::class_data.frame),
      default = NULL,
      setter = function(self, value) {
        attr(self, "col_data") <- standardize_col_data(value)
        self
      }
    ),
    dea = S7::new_property(class = S7::class_list, default = list()),
    gsea = S7::new_property(class = S7::class_list, default = list()),
    metadata = S7::new_property(class = S7::class_list, default = list()),
    annot_metadata = S7::new_property(
      class = S7::class_list,
      default = list(),
      getter = function(self) {
        stored <- attr(self, "annot_metadata", exact = TRUE)
        if (length(stored) > 0L || is.null(prop_raw(self, "annot_path"))) {
          return(stored)
        }
        kv_to_list(get_annot_table(self, "metadata"))
      }
    ),
    tx_type_palette = S7::new_property(
      class = S7::class_character,
      default = quote(tx_type_palette()),
      setter = function(self, value) {
        if (!length(value)) {
          value <- tx_type_palette()
        }
        if (
          !is.character(value) ||
            is.null(names(value)) ||
            !all(nzchar(names(value)))
        ) {
          cli::cli_abort(
            message = c(
              x = "{.var tx_type_palette} must be a named character vector."
            ),
            class = "isoformic_invalid_tx_type_palette"
          )
        }
        attr(self, "tx_type_palette") <- value
        self
      }
    ),
    annot_data_transcripts = S7::new_property(
      class = S7::class_data.frame,
      getter = function(self) get_annot_table(self, "transcripts")
    ),
    annot_data_genes = S7::new_property(
      class = S7::class_data.frame,
      getter = function(self) get_annot_table(self, "genes")
    ),
    annot_data_exons = S7::new_property(
      class = S7::class_data.frame,
      getter = function(self) get_annot_table(self, "exons")
    ),
    .cache = S7::new_property(
      class = S7::class_environment,
      default = quote(new.env(parent = emptyenv()))
    )
  ),
  constructor = function(
    experiment_name = NULL,
    annot_path = NULL,
    assays = list(),
    col_data = NULL,
    gene_assays = list(),
    dea = list(),
    gsea = list(),
    metadata = list(),
    data_path = NULL,
    tx_type_palette = NULL
  ) {
    if (is.null(experiment_name) || isTRUE(is.na(experiment_name))) {
      experiment_name <- random_experiment_name()
    }
    check_string(experiment_name)
    assays <- standardize_assays(assays)
    if (is.null(col_data) && length(assays)) {
      col_data <- tibble::tibble(sample_id = colnames(assays[[1L]]))
    }
    col_data <- standardize_col_data(col_data)
    if (!is.null(annot_path)) {
      annot_path <- resolve_annot_path(annot_path)
    }
    if (is.null(tx_type_palette)) {
      tx_type_palette <- tx_type_palette()
    }
    dea <- standardize_dea(dea, annot_path)
    obj <- S7::new_object(
      S7::S7_object(),
      experiment_name = experiment_name,
      annot_path = annot_path,
      data_path = data_path,
      assays = assays,
      gene_assays = standardize_assays(gene_assays, level = "gene"),
      col_data = col_data,
      dea = dea,
      gsea = gsea,
      metadata = metadata,
      tx_type_palette = tx_type_palette,
      .cache = new.env(parent = emptyenv())
    )
    obj
  },
  validator = function(self) {
    validate_isoformic(self)
  }
)

# Construction helpers --------------------------------------------------------

#' Generate a Random Experiment Name
#' @keywords internal
#' @noRd
random_experiment_name <- function() {
  seed <- list(Sys.time(), Sys.getpid(), stats::runif(1L))
  paste0("isoformic_", substr(rlang::hash(seed), 1L, 10L))
}

#' Coerce an assay list into a named list of numeric matrices
#'
#' Data frames with a `transcript_id`/`gene_id`/`feature_id` column (or row
#' names) are converted to matrices. All matrices must share dimnames.
#'
#' @keywords internal
#' @noRd
standardize_assays <- function(assays, level = c("transcript", "gene")) {
  level <- rlang::arg_match(level)
  if (is.null(assays)) {
    return(list())
  }
  if (is.matrix(assays) || is.data.frame(assays)) {
    assays <- list(abundance = assays)
  }
  if (!is.list(assays)) {
    cli::cli_abort(
      message = c(x = "{.arg assays} must be a named list of matrices."),
      class = "isoformic_invalid_assays"
    )
  }
  if (!length(assays)) {
    return(list())
  }
  if (is.null(names(assays)) || !all(nzchar(names(assays)))) {
    cli::cli_abort(
      message = c(
        x = "Every element of {.arg assays} must be named.",
        i = "Typical names are {.val abundance}, {.val counts} and {.val length}."
      ),
      class = "isoformic_invalid_assays"
    )
  }
  id_cols <- c(
    paste0(level, "_id"),
    "feature_id",
    "transcript_id",
    "gene_id",
    "id"
  )
  assays <- lapply(assays, function(m) {
    if (is.data.frame(m)) {
      id_col <- intersect(id_cols, colnames(m))[1L]
      if (!is.na(id_col)) {
        ids <- as.character(m[[id_col]])
        m <- m[, setdiff(colnames(m), id_col), drop = FALSE]
        m <- as.matrix(m)
        rownames(m) <- ids
      } else {
        m <- as.matrix(m)
      }
    }
    if (!is.matrix(m) || !is.numeric(m)) {
      cli::cli_abort(
        message = c(x = "Each assay must be a numeric matrix or a data frame."),
        class = "isoformic_invalid_assays"
      )
    }
    if (is.null(rownames(m)) || is.null(colnames(m))) {
      cli::cli_abort(
        message = c(
          x = "Each assay needs {level} IDs as row names and sample IDs as column names."
        ),
        class = "isoformic_invalid_assays"
      )
    }
    m
  })
  ref <- dimnames(assays[[1L]])
  same <- vapply(assays, function(m) identical(dimnames(m), ref), logical(1L))
  if (!all(same)) {
    cli::cli_abort(
      message = c(
        x = "All assays must have identical row and column names.",
        i = "Assays differing from {.val {names(assays)[1]}}: {.val {names(assays)[!same]}}."
      ),
      class = "isoformic_invalid_assays"
    )
  }
  if (anyDuplicated(ref[[1L]])) {
    cli::cli_abort(
      message = c(x = "Assay row names ({level} IDs) must be unique."),
      class = "isoformic_invalid_assays"
    )
  }
  assays
}

#' Coerce sample metadata into a tibble with a leading `sample_id` column
#' @keywords internal
#' @noRd
standardize_col_data <- function(col_data) {
  if (is.null(col_data)) {
    return(NULL)
  }
  if (!is.data.frame(col_data)) {
    col_data <- rlang::try_fetch(
      as.data.frame(col_data, stringsAsFactors = FALSE, optional = TRUE),
      error = function(cnd) {
        cli::cli_abort(
          message = c(x = "{.arg col_data} must be a data frame."),
          class = "isoformic_invalid_col_data",
          parent = cnd
        )
      }
    )
  }
  if (!"sample_id" %in% colnames(col_data)) {
    has_rownames <- !is.null(rownames(col_data)) &&
      !identical(rownames(col_data), as.character(seq_len(nrow(col_data))))
    if (has_rownames) {
      col_data <- cbind(
        sample_id = rownames(col_data),
        as.data.frame(col_data, stringsAsFactors = FALSE)
      )
    } else if (ncol(col_data) >= 1L) {
      colnames(col_data)[1L] <- "sample_id"
    } else {
      cli::cli_abort(
        message = c(
          x = "{.arg col_data} must contain a {.field sample_id} column."
        ),
        class = "isoformic_invalid_col_data"
      )
    }
  }
  col_data <- tibble::as_tibble(as.data.frame(
    col_data,
    stringsAsFactors = FALSE
  ))
  col_data$sample_id <- as.character(col_data$sample_id)
  if (anyDuplicated(col_data$sample_id) || anyNA(col_data$sample_id)) {
    cli::cli_abort(
      message = c(
        x = "{.field sample_id} in {.arg col_data} must be unique and non-missing."
      ),
      class = "isoformic_invalid_col_data"
    )
  }
  dplyr::relocate(col_data, "sample_id")
}

# Validation ------------------------------------------------------------------

#' Validate the stored state of an IsoformicExperiment
#' @keywords internal
#' @noRd
validate_isoformic <- function(self) {
  name <- prop_raw(self, "experiment_name")
  if (!is.character(name) || length(name) != 1L) {
    cli::cli_abort(
      message = c(
        x = "{.var experiment_name} must be a single character string."
      ),
      class = "isoformic_invalid_experiment_name"
    )
  }
  annot_path <- prop_raw(self, "annot_path")
  if (!is.null(annot_path) && !is_annotation_dir(annot_path)) {
    cli::cli_abort(
      message = c(
        x = "The annotation path {.path {annot_path}} is not a complete annotation directory."
      ),
      class = "isoformic_annot_path_dont_exist"
    )
  }
  data_path <- prop_raw(self, "data_path")
  if (!is.null(data_path) && !rlang::is_string(data_path)) {
    cli::cli_abort(
      message = c(
        x = "{.var data_path} must be a single character string or NULL."
      ),
      class = "isoformic_invalid_data_path"
    )
  }
  assays <- prop_raw(self, "assays")
  col_data <- prop_raw(self, "col_data")
  if (length(assays) > 0L) {
    if (is.null(annot_path)) {
      cli::cli_abort(
        message = c(
          x = "Assays require an annotation: set {.var annot_path} first."
        ),
        class = "isoformic_missing_annotation"
      )
    }
    validate_assay_ids(
      ids = rownames(assays[[1L]]),
      universe = annot_ids(self, "transcripts"),
      what = "transcript",
      where = "annotation"
    )
    validate_assay_ids(
      ids = colnames(assays[[1L]]),
      universe = col_data$sample_id,
      what = "sample",
      where = "col_data"
    )
  }
  gene_assays <- prop_raw(self, "gene_assays")
  if (length(gene_assays) > 0L) {
    validate_assay_ids(
      ids = rownames(gene_assays[[1L]]),
      universe = annot_ids(self, "genes"),
      what = "gene",
      where = "annotation"
    )
    validate_assay_ids(
      ids = colnames(gene_assays[[1L]]),
      universe = col_data$sample_id,
      what = "sample",
      where = "col_data"
    )
  }
  dea <- prop_raw(self, "dea")
  if (length(dea) > 0L) {
    validate_dea(self, dea)
  }
  NULL
}

validate_assay_ids <- function(ids, universe, what, where) {
  missing <- ids[!ids %in% universe]
  if (length(missing)) {
    cli::cli_abort(
      message = c(
        x = "{length(missing)} {what} ID{?s} in the assays {?is/are} not present in the {where}.",
        i = "First missing: {.val {utils::head(missing, 5)}}."
      ),
      class = paste0("isoformic_invalid_assay_", what, "_ids")
    )
  }
  invisible(TRUE)
}

validate_dea <- function(self, dea) {
  if (!is.null(dea[["det"]])) {
    det <- dea[["det"]]
    check_de_columns(det, "transcript_id", "det")
    validate_assay_ids(
      ids = unique(det$transcript_id),
      universe = annot_ids(self, "transcripts"),
      what = "transcript",
      where = "annotation (in the DET table)"
    )
  }
  if (!is.null(dea[["deg"]])) {
    deg <- dea[["deg"]]
    check_de_columns(deg, "gene_id", "deg")
    validate_assay_ids(
      ids = unique(deg$gene_id),
      universe = annot_ids(self, "genes"),
      what = "gene",
      where = "annotation (in the DEG table)"
    )
  }
  invisible(TRUE)
}

check_de_columns <- function(x, id_col, what) {
  required <- c(id_col, "log2FC", "pvalue")
  missing_cols <- setdiff(required, colnames(x))
  if (length(missing_cols)) {
    cli::cli_abort(
      message = c(
        x = "The {toupper(what)} table is missing required column{?s}: {.field {missing_cols}}."
      ),
      class = paste0("isoformic_dea_invalid_", what)
    )
  }
  invisible(TRUE)
}

# Annotation access and caching ----------------------------------------------

# Package-level cache of annotation tables, keyed by file signature so that
# objects sharing an annotation share the in-memory tables.
annotation_cache <- new.env(parent = emptyenv())

#' Fetch a cached value or compute and store it
#' @keywords internal
#' @noRd
cache_fetch <- function(cache, key, signature, compute) {
  hit <- cache[[key]]
  if (!is.null(hit) && identical(hit$signature, signature)) {
    return(hit$value)
  }
  value <- compute()
  assign(key, list(signature = signature, value = value), envir = cache)
  value
}

#' Read (cached) one annotation table for an object
#' @keywords internal
#' @noRd
get_annot_table <- function(self, type) {
  annot_path <- prop_raw(self, "annot_path")
  if (is.null(annot_path)) {
    return(empty_annotation_table(type))
  }
  signature <- annot_signature(annot_path, type)
  cache_fetch(
    cache = annotation_cache,
    key = paste0(type, "::", annot_path),
    signature = signature,
    compute = function() read_annotation_table(annot_path, type)
  )
}

#' Is an annotation table of an object in the annotation cache, up to date?
#' @keywords internal
#' @noRd
annot_table_cached <- function(self, type) {
  annot_path <- prop_raw(self, "annot_path")
  if (is.null(annot_path)) {
    return(FALSE)
  }
  hit <- annotation_cache[[paste0(type, "::", annot_path)]]
  !is.null(hit) && identical(hit$signature, annot_signature(annot_path, type))
}

#' Row count of an annotation table without loading it when not cached
#' @keywords internal
#' @noRd
annot_nrow <- function(self, type) {
  annot_path <- prop_raw(self, "annot_path")
  if (is.null(annot_path)) {
    return(0L)
  }
  hit <- annotation_cache[[paste0(type, "::", annot_path)]]
  if (
    !is.null(hit) && identical(hit$signature, annot_signature(annot_path, type))
  ) {
    return(nrow(hit$value))
  }
  parquet_nrow(annot_file_path(annot_path, type))
}

#' Identifier vector of an annotation table (cached)
#' @keywords internal
#' @noRd
annot_ids <- function(self, type = c("transcripts", "genes")) {
  type <- rlang::arg_match(type)
  tbl <- get_annot_table(self, type)
  if (identical(type, "transcripts")) tbl$transcript_id else tbl$gene_id
}

#' Signature of the annotation directory as a whole
#' @keywords internal
#' @noRd
annot_full_signature <- function(self) {
  annot_path <- prop_raw(self, "annot_path")
  if (is.null(annot_path)) {
    return("none")
  }
  paste(
    annot_signature(annot_path, "genes"),
    annot_signature(annot_path, "transcripts"),
    annot_signature(annot_path, "exons"),
    sep = ";"
  )
}

#' Signature of the transcript assays (dimnames and values), hashed when set
#' @keywords internal
#' @noRd
assays_signature <- function(self) {
  prop_raw(self, "assays_hash") %||%
    rlang::hash(prop_raw(self, "assays"))
}

#' Signature of the stored (explicitly set) gene assays
#' @keywords internal
#' @noRd
gene_assays_signature <- function(self) {
  stored <- prop_raw(self, "gene_assays")
  if (length(stored) > 0L) rlang::hash(stored) else "derived"
}

#' Fetch a derived value from the per-object cache
#' @keywords internal
#' @noRd
object_cache_fetch <- function(self, key, signature, compute) {
  cache_fetch(prop_raw(self, ".cache"), key, signature, compute)
}

#' A derived value of the per-object cache, or `NULL` when absent or stale
#' @keywords internal
#' @noRd
object_cache_peek <- function(self, key, signature) {
  cache <- prop_raw(self, ".cache")
  hit <- if (is.null(cache)) NULL else cache[[key]]
  if (!is.null(hit) && identical(hit$signature, signature)) hit$value else NULL
}

# Gene-level derivation -------------------------------------------------------

#' Derive gene-level assays by summarising transcript assays
#'
#' Sums isoform values per gene (`rowsum()`, C-level). An assay named
#' `length` is summarised as the abundance-weighted mean of isoform lengths
#' (as `tximport` does), falling back to a plain mean when a gene has zero
#' abundance.
#'
#' @keywords internal
#' @noRd
derive_gene_assays <- function(self) {
  assays <- prop_raw(self, "assays")
  if (!length(assays)) {
    return(list())
  }
  signature <- paste(assays_signature(self), annot_full_signature(self))
  object_cache_fetch(self, "gene_assays", signature, function() {
    tx2gene <- get_annot_table(self, "transcripts")
    gene_ids <- tx2gene$gene_id[
      match(rownames(assays[[1L]]), tx2gene$transcript_id)
    ]
    summarize_assays_to_gene(assays, gene_ids)
  })
}

#' Summarise a list of transcript matrices to gene matrices
#'
#' @param assays Named list of transcript-level matrices.
#' @param gene_ids Character vector, the gene of each row of the matrices.
#'
#' @keywords internal
#' @noRd
summarize_assays_to_gene <- function(assays, gene_ids) {
  keep <- !is.na(gene_ids)
  all_kept <- all(keep)
  groups <- if (all_kept) gene_ids else gene_ids[keep]
  # No copy of an assay when every transcript has a gene, the usual case
  kept_rows <- function(m) if (all_kept) m else m[keep, , drop = FALSE]
  abundance <- assays[["abundance"]]
  out <- lapply(names(assays), function(nm) {
    m <- kept_rows(assays[[nm]])
    if (identical(tolower(nm), "length") && !is.null(abundance)) {
      w <- kept_rows(abundance)
      num <- rowsum(m * w, groups, reorder = TRUE)
      den <- rowsum(w, groups, reorder = TRUE)
      n <- rowsum(rep(1, nrow(m)), groups, reorder = TRUE)
      simple <- rowsum(m, groups, reorder = TRUE) / as.vector(n)
      res <- num / den
      zero <- den == 0
      res[zero] <- simple[zero]
      return(res)
    }
    rowsum(m, groups, reorder = TRUE)
  })
  names(out) <- names(assays)
  out
}

#' Summarize Transcript-Level Matrices to Gene Level
#'
#' Sums the values of all isoforms of a gene. A `length` assay is combined
#' as the abundance-weighted mean when an `abundance` assay is supplied in
#' the same list (the `tximport` convention). This is the function used by
#' [gene_assays()] and can be called directly on matrices.
#'
#' @param x A numeric matrix (transcripts in rows) or a named list of such
#'   matrices sharing row names.
#' @param tx_to_gene Data frame with `transcript_id` and `gene_id` columns.
#'
#' @returns A matrix, or a named list of matrices, with genes in rows.
#'
#' @examples
#' iso <- example_isoformic()
#' g <- summarize_to_gene(tx_assay(iso, "abundance"), tx_to_gene(iso))
#' head(g)
#'
#' @export
summarize_to_gene <- function(x, tx_to_gene) {
  single <- is.matrix(x)
  assays <- standardize_assays(if (single) list(abundance = x) else x)
  gene_ids <- tx_to_gene$gene_id[
    match(rownames(assays[[1L]]), tx_to_gene$transcript_id)
  ]
  res <- summarize_assays_to_gene(assays, gene_ids)
  if (single) res[[1L]] else res
}

# Base generics ---------------------------------------------------------------

S7::method(print, IsoformicExperiment) <- function(x, ...) {
  assays <- prop_raw(x, "assays")
  col_data <- prop_raw(x, "col_data")
  annot_path <- prop_raw(x, "annot_path")
  fmt <- function(n) format(n, big.mark = ",")

  cat("<IsoformicExperiment> ", x@experiment_name, "\n", sep = "")
  if (is.null(annot_path)) {
    cat("  annotation : <none>\n")
  } else {
    meta <- x@annot_metadata
    annot_label <- meta[["annot_name"]] %||% fs::path_file(annot_path)
    cat(
      "  annotation : ",
      annot_label,
      " (",
      fmt(annot_nrow(x, "genes")),
      " genes, ",
      fmt(annot_nrow(x, "transcripts")),
      " transcripts, ",
      fmt(annot_nrow(x, "exons")),
      " exons)\n",
      sep = ""
    )
  }
  n_samples <- if (is.null(col_data)) 0L else nrow(col_data)
  cd_cols <- if (is.null(col_data)) {
    character()
  } else {
    setdiff(colnames(col_data), "sample_id")
  }
  cat(
    "  samples    : ",
    fmt(n_samples),
    if (length(cd_cols)) {
      paste0(" [", paste(cd_cols, collapse = ", "), "]")
    } else {
      ""
    },
    "\n",
    sep = ""
  )
  if (length(assays)) {
    n_genes <- length(unique(gene_ids_of(x, rownames(assays[[1L]]))))
    cat(
      "  transcripts: ",
      fmt(nrow(assays[[1L]])),
      " assayed from ",
      fmt(n_genes),
      " genes\n",
      sep = ""
    )
    cat(
      "  assays     : ",
      paste(names(assays), collapse = ", "),
      "\n",
      sep = ""
    )
    stored_gene <- prop_raw(x, "gene_assays")
    cat(
      "  gene assays: ",
      if (length(stored_gene)) {
        paste(names(stored_gene), collapse = ", ")
      } else {
        "derived from transcript assays on demand"
      },
      "\n",
      sep = ""
    )
  } else {
    cat("  assays     : <none>\n")
  }
  dea <- prop_raw(x, "dea")
  de_summary <- character()
  if (!is.null(dea[["det"]])) {
    de_summary <- c(
      de_summary,
      paste0("det (", fmt(nrow(dea$det)), " transcripts)")
    )
  }
  if (!is.null(dea[["deg"]])) {
    de_summary <- c(de_summary, paste0("deg (", fmt(nrow(dea$deg)), " genes)"))
  }
  cat(
    "  DE results : ",
    if (length(de_summary)) paste(de_summary, collapse = ", ") else "<none>",
    "\n",
    sep = ""
  )
  gsea <- prop_raw(x, "gsea")
  if (length(gsea)) {
    cat("  enrichment : ", paste(names(gsea), collapse = ", "), "\n", sep = "")
  }
  data_path <- prop_raw(x, "data_path")
  if (!is.null(data_path)) {
    cat("  on disk    : ", display_data_path(data_path), "\n", sep = "")
  }
  invisible(x)
}

# The printed form of `data_path`: relative to the working directory when the
# directory lies under it, absolute otherwise. `data_path` itself stays
# absolute, since the lazy readers and write_isoformic() rely on it.
display_data_path <- function(data_path) {
  real_or_abs <- function(path) {
    out <- tryCatch(fs::path_real(path), error = function(e) NULL)
    as.character(out %||% fs::path_abs(path))
  }
  path <- real_or_abs(data_path)
  wd <- real_or_abs(getwd())
  if (!identical(path, wd) && fs::path_has_parent(path, wd)) {
    return(as.character(fs::path_rel(path, start = wd)))
  }
  as.character(data_path)
}

S7::method(dim, IsoformicExperiment) <- function(x) {
  assays <- prop_raw(x, "assays")
  if (!length(assays)) {
    return(c(0L, 0L))
  }
  dim(assays[[1L]])
}

S7::method(dimnames, IsoformicExperiment) <- function(x) {
  assays <- prop_raw(x, "assays")
  if (!length(assays)) {
    return(list(NULL, NULL))
  }
  dimnames(assays[[1L]])
}

S7::method(`[`, IsoformicExperiment) <- function(x, i, j, ..., drop = FALSE) {
  assays <- prop_raw(x, "assays")
  if (!length(assays)) {
    cli::cli_abort(
      message = c(
        x = "Cannot subset an {.cls IsoformicExperiment} without assays."
      ),
      class = "isoformic_no_assays"
    )
  }
  has_i <- !missing(i)
  has_j <- !missing(j)
  all_tx <- rownames(assays[[1L]])
  all_samples <- colnames(assays[[1L]])
  # Everything below is subset by ID, never by position
  tx_ids <- if (has_i) subset_ids(i, all_tx, "transcript") else all_tx
  sample_ids <- if (has_j) subset_ids(j, all_samples, "sample") else all_samples
  # The assays share their dimnames: the positions are matched once
  tx_rows <- match(tx_ids, all_tx)
  sample_cols <- match(sample_ids, all_samples)
  new_assays <- lapply(assays, function(m) {
    m[tx_rows, sample_cols, drop = FALSE]
  })
  col_data <- prop_raw(x, "col_data")
  new_col_data <- col_data[
    match(sample_ids, col_data$sample_id),
    ,
    drop = FALSE
  ]
  stored_gene <- prop_raw(x, "gene_assays")
  if (length(stored_gene)) {
    gene_ids <- rownames(stored_gene[[1L]])
    if (has_i) {
      gene_ids <- intersect(gene_ids, gene_ids_of(x, tx_ids))
    }
    # Stored gene assays may hold some of the samples only
    stored_gene <- lapply(stored_gene, function(m) {
      m[
        gene_ids[gene_ids %in% rownames(m)],
        sample_ids[sample_ids %in% colnames(m)],
        drop = FALSE
      ]
    })
    if (any(vapply(stored_gene, function(m) min(dim(m)), integer(1L)) == 0L)) {
      # Nothing stored is left: the gene assays are derived again
      stored_gene <- list()
    }
  }
  metadata <- prop_raw(x, "metadata")
  if (!is.null(metadata$inf_reps)) {
    metadata$inf_reps <- subset_inf_reps(metadata$inf_reps, tx_ids, sample_ids)
  }
  dea <- prop_raw(x, "dea")
  has_tables <- any(uncertainty_table_names() %in% names(dea))
  if (has_tables && has_i) {
    dea <- subset_uncertainty_tables(dea, tx_ids)
    if (!is.null(dea[["uncertainty"]])) {
      metadata$uncertainty_n_assessed <- sum(dea[["uncertainty"]]$assessed)
    }
  }
  if (has_tables && !setequal(sample_ids, all_samples)) {
    cli::cli_inform(
      message = c(
        i = "The stored uncertainty tables are not recomputed for this subset of the samples.",
        i = "Recompute them with {.code add_uncertainty(x)} (and {.arg group}) if they must follow the subset."
      ),
      class = "isoformic_stale_uncertainty"
    )
  }
  if (!identical(tx_ids, all_tx)) {
    # Other transcripts: the subset gets its own cache, so that it does not
    # replace the entries of `x` (keyed by the transcripts) nor keep them
    # alive. The annotation join does not depend on the transcripts.
    old_cache <- prop_raw(x, ".cache")
    new_cache <- new.env(parent = emptyenv())
    shared <- old_cache[["tx_gene_annotation"]]
    if (!is.null(shared)) {
      assign("tx_gene_annotation", shared, envir = new_cache)
    }
    attr(x, ".cache") <- new_cache
  }
  # Everything is set before the one validated assignment (`col_data`), which
  # checks the new assays and the stored gene assays against the new samples.
  # The new assays are parts of validated assays: only their hash is renewed.
  attr(x, "gene_assays") <- stored_gene
  attr(x, "metadata") <- metadata
  attr(x, "dea") <- dea
  attr(x, "assays") <- new_assays
  attr(x, "assays_hash") <- rlang::hash(new_assays)
  x@col_data <- new_col_data
  x
}

#' IDs selected by a subscript of `[`
#'
#' Accepts what a matrix accepts (positions, negative positions, logical,
#' IDs); a factor is read as IDs, not as its integer codes.
#' @keywords internal
#' @noRd
subset_ids <- function(index, ids, what) {
  if (is.factor(index)) {
    index <- as.character(index)
  }
  if (is.logical(index) && length(index) > length(ids)) {
    cli::cli_abort(
      message = c(
        x = "The logical {what} subscript has {length(index)} values, but the object has {length(ids)} {what}{?s}."
      ),
      class = "isoformic_invalid_subset"
    )
  }
  positions <- stats::setNames(seq_along(ids), ids)
  picked <- tryCatch(positions[index], error = function(cnd) NA_integer_)
  if (anyNA(picked)) {
    cli::cli_abort(
      message = c(
        x = "The {what} subscript selects {what}s that are not in the object.",
        i = "Use {what} IDs, positions or a logical vector without missing values."
      ),
      class = "isoformic_invalid_subset"
    )
  }
  if (anyDuplicated(picked)) {
    cli::cli_abort(
      message = c(x = "The {what} subscript selects a {what} more than once."),
      class = "isoformic_invalid_subset"
    )
  }
  if (!length(picked)) {
    cli::cli_abort(
      message = c(
        x = "The {what} subscript selects no {what}.",
        i = "An {.cls IsoformicExperiment} needs at least one transcript and one sample."
      ),
      class = "isoformic_invalid_subset"
    )
  }
  ids[unname(picked)]
}

#' Inferential replicates of the given transcripts and samples, as a view
#'
#' The draws are not copied: the list keeps the matrices of `reps` and gets
#' a `selection` attribute with the IDs of the kept transcripts and samples
#' (the list itself, a vector of pointers, is the only new object). Every
#' consumer reads the draws through `inf_reps_reader()` (one selected draw at
#' a time), `inf_reps_draw_index()` (the positions of the selection in a
#' draw) or `inf_reps_moments()`, which apply the selection. A subset of a
#' subset replaces the selection: its transcripts and samples are among
#' those of the first one, so selecting them in the parent draws gives the
#' same rows and columns as selecting them in the selected draws.
#'
#' The draws may hold fewer transcripts than the assays ([add_uncertainty()]
#' keeps the expressed ones), so rows are matched by ID. Draws that do not
#' have every sample are left as they are; their consumers report them.
#' When no kept transcript has draws, the result is `NULL`, unless the
#' `moments` attribute covers the kept samples: the subset then keeps the
#' moments (which cover every transcript) with zero-row draws, as plain
#' copies, so that it does not hold the draws of the parent. Its consumers
#' treat it as an object without draws (`draws_view()` aborts with
#' `isoformic_no_inferential_replicates`, `resolve_inf_reps()` gives
#' `NULL`), but its overdispersion comes from the moments.
#'
#' The `selection` attribute is lost by list operations (`lapply()`, `[`);
#' the consumers therefore read the stored draws through
#' `object_inf_reps()`, which restores the selection from the transcripts
#' and samples of the object.
#' @keywords internal
#' @noRd
subset_inf_reps <- function(reps, tx_ids, sample_ids) {
  if (!is.list(reps) || is.data.frame(reps) || !length(reps)) {
    return(reps)
  }
  attr(reps, "selection") <- list(transcripts = tx_ids, samples = sample_ids)
  if (!inf_reps_has_rows(reps)) {
    # No kept transcript has draws. Without moments of the kept samples the
    # object has no replicates left; with them, only the moments are kept
    # (zero-row copies of the draws), not the draws of the parent
    mom <- attr(reps, "moments", exact = TRUE)
    if (is.null(mom) || !all(sample_ids %in% colnames(mom$mean))) {
      return(NULL)
    }
    return(inf_reps_materialize(reps))
  }
  # The moments keep every transcript (the overdispersion prior uses them
  # all); inf_reps_moments() gives the kept samples
  reps
}

#' The stored inferential replicates of an object, with its selection
#'
#' `x@metadata$inf_reps` read as the draws of the transcripts and samples of
#' the object. The `selection` that `x[i, j]` records is an attribute, which
#' list operations on the draws drop (`lapply()`, `[`); without it, the
#' matrices of a subset would still be those of its parent. When the draws
#' have no selection but hold transcripts or samples that are not in the
#' object, the selection is rebuilt from the object (its transcripts and
#' samples, matched by ID), so that the results are those of the kept
#' draws. Draws that only cover the object are returned as they are.
#' @keywords internal
#' @noRd
object_inf_reps <- function(self, reps = self@metadata$inf_reps) {
  if (
    !is.list(reps) ||
      is.data.frame(reps) ||
      !length(reps) ||
      !is.null(inf_reps_selection(reps))
  ) {
    return(reps)
  }
  dn <- dimnames(reps[[1L]])
  if (length(dn) != 2L || is.null(dn[[1L]]) || is.null(dn[[2L]])) {
    return(reps)
  }
  obj_dn <- dimnames(self)
  if (all(dn[[1L]] %in% obj_dn[[1L]]) && all(dn[[2L]] %in% obj_dn[[2L]])) {
    return(reps)
  }
  attr(reps, "selection") <- list(
    transcripts = obj_dn[[1L]],
    samples = obj_dn[[2L]]
  )
  reps
}

#' The selection of a view of inferential replicates, `NULL` when none
#' @keywords internal
#' @noRd
inf_reps_selection <- function(reps) {
  attr(reps, "selection", exact = TRUE)
}

#' Positions of the selection in one draw
#'
#' `NULL` when `reps` has no selection or the draw is left as it is (not a
#' matrix with transcript names, or without every selected sample);
#' otherwise the index of `inf_reps_index()`. `cache` (an environment) keeps
#' the index of the last dimnames seen, so that draws that share their
#' dimnames are matched once.
#' @keywords internal
#' @noRd
inf_reps_draw_index <- function(m, selection, cache = NULL) {
  if (is.null(selection)) {
    return(NULL)
  }
  dn <- dimnames(m)
  if (
    is.null(dim(m)) ||
      length(dn) != 2L ||
      is.null(dn[[1L]]) ||
      !all(selection$samples %in% dn[[2L]])
  ) {
    return(NULL)
  }
  if (!is.null(cache) && identical(dn, cache$dimnames)) {
    return(cache$index)
  }
  index <- inf_reps_index(dn, selection$transcripts, selection$samples)
  if (!is.null(cache)) {
    cache$dimnames <- dn
    cache$index <- index
  }
  index
}

#' Read the draws of inferential replicates with their selection applied
#'
#' Returns a function of `k` that gives draw `k` restricted to the selected
#' transcripts and samples (see `subset_inf_reps()`). Replicates without a
#' selection, and a draw that already is the selection, are returned without
#' a copy; otherwise only the one draw read is copied.
#' @keywords internal
#' @noRd
inf_reps_reader <- function(reps) {
  selection <- inf_reps_selection(reps)
  cache <- new.env(parent = emptyenv())
  function(k) {
    m <- reps[[k]]
    index <- inf_reps_draw_index(m, selection, cache)
    if (is.null(index)) {
      return(m)
    }
    subset_one_draw(m, index)
  }
}

#' Moments of inferential replicates with their selection applied
#'
#' The moments keep every transcript; with a selection that they cover, they
#' are restricted to the selected samples.
#' @keywords internal
#' @noRd
inf_reps_moments <- function(reps) {
  mom <- attr(reps, "moments", exact = TRUE)
  samples <- inf_reps_selection(reps)$samples
  if (
    is.null(mom) || is.null(samples) || !all(samples %in% colnames(mom$mean))
  ) {
    return(mom)
  }
  if (identical(samples, colnames(mom$mean))) {
    return(mom)
  }
  cols <- match(samples, colnames(mom$mean))
  mom$mean <- mom$mean[, cols, drop = FALSE]
  mom$var <- mom$var[, cols, drop = FALSE]
  mom
}

#' Inferential replicates with their selection applied, as a plain list
#'
#' A copy of the selected draws: what `x[i, j]` stored before the draws
#' became a view. Replicates without a selection are returned as they are.
#' @keywords internal
#' @noRd
inf_reps_materialize <- function(reps) {
  if (is.null(inf_reps_selection(reps))) {
    return(reps)
  }
  read <- inf_reps_reader(reps)
  out <- lapply(seq_along(reps), read)
  extra <- attributes(reps)
  extra$selection <- NULL
  if (!is.null(extra$moments)) {
    extra$moments <- inf_reps_moments(reps)
  }
  attributes(out) <- extra
  out
}

#' Does any draw have rows once its selection is applied?
#'
#' Draws left as they are count by their own rows; draws that are not
#' matrices count as having rows.
#' @keywords internal
#' @noRd
inf_reps_has_rows <- function(reps) {
  selection <- inf_reps_selection(reps)
  cache <- new.env(parent = emptyenv())
  for (k in seq_along(reps)) {
    m <- reps[[k]]
    if (is.null(dim(m))) {
      return(TRUE)
    }
    index <- inf_reps_draw_index(m, selection, cache)
    n <- if (is.null(index)) nrow(m) else length(index$rows)
    if (n > 0L) {
      return(TRUE)
    }
  }
  FALSE
}

#' Positions and new dimnames of the kept transcripts and samples of a draw
#' @keywords internal
#' @noRd
inf_reps_index <- function(dn, tx_ids, sample_ids) {
  rows <- match(tx_ids, dn[[1L]])
  rows <- rows[!is.na(rows)]
  cols <- match(sample_ids, dn[[2L]])
  new_dn <- dn
  new_dn[[1L]] <- dn[[1L]][rows]
  new_dn[[2L]] <- dn[[2L]][cols]
  list(
    rows = rows,
    cols = cols,
    dimnames = new_dn,
    unchanged = identical(rows, seq_along(dn[[1L]])) &&
      identical(cols, seq_along(dn[[2L]]))
  )
}

#' One draw restricted to an index of `inf_reps_index()`
#'
#' A plain matrix that already is the selection is returned without a copy.
#' @keywords internal
#' @noRd
subset_one_draw <- function(m, index) {
  if (index$unchanged && is.matrix(m) && length(attributes(m)) == 2L) {
    return(m)
  }
  res <- m[index$rows, index$cols, drop = FALSE]
  if (is.matrix(res)) {
    dimnames(res) <- index$dimnames
  }
  res
}

#' Uncertainty tables of the given transcripts
#' @keywords internal
#' @noRd
subset_uncertainty_tables <- function(dea, tx_ids) {
  unc <- dea[["uncertainty"]]
  if (is.data.frame(unc) && "transcript_id" %in% colnames(unc)) {
    dea[["uncertainty"]] <- unc[unc$transcript_id %in% tx_ids, , drop = FALSE]
  }
  pairs <- dea[["ambiguity"]]
  if (
    is.data.frame(pairs) &&
      all(c("transcript_id_1", "transcript_id_2") %in% colnames(pairs))
  ) {
    dea[["ambiguity"]] <- pairs[
      pairs$transcript_id_1 %in% tx_ids & pairs$transcript_id_2 %in% tx_ids,
      ,
      drop = FALSE
    ]
  }
  conf <- dea[["confidence"]]
  if (
    is.data.frame(conf) &&
      all(c("feature_id", "feature_type") %in% colnames(conf))
  ) {
    # A collapsed ambiguity group is kept only when all its members are
    keep <- conf$feature_id %in% tx_ids
    is_group <- conf$feature_type %in% "group"
    if (any(is_group) && "transcript_ids" %in% colnames(conf)) {
      members <- strsplit(
        as.character(conf$transcript_ids[is_group]),
        ";",
        fixed = TRUE
      )
      keep[is_group] <- vapply(
        members,
        function(ids) all(ids %in% tx_ids),
        logical(1L)
      )
    } else {
      keep[is_group] <- TRUE
    }
    dea[["confidence"]] <- conf[keep, , drop = FALSE]
  }
  dea
}

#' Map transcript IDs to gene IDs using the cached annotation
#' @keywords internal
#' @noRd
gene_ids_of <- function(self, transcript_ids) {
  tx <- get_annot_table(self, "transcripts")
  tx$gene_id[match(transcript_ids, tx$transcript_id)]
}

# Read a property as stored, without its getter. S7 >= 1.0.0 stores a NULL
# property as an internal marker symbol; older versions drop the attribute.
# No property of the class holds a symbol, so any symbol means NULL.
prop_raw <- function(x, name) {
  value <- attr(x, name, exact = TRUE)
  if (is.symbol(value)) {
    return(NULL)
  }
  value
}
