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
#' Derived results (gene assays, row data, isoform fractions) are cached per
#' object in the `.cache` environment, keyed by a hash of their inputs, and
#' therefore survive copies but never go stale. See [clear_cache()].
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
#'     (transcripts), `deg` (genes) and `cutoffs`. Use [de_tx<-] and
#'     [de_gene<-] to set them.}
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
#' @param experiment_name Single string naming the experiment.
#' @param annot_path Annotation directory or GFF3/GTF file.
#' @param assays Named list of transcript-level matrices (or data frames
#'   with a `transcript_id` column).
#' @param col_data Data frame of sample metadata with a `sample_id` column
#'   (or sample IDs as row names / first column).
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
#' # Subset like a SummarizedExperiment
#' iso[1:10, c("SRR11498039", "SRR11498063")]
#'
#' # Build one from matrices
#' annot <- system.file("extdata", "example_experiment", "annotation", package = "isoformic")
#' tpm <- assay(iso, "abundance")[1:5, 1:3]
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
        stored <- attr(self, "gene_assays", exact = TRUE)
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
        if (length(stored) > 0L || is.null(attr(self, "annot_path"))) {
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
  name <- attr(self, "experiment_name", exact = TRUE)
  if (!is.character(name) || length(name) != 1L) {
    cli::cli_abort(
      message = c(
        x = "{.var experiment_name} must be a single character string."
      ),
      class = "isoformic_invalid_experiment_name"
    )
  }
  annot_path <- attr(self, "annot_path", exact = TRUE)
  if (!is.null(annot_path) && !is_annotation_dir(annot_path)) {
    cli::cli_abort(
      message = c(
        x = "The annotation path {.path {annot_path}} is not a complete annotation directory."
      ),
      class = "isoformic_annot_path_dont_exist"
    )
  }
  data_path <- attr(self, "data_path", exact = TRUE)
  if (!is.null(data_path) && !rlang::is_string(data_path)) {
    cli::cli_abort(
      message = c(
        x = "{.var data_path} must be a single character string or NULL."
      ),
      class = "isoformic_invalid_data_path"
    )
  }
  assays <- attr(self, "assays", exact = TRUE)
  col_data <- attr(self, "col_data", exact = TRUE)
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
  gene_assays <- attr(self, "gene_assays", exact = TRUE)
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
  dea <- attr(self, "dea", exact = TRUE)
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
  annot_path <- attr(self, "annot_path", exact = TRUE)
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

#' Row count of an annotation table without loading it when not cached
#' @keywords internal
#' @noRd
annot_nrow <- function(self, type) {
  annot_path <- attr(self, "annot_path", exact = TRUE)
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
  annot_path <- attr(self, "annot_path", exact = TRUE)
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
  attr(self, "assays_hash", exact = TRUE) %||%
    rlang::hash(attr(self, "assays", exact = TRUE))
}

#' Fetch a derived value from the per-object cache
#' @keywords internal
#' @noRd
object_cache_fetch <- function(self, key, signature, compute) {
  cache_fetch(attr(self, ".cache", exact = TRUE), key, signature, compute)
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
  assays <- attr(self, "assays", exact = TRUE)
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
  groups <- gene_ids[keep]
  abundance <- assays[["abundance"]]
  out <- lapply(names(assays), function(nm) {
    m <- assays[[nm]][keep, , drop = FALSE]
    if (identical(tolower(nm), "length") && !is.null(abundance)) {
      w <- abundance[keep, , drop = FALSE]
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
#' g <- summarize_to_gene(assay(iso, "abundance"), tx_to_gene(iso))
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
  assays <- attr(x, "assays", exact = TRUE)
  col_data <- attr(x, "col_data", exact = TRUE)
  annot_path <- attr(x, "annot_path", exact = TRUE)
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
    stored_gene <- attr(x, "gene_assays", exact = TRUE)
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
  dea <- attr(x, "dea", exact = TRUE)
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
  gsea <- attr(x, "gsea", exact = TRUE)
  if (length(gsea)) {
    cat("  enrichment : ", paste(names(gsea), collapse = ", "), "\n", sep = "")
  }
  data_path <- attr(x, "data_path", exact = TRUE)
  if (!is.null(data_path)) {
    cat("  on disk    : ", data_path, "\n", sep = "")
  }
  invisible(x)
}

S7::method(dim, IsoformicExperiment) <- function(x) {
  assays <- attr(x, "assays", exact = TRUE)
  if (!length(assays)) {
    return(c(0L, 0L))
  }
  dim(assays[[1L]])
}

S7::method(dimnames, IsoformicExperiment) <- function(x) {
  assays <- attr(x, "assays", exact = TRUE)
  if (!length(assays)) {
    return(list(NULL, NULL))
  }
  dimnames(assays[[1L]])
}

S7::method(`[`, IsoformicExperiment) <- function(x, i, j, ..., drop = FALSE) {
  assays <- attr(x, "assays", exact = TRUE)
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
  subset_mat <- function(m) {
    if (has_i && has_j) {
      m[i, j, drop = FALSE]
    } else if (has_i) {
      m[i, , drop = FALSE]
    } else if (has_j) {
      m[, j, drop = FALSE]
    } else {
      m
    }
  }
  new_assays <- lapply(assays, subset_mat)
  sample_ids <- colnames(new_assays[[1L]])
  col_data <- attr(x, "col_data", exact = TRUE)
  new_col_data <- col_data[
    match(sample_ids, col_data$sample_id),
    ,
    drop = FALSE
  ]
  stored_gene <- attr(x, "gene_assays", exact = TRUE)
  if (length(stored_gene) && has_j) {
    stored_gene <- lapply(stored_gene, function(m) {
      m[, sample_ids, drop = FALSE]
    })
  }
  x@assays <- new_assays
  x@col_data <- new_col_data
  attr(x, "gene_assays") <- stored_gene
  x
}

#' Map transcript IDs to gene IDs using the cached annotation
#' @keywords internal
#' @noRd
gene_ids_of <- function(self, transcript_ids) {
  tx <- get_annot_table(self, "transcripts")
  tx$gene_id[match(transcript_ids, tx$transcript_id)]
}
