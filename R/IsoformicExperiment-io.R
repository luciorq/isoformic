# Persistence and conversion --------------------------------------------------

#' Save and Restore IsoformicExperiment Objects as Parquet
#'
#' @description
#' `write_isoformic()` writes an [IsoformicExperiment] to a directory of
#' Parquet files and `read_isoformic()` restores it. The layout is:
#'
#' ```
#' <path>/
#'   experiment.parquet        key/value: experiment_name, annot_path, ...
#'   col_data.parquet          sample metadata
#'   metadata.parquet          key/value: scalar entries of `metadata`
#'   tx_type_palette.parquet   transcript_type, color
#'   assays/<name>.parquet     transcript_id + one column per sample
#'   gene_assays/<name>.parquet  only explicitly set gene assays
#'   dea/det.parquet, dea/deg.parquet, dea/cutoffs.parquet
#'   gsea/<name>.parquet
#'   annotation/               copy of the annotation directory (optional)
#' ```
#'
#' With `copy_annotation = TRUE` the annotation directory is copied inside
#' the experiment directory so that the result is self-contained and
#' portable. Otherwise `annot_path` is stored as an absolute path and must
#' still exist when the object is read back (it can be overridden with the
#' `annot_path` argument of `read_isoformic()`).
#'
#' @param self An [IsoformicExperiment].
#' @param path Directory to write to or read from. For `write_isoformic()`
#'   defaults to the object's `data_path`, or to
#'   `get_isoformic_cache("experiments", experiment_name)`.
#' @param overwrite Logical, replace an existing directory.
#' @param copy_annotation Logical, copy the annotation directory inside
#'   `path`.
#' @param annot_path Optional annotation directory overriding the stored
#'   one when reading.
#' @param ... Passed on to methods.
#'
#' @returns `write_isoformic()` returns the object with `data_path` set,
#'   invisibly. `read_isoformic()` returns an [IsoformicExperiment].
#'
#' @examples
#' iso <- example_isoformic()
#' dir <- file.path(tempdir(), "iso_example")
#' iso <- write_isoformic(iso, dir, overwrite = TRUE, copy_annotation = TRUE)
#' list.files(dir, recursive = TRUE)
#' iso2 <- read_isoformic(dir)
#' identical(assay(iso), assay(iso2))
#'
#' @name isoformic-io
#' @rdname isoformic-io
NULL

#' @rdname isoformic-io
#' @export
write_isoformic <- S7::new_generic(
  "write_isoformic",
  "self",
  function(self, path = NULL, overwrite = FALSE, copy_annotation = FALSE, ...) {
    S7::S7_dispatch()
  }
)

S7::method(write_isoformic, IsoformicExperiment) <- function(
  self,
  path = NULL,
  overwrite = FALSE,
  copy_annotation = FALSE,
  ...
) {
  if (is.null(path)) {
    path <- self@data_path %||%
      get_isoformic_cache("experiments", self@experiment_name)
  }
  check_string(path)
  path <- as.character(fs::path_abs(path))
  if (fs::dir_exists(path) && length(fs::dir_ls(path))) {
    if (!isTRUE(overwrite)) {
      cli::cli_abort(
        message = c(
          x = "Directory {.path {path}} already exists and is not empty.",
          i = "Use {.code overwrite = TRUE} to replace it."
        ),
        class = "isoformic_dir_exists"
      )
    }
    if (!fs::file_exists(fs::path(path, "experiment.parquet"))) {
      cli::cli_abort(
        message = c(
          x = "{.path {path}} is not an isoformic experiment directory.",
          i = "Refusing to delete it; remove it yourself or choose another {.arg path}."
        ),
        class = "isoformic_not_experiment_dir"
      )
    }
  }
  # Write into a sibling directory and swap it in at the end, so that the
  # data being replaced (including an annotation copied inside it) is still
  # there while the new directory is built, and survives a failed write.
  final_path <- path
  path <- paste0(final_path, ".writing-", Sys.getpid())
  if (fs::dir_exists(path)) {
    fs::dir_delete(path)
  }
  on.exit(
    if (fs::dir_exists(path)) fs::dir_delete(path),
    add = TRUE
  )
  fs::dir_create(path)

  annot_path <- self@annot_path
  # An annotation stored inside the directory being replaced would be lost
  annot_inside <- !is.null(annot_path) &&
    startsWith(
      as.character(fs::path_abs(annot_path)),
      paste0(final_path, .Platform$file.sep)
    )
  if ((isTRUE(copy_annotation) || annot_inside) && !is.null(annot_path)) {
    fs::dir_copy(annot_path, fs::path(path, "annotation"))
    annot_path_stored <- "annotation"
  } else {
    annot_path_stored <- annot_path %||% ""
  }
  write_parquet(
    list_to_kv(list(
      experiment_name = self@experiment_name,
      annot_path = annot_path_stored,
      isoformic_version = as.character(utils::packageVersion("isoformic")),
      written_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
    )),
    fs::path(path, "experiment.parquet")
  )
  if (!is.null(self@col_data)) {
    write_parquet(self@col_data, fs::path(path, "col_data.parquet"))
  }
  # Inferential replicates stay in memory: they are large and are re-read from
  # the Salmon output; the tables derived from them are written below.
  meta <- self@metadata
  meta[["inf_reps"]] <- NULL
  if (length(meta)) {
    write_parquet(list_to_kv(meta), fs::path(path, "metadata.parquet"))
  }
  pal <- self@tx_type_palette
  write_parquet(
    tibble::tibble(transcript_type = names(pal), color = unname(pal)),
    fs::path(path, "tx_type_palette.parquet")
  )
  write_assay_dir(self@assays, fs::path(path, "assays"), "transcript_id")
  write_assay_dir(
    attr(self, "gene_assays", exact = TRUE),
    fs::path(path, "gene_assays"),
    "gene_id"
  )
  dea <- self@dea
  for (nm in intersect(
    c("det", "deg", uncertainty_table_names()),
    names(dea)
  )) {
    write_parquet(dea[[nm]], fs::path(path, "dea", nm, ext = "parquet"))
  }
  if (!is.null(dea[["cutoffs"]])) {
    write_parquet(
      list_to_kv(dea[["cutoffs"]]),
      fs::path(path, "dea", "cutoffs.parquet")
    )
  }
  for (nm in names(self@gsea)) {
    write_parquet(
      flatten_list_columns(self@gsea[[nm]]),
      fs::path(path, "gsea", nm, ext = "parquet")
    )
  }
  if (fs::dir_exists(final_path)) {
    fs::dir_delete(final_path)
  }
  fs::file_move(path, final_path)
  self@data_path <- final_path
  invisible(self)
}

write_assay_dir <- function(assays, dir, id_col) {
  if (!length(assays)) {
    return(invisible(NULL))
  }
  fs::dir_create(dir)
  for (nm in names(assays)) {
    write_parquet(
      matrix_to_table(assays[[nm]], id_col),
      fs::path(dir, nm, ext = "parquet")
    )
  }
  invisible(NULL)
}

matrix_to_table <- function(m, id_col) {
  df <- as.data.frame(m, stringsAsFactors = FALSE, optional = TRUE)
  df <- cbind(
    stats::setNames(data.frame(rownames(m), stringsAsFactors = FALSE), id_col),
    df
  )
  rownames(df) <- NULL
  df
}

table_to_matrix <- function(df, id_col) {
  ids <- as.character(df[[id_col]])
  m <- as.matrix(df[, setdiff(colnames(df), id_col), drop = FALSE])
  storage.mode(m) <- "double"
  rownames(m) <- ids
  m
}

read_assay_dir <- function(dir, id_col) {
  if (!fs::dir_exists(dir)) {
    return(list())
  }
  files <- fs::dir_ls(dir, glob = "*.parquet")
  out <- lapply(files, function(f) table_to_matrix(read_parquet(f), id_col))
  names(out) <- fs::path_ext_remove(fs::path_file(files))
  out
}

flatten_list_columns <- function(df) {
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  is_list <- vapply(df, is.list, logical(1L))
  df[is_list] <- lapply(df[is_list], function(col) {
    vapply(
      col,
      function(v) paste(as.character(v), collapse = ";"),
      character(1L)
    )
  })
  df
}

#' @rdname isoformic-io
#' @export
read_isoformic <- function(path, annot_path = NULL) {
  assert_dir_exists(path)
  path <- as.character(fs::path_abs(path))
  exp_file <- fs::path(path, "experiment.parquet")
  if (!fs::file_exists(exp_file)) {
    cli::cli_abort(
      message = c(
        x = "{.path {path}} is not an IsoformicExperiment directory (no {.file experiment.parquet})."
      ),
      class = "isoformic_invalid_experiment_dir"
    )
  }
  info <- kv_to_list(read_parquet(exp_file))
  if (is.null(annot_path)) {
    stored <- info[["annot_path"]] %||% ""
    if (
      identical(stored, "annotation") ||
        is_annotation_dir(fs::path(path, "annotation"))
    ) {
      annot_path <- as.character(fs::path(path, "annotation"))
    } else if (nzchar(stored)) {
      annot_path <- stored
    }
  }
  col_data_file <- fs::path(path, "col_data.parquet")
  col_data <- if (fs::file_exists(col_data_file)) {
    read_parquet(col_data_file)
  } else {
    NULL
  }
  meta_file <- fs::path(path, "metadata.parquet")
  metadata <- if (fs::file_exists(meta_file)) {
    kv_to_list(read_parquet(meta_file))
  } else {
    list()
  }
  pal_file <- fs::path(path, "tx_type_palette.parquet")
  palette <- NULL
  if (fs::file_exists(pal_file)) {
    pal <- read_parquet(pal_file)
    palette <- stats::setNames(pal$color, pal$transcript_type)
  }
  dea <- list()
  for (nm in c("det", "deg", uncertainty_table_names())) {
    f <- fs::path(path, "dea", nm, ext = "parquet")
    if (fs::file_exists(f)) {
      dea[[nm]] <- read_parquet(f)
    }
  }
  cut_file <- fs::path(path, "dea", "cutoffs.parquet")
  if (fs::file_exists(cut_file)) {
    raw <- kv_to_list(read_parquet(cut_file))
    dea[["cutoffs"]] <- list(
      log2fc = as.numeric(raw$log2fc),
      pvalue = as.numeric(raw$pvalue),
      use_fdr = as.logical(raw$use_fdr)
    )
  }
  gsea <- list()
  gsea_dir <- fs::path(path, "gsea")
  if (fs::dir_exists(gsea_dir)) {
    files <- fs::dir_ls(gsea_dir, glob = "*.parquet")
    gsea <- lapply(files, read_parquet)
    names(gsea) <- fs::path_ext_remove(fs::path_file(files))
  }
  IsoformicExperiment(
    experiment_name = info[["experiment_name"]],
    annot_path = annot_path,
    assays = read_assay_dir(fs::path(path, "assays"), "transcript_id"),
    col_data = col_data,
    gene_assays = read_assay_dir(fs::path(path, "gene_assays"), "gene_id"),
    dea = dea,
    gsea = gsea,
    metadata = metadata,
    data_path = path,
    tx_type_palette = palette
  )
}

# Conversion ------------------------------------------------------------------

#' Convert Between IsoformicExperiment and SummarizedExperiment
#'
#' @description
#' `as_isoformic()` builds an [IsoformicExperiment] from a
#' `SummarizedExperiment` (for example the output of `tximeta`), from a
#' `tximport` list (`abundance`, `counts`, `length`) or from a named list of
#' matrices. `as_summarized_experiment()` does the reverse at the transcript
#' or gene level, which is convenient for running Bioconductor tools such as
#' `DESeq2` or `fishpond` on the stored assays (a `counts` assay is placed
#' first, as those tools expect).
#'
#' @param x A `SummarizedExperiment`, a `tximport` list or a named list of
#'   matrices.
#' @param annot_path Annotation directory or GFF3/GTF file.
#' @param col_data Sample metadata; taken from `colData(x)` when `x` is a
#'   `SummarizedExperiment`.
#' @param experiment_name Optional experiment name.
#' @param ignore_after_bar Logical, truncate GENCODE `|`-separated
#'   transcript IDs.
#' @param drop_unannotated Logical, drop transcripts absent from the
#'   annotation instead of failing.
#' @param self An [IsoformicExperiment].
#' @param level `"transcript"` (default) or `"gene"`.
#' @param ... Passed on to methods.
#'
#' @returns `as_isoformic()` returns an [IsoformicExperiment];
#'   `as_summarized_experiment()` returns a `SummarizedExperiment`.
#'
#' @examples
#' iso <- example_isoformic()
#' txi <- list(abundance = assay(iso, "abundance"), counts = assay(iso, "counts"))
#' annot <- system.file("extdata", "example_experiment", "annotation", package = "isoformic")
#' iso2 <- as_isoformic(txi, annot_path = annot, col_data = col_data(iso))
#' iso2
#'
#' if (requireNamespace("SummarizedExperiment", quietly = TRUE)) {
#'   se <- as_summarized_experiment(iso, level = "gene")
#'   se
#' }
#'
#' @name isoformic-conversion
#' @rdname isoformic-conversion
NULL

#' @rdname isoformic-conversion
#' @export
as_isoformic <- function(
  x,
  annot_path,
  col_data = NULL,
  experiment_name = NULL,
  ignore_after_bar = TRUE,
  drop_unannotated = FALSE
) {
  rlang::check_required(x)
  if (S7::S7_inherits(x, IsoformicExperiment)) {
    return(x)
  }
  rlang::check_required(annot_path)
  if (inherits(x, "SummarizedExperiment")) {
    rlang::check_installed("SummarizedExperiment")
    assays <- as.list(SummarizedExperiment::assays(x))
    assays <- lapply(assays, function(m) {
      m <- as.matrix(m)
      if (is.null(rownames(m))) {
        rownames(m) <- rownames(x)
      }
      if (is.null(colnames(m))) {
        colnames(m) <- colnames(x)
      }
      m
    })
    if (is.null(col_data)) {
      col_data <- as.data.frame(SummarizedExperiment::colData(x))
      if (!"sample_id" %in% colnames(col_data)) {
        col_data <- cbind(sample_id = colnames(x), col_data)
      }
      col_data$files <- NULL
    }
  } else if (is.list(x) && !is.data.frame(x)) {
    assays <- x[vapply(
      x,
      function(v) is.matrix(v) || is.data.frame(v),
      logical(1L)
    )]
    if (!length(assays)) {
      cli::cli_abort(
        message = c(
          x = "{.arg x} must contain at least one matrix (e.g. {.field abundance})."
        ),
        class = "isoformic_invalid_assays"
      )
    }
  } else {
    cli::cli_abort(
      message = c(
        x = "Cannot convert an object of class {.cls {class(x)}} to {.cls IsoformicExperiment}."
      ),
      class = "isoformic_convert_error"
    )
  }
  assays <- standardize_assays(assays)
  if (isTRUE(ignore_after_bar)) {
    assays <- lapply(assays, function(m) {
      rownames(m) <- sub("\\|.*$", "", rownames(m))
      m
    })
  }
  annot_dir <- resolve_annot_path(annot_path)
  if (isTRUE(drop_unannotated)) {
    assays <- drop_unannotated_rows(assays, annot_dir)
  }
  IsoformicExperiment(
    experiment_name = experiment_name,
    annot_path = annot_dir,
    assays = assays,
    col_data = col_data
  )
}

#' @rdname isoformic-conversion
#' @export
as_summarized_experiment <- S7::new_generic(
  "as_summarized_experiment",
  "self",
  function(self, level = c("transcript", "gene"), ...) {
    S7::S7_dispatch()
  }
)

S7::method(as_summarized_experiment, IsoformicExperiment) <- function(
  self,
  level = c("transcript", "gene"),
  ...
) {
  rlang::check_installed("SummarizedExperiment")
  level <- rlang::arg_match(level)
  if (identical(level, "transcript")) {
    assays <- self@assays
    row_df <- as.data.frame(row_data(self))
    rownames(row_df) <- row_df$transcript_id
  } else {
    assays <- self@gene_assays
    row_df <- as.data.frame(gene_data(self))
    rownames(row_df) <- row_df$gene_id
    # `strand` is a reserved column name in Bioconductor row ranges
    names(row_df)[names(row_df) == "strand"] <- "gene_strand"
  }
  col_df <- as.data.frame(self@col_data)
  rownames(col_df) <- col_df$sample_id
  col_df <- col_df[colnames(assays[[1L]]), , drop = FALSE]
  # Bioconductor tools (DESeq2) expect `counts` to be the first assay
  assays <- assays[order(names(assays) != "counts")]
  SummarizedExperiment::SummarizedExperiment(
    assays = assays,
    rowData = row_df,
    colData = col_df,
    metadata = c(
      list(
        experiment_name = self@experiment_name,
        annot_path = self@annot_path
      ),
      self@metadata
    )
  )
}

#' Example IsoformicExperiment
#'
#' Loads a small [IsoformicExperiment] bundled with the package: transcript
#' abundances (TPM) for the isoforms of a handful of genes from a
#' preeclampsia placenta RNA-seq study (25 samples), GENCODE v33
#' annotation for those genes, and transcript- and gene-level differential
#' expression results. It is used throughout the documentation and is
#' suitable for trying every method of the class in well under a second.
#'
#' @returns An [IsoformicExperiment].
#'
#' @examples
#' iso <- example_isoformic()
#' iso
#'
#' @export
example_isoformic <- function() {
  path <- system.file("extdata", "example_experiment", package = "isoformic")
  if (!nzchar(path)) {
    cli::cli_abort(
      message = c(
        x = "The example experiment is not available in this installation."
      ),
      class = "isoformic_missing_example"
    )
  }
  iso <- read_isoformic(path)
  iso@data_path <- NULL
  iso
}
