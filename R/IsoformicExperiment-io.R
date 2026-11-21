# Persistence and conversion --------------------------------------------------

#' Save and Restore IsoformicExperiment Objects as Parquet
#'
#' @description
#' `write_isoformic()` writes an [IsoformicExperiment] to a directory of
#' Parquet files and `read_isoformic()` restores it. The layout is:
#'
#' ```
#' <path>/
#'   experiment.parquet        key/value: experiment_name, annot_path,
#'                             isoformic_version, written_at
#'   col_data.parquet          sample metadata
#'   metadata.parquet          key/value: scalar entries of `metadata`
#'   tx_type_palette.parquet   transcript_type, color
#'   assays/<name>.parquet     transcript_id + one column per sample
#'   gene_assays/<name>.parquet  only explicitly set gene assays
#'   dea/det.parquet, dea/deg.parquet, dea/cutoffs.parquet
#'   dea/uncertainty.parquet, dea/ambiguity.parquet, dea/confidence.parquet
#'   gsea/<name>.parquet       list columns joined with ";"
#'   inf_reps/draw_0001.parquet, ...  one draw: transcript_id + samples
#'   inf_reps/moments/mean.parquet, var.parquet  replicate moments
#'   annotation/               copy of the annotation directory (optional)
#' ```
#'
#' Each file is written only when the object holds that part. Entries of
#' `metadata` that are not single values are dropped with a warning.
#'
#' The inferential replicates (`self@metadata$inf_reps`, see
#' [read_inferential_replicates()] and [add_uncertainty()]) are written to
#' `inf_reps/`: one Parquet file per draw (`draw_0001.parquet`, ...) with
#' `transcript_id` and one column per sample, and, when the replicates have
#' a `moments` attribute, its `mean` and `var` matrices in
#' `inf_reps/moments/`. Their number, their type (the `samp_type`
#' attribute) and the number of draws behind the moments are recorded in
#' `metadata.parquet` under the keys `inf_reps_n_draws`,
#' `inf_reps_samp_type` and `inf_reps_moments_n`; when replicates are
#' written, these keys replace entries of `metadata` with those names, with
#' a warning. One draw is written, and read, at a time, so memory holds the
#' replicates and a single draw on top of them. The replicates of a subset
#' (`x[i, j]` shares the draws of `x` and records a `selection`, see
#' [IsoformicExperiment]) are written as the selected transcripts and
#' samples only, and the moments as the selected samples. `read_isoformic()`
#' restores the list of matrices with its `samp_type` and `moments`
#' attributes, as they were (including the storage mode of the matrices and
#' of `n`); for a subset, as plain matrices of the selected draws, without
#' the `selection`. Only
#' these are kept: names of the draws, other attributes of the list or of
#' the matrices and names of their dimnames are dropped, with a warning.
#' Replicates that are not numeric matrices with transcript and sample
#' names are not written, with a warning. Directories written by earlier
#' versions have no `inf_reps/` and are read without replicates. Use
#' `add_uncertainty(keep_draws = FALSE)`, or remove `metadata$inf_reps`, to
#' save an object without its replicates.
#'
#' With `copy_annotation = TRUE` the annotation directory is copied inside
#' the experiment directory so that the result is self-contained and
#' portable. The cached ORF table (`annot_data_orfs.parquet`) is left out
#' of the copy and computed again when needed. An annotation that lies
#' inside the directory being replaced is always copied, and so is an
#' annotation directory under [tempdir()] (such as the private copy made by
#' [example_isoformic()]) when `path` is not under [tempdir()] itself: the
#' temporary directory is removed when the R session ends. A message says
#' so. Otherwise
#' `annot_path` is stored as an absolute path and must still exist when the
#' object is read back (it can be overridden with the `annot_path` argument
#' of `read_isoformic()`). Without that argument, `read_isoformic()` uses
#' an `annotation/` directory inside `path` whenever there is one.
#'
#' The new directory is built next to `path` and moved in place at the end,
#' so a failed write leaves an existing experiment as it was. With
#' `overwrite = TRUE`, a non-empty `path` is replaced only when it holds an
#' `experiment.parquet` file.
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
#' @param ... Not used; reserved for methods.
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
#' identical(tx_assay(iso), tx_assay(iso2))
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

  annot_path_stored <- store_annotation(
    self@annot_path,
    final_path,
    path,
    copy_annotation
  )
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
  # The inferential replicates go to inf_reps/, one file per draw; their
  # type, number and moments count go to metadata.parquet
  meta <- self@metadata
  # Read as the draws of the object: the selection of a subset applies
  # even when list operations dropped its attribute
  reps <- object_inf_reps(self, meta[["inf_reps"]])
  meta[["inf_reps"]] <- NULL
  reps_meta <- write_inf_reps(reps, fs::path(path, "inf_reps"))
  rm(reps)
  if (length(reps_meta)) {
    replaced <- intersect(names(meta), inf_reps_keys())
    if (length(replaced)) {
      cli::cli_warn(
        message = c(
          "!" = paste0(
            "The {.code metadata} {cli::qty(length(replaced))}entr{?y/ies} ",
            "{.val {replaced}} {cli::qty(length(replaced))}",
            "{?is/are} replaced by the description of the inferential ",
            "replicates."
          )
        ),
        class = "isoformic_metadata_replaced"
      )
    }
    meta[inf_reps_keys()] <- NULL
    meta <- c(meta, reps_meta)
  }
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
    prop_raw(self, "gene_assays"),
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

#' Copy the annotation into an experiment directory when it must be
#'
#' The annotation is copied with `copy_annotation = TRUE`, and always when it
#' lies inside the directory being replaced (it would be lost). Returns the
#' annotation path to record: `"annotation"` for a copy, the path otherwise.
#' @keywords internal
#' @noRd
store_annotation <- function(annot_path, final_path, path, copy_annotation) {
  if (is.null(annot_path)) {
    return("")
  }
  annot_inside <- startsWith(
    as.character(fs::path_abs(annot_path)),
    paste0(final_path, .Platform$file.sep)
  )
  # An annotation in the session temporary directory is removed when the
  # session ends: copy it unless the experiment is written there as well
  annot_temporary <- !isTRUE(copy_annotation) &&
    !annot_inside &&
    fs::dir_exists(annot_path) &&
    in_session_tempdir(annot_path) &&
    !in_session_tempdir(final_path)
  if (!isTRUE(copy_annotation) && !annot_inside && !annot_temporary) {
    return(annot_path)
  }
  if (annot_temporary) {
    cli::cli_inform(
      message = c(
        i = paste0(
          "The annotation directory is in the temporary directory of this ",
          "R session; it is copied into the experiment directory so that ",
          "the experiment can be read in a later session."
        )
      ),
      class = "isoformic_annotation_copied"
    )
  }
  fs::dir_copy(annot_path, fs::path(path, "annotation"))
  # The stored ORF table is not valid for the copy (it is keyed by the file
  # times of the annotation): it is computed again when it is needed
  copied_orfs <- fs::path(path, "annotation", "annot_data_orfs.parquet")
  if (fs::file_exists(copied_orfs)) {
    fs::file_delete(copied_orfs)
  }
  "annotation"
}

#' Is a path under the temporary directory of the R session?
#'
#' [tempdir()] is removed when the R session ends, so an annotation stored
#' there cannot be referenced by an experiment written elsewhere. The path
#' need not exist: symbolic links are resolved on its deepest existing
#' parent.
#' @keywords internal
#' @noRd
in_session_tempdir <- function(path) {
  tmp <- as.character(fs::path_real(session_tempdir()))
  path <- as.character(fs::path_abs(path))
  rest <- character()
  while (!fs::file_exists(path) && !identical(fs::path_dir(path), path)) {
    rest <- c(fs::path_file(path), rest)
    path <- as.character(fs::path_dir(path))
  }
  resolved <- as.character(fs::path_real(path))
  if (length(rest) > 0L) {
    resolved <- as.character(do.call(fs::path, as.list(c(resolved, rest))))
  }
  startsWith(resolved, paste0(tmp, .Platform$file.sep))
}

#' The session temporary directory (a function so tests can replace it)
#' @keywords internal
#' @noRd
session_tempdir <- function() {
  tempdir()
}

write_assay_dir <- function(assays, dir, id_col, collect = gc_pacer()) {
  if (!length(assays)) {
    return(invisible(NULL))
  }
  fs::dir_create(dir)
  for (nm in names(assays)) {
    m <- assays[[nm]]
    write_parquet(
      matrix_to_table(m, id_col),
      fs::path(dir, nm, ext = "parquet")
    )
    # Free the tables of the matrices written so far once they add up: the
    # tables of many matrices (assays, or the draws of the inferential
    # replicates) otherwise pile up before R collects them, about one copy
    # of the draws at 295,000 transcripts x 12 samples x 25 draws
    collect(matrix_cells(m))
    rm(m)
  }
  invisible(NULL)
}

#' Number of cells of a matrix and of its row names
#' @keywords internal
#' @noRd
matrix_cells <- function(m) {
  d <- dim(m)
  if (length(d) != 2L) {
    return(as.double(length(m)))
  }
  as.double(d[[1L]]) * (as.double(d[[2L]]) + 1)
}

#' Run a full garbage collection once enough data has gone through
#'
#' Returns a function that takes the number of cells just written or read
#' and calls `gc()` when the cells since the last collection reach
#' `threshold` (2 million cells, about 16 MB of doubles). A large matrix is
#' followed by a collection, as before, so that the tables of large draws
#' do not pile up, while small matrices, such as hundreds of small draws,
#' do not pay for a full collection each (about 0.1 s).
#' @keywords internal
#' @noRd
gc_pacer <- function(threshold = 2e6) {
  state <- new.env(parent = emptyenv())
  state$pending <- 0
  function(cells) {
    state$pending <- state$pending + cells
    if (state$pending >= threshold) {
      state$pending <- 0
      full_gc()
    }
    invisible(NULL)
  }
}

#' A full garbage collection (a function so tests can count the calls)
#' @keywords internal
#' @noRd
full_gc <- function() {
  invisible(gc())
}

matrix_to_table <- function(m, id_col) {
  # A matrix without rows has no row names: R drops empty dimnames
  if (
    !is.matrix(m) ||
      (is.null(rownames(m)) && nrow(m) > 0L) ||
      is.null(colnames(m))
  ) {
    df <- as.data.frame(m, stringsAsFactors = FALSE, optional = TRUE)
    df <- cbind(
      stats::setNames(
        data.frame(rownames(m), stringsAsFactors = FALSE),
        id_col
      ),
      df
    )
    rownames(df) <- NULL
    return(df)
  }
  # One vector per column and nothing else: `as.data.frame()` allocates four
  # per column. The offsets are doubles, so they do not overflow.
  n <- nrow(m)
  rows <- seq_len(n)
  cols <- lapply(seq_len(ncol(m)), function(j) m[(j - 1) * n + rows])
  names(cols) <- colnames(m)
  df <- c(stats::setNames(list(rownames(m) %||% character()), id_col), cols)
  attr(df, "row.names") <- .set_row_names(n)
  class(df) <- "data.frame"
  df
}

table_to_matrix <- function(df, id_col, keep_type = FALSE) {
  ids <- as.character(df[[id_col]])
  values <- df[, setdiff(colnames(df), id_col), drop = FALSE]
  # Integer columns stay integer with `keep_type`; a table without rows
  # would otherwise give a logical matrix
  integer_cols <- isTRUE(keep_type) &&
    length(values) > 0L &&
    all(vapply(values, is.integer, logical(1L)))
  m <- as.matrix(values)
  rm(values)
  storage.mode(m) <- if (integer_cols) "integer" else "double"
  rownames(m) <- ids
  m
}

read_assay_dir <- function(dir, id_col, keep_type = FALSE) {
  if (!fs::dir_exists(dir)) {
    return(list())
  }
  files <- fs::dir_ls(dir, glob = "*.parquet")
  out <- lapply(files, function(f) {
    table_to_matrix(read_parquet(f), id_col, keep_type = keep_type)
  })
  names(out) <- fs::path_ext_remove(fs::path_file(files))
  out
}

#' Keys of metadata.parquet that describe the stored inferential replicates
#' @keywords internal
#' @noRd
inf_reps_keys <- function() {
  c("inf_reps_n_draws", "inf_reps_samp_type", "inf_reps_moments_n")
}

#' File name of draw `k` of `n` (zero-padded so that the files sort in order)
#' @keywords internal
#' @noRd
inf_reps_file <- function(k, n) {
  width <- max(4L, nchar(as.character(n)))
  sprintf(paste0("draw_%0", width, "d.parquet"), k)
}

#' Can the inferential replicates be written as wide Parquet files?
#'
#' Every draw must be a numeric matrix with sample names (and transcript
#' names when it has rows); the moments, when present, a list with the
#' `mean` and `var` matrices, named the same way, and a single `n`.
#' @keywords internal
#' @noRd
inf_reps_writable <- function(reps) {
  named_matrix <- function(m) {
    is.matrix(m) &&
      is.numeric(m) &&
      !is.null(colnames(m)) &&
      (nrow(m) == 0L || !is.null(rownames(m)))
  }
  if (!is.list(reps) || !all(vapply(reps, named_matrix, logical(1L)))) {
    return(FALSE)
  }
  samp_type <- attr(reps, "samp_type", exact = TRUE)
  if (
    !is.null(samp_type) && !(is.character(samp_type) && length(samp_type) == 1L)
  ) {
    return(FALSE)
  }
  mom <- attr(reps, "moments", exact = TRUE)
  is.null(mom) ||
    (is.list(mom) &&
      identical(names(mom), c("mean", "var", "n")) &&
      named_matrix(mom$mean) &&
      named_matrix(mom$var) &&
      is.numeric(mom$n) &&
      length(mom$n) == 1L)
}

#' Write the inferential replicates as one wide Parquet file per draw
#'
#' `dir/draw_0001.parquet`, ... hold `transcript_id` and one column per
#' sample; `dir/moments/mean.parquet` and `dir/moments/var.parquet` the
#' moments, when the replicates have them. One draw is converted at a time,
#' and the long table of all the draws is never built. Returns the entries
#' to record in metadata.parquet (none when there is nothing to write).
#' @keywords internal
#' @noRd
write_inf_reps <- function(reps, dir) {
  if (is.null(reps)) {
    return(list())
  }
  if (!inf_reps_writable(reps)) {
    cli::cli_warn(
      message = c(
        "!" = "The inferential replicates are not written.",
        i = paste0(
          "{.code metadata$inf_reps} must be a list of numeric matrices ",
          "with transcript and sample names, as returned by ",
          "{.fn read_inferential_replicates}."
        )
      ),
      class = "isoformic_inf_reps_not_written"
    )
    return(list())
  }
  lost <- inf_reps_lost(reps)
  if (length(lost)) {
    cli::cli_warn(
      message = c(
        "!" = "Some properties of the inferential replicates are not written.",
        i = "Not kept: {lost}.",
        i = paste0(
          "Only the draws (values and dimnames), the {.field samp_type} ",
          "attribute and the {.field moments} attribute are written."
        )
      ),
      class = "isoformic_inf_reps_changed"
    )
  }
  fs::dir_create(dir)
  n <- length(reps)
  # One draw at a time, with the selection of a subset applied (a copy of
  # the selected rows and columns of that draw only; none without one)
  read_draw <- inf_reps_reader(reps)
  files <- fs::path_ext_remove(inf_reps_file(seq_len(n), n))
  collect <- gc_pacer()
  for (k in seq_len(n)) {
    write_assay_dir(
      stats::setNames(list(read_draw(k)), files[[k]]),
      dir,
      "transcript_id",
      collect = collect
    )
  }
  out <- list(inf_reps_n_draws = n)
  samp_type <- attr(reps, "samp_type", exact = TRUE)
  if (!is.null(samp_type)) {
    out$inf_reps_samp_type <- samp_type
  }
  mom <- inf_reps_moments(reps)
  if (!is.null(mom)) {
    write_assay_dir(
      list(mean = mom$mean, var = mom$var),
      fs::path(dir, "moments"),
      "transcript_id"
    )
    out$inf_reps_moments_n <- encode_moments_n(mom$n)
  }
  out
}

#' Properties of the inferential replicates that write_inf_reps() drops
#'
#' Returns a description of each one present: names of the draws, other
#' attributes of the list, attributes of the draws or the moments matrices
#' other than `dim` and `dimnames`, and names of their dimnames.
#' @keywords internal
#' @noRd
inf_reps_lost <- function(reps) {
  lost <- character()
  if (!is.null(names(reps))) {
    lost <- c(lost, "the names of the draws")
  }
  # The selection of a subset is applied when the draws are written
  extra <- setdiff(
    names(attributes(reps)),
    c("names", "samp_type", "moments", "selection")
  )
  if (length(extra)) {
    lost <- c(lost, paste0("the list attributes ", toString(extra)))
  }
  mom <- attr(reps, "moments", exact = TRUE)
  mats <- c(as.list(reps), if (!is.null(mom)) list(mom$mean, mom$var))
  matrix_extra <- function(m) {
    length(setdiff(names(attributes(m)), c("dim", "dimnames"))) > 0L ||
      !is.null(names(dimnames(m)))
  }
  if (any(vapply(mats, matrix_extra, logical(1L)))) {
    lost <- c(
      lost,
      "attributes of the matrices other than dim and dimnames (or names of their dimnames)"
    )
  }
  if (!is.null(mom) && length(setdiff(names(attributes(mom)), "names"))) {
    lost <- c(lost, "attributes of the moments list")
  }
  lost
}

#' Encode and decode the number of draws behind the moments
#'
#' metadata.parquet stores text; a double `n` gets a decimal point (or is
#' written in full precision) so that it reads back as a double, an integer
#' `n` reads back as an integer.
#' @keywords internal
#' @noRd
encode_moments_n <- function(n) {
  if (is.integer(n)) {
    return(as.character(n))
  }
  out <- sprintf("%.17g", as.double(n))
  if (grepl("^-?[0-9]+$", out)) {
    out <- paste0(out, ".0")
  }
  out
}

decode_moments_n <- function(x) {
  x <- as.character(x)
  if (is.na(x) || grepl("^-?[0-9]+$", x)) {
    return(as.integer(x))
  }
  suppressWarnings(as.double(x))
}

#' Read the inferential replicates written by `write_inf_reps()`
#'
#' One draw is read at a time. Draws with the same dimnames as the previous
#' one share its dimnames, as the draws of `read_inferential_replicates()`
#' do.
#' @keywords internal
#' @noRd
read_inf_reps <- function(dir, metadata) {
  n <- suppressWarnings(as.integer(metadata[["inf_reps_n_draws"]]))
  if (length(n) != 1L || is.na(n) || n < 0L) {
    cli::cli_abort(
      message = c(
        x = paste0(
          "{.path {dir}} holds inferential replicates but ",
          "{.file metadata.parquet} has no valid {.field inf_reps_n_draws}."
        )
      ),
      class = "isoformic_invalid_experiment_dir"
    )
  }
  reps <- vector("list", n)
  last_dn <- NULL
  collect <- gc_pacer()
  for (k in seq_len(n)) {
    f <- fs::path(dir, inf_reps_file(k, n))
    if (!fs::file_exists(f)) {
      cli::cli_abort(
        message = c(x = "Inferential replicate file {.path {f}} is missing."),
        class = "isoformic_invalid_experiment_dir"
      )
    }
    m <- table_to_matrix(read_parquet(f), "transcript_id", keep_type = TRUE)
    if (identical(dimnames(m), last_dn)) {
      dimnames(m) <- last_dn
    } else {
      last_dn <- dimnames(m)
    }
    reps[[k]] <- m
    collect(matrix_cells(m))
    rm(m)
  }
  samp_type <- metadata[["inf_reps_samp_type"]]
  if (!is.null(samp_type)) {
    attr(reps, "samp_type") <- samp_type
  }
  mom_n <- metadata[["inf_reps_moments_n"]]
  if (!is.null(mom_n)) {
    mom <- read_assay_dir(
      fs::path(dir, "moments"),
      "transcript_id",
      keep_type = TRUE
    )
    if (!all(c("mean", "var") %in% names(mom))) {
      cli::cli_abort(
        message = c(
          x = "The moments of the inferential replicates in {.path {dir}} are incomplete."
        ),
        class = "isoformic_invalid_experiment_dir"
      )
    }
    attr(reps, "moments") <- list(
      mean = mom$mean,
      var = mom$var,
      n = decode_moments_n(mom_n)
    )
  }
  reps
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
  # Directories written before the replicates were saved have no inf_reps/
  reps_dir <- fs::path(path, "inf_reps")
  if (fs::dir_exists(reps_dir)) {
    reps <- read_inf_reps(reps_dir, metadata)
    metadata[inf_reps_keys()] <- NULL
    metadata["inf_reps"] <- list(reps)
    rm(reps)
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
  # Fold-change tables written by earlier development versions carry the
  # removed confidence class: derive the current columns; nothing is written
  conf <- dea[["confidence"]]
  if (!is.null(conf) && !"unstable" %in% colnames(conf)) {
    min_agreement <- suppressWarnings(
      as.numeric(metadata[["uncertainty_min_agreement"]])
    )
    if (length(min_agreement) != 1L || is.na(min_agreement)) {
      min_agreement <- uncertainty_defaults()$min_agreement
    }
    dea[["confidence"]] <- upgrade_confidence_table(conf, min_agreement)
    cli::cli_inform(
      message = c(
        i = paste0(
          "The stored fold-change table was written by an earlier ",
          "development version; its confidence class was removed and ",
          "{.field unstable} derived from {.field sign_agreement}."
        )
      ),
      class = "isoformic_legacy_confidence"
    )
  }
  rm(conf)
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
#' first, as those tools expect). The `metadata` of the result holds
#' `experiment_name`, `annot_path` and the object's `metadata` without the
#' inferential replicates.
#'
#' @param x A `SummarizedExperiment`, a `tximport` list or a named list of
#'   matrices. An [IsoformicExperiment] is returned as it is.
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
#' @param ... Not used; reserved for methods.
#'
#' @returns `as_isoformic()` returns an [IsoformicExperiment];
#'   `as_summarized_experiment()` returns a `SummarizedExperiment`.
#'
#' @examples
#' iso <- example_isoformic()
#' txi <- list(
#'   abundance = tx_assay(iso, "abundance"),
#'   counts = tx_assay(iso, "counts")
#' )
#' annot <- system.file("extdata", "example_experiment", "annotation", package = "isoformic")
#' iso2 <- as_isoformic(txi, annot_path = annot, col_data = col_data(iso))
#' iso2
#'
#' # Loading the SummarizedExperiment namespace takes several seconds
#' if (interactive() && requireNamespace("SummarizedExperiment", quietly = TRUE)) {
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
  if (isTRUE(ignore_after_bar) && length(assays)) {
    # The assays share their row names; they are renamed, which copies them,
    # only when an ID has a bar
    ids <- rownames(assays[[1L]])
    if (any(grepl("|", ids, fixed = TRUE))) {
      ids <- sub("\\|.*$", "", ids)
      assays <- lapply(assays, function(m) {
        rownames(m) <- ids
        m
      })
    }
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
  # The inferential replicates stay in the object: in the metadata of the
  # result they would be written with every saved or serialised copy of it
  meta <- self@metadata
  meta[["inf_reps"]] <- NULL
  SummarizedExperiment::SummarizedExperiment(
    assays = assays,
    rowData = row_df,
    colData = col_df,
    metadata = c(
      list(
        experiment_name = self@experiment_name,
        annot_path = self@annot_path
      ),
      meta
    )
  )
}

#' Example IsoformicExperiment
#'
#' Loads a small [IsoformicExperiment] bundled with the package: transcript
#' abundances (TPM) and demonstration counts (TPM x 20, rounded) for the
#' isoforms of 69 genes from a preeclampsia placenta RNA-seq study (25
#' samples), GENCODE v50 annotation and transcript sequences for those
#' genes, and transcript- and gene-level differential expression results.
#' It is used throughout the documentation and is suitable for trying every
#' method of the class in well under a second.
#'
#' @details
#' The annotation directory of the bundled experiment is copied on every call
#' into a new directory under [tempdir()], and the returned object points
#' there (`iso@annot_path`). Functions that store tables in the
#' annotation directory, such as [extract_transcript_sequences()],
#' [transcript_orfs()], [add_appris()], [add_riboseq_orfs()] or
#' [add_nmd_exons()], therefore change only that private copy: the installed
#' package is never modified, two calls never see each other's changes, and
#' the copy is removed with the temporary directory when the R session ends.
#' [write_isoformic()] copies such a temporary annotation into an experiment
#' directory outside [tempdir()], so an object saved there can be read again
#' in a later session.
#'
#' @returns An [IsoformicExperiment] whose annotation directory is a private
#'   temporary copy of the bundled one.
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
  iso <- read_isoformic(
    path,
    annot_path = private_annotation_copy(fs::path(path, "annotation"))
  )
  iso@data_path <- NULL
  iso
}

#' Copy a bundled annotation directory to a new writable temporary directory
#'
#' Bundled data lives inside the installed package, which can be read-only
#' and must never be modified: the copy gives the caller a directory that
#' functions storing tables in the annotation can write to.
#' @keywords internal
#' @noRd
private_annotation_copy <- function(source) {
  dest <- fs::path(
    fs::file_temp(pattern = "isoformic_example_", tmp_dir = tempdir()),
    "annotation"
  )
  fs::dir_copy(source, dest)
  # Installed files can be read-only; the new directories are not
  files <- fs::dir_ls(dest, recurse = TRUE, type = "file")
  if (length(files) > 0L) {
    fs::file_chmod(files, "u+w")
  }
  as.character(fs::path_real(dest))
}
