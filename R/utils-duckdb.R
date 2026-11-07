#' Open a DuckDB Connection That Closes Itself
#'
#' Opens an in-memory (default) or file-backed DuckDB connection and
#' registers its disconnection in the calling frame, so callers do not need
#' to manage the connection lifetime.
#'
#' @param db_file_path Optional path to a DuckDB database file.
#'   `NULL` (default) opens an in-memory database.
#' @param read_only Logical, open the database in read-only mode.
#' @param envir Environment whose exit triggers the disconnection.
#'
#' @returns A `DBIConnection` to DuckDB.
#'
#' @keywords internal
#' @noRd
duckdb_connect <- function(
  db_file_path = NULL,
  read_only = FALSE,
  envir = parent.frame()
) {
  if (is.null(db_file_path) || identical(db_file_path, "")) {
    db_dir <- ":memory:"
  } else {
    db_dir <- as.character(db_file_path)
  }
  con <- DBI::dbConnect(
    drv = duckdb::duckdb(config = duckdb_storage_config()),
    dbdir = db_dir,
    read_only = isTRUE(read_only)
  )
  withr::defer(DBI::dbDisconnect(con, shutdown = TRUE), envir = envir)
  DBI::dbExecute(con, "SET enable_progress_bar = false;")
  return(con)
}

#' Storage Locations for the DuckDB Connections
#'
#' isoformic uses no DuckDB extension or secret, so extensions, secrets and
#' spill files go to a per-session directory under [tempdir()]. Setting them
#' explicitly keeps duckdb from announcing its storage location on every new
#' session and from offering to create `~/.duckdb` in interactive sessions.
#' A `duckdb.home` option or `DUCKDB_R_HOME` variable set by the user wins.
#'
#' @returns A named list for the `config` argument of [duckdb::duckdb()].
#'
#' @keywords internal
#' @noRd
duckdb_storage_config <- function() {
  if (
    !is.null(getOption("duckdb.home")) || nzchar(Sys.getenv("DUCKDB_R_HOME"))
  ) {
    return(list())
  }
  root <- fs::path(tempdir(), "isoformic-duckdb")
  list(
    extension_directory = as.character(fs::path(root, "extensions")),
    secret_directory = as.character(fs::path(root, "stored_secrets")),
    temp_directory = as.character(fs::path(root, "tmp"))
  )
}

#' Run One or More DuckDB SQL Statements
#'
#' Interpolates `sql_string` with [glue::glue_sql()] (so paths and values
#' are quoted safely) and executes it in a short-lived DuckDB session.
#'
#' @param sql_string SQL text, possibly containing `{}` glue placeholders
#'   resolved from the named values in `...`.
#' @param ... Named values interpolated into `sql_string`.
#' @param db_file_path Optional DuckDB database file. `NULL` runs in-memory.
#' @param read_only Logical, open the database read-only.
#'
#' @returns `NULL`, invisibly.
#'
#' @keywords internal
#' @noRd
duckdb_run <- function(
  sql_string,
  ...,
  db_file_path = NULL,
  read_only = FALSE
) {
  if (!rlang::is_scalar_logical(read_only)) {
    cli::cli_abort(
      message = c(
        x = "{.arg read_only} must be a single {.cls logical} value."
      ),
      class = "isoformic_read_only_not_logical"
    )
  }
  con <- duckdb_connect(db_file_path = db_file_path, read_only = read_only)
  values <- rlang::list2(...)
  if (
    length(values) && (is.null(names(values)) || !all(nzchar(names(values))))
  ) {
    cli::cli_abort(
      message = c(
        x = "All interpolation values passed to {.fn duckdb_run} must be named."
      ),
      class = "isoformic_duckdb_run_unnamed_values"
    )
  }
  sql_str <- glue::glue_sql(
    sql_string,
    .con = con,
    .envir = list2env(values, parent = emptyenv()),
    .null = DBI::SQL("NULL"),
    .na = DBI::SQL("NULL")
  )
  DBI::dbExecute(con, sql_str)
  return(invisible(NULL))
}

#' Read and Write Parquet Files
#'
#' Thin, dependency-light wrappers around DuckDB for reading tabular files
#' into tibbles and for writing data frames as Parquet. These are the readers
#' used throughout `isoformic`: annotation tables, experiment caches and the
#' bundled example data are all stored as Parquet.
#'
#' `read_parquet()` reads a Parquet file. `read_tabular()` reads Parquet,
#' CSV or TSV files (`.csv`, `.tsv`, `.tab`, `.txt`, Salmon `.sf`; optionally
#' gzip-compressed), auto-detecting delimiter and types. `write_parquet()` writes any data frame to a
#' Parquet file, creating parent directories as needed.
#'
#' @param path Path to the file to read or write.
#' @param x A data frame (or tibble) to write.
#' @param columns Optional character vector of columns to read.
#'   `NULL` reads all columns.
#' @param n_max Maximum number of rows to read. `Inf` (default) reads all.
#' @param compression Parquet compression codec.
#' @param ... Reserved for future use; must be empty.
#'
#' @returns `read_parquet()` and `read_tabular()` return a tibble.
#'   `write_parquet()` returns `path`, invisibly.
#'
#' @details
#' Factor columns are stored as Parquet strings and come back as character.
#' Row names are not stored; keep identifiers in a regular column.
#'
#' @examples
#' tmp <- tempfile(fileext = ".parquet")
#' write_parquet(data.frame(id = c("a", "b"), value = c(1.5, 2)), tmp)
#' read_parquet(tmp)
#' read_parquet(tmp, columns = "id", n_max = 1)
#' read_tabular(tmp)
#'
#' @name parquet
#' @export
read_parquet <- function(path, columns = NULL, n_max = Inf, ...) {
  rlang::check_dots_empty()
  assert_file_exists(path)
  path <- as.character(fs::path_abs(path))
  con <- duckdb_connect()
  sql <- glue::glue_sql(
    "SELECT {cols*} FROM read_parquet({path}){limit_sql}",
    cols = select_cols_sql(con, columns),
    limit_sql = limit_sql(con, n_max),
    .con = con
  )
  res <- DBI::dbGetQuery(con, sql)
  return(tibble::as_tibble(res))
}

#' @rdname parquet
#' @export
read_tabular <- function(path, columns = NULL, n_max = Inf, ...) {
  rlang::check_dots_empty()
  assert_file_exists(path)
  path <- as.character(fs::path_abs(path))
  ext <- tolower(fs::path_ext(sub("\\.(gz|zst|bz2)$", "", path)))
  con <- duckdb_connect()
  reader <- switch(
    ext,
    parquet = "read_parquet",
    csv = "read_csv",
    tsv = ,
    tab = ,
    sf = "read_csv",
    txt = "read_csv",
    cli::cli_abort(
      message = c(
        x = "Unsupported file extension {.val {ext}} in {.path {path}}.",
        i = "Supported: parquet, csv, tsv, tab, txt, sf (optionally gzip-compressed)."
      ),
      class = "isoformic_unsupported_file_type"
    )
  )
  delim_sql <- switch(
    ext,
    csv = ", delim = ','",
    tsv = ,
    tab = ,
    sf = ", delim = '\t'",
    ""
  )
  sql <- glue::glue_sql(
    "SELECT {cols*} FROM {DBI::SQL(reader)}({path}{DBI::SQL(delim_sql)}){limit_sql}",
    cols = select_cols_sql(con, columns),
    limit_sql = limit_sql(con, n_max),
    .con = con
  )
  res <- DBI::dbGetQuery(con, sql)
  return(tibble::as_tibble(res))
}

#' @rdname parquet
#' @export
write_parquet <- function(
  x,
  path,
  compression = c("zstd", "snappy", "gzip", "uncompressed"),
  ...
) {
  rlang::check_dots_empty()
  compression <- rlang::arg_match(compression)
  check_string(path)
  if (!is.data.frame(x)) {
    cli::cli_abort(
      message = c(x = "{.arg x} must be a data frame."),
      class = "isoformic_not_data_frame"
    )
  }
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  # Drop a hidden rownames column: identifiers must be regular columns
  rownames(x) <- NULL
  x[] <- lapply(x, function(col) if (is.factor(col)) as.character(col) else col)
  fs::dir_create(fs::path_dir(path))
  path <- as.character(fs::path_abs(path))
  con <- duckdb_connect()
  duckdb::duckdb_register(con, "isoformic_write_tbl", x)
  sql <- glue::glue_sql(
    "COPY isoformic_write_tbl TO {path} (FORMAT parquet, COMPRESSION {DBI::SQL(compression)})",
    .con = con
  )
  DBI::dbExecute(con, sql)
  return(invisible(path))
}

#' Number of Rows in a Parquet File Without Reading It
#' @keywords internal
#' @noRd
parquet_nrow <- function(path) {
  assert_file_exists(path)
  path <- as.character(fs::path_abs(path))
  con <- duckdb_connect()
  sql <- glue::glue_sql(
    "SELECT count(*) AS n FROM read_parquet({path})",
    .con = con
  )
  return(as.integer(DBI::dbGetQuery(con, sql)$n))
}

select_cols_sql <- function(con, columns) {
  if (is.null(columns)) {
    return(DBI::SQL("*"))
  }
  return(DBI::dbQuoteIdentifier(con, as.character(columns)))
}

limit_sql <- function(con, n_max) {
  if (is.infinite(n_max)) {
    return(DBI::SQL(""))
  }
  return(glue::glue_sql(" LIMIT {as.integer(n_max)}", .con = con))
}

#' Read a delimited text table with DuckDB, with explicit parsing options
#'
#' For source files of the frozen datasets, whose quoting and missing-value
#' conventions differ from the defaults of [read_tabular()]. Text is never
#' unquoted (`quote = ""`), the column names are kept as written, and the types
#' are inferred from the whole file.
#' @keywords internal
#' @noRd
read_delimited <- function(
  path,
  delim = "\t",
  null_strings = c("", "NA"),
  envir = parent.frame()
) {
  assert_file_exists(path)
  path <- as.character(fs::path_abs(path))
  con <- duckdb_connect(envir = envir)
  nulls <- paste0(
    "[",
    paste0(
      "'",
      gsub("'", "''", null_strings, fixed = TRUE),
      "'",
      collapse = ", "
    ),
    "]"
  )
  sql <- glue::glue_sql(
    "SELECT * FROM read_csv({path}, delim = {delim}, quote = '', header = true, ",
    "nullstr = {DBI::SQL(nulls)}, sample_size = -1)",
    .con = con
  )
  tibble::as_tibble(DBI::dbGetQuery(con, sql))
}
