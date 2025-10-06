#' Run DuckDB queries
#' @keywords internal
#' @noRd
duckdb_run <- function(
    sql_string,
    db_type = c("duckdb_memory", "duckdb_tempfile", "duckdb_file"),
    db_file_path = NULL,
    read_only = FALSE,
    envir = parent.frame()) {
  db_type <- rlang::arg_match(db_type)
  if (!isTRUE(rlang::is_scalar_logical(read_only))) {
    cli::cli_abort(
      message = c(
        x = "{.arg read_only} must be a single {.cls logical} value."
      ),
      class = "isoformic_read_only_not_logical"
    )
  }
  if (isTRUE(rlang::is_null(db_file_path) || identical(db_file_path, ""))) {
    db_type <- "duckdb_memory"
    db_file_path <- "none"
  } else {
    db_type <- rlang::arg_match(db_type)
  }
  if (identical(db_type, "duckdb_memory")) {
    db_storage_str <- ":memory:"
  } else if (identical(db_type, "duckdb_tempfile")) {
    db_storage_str <- fs::file_temp(ext = "duckdb")
  } else if (
    identical(db_type, "duckdb_file") && !identical(db_file_path, "none")
  ) {
    db_storage_str <- fs::path(
      fs::path_ext_remove(db_file_path),
      ext = "duckdb"
    )
  } else {
    db_storage_str <- fs::path("_isoformic_annot_db", ext = "duckdb")
  }
  conn_obj <- withr::local_db_connection(
    DBI::dbConnect(
      drv = duckdb::duckdb(),
      dbdir = db_storage_str,
      read_only = TRUE
    )
  )
  sql_str <- glue::glue_sql(sql_string, .con = conn_obj, .envir = envir)
  DBI::dbExecute(conn_obj, sql_str)

  return(invisible(NULL))
}
