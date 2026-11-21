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
  # Built before dbConnect(), so an invalid option is reported as is and not
  # wrapped in an S4 dispatch error
  config <- c(duckdb_storage_config(), duckdb_limits_config())
  con <- DBI::dbConnect(
    drv = duckdb::duckdb(config = config),
    dbdir = db_dir,
    read_only = isTRUE(read_only)
  )
  withr::defer(DBI::dbDisconnect(con, shutdown = TRUE), envir = envir)
  DBI::dbExecute(con, "SET enable_progress_bar = false;")
  return(con)
}

#' Query Data Frames Registered in a DuckDB Connection
#'
#' Registers each data frame of `tables` under its name and runs `sql`. With
#' `con = NULL` a connection is opened for the call and closed on exit;
#' otherwise the tables are unregistered on exit, so that one connection
#' serves several calls (opening one takes tens of milliseconds).
#' @keywords internal
#' @noRd
query_registered <- function(con, tables, sql) {
  if (is.null(con)) {
    con <- duckdb_connect()
  } else {
    on.exit(
      for (name in names(tables)) {
        duckdb::duckdb_unregister(con, name)
      },
      add = TRUE
    )
  }
  for (name in names(tables)) {
    duckdb::duckdb_register(con, name, tables[[name]])
  }
  DBI::dbGetQuery(con, sql)
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
  # DuckDB creates only the last component of these paths, so a query that
  # spills fails unless the root exists
  if (!fs::dir_exists(root)) {
    fs::dir_create(root)
  }
  list(
    extension_directory = as.character(fs::path(root, "extensions")),
    secret_directory = as.character(fs::path(root, "stored_secrets")),
    temp_directory = as.character(fs::path(root, "tmp"))
  )
}

#' Thread and Memory Limits for the DuckDB Connections
#'
#' DuckDB 1.5 sizes its thread pool from the CPUs of the machine and the
#' cgroup CPU quota, and its `memory_limit` as 80% of the RAM or of the
#' cgroup memory limit. It ignores the CPU affinity mask (`taskset`, Slurm
#' and other schedulers that pin cores without a cgroup quota) and the Slurm
#' memory allocation when Slurm does not enforce it with a cgroup.
#'
#' `threads`: the `isoformic.duckdb_threads` option; else, when
#' `SLURM_CPUS_PER_TASK` is set or the affinity mask has fewer CPUs than the
#' machine, the smallest of those and the cgroup v2 CPU quota; else unset.
#' `memory_limit`: the `isoformic.duckdb_memory_limit` option; else 80% of
#' `SLURM_MEM_PER_NODE` or of `SLURM_MEM_PER_CPU` times the CPUs of the task,
#' or of the cgroup v2 memory limit when it is smaller (at least 64 MiB);
#' else unset. Unset settings keep the DuckDB defaults.
#'
#' The config applies when a database is opened: a second connection to a
#' database file that is already open in the session shares the open
#' instance and keeps its settings.
#'
#' @returns A named list for the `config` argument of [duckdb::duckdb()].
#'
#' @keywords internal
#' @noRd
duckdb_limits_config <- function() {
  threads <- check_duckdb_threads(getOption("isoformic.duckdb_threads"))
  threads <- threads %||% detect_duckdb_threads()
  memory <- check_duckdb_memory_limit(
    getOption("isoformic.duckdb_memory_limit")
  )
  memory <- memory %||% detect_duckdb_memory_limit()
  config <- list()
  if (!is.null(threads)) {
    config$threads <- as.character(threads)
  }
  if (!is.null(memory)) {
    config$memory_limit <- memory
  }
  return(config)
}

check_duckdb_threads <- function(x) {
  if (is.null(x)) {
    return(NULL)
  }
  ok <- is.numeric(x) &&
    length(x) == 1L &&
    !is.na(x) &&
    is.finite(x) &&
    x >= 1 &&
    x <= .Machine$integer.max &&
    x == trunc(x)
  if (!ok) {
    what <- if (is.numeric(x) && length(x) == 1L) {
      "It is {.val {x}}."
    } else {
      "It is {.obj_type_friendly {x}}."
    }
    cli::cli_abort(
      message = c(
        x = "Option {.code isoformic.duckdb_threads} must be a single positive whole number.",
        i = what,
        i = "Use {.code options(isoformic.duckdb_threads = NULL)} for the default."
      ),
      class = "isoformic_invalid_option"
    )
  }
  return(as.integer(x))
}

check_duckdb_memory_limit <- function(x) {
  if (is.null(x)) {
    return(NULL)
  }
  pattern <- "^[0-9]+(\\.[0-9]+)?\\s*([KMGT]i?B)$"
  ok <- is.character(x) &&
    length(x) == 1L &&
    !is.na(x) &&
    grepl(pattern, trimws(x), ignore.case = TRUE) &&
    memory_size_bytes(x) >= 1024^2
  if (!ok) {
    cli::cli_abort(
      message = c(
        x = "Option {.code isoformic.duckdb_memory_limit} must be a single string with a size and a unit.",
        i = "The size must be at least 1 MiB.",
        i = "For example {.val 8GB}, {.val 512MB} or {.val 16GiB} (units KB, MB, GB, TB, KiB, MiB, GiB, TiB).",
        i = "Use {.code options(isoformic.duckdb_memory_limit = NULL)} for the default."
      ),
      class = "isoformic_invalid_option"
    )
  }
  return(trimws(x))
}

# Bytes of a valid memory size such as "8GB" or "512 MiB" (KB = 1000 bytes,
# KiB = 1024 bytes, as in DuckDB).
memory_size_bytes <- function(x) {
  x <- toupper(trimws(x))
  unit <- sub("^[0-9.]+\\s*", "", x)
  base <- if (grepl("I", unit, fixed = TRUE)) 1024 else 1000
  power <- match(substr(unit, 1L, 1L), c("K", "M", "G", "T"))
  return(as.numeric(sub("\\s*[KMGT]I?B$", "", x)) * base^power)
}

# Positive integer from an environment variable such as `SLURM_CPUS_PER_TASK`
# or `SLURM_MEM_PER_NODE` (in megabytes, optionally with an M/G/T suffix),
# NULL when unset, empty, zero or malformed.
env_positive_number <- function(name, size = FALSE) {
  value <- trimws(Sys.getenv(name))
  pattern <- if (size) "^[0-9]+[KMGT]?$" else "^[0-9]+$"
  if (!grepl(pattern, value, ignore.case = TRUE)) {
    return(NULL)
  }
  if (!size) {
    n <- as.numeric(value)
    return(if (n > 0) n else NULL)
  }
  unit <- toupper(sub("^[0-9]+", "", value))
  scale <- switch(unit, K = 1 / 1024, G = 1024, T = 1024^2, 1)
  n <- as.numeric(sub("[A-Za-z]$", "", value)) * scale
  return(if (n > 0) n else NULL)
}

detect_duckdb_threads <- function() {
  slurm <- env_positive_number("SLURM_CPUS_PER_TASK")
  affinity <- affinity_cpus()
  if (!is.null(affinity) && affinity >= total_cpus()) {
    affinity <- NULL
  }
  if (is.null(slurm) && is.null(affinity)) {
    return(NULL)
  }
  n <- min(c(slurm, affinity, cgroup_cpu_quota()))
  return(max(1L, as.integer(n)))
}

detect_duckdb_memory_limit <- function() {
  mb <- env_positive_number("SLURM_MEM_PER_NODE", size = TRUE)
  if (is.null(mb)) {
    per_cpu <- env_positive_number("SLURM_MEM_PER_CPU", size = TRUE)
    if (!is.null(per_cpu)) {
      cpus <- env_positive_number("SLURM_CPUS_PER_TASK") %||% 1
      mb <- per_cpu * cpus
    }
  }
  if (is.null(mb)) {
    return(NULL)
  }
  # A cgroup smaller than the Slurm allocation (a job step, a container or
  # systemd scope inside the job) wins, as DuckDB would have applied it.
  mb <- min(mb, cgroup_memory_max())
  # At least 64 MiB, so a tiny value does not make every query fail
  return(paste0(max(64, floor(0.8 * mb)), "MiB"))
}

# CPUs of the machine. `parallel::detectCores()` runs a shell command on
# Linux, so the count is computed once per session.
duckdb_limits_cache <- new.env(parent = emptyenv())

total_cpus <- function() {
  if (is.null(duckdb_limits_cache$total_cpus)) {
    n <- suppressWarnings(parallel::detectCores(logical = TRUE))
    duckdb_limits_cache$total_cpus <- if (is.na(n)) Inf else n
  }
  return(duckdb_limits_cache$total_cpus)
}

# Number of CPUs in the affinity mask of this process (Linux only).
affinity_cpus <- function(status_file = "/proc/self/status") {
  if (!file.exists(status_file)) {
    return(NULL)
  }
  lines <- tryCatch(readLines(status_file, warn = FALSE), error = function(e) {
    character()
  })
  line <- grep("^Cpus_allowed_list:", lines, value = TRUE)
  if (length(line) != 1L) {
    return(NULL)
  }
  return(count_cpu_list(sub("^Cpus_allowed_list:\\s*", "", line)))
}

# Count the CPUs of a Linux CPU list such as "0-3,8,10-11".
count_cpu_list <- function(x) {
  parts <- strsplit(trimws(x), ",", fixed = TRUE)[[1L]]
  parts <- parts[nzchar(parts)]
  if (!length(parts) || !all(grepl("^[0-9]+(-[0-9]+)?$", parts))) {
    return(NULL)
  }
  n <- vapply(
    strsplit(parts, "-", fixed = TRUE),
    function(r) {
      r <- as.numeric(r)
      if (length(r) == 1L) 1 else r[2L] - r[1L] + 1
    },
    numeric(1L)
  )
  return(as.integer(sum(n)))
}

# Smallest value of a cgroup v2 interface file (such as `cpu.max` or
# `memory.max`) along the cgroup of this process and its ancestors, as
# parsed by `parse` (which returns NULL for "max" or a malformed value);
# NULL when no level sets a limit.
cgroup_min_limit <- function(
  file_name,
  parse,
  cgroup_file = "/proc/self/cgroup",
  root = "/sys/fs/cgroup"
) {
  if (!file.exists(cgroup_file)) {
    return(NULL)
  }
  lines <- tryCatch(readLines(cgroup_file, warn = FALSE), error = function(e) {
    character()
  })
  line <- grep("^0::", lines, value = TRUE)
  if (length(line) != 1L) {
    return(NULL)
  }
  path <- sub("^0::", "", line)
  limit <- Inf
  repeat {
    limit_file <- file.path(root, sub("^/", "", path), file_name)
    if (file.exists(limit_file)) {
      value <- tryCatch(
        readLines(limit_file, n = 1L, warn = FALSE),
        error = function(e) character()
      )
      if (length(value) == 1L) {
        limit <- min(limit, parse(value) %||% Inf)
      }
    }
    if (path %in% c("/", "")) {
      break
    }
    path <- dirname(path)
  }
  if (is.infinite(limit)) {
    return(NULL)
  }
  return(limit)
}

# Smallest cgroup v2 CPU quota (`cpu.max`) along the cgroup of this process
# and its ancestors, rounded up to whole CPUs; NULL without a quota.
cgroup_cpu_quota <- function(
  cgroup_file = "/proc/self/cgroup",
  root = "/sys/fs/cgroup"
) {
  quota <- cgroup_min_limit(
    "cpu.max",
    parse = function(value) {
      value <- strsplit(value, " ", fixed = TRUE)[[1L]]
      if (length(value) != 2L || !grepl("^[0-9]+$", value[1L])) {
        return(NULL)
      }
      return(ceiling(as.numeric(value[1L]) / as.numeric(value[2L])))
    },
    cgroup_file = cgroup_file,
    root = root
  )
  if (is.null(quota)) {
    return(NULL)
  }
  return(max(1L, as.integer(quota)))
}

# Smallest cgroup v2 memory limit (`memory.max`) along the cgroup of this
# process and its ancestors, in MiB; NULL without a limit.
cgroup_memory_max <- function(
  cgroup_file = "/proc/self/cgroup",
  root = "/sys/fs/cgroup"
) {
  bytes <- cgroup_min_limit(
    "memory.max",
    parse = function(value) {
      value <- trimws(value)
      if (!grepl("^[0-9]+$", value)) {
        return(NULL)
      }
      return(as.numeric(value))
    },
    cgroup_file = cgroup_file,
    root = root
  )
  if (is.null(bytes)) {
    return(NULL)
  }
  return(bytes / 1024^2)
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
#' gzip-compressed) and auto-detects the column types (and the delimiter of
#' `.txt` files; `.tsv`, `.tab` and `.sf` are read as tab-separated).
#' `write_parquet()` writes any data frame to a Parquet file, creating
#' parent directories as needed.
#'
#' @param path Path to the file to read or write.
#' @param x A data frame (or tibble) to write.
#' @param columns Optional character vector of columns to read.
#'   `NULL` reads all columns.
#' @param n_max Maximum number of rows to read. `Inf` (default) reads all.
#' @param compression Parquet compression codec: `"zstd"` (default),
#'   `"snappy"`, `"gzip"` or `"uncompressed"`.
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

#' Read the Rows of a Parquet File Whose Key Is in a Set of Values
#'
#' Same result as `read_parquet(path, columns)` followed by
#' `res[res[[key]] %in% values, ]`, with the selection done by DuckDB: the
#' values are registered as a table and joined with a semi join, so only the
#' selected rows (and only `columns`) reach R. The rows keep the order of the
#' file; it is restored in R from the file row number, which costs less than
#' sorting the rows (possibly long sequences) in DuckDB.
#'
#' With several key columns, `values` is a data frame with those columns and
#' a row is kept when all its key columns match one row of `values`.
#'
#' @param path Path to a Parquet file.
#' @param columns Columns to read; `NULL` reads all.
#' @param key Name of the column to match, or names of the columns.
#' @param values Values of `key` to keep (missing values match missing
#'   values, as with `%in%`); a data frame with the `key` columns when there
#'   are several.
#' @returns A tibble.
#' @keywords internal
#' @noRd
read_parquet_matching <- function(path, columns = NULL, key, values) {
  assert_file_exists(path)
  path <- as.character(fs::path_abs(path))
  con <- duckdb_connect()
  if (length(key) > 1L) {
    wanted <- unique(as.data.frame(values)[key])
  } else {
    wanted <- data.frame(
      key = unique(as.character(values)),
      stringsAsFactors = FALSE
    )
  }
  names(wanted) <- paste0("key", seq_along(key))
  duckdb::duckdb_register(con, "isoformic_wanted_keys", wanted)
  cols <- if (is.null(columns)) {
    DBI::SQL("t.* EXCLUDE (file_row_number)")
  } else {
    DBI::SQL(paste0(
      "t.",
      DBI::dbQuoteIdentifier(con, as.character(columns)),
      collapse = ", "
    ))
  }
  on <- DBI::SQL(paste0(
    "t.",
    DBI::dbQuoteIdentifier(con, key),
    " IS NOT DISTINCT FROM w.",
    names(wanted),
    collapse = " AND "
  ))
  sql <- glue::glue_sql(
    "SELECT {cols}, t.file_row_number AS isoformic_file_row
     FROM read_parquet({path}, file_row_number = true) t
     SEMI JOIN isoformic_wanted_keys w ON {on}",
    .con = con
  )
  res <- tibble::as_tibble(DBI::dbGetQuery(con, sql))
  if (is.unsorted(res$isoformic_file_row)) {
    res <- res[order(res$isoformic_file_row), ]
  }
  res$isoformic_file_row <- NULL
  return(res)
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
