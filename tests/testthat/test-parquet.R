testthat::test_that("write_parquet and read_parquet round-trip a data frame", {
  path <- withr::local_tempfile(fileext = ".parquet")
  df <- data.frame(
    id = c("a", "b", "c"),
    value = c(1.5, 2, NA),
    flag = c(TRUE, FALSE, TRUE),
    f = factor(c("x", "y", "x")),
    stringsAsFactors = FALSE
  )
  testthat::expect_identical(
    write_parquet(df, path),
    as.character(fs::path_abs(path))
  )
  res <- read_parquet(path)
  testthat::expect_s3_class(res, "tbl_df")
  testthat::expect_equal(res$id, df$id)
  testthat::expect_equal(res$value, df$value)
  testthat::expect_equal(res$flag, df$flag)
  testthat::expect_type(res$f, "character")
})

testthat::test_that("read_parquet supports column selection and n_max", {
  path <- withr::local_tempfile(fileext = ".parquet")
  write_parquet(data.frame(a = 1:5, b = letters[1:5]), path)
  res <- read_parquet(path, columns = "b", n_max = 2)
  testthat::expect_named(res, "b")
  testthat::expect_equal(nrow(res), 2L)
  testthat::expect_equal(parquet_nrow(path), 5L)
})

testthat::test_that("write_parquet creates parent directories and rejects non data frames", {
  dir <- withr::local_tempdir()
  path <- fs::path(dir, "nested", "deeper", "x.parquet")
  write_parquet(data.frame(a = 1), path)
  testthat::expect_true(fs::file_exists(path))
  testthat::expect_error(
    write_parquet(1:3, path),
    class = "isoformic_not_data_frame"
  )
})

testthat::test_that("read_tabular reads csv, tsv and parquet", {
  csv <- withr::local_tempfile(fileext = ".csv")
  tsv <- withr::local_tempfile(fileext = ".tsv")
  pq <- withr::local_tempfile(fileext = ".parquet")
  df <- data.frame(id = c("x", "y"), v = c(1, 2))
  utils::write.csv(df, csv, row.names = FALSE)
  utils::write.table(df, tsv, sep = "\t", row.names = FALSE, quote = FALSE)
  write_parquet(df, pq)
  testthat::expect_equal(as.data.frame(read_tabular(csv)), df)
  testthat::expect_equal(as.data.frame(read_tabular(tsv)), df)
  testthat::expect_equal(as.data.frame(read_tabular(pq)), df)
  testthat::expect_error(
    read_tabular("/no/such/file.csv"),
    class = "isoformic_annot_file_dont_exist"
  )
  bad <- withr::local_tempfile(fileext = ".xyz")
  writeLines("a", bad)
  testthat::expect_error(
    read_tabular(bad),
    class = "isoformic_unsupported_file_type"
  )
})

testthat::test_that("duckdb_run validates read_only and executes statements", {
  testthat::expect_error(
    duckdb_run("SELECT 1", read_only = "yes"),
    class = "isoformic_read_only_not_logical"
  )
  testthat::expect_null(duckdb_run("SELECT 1 AS x;"))
  db <- withr::local_tempfile(fileext = ".duckdb")
  testthat::expect_null(duckdb_run(
    "CREATE TABLE t AS SELECT 1 AS x;",
    db_file_path = db
  ))
  testthat::expect_true(fs::file_exists(db))
})

testthat::test_that("a query that spills to disk succeeds", {
  root <- fs::path(tempdir(), "isoformic-duckdb")
  if (fs::dir_exists(root)) {
    fs::dir_delete(root)
  }
  withr::local_options(
    isoformic.duckdb_threads = 1,
    isoformic.duckdb_memory_limit = "100MB"
  )
  con <- duckdb_connect()
  res <- DBI::dbGetQuery(
    con,
    paste(
      "SELECT count(*) AS n, count(DISTINCT k) AS k FROM (",
      "SELECT md5(CAST(i AS VARCHAR)) AS k FROM range(3000000) t(i)",
      "ORDER BY k)"
    )
  )
  testthat::expect_equal(res$n, 3e6)
  testthat::expect_equal(res$k, 3e6)
  testthat::expect_true(fs::dir_exists(root))
})

duckdb_settings <- function() {
  con <- duckdb_connect()
  DBI::dbGetQuery(
    con,
    paste(
      "SELECT current_setting('threads') AS threads,",
      "current_setting('memory_limit') AS memory_limit"
    )
  )
}

testthat::test_that("duckdb options set threads and memory_limit", {
  withr::local_options(
    isoformic.duckdb_threads = 2,
    isoformic.duckdb_memory_limit = "1GB"
  )
  res <- duckdb_settings()
  testthat::expect_equal(as.integer(res$threads), 2L)
  # DuckDB reports 1GB (10^9 bytes) as 953.6 MiB
  testthat::expect_match(res$memory_limit, "^953\\.[0-9] MiB$")
  testthat::expect_identical(
    duckdb_limits_config(),
    list(threads = "2", memory_limit = "1GB")
  )
  withr::local_options(isoformic.duckdb_memory_limit = " 512 MiB ")
  testthat::expect_identical(duckdb_settings()$memory_limit, "512.0 MiB")
})

testthat::test_that("invalid duckdb options raise clear errors", {
  for (bad in list(0, -1, 1.5, 3e9, 1e10, NA_real_, Inf, "2", c(1, 2), TRUE)) {
    withr::with_options(
      list(isoformic.duckdb_threads = bad),
      testthat::expect_error(
        duckdb_limits_config(),
        regexp = "isoformic.duckdb_threads",
        class = "isoformic_invalid_option"
      )
    )
  }
  for (bad in list(
    "1",
    "8 gigabytes",
    "80%",
    1e9,
    "0GB",
    "0.0MiB",
    "0.5KB",
    "1kb",
    "1000KB",
    NA_character_,
    c("1GB", "2GB")
  )) {
    withr::with_options(
      list(isoformic.duckdb_memory_limit = bad),
      testthat::expect_error(
        read_parquet(system.file("DESCRIPTION", package = "isoformic")),
        regexp = "isoformic.duckdb_memory_limit",
        class = "isoformic_invalid_option"
      )
    )
  }
  # The smallest accepted size
  withr::local_options(isoformic.duckdb_memory_limit = "1MiB")
  testthat::expect_identical(duckdb_limits_config()$memory_limit, "1MiB")
  # A numeric value is shown in the error
  withr::local_options(isoformic.duckdb_threads = 0)
  testthat::expect_error(duckdb_limits_config(), regexp = "It is 0")
})

testthat::test_that("Slurm variables set the default threads and memory_limit", {
  withr::local_options(
    isoformic.duckdb_threads = NULL,
    isoformic.duckdb_memory_limit = NULL
  )
  withr::local_envvar(
    SLURM_CPUS_PER_TASK = "2",
    SLURM_MEM_PER_NODE = NA,
    SLURM_MEM_PER_CPU = "1000"
  )
  testthat::local_mocked_bindings(cgroup_memory_max = function(...) NULL)
  config <- duckdb_limits_config()
  testthat::expect_lte(as.integer(config$threads), 2L)
  testthat::expect_identical(config$memory_limit, "1600MiB")
  res <- duckdb_settings()
  testthat::expect_equal(as.integer(res$threads), as.integer(config$threads))
  testthat::expect_identical(res$memory_limit, "1.5 GiB")

  withr::local_envvar(SLURM_MEM_PER_NODE = "4000")
  testthat::expect_identical(detect_duckdb_memory_limit(), "3200MiB")
  withr::local_envvar(SLURM_MEM_PER_NODE = "2G")
  testthat::expect_identical(detect_duckdb_memory_limit(), "1638MiB")
  # 0 means the whole node in Slurm, malformed values are ignored
  # A tiny allocation still leaves DuckDB 64 MiB
  withr::local_envvar(SLURM_MEM_PER_NODE = "1")
  testthat::expect_identical(detect_duckdb_memory_limit(), "64MiB")
  # A cgroup memory limit below the Slurm allocation wins
  withr::local_envvar(SLURM_MEM_PER_NODE = "20000")
  testthat::local_mocked_bindings(cgroup_memory_max = function(...) 4096)
  testthat::expect_identical(detect_duckdb_memory_limit(), "3276MiB")
  testthat::local_mocked_bindings(cgroup_memory_max = function(...) 1e6)
  testthat::expect_identical(detect_duckdb_memory_limit(), "16000MiB")
  # 0 means the whole node in Slurm, malformed values are ignored
  withr::local_envvar(SLURM_MEM_PER_NODE = "0", SLURM_MEM_PER_CPU = "abc")
  testthat::expect_null(detect_duckdb_memory_limit())
})

testthat::test_that("without options or Slurm, the DuckDB defaults are kept", {
  withr::local_options(
    isoformic.duckdb_threads = NULL,
    isoformic.duckdb_memory_limit = NULL
  )
  withr::local_envvar(
    SLURM_CPUS_PER_TASK = NA,
    SLURM_MEM_PER_NODE = NA,
    SLURM_MEM_PER_CPU = NA
  )
  testthat::local_mocked_bindings(affinity_cpus = function(...) NULL)
  testthat::expect_identical(duckdb_limits_config(), list())
})

testthat::test_that("the affinity mask and the cgroup CPU quota are read", {
  testthat::expect_identical(count_cpu_list("0-3,8,10-11"), 7L)
  testthat::expect_identical(count_cpu_list("5"), 1L)
  testthat::expect_null(count_cpu_list("garbage"))

  dir <- withr::local_tempdir()
  status <- fs::path(dir, "status")
  writeLines(c("Name:\tR", "Cpus_allowed_list:\t0-1,4"), status)
  testthat::expect_identical(affinity_cpus(status), 3L)
  testthat::expect_null(affinity_cpus(fs::path(dir, "missing")))

  cgroup <- fs::path(dir, "cgroup")
  writeLines("0::/a.slice/b.scope", cgroup)
  root <- fs::path(dir, "sys")
  fs::dir_create(fs::path(root, "a.slice", "b.scope"))
  writeLines("max 100000", fs::path(root, "a.slice", "b.scope", "cpu.max"))
  testthat::expect_null(cgroup_cpu_quota(cgroup, root))
  writeLines("250000 100000", fs::path(root, "a.slice", "cpu.max"))
  testthat::expect_identical(cgroup_cpu_quota(cgroup, root), 3L)

  testthat::expect_null(cgroup_memory_max(cgroup, root))
  writeLines("max", fs::path(root, "a.slice", "b.scope", "memory.max"))
  testthat::expect_null(cgroup_memory_max(cgroup, root))
  writeLines(
    as.character(4 * 1024^3),
    fs::path(root, "a.slice", "memory.max")
  )
  writeLines(as.character(8 * 1024^3), fs::path(root, "memory.max"))
  testthat::expect_identical(cgroup_memory_max(cgroup, root), 4096)
  testthat::expect_null(cgroup_memory_max(fs::path(dir, "missing"), root))

  testthat::local_mocked_bindings(
    affinity_cpus = function(...) 3L,
    total_cpus = function() 16,
    cgroup_cpu_quota = function(...) 2L
  )
  withr::local_envvar(SLURM_CPUS_PER_TASK = NA)
  testthat::expect_identical(detect_duckdb_threads(), 2L)
  testthat::local_mocked_bindings(affinity_cpus = function(...) 16L)
  testthat::expect_null(detect_duckdb_threads())
})
