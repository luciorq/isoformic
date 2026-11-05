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
