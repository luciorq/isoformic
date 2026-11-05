testthat::test_that("get_isoformic_cache builds paths under the user cache directory", {
  base <- get_isoformic_cache()
  testthat::expect_match(base, "isoformic")
  testthat::expect_equal(
    get_isoformic_cache("a", "b"),
    fs::path(base, "a", "b")
  )
  testthat::expect_equal(
    get_isoformic_cache("x", ext = ".RDS"),
    fs::path(base, "x", ext = "rds")
  )
})

testthat::test_that("random experiment names are unique strings", {
  n1 <- random_experiment_name()
  n2 <- random_experiment_name()
  testthat::expect_match(n1, "^isoformic_[0-9a-f]{10}$")
  testthat::expect_false(identical(n1, n2))
})

testthat::test_that("kv_to_list and list_to_kv round-trip scalar metadata", {
  x <- list(a = "1", b = "two")
  testthat::expect_equal(kv_to_list(list_to_kv(x)), x)
  testthat::expect_warning(kv <- list_to_kv(list(a = 1, v = 1:3)), "non-scalar")
  testthat::expect_equal(kv$key, "a")
  testthat::expect_equal(kv_to_list(NULL), list())
})
