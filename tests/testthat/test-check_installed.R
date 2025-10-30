testthat::test_that("check_installed interactive", {
  rlang::local_interactive(TRUE)

  res <- check_installed("testthat")
  testthat::expect_equal(res, TRUE)

  res <- check_installed(pkgs = c("testthat", "rlang"))
  testthat::expect_equal(res, TRUE)
})

testthat::test_that("check_installed non interactive", {
  rlang::local_interactive(FALSE)

  res <- check_installed("testthat")
  testthat::expect_equal(res, TRUE)

  res_missing <- check_installed("someNonExistentPackage12345")
  testthat::expect_equal(res_missing, FALSE)
})
