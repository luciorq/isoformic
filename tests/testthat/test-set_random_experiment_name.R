testthat::test_that("set_random_experiment_name - generates a valid name", {
  name <- set_random_experiment_name()
  testthat::expect_true(is.character(name))
  testthat::expect_equal(length(name), 1)
})

testthat::test_that("set_random_experiment_name - creates data_path directory", {
  data_path <- withr::local_tempdir("isoformic_data")
  name <- set_random_experiment_name(data_path = data_path)
  dir_path <- fs::path(data_path, name)
  testthat::expect_true(fs::dir_exists(dir_path))
})

testthat::test_that("set_random_experiment_name - uses provided experiment_name", {
  custom_name <- "custom_experiment"
  data_path <- withr::local_tempdir(custom_name)
  name <- set_random_experiment_name(
    data_path = data_path,
    experiment_name = custom_name
  )
  testthat::expect_equal(name, custom_name)
  dir_path <- fs::path(data_path, name)
  testthat::expect_true(fs::dir_exists(fs::path_dir(dir_path)))
  testthat::expect_false(fs::dir_exists(dir_path))
})

testthat::test_that("set_random_experiment_name - handles invalid data_path", {
  invalid_path <- "/invalid_path_that_does_not_exist"
  testthat::expect_error(
    set_random_experiment_name(data_path = invalid_path),
    class = "isoformic_data_path_access_error"
  )
})
