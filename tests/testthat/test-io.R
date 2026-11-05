testthat::test_that("write_isoformic can overwrite the directory it was read from", {
  iso <- example_isoformic()
  dir <- withr::local_tempdir()
  out <- file.path(dir, "exp")
  iso <- write_isoformic(iso, out, copy_annotation = TRUE)
  again <- read_isoformic(out)
  # Same directory, with and without copying the annotation
  testthat::expect_no_error(
    again <- write_isoformic(again, overwrite = TRUE, copy_annotation = TRUE)
  )
  testthat::expect_no_error(
    again <- write_isoformic(again, overwrite = TRUE)
  )
  back <- read_isoformic(out)
  testthat::expect_identical(dim(back), dim(iso))
  testthat::expect_true(fs::dir_exists(fs::path(out, "annotation")))
  testthat::expect_length(fs::dir_ls(dir, glob = "*.writing-*"), 0L)
})

testthat::test_that("write_isoformic refuses to delete other directories", {
  iso <- example_isoformic()
  dir <- withr::local_tempdir()
  writeLines("keep me", file.path(dir, "notes.txt"))
  testthat::expect_error(
    write_isoformic(iso, dir, overwrite = TRUE),
    class = "isoformic_not_experiment_dir"
  )
  testthat::expect_true(file.exists(file.path(dir, "notes.txt")))
})
