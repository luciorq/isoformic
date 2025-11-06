testthat::test_that("IsoformicExperiment - class is registered", {
  testthat::expect_true(S7::S7_inherits(IsoformicExperiment))

  testthat::expect_true(inherits(x = IsoformicExperiment, what = "S7_class"))
})


testthat::test_that("IsoformicExperiment - object is instantiated", {
  data_path <- withr::local_tempdir("isoformic_data")
  annot_path <- withr::local_tempdir("isoformic_annot")

  iso <- IsoformicExperiment()

  testthat::expect_equal(iso@experiment_name, NA_character_)

  testthat::expect_s7_class(iso, IsoformicExperiment)

  iso@experiment_name
  testthat::expect_equal(iso@experiment_name, NA_character_)

  iso@experiment_name <- "123"

  testthat::expect_equal(iso@experiment_name, "123")

  testthat::expect_error(
    object = {
      iso@experiment_name <- NULL
    }
  )

  testthat::expect_error(
    object = {
      iso@experiment_name <- 123
    }
  )

  testthat::expect_error(
    object = {
      iso@experiment_name <- list()
    }
  )

  testthat::expect_error(
    object = {
      iso@experiment_name <- c("name1", "name2")
    }
  )

  iso@experiment_name <- "new_name"

  testthat::expect_equal(iso@experiment_name, "new_name")

  iso@data_path <- NULL

  testthat::expect_equal(iso@data_path, NULL)

  iso@data_path <- data_path

  testthat::expect_equal(iso@data_path, data_path)

  testthat::expect_true(fs::dir_exists(iso@data_path))

  iso@annot_path <- NULL

  testthat::expect_equal(iso@annot_path, NULL)

  testthat::expect_equal(rownames(iso), NULL)

  testthat::expect_equal(colnames(iso), NULL)

  testthat::expect_equal(dimnames(iso), list(NULL, NULL))

  testthat::expect_equal(dim(iso), c(0L, 0L))

  iso@annot_path <- annot_path

  testthat::expect_equal(iso@annot_path, annot_path)

  testthat::expect_true(fs::file_exists(iso@annot_path))

  iso@annot_path <- NULL

  testthat::expect_equal(iso@annot_path, NULL)

  testthat::expect_equal(iso@annot_metadata, NULL)

  iso@annot_metadata <- NULL

  iso@annot_metadata <- list()

  iso@annot_metadata[["test"]] <- "value"

  testthat::expect_true(inherits(iso@annot_metadata, "list"))

  # attributes(iso@annot_metadata)

  iso@annot_metadata <- NULL

  # TODO: @luciorq - This should be an error in a newer version
  iso@annot_metadata <- data.frame(
    key = c("a", "b"),
    value = c("1", "2"),
    stringsAsFactors = FALSE
  )

  iso@annot_metadata <- NULL

  testthat::expect_equal(iso@annot_metadata, NULL)
})


testthat::test_that("IsoformicExperiment - generics are defined", {
  data_path <- withr::local_tempdir("isoformic_data")
  annot_path <- withr::local_tempdir("isoformic_annot")

  iso <- IsoformicExperiment(
    experiment_name = "test_experiment",
    data_path = data_path,
    annot_path = annot_path
  )

  # S7 Public API
  testthat::expect_true(S7::S7_inherits(col_data))
  testthat::expect_true(S7::S7_inherits(row_data))
  testthat::expect_true(S7::S7_inherits(annot_data))
  testthat::expect_true(S7::S7_inherits(annot_data_transcripts))
  testthat::expect_true(S7::S7_inherits(annot_data_genes))
  testthat::expect_true(S7::S7_inherits(annot_data_exons))
  testthat::expect_true(S7::S7_inherits(tx_to_gene))
  testthat::expect_true(S7::S7_inherits(de_gene))
  testthat::expect_true(S7::S7_inherits(de_tx))
  testthat::expect_true(S7::S7_inherits(tx_annot))

  testthat::expect_s3_class(row_data(iso), "data.frame")
  testthat::expect_s3_class(row_data(iso), "tbl_df")

  testthat::expect_equal(iso@experiment_name, "test_experiment")
  testthat::expect_equal(iso@data_path, data_path)
  testthat::expect_equal(iso@annot_path, annot_path)
})
