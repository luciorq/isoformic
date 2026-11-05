read_knockdown_example <- function(series) {
  dir <- fs::path_package("isoformic", "extdata", "knockdown", series)
  iso <- IsoformicExperiment(
    experiment_name = series,
    annot_path = fs::path_package(
      "isoformic",
      "extdata",
      "knockdown_annotation"
    ),
    assays = list(
      abundance = read_parquet(fs::path(dir, "abundance.parquet")),
      counts = read_parquet(fs::path(dir, "counts.parquet"))
    ),
    col_data = read_parquet(fs::path(dir, "col_data.parquet"))
  )
  de_tx(iso) <- read_parquet(fs::path(dir, "det.parquet"))
  iso@metadata$inf_reps <- as_inferential_replicates(
    read_parquet(fs::path(dir, "draws.parquet"))
  )
  list(
    iso = iso,
    equivalence = read_parquet(fs::path(dir, "equivalence.parquet"))
  )
}

testthat::test_that("the bundled knockdown examples carry draws for every sample", {
  for (series in c("GSE129540", "GSE75489")) {
    x <- read_knockdown_example(series)
    draws <- x$iso@metadata$inf_reps
    testthat::expect_length(draws, 30L)
    testthat::expect_identical(
      colnames(draws[[1L]]),
      colnames(assay(x$iso, "counts"))
    )
    testthat::expect_true(all(rownames(draws[[1L]]) %in% rownames(x$iso)))
    testthat::expect_identical(
      col_data(x$iso)$num_replicates,
      rep(30, 6L)
    )
    testthat::expect_true(all(col_data(x$iso)$library_size > 1e6))
  }
})

testthat::test_that("the ESRP knockdown reproduces the expected isoform calls", {
  prostate <- read_knockdown_example("GSE129540")
  iso <- add_uncertainty(
    prostate$iso,
    group = "condition",
    levels = c("control", "knockdown"),
    equivalence = prostate$equivalence
  )
  conf <- de_confidence(iso, gene_name = "SCRIB")
  scrib <- conf[conf$transcript_name == "SCRIB-207", ]
  testthat::expect_identical(as.character(scrib$confidence), "high")
  testthat::expect_lt(scrib$lfc_hi, -1)
  testthat::expect_lt(scrib$mean_inf_rv, 0.5)
  amb <- isoform_ambiguity(iso, gene_name = "SCRIB", ambiguous_only = TRUE)
  testthat::expect_setequal(
    c(amb$transcript_name_1, amb$transcript_name_2),
    c("SCRIB-201", "SCRIB-211")
  )

  lung <- read_knockdown_example("GSE75489")
  iso_lung <- add_uncertainty(
    lung$iso,
    group = "condition",
    levels = c("control", "knockdown"),
    equivalence = lung$equivalence
  )
  cd44 <- de_confidence(iso_lung, gene_name = "CD44")
  short <- cd44[cd44$transcript_name == "CD44-201", ]
  testthat::expect_gt(short$lfc_lo, 1)
  testthat::expect_gte(short$sign_agreement, 0.95)
  # The reads cannot separate the CD44 transcripts
  testthat::expect_identical(as.character(short$confidence), "moderate")
})

testthat::test_that("library_size keeps the scaling of a transcript subset", {
  x <- read_knockdown_example("GSE129540")
  with_size <- add_uncertainty(
    x$iso,
    group = "condition",
    levels = c("control", "knockdown")
  )
  sample_table <- col_data(x$iso)
  sample_table$library_size <- NULL
  col_data(x$iso) <- sample_table
  without_size <- add_uncertainty(
    x$iso,
    group = "condition",
    levels = c("control", "knockdown")
  )
  a <- de_confidence(with_size)
  b <- de_confidence(without_size)
  both <- merge(a, b, by = "feature_id", suffixes = c("_full", "_subset"))
  # The subset has other library sizes than the complete quantification
  testthat::expect_gt(
    max(abs(both$lfc_median_full - both$lfc_median_subset)),
    0.01
  )
})
