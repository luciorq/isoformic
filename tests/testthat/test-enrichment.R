testthat::test_that("read_gmt parses gene sets", {
  gmt <- withr::local_tempfile(fileext = ".gmt")
  writeLines(c("SET1\tdesc\tA\tB\tC", "SET2\tna\tB\tD", ""), gmt)
  sets <- read_gmt(gmt)
  testthat::expect_named(sets, c("SET1", "SET2"))
  testthat::expect_equal(sets$SET1, c("A", "B", "C"))
  testthat::expect_equal(sets$SET2, c("B", "D"))
  bundled <- read_gmt(system.file(
    "extdata",
    "reactome_subset.gmt",
    package = "isoformic"
  ))
  testthat::expect_true(length(bundled) >= 10L)
})

testthat::test_that("default_enrichment_types picks the processed_transcript label", {
  old <- default_enrichment_types(c("protein_coding", "processed_transcript"))
  testthat::expect_true("processed_transcript" %in% names(old))
  new <- default_enrichment_types(c(
    "protein_coding",
    "protein_coding_CDS_not_defined"
  ))
  testthat::expect_true("protein_coding_CDS_not_defined" %in% names(new))
  testthat::expect_length(new$unproductive, 3L)
})

testthat::test_that("run_enrichment runs fgsea per transcript category", {
  testthat::skip_if_not_installed("fgsea")
  iso <- example_isoformic()
  sets <- read_gmt(system.file(
    "extdata",
    "reactome_subset.gmt",
    package = "isoformic"
  ))
  set.seed(1)
  res <- suppressWarnings(run_enrichment(
    iso,
    sets,
    pval_cutoff = NULL,
    min_size = 3,
    max_size = 500
  ))
  testthat::expect_s3_class(res, "tbl_df")
  testthat::expect_true(all(
    c("pathway", "NES", "padj", "size", "experiment") %in% colnames(res)
  ))
  testthat::expect_true("protein_coding" %in% res$experiment)

  det <- dplyr::left_join(de_tx(iso), tx_annot(iso), by = "transcript_id")
  res_tbl <- suppressWarnings(run_enrichment(
    det,
    sets,
    tx_to_gene = tx_annot(iso),
    pval_cutoff = NULL,
    min_size = 3
  ))
  testthat::expect_true(nrow(res_tbl) > 0L)
  testthat::expect_error(
    run_enrichment(iso, list(1, 2)),
    class = "isoformic_invalid_genesets"
  )
  testthat::expect_error(
    run_enrichment(
      det[, "log2FC", drop = FALSE],
      sets,
      tx_to_gene = tx_annot(iso)
    ),
    class = "isoformic_invalid_de_table"
  )
})
