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

testthat::test_that("run_enrichment ranks by fold-change, signed p-value or replicate z-score", {
  det <- tibble::tibble(
    transcript_name = paste0("t", 1:6),
    transcript_type = "protein_coding",
    log2FC = c(2, -1, 0.5, -3, 1, NA),
    pvalue = c(0.01, 0.1, 1, 0, 0.001, 0.5),
    lfc_median = c(2, -1, 0.5, -3, 0, 1),
    lfc_sd = c(1, 0.5, NA, 0, 0, 2)
  )
  map <- tibble::tibble(
    transcript_name = paste0("t", 1:6),
    gene_name = paste0("g", 1:6)
  )
  sets <- list(set = paste0("g", 1:6))
  types <- list(protein_coding = "protein_coding")
  calls <- record_fgsea()
  ranks_of <- function(...) {
    calls$stats <- list()
    res <- run_enrichment(
      det,
      sets,
      tx_to_gene = map,
      tx_types = types,
      pval_cutoff = NULL,
      ...
    )
    list(res = res, ranks = calls$stats[[1L]])
  }
  by_lfc <- ranks_of()
  testthat::expect_identical(unique(by_lfc$res$rank_by), "log2FC")
  testthat::expect_identical(
    by_lfc$ranks,
    sort(c(t1 = 2, t2 = -1, t3 = 0.5, t4 = -3, t5 = 1))
  )
  # -log10(p) with the sign of the change; p = 0 is the largest of the table
  by_p <- ranks_of(rank_by = "signed_p")
  testthat::expect_identical(unique(by_p$res$rank_by), "signed_p")
  testthat::expect_equal(
    by_p$ranks[paste0("t", 1:5)],
    c(t1 = 2, t2 = -1, t3 = 0, t4 = -3, t5 = 3)
  )
  # Median over standard deviation; no deviation is the most certain change,
  # no change at all is zero, and a missing deviation is not ranked
  by_z <- ranks_of(rank_by = "lfc_z")
  testthat::expect_identical(unique(by_z$res$rank_by), "lfc_z")
  testthat::expect_equal(
    by_z$ranks[c("t1", "t2", "t4", "t5", "t6")],
    c(t1 = 2, t2 = -2, t4 = -2, t5 = 0, t6 = 0.5)
  )
  testthat::expect_false("t3" %in% names(by_z$ranks))

  testthat::expect_error(
    ranks_of(rank_by = "other"),
    class = "rlang_error"
  )
  det$pvalue <- NULL
  testthat::expect_error(
    ranks_of(rank_by = "signed_p"),
    class = "isoformic_invalid_de_table"
  )
  det$lfc_sd <- NULL
  testthat::expect_error(
    ranks_of(rank_by = "lfc_z"),
    class = "isoformic_no_inferential_replicates"
  )
})

testthat::test_that("run_enrichment takes the z-score of the object from its confidence table", {
  x <- uncertainty_experiment()
  iso <- add_uncertainty(x$iso, paths = x$fixture$dirs)
  sets <- list(set = unique(tx_annot(iso)$gene_name))
  types <- list(all = unique(tx_annot(iso)$transcript_type))
  calls <- record_fgsea()
  testthat::expect_error(
    run_enrichment(iso, sets, tx_types = types, rank_by = "lfc_z"),
    class = "isoformic_no_inferential_replicates"
  )
  iso <- add_uncertainty(iso, group = "condition")
  conf <- de_confidence(iso)
  res <- run_enrichment(
    iso,
    sets,
    tx_types = types,
    pval_cutoff = NULL,
    rank_by = "lfc_z"
  )
  testthat::expect_identical(unique(res$rank_by), "lfc_z")
  ranks <- calls$stats[[1L]]
  conf <- conf[conf$feature_type == "transcript", ]
  expected <- certainty_z(conf$lfc_median, conf$lfc_sd)
  names(expected) <- conf$feature_id
  expected <- expected[
    !is.na(expected) & names(expected) %in% de_tx(iso)$transcript_id
  ]
  testthat::expect_gt(length(ranks), 0L)
  testthat::expect_setequal(names(ranks), names(expected))
  testthat::expect_equal(ranks, sort(expected)[names(ranks)])
})

testthat::test_that("several collections of gene sets are tested one by one", {
  det <- tibble::tibble(
    transcript_name = paste0("t", 1:4),
    transcript_type = "protein_coding",
    log2FC = c(2, -1, 0.5, -3)
  )
  map <- tibble::tibble(
    transcript_name = paste0("t", 1:4),
    gene_name = paste0("g", 1:4)
  )
  types <- list(protein_coding = "protein_coding")
  first <- list(a = c("g1", "g2"), b = c("g3", "g4"))
  second <- list(c = c("g1", "g4"))
  calls <- record_fgsea()
  flat <- run_enrichment(
    det,
    first,
    tx_to_gene = map,
    tx_types = types,
    pval_cutoff = NULL
  )
  testthat::expect_false("collection" %in% colnames(flat))
  calls$pathways <- list()
  nested <- run_enrichment(
    det,
    list(first = first, second = second),
    tx_to_gene = map,
    tx_types = types,
    pval_cutoff = NULL
  )
  testthat::expect_identical(nested$collection, c("first", "first", "second"))
  testthat::expect_identical(nested$pathway, c("a", "b", "c"))
  # One call per collection, each with its own sets only
  testthat::expect_identical(
    lapply(calls$pathways, names),
    list(c("a", "b"), "c")
  )
  testthat::expect_identical(
    as.data.frame(nested[nested$collection == "first", colnames(flat)]),
    as.data.frame(flat)
  )
  for (bad in list(
    list(first = first, second = "g1"),
    list(first, second),
    list(first = first, first = second),
    list(first = list(a = list("g1"))),
    "g1"
  )) {
    testthat::expect_error(
      run_enrichment(det, bad, tx_to_gene = map, tx_types = types),
      class = "isoformic_invalid_genesets"
    )
  }
})
