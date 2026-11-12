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
    list(first = list(a = list(c("g1", "g2")))),
    stats::setNames(list(first, second), c("first", NA)),
    stats::setNames(first, c("a", NA)),
    stats::setNames(first, c("a", "")),
    "g1"
  )) {
    testthat::expect_error(
      run_enrichment(det, bad, tx_to_gene = map, tx_types = types),
      class = "isoformic_invalid_genesets"
    )
  }
})

testthat::test_that("gene sets may be empty or lists of symbols, and no result keeps its columns", {
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
  calls <- record_fgsea()
  plain <- run_enrichment(
    det,
    list(a = c("g1", "g2"), b = c("g3", "g4")),
    tx_to_gene = map,
    tx_types = types,
    pval_cutoff = NULL
  )
  odd <- run_enrichment(
    det,
    list(a = list("g1", "g2"), b = factor(c("g3", "g4")), none = NULL),
    tx_to_gene = map,
    tx_types = types,
    pval_cutoff = NULL
  )
  testthat::expect_identical(odd, plain)
  testthat::expect_identical(
    geneset_collections(list(a = list("g1", "g2"), none = NULL)),
    list(list(a = c("g1", "g2"), none = character()))
  )
  # No category with two transcripts: the columns of fgsea, no rows
  for (sets in list(
    list(a = c("g1", "g2")),
    list(first = list(a = c("g1", "g2")), second = list(b = "g3"))
  )) {
    none <- run_enrichment(
      det,
      sets,
      tx_to_gene = map,
      tx_types = list(other = "lncRNA")
    )
    testthat::expect_identical(nrow(none), 0L)
    testthat::expect_true(all(
      c("pathway", "padj", "NES", "size", "experiment", "rank_by") %in%
        colnames(none)
    ))
  }
  testthat::expect_true("collection" %in% colnames(none))
})

testthat::test_that("lfc_z takes the de_confidence() table and says when it cannot rank", {
  x <- uncertainty_experiment()
  iso <- add_uncertainty(x$iso, paths = x$fixture$dirs)
  iso <- add_uncertainty(iso, group = "condition")
  conf <- de_confidence(iso)
  conf <- conf[conf$feature_type == "transcript", ]
  annot <- tx_annot(iso)
  sets <- list(set = unique(annot$gene_name))
  types <- list(all = unique(annot$transcript_type))
  calls <- record_fgsea()
  from_object <- run_enrichment(
    iso,
    sets,
    tx_types = types,
    pval_cutoff = NULL,
    rank_by = "lfc_z"
  )
  from_table <- run_enrichment(
    conf,
    sets,
    tx_to_gene = annot,
    tx_types = types,
    pval_cutoff = NULL,
    rank_by = "lfc_z"
  )
  # The object ranks the transcripts of its results, the table all its rows
  by_name <- stats::setNames(annot$transcript_name, annot$transcript_id)
  of_object <- calls$stats[[1L]]
  names(of_object) <- by_name[names(of_object)]
  of_table <- calls$stats[[2L]]
  testthat::expect_gt(length(of_object), 0L)
  testthat::expect_true(all(names(of_object) %in% names(of_table)))
  testthat::expect_equal(of_table[names(of_object)], of_object)
  testthat::expect_identical(unique(from_table$rank_by), "lfc_z")
  testthat::expect_identical(unique(from_object$rank_by), "lfc_z")
  # One replicate has no standard deviation
  conf$lfc_sd <- NA_real_
  testthat::expect_error(
    run_enrichment(
      conf,
      sets,
      tx_to_gene = annot,
      tx_types = types,
      rank_by = "lfc_z"
    ),
    class = "isoformic_no_inferential_replicates"
  )
  # The confidence of the opposite contrast is not the one of the results
  reversed <- add_uncertainty(
    iso,
    group = "condition",
    levels = rev(levels(factor(col_data(iso)$condition)))
  )
  agree <- sign(de_confidence(iso)$lfc_median) ==
    sign(de_confidence(iso)$det_log2FC)
  expect_mismatch <- if (mean(agree, na.rm = TRUE) >= 0.5) reversed else iso
  testthat::expect_warning(
    run_enrichment(
      expect_mismatch,
      sets,
      tx_types = types,
      pval_cutoff = NULL,
      rank_by = "lfc_z"
    ),
    class = "isoformic_contrast_mismatch"
  )
  expect_match <- if (mean(agree, na.rm = TRUE) >= 0.5) iso else reversed
  testthat::expect_no_warning(
    run_enrichment(
      expect_match,
      sets,
      tx_types = types,
      pval_cutoff = NULL,
      rank_by = "lfc_z"
    )
  )
})

testthat::test_that("plot_enrichment draws one panel per collection", {
  res <- tibble::tibble(
    pathway = c("a", "b", "a"),
    NES = c(1, -1, 2),
    padj = c(0.01, 0.02, 0.03),
    size = c(10L, 20L, 10L),
    experiment = "protein_coding",
    collection = c("first", "first", "second")
  )
  p <- plot_enrichment(res)
  built <- ggplot2::ggplot_build(p)
  testthat::expect_identical(length(unique(built$data[[1L]]$PANEL)), 2L)
  single <- plot_enrichment(res[
    res$collection == "first",
    colnames(res) != "collection"
  ])
  testthat::expect_identical(
    length(unique(ggplot2::ggplot_build(single)$data[[1L]]$PANEL)),
    1L
  )
})
