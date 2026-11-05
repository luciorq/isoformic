testthat::test_that("test_usage finds the isoform switches with every backend", {
  x <- usage_experiment()
  for (method in c("diffsplice", "satuRn", "DEXSeq")) {
    testthat::skip_if_not_installed(switch(
      method,
      diffsplice = "edgeR",
      satuRn = "satuRn",
      DEXSeq = "DEXSeq"
    ))
    res <- suppressWarnings(
      test_usage(x$iso, group = "condition", method = method)
    )
    testthat::expect_identical(
      colnames(res),
      c(
        "gene_id",
        "gene_name",
        "transcript_id",
        "transcript_name",
        "transcript_type",
        "contrast",
        "usage_log2FC",
        "pvalue",
        "qvalue",
        "gene_pvalue",
        "gene_qvalue",
        "method",
        "scaled"
      )
    )
    testthat::expect_identical(unique(res$method), method)
    testthat::expect_identical(unique(res$contrast), "b vs a")
    testthat::expect_false(unique(res$scaled))
    testthat::expect_true(all(res$pvalue >= 0 & res$pvalue <= 1))
    # The switching genes are found, the others mostly not
    genes <- unique(res[, c("gene_id", "gene_qvalue")])
    found <- genes$gene_id[genes$gene_qvalue < 0.05]
    testthat::expect_gte(sum(x$switch_genes %in% found), 6L)
    testthat::expect_lte(sum(!found %in% x$switch_genes), 4L)
    # The second isoform loses usage and the third gains it
    t2 <- res[res$transcript_id == "g1t2", ]
    t3 <- res[res$transcript_id == "g1t3", ]
    testthat::expect_lt(t2$usage_log2FC, 0)
    testthat::expect_gt(t3$usage_log2FC, 0)
  }
})

testthat::test_that("test_usage handles groups, filters, scaling and missing backends", {
  testthat::skip_if_not_installed("edgeR")
  x <- usage_experiment()
  iso <- x$iso
  sample_table <- col_data(iso)
  sample_table$condition <- c("a", "a", "a", "b", "b", "other")
  col_data(iso) <- sample_table
  # A third group, even of one sample, becomes a second contrast
  all_groups <- test_usage(iso, group = "condition")
  testthat::expect_identical(
    unique(all_groups$contrast),
    c("b vs a", "other vs a")
  )
  res <- test_usage(iso, group = "condition", levels = c("a", "b"))
  testthat::expect_gt(nrow(res), 0L)
  testthat::expect_identical(unique(res$contrast), "b vs a")
  one_gene <- test_usage(x$iso, group = "condition", gene_name = "GENE1")
  testthat::expect_identical(unique(one_gene$gene_id), "g1")
  testthat::expect_error(
    test_usage(x$iso, group = "condition", min_count = 1e6),
    class = "isoformic_no_usage_features"
  )
  testthat::expect_error(
    test_usage(x$iso, group = "condition", assay = "nope"),
    class = "isoformic_assay_not_found"
  )
  testthat::expect_error(
    test_usage(x$iso, group = "condition", method = "other")
  )
  # With replicates in the object the counts are divided by the overdispersion
  draws <- lapply(1:10, function(k) {
    m <- assay(x$iso, "counts")
    m[] <- stats::rpois(length(m), as.vector(m) + 1)
    m
  })
  attr(draws, "samp_type") <- "gibbs"
  x$iso@metadata$inf_reps <- draws
  scaled <- test_usage(x$iso, group = "condition")
  testthat::expect_true(unique(scaled$scaled))
  unscaled <- test_usage(x$iso, group = "condition", scale = FALSE)
  testthat::expect_false(unique(unscaled$scaled))
})

testthat::test_that("filter_usage keeps genes with two expressed transcripts", {
  counts <- rbind(
    a1 = c(50, 60, 55, 40),
    a2 = c(30, 25, 40, 35),
    b1 = c(100, 90, 80, 95),
    b2 = c(0, 1, 0, 1),
    c1 = c(40, 50, 60, 70)
  )
  gene_of <- c("A", "A", "B", "B", "C")
  keep <- filter_usage(counts, gene_of, factor(c("x", "x", "y", "y")), 10, 0.1)
  testthat::expect_identical(unname(keep), c(TRUE, TRUE, FALSE, FALSE, FALSE))
  testthat::expect_equal(simes_pvalue(c(0.01, 0.04, 0.5)), 0.03)
  testthat::expect_true(is.na(simes_pvalue(NA_real_)))
})

testthat::test_that("test_usage tests every contrast of several groups with every backend", {
  x <- usage_experiment(n_groups = 3L)
  found <- function(res, contrast) {
    genes <- unique(res[res$contrast == contrast, c("gene_id", "gene_qvalue")])
    genes$gene_id[genes$gene_qvalue < 0.05]
  }
  for (method in c("diffsplice", "satuRn", "DEXSeq")) {
    testthat::skip_if_not_installed(switch(
      method,
      diffsplice = "edgeR",
      satuRn = "satuRn",
      DEXSeq = "DEXSeq"
    ))
    res <- suppressWarnings(
      test_usage(x$iso, group = "condition", method = method)
    )
    # Default: every level against the first
    testthat::expect_identical(unique(res$contrast), c("b vs a", "c vs a"))
    b_vs_a <- found(res, "b vs a")
    c_vs_a <- found(res, "c vs a")
    testthat::expect_gte(sum(x$switch_genes %in% b_vs_a), 6L)
    testthat::expect_lte(sum(x$c_genes %in% b_vs_a), 2L)
    testthat::expect_gte(sum(x$c_genes %in% c_vs_a), 6L)
    testthat::expect_lte(sum(x$switch_genes %in% c_vs_a), 2L)
    # Directions within each contrast
    lfc <- function(contrast, id) {
      res$usage_log2FC[res$contrast == contrast & res$transcript_id == id]
    }
    testthat::expect_lt(lfc("b vs a", "g1t2"), 0)
    testthat::expect_lt(lfc("c vs a", "g9t1"), 0)
    testthat::expect_gt(lfc("c vs a", "g9t2"), 0)
    # Adjustments are within each contrast
    for (nm in unique(res$contrast)) {
      part <- res[res$contrast == nm, ]
      testthat::expect_equal(
        part$qvalue,
        stats::p.adjust(part$pvalue, method = "BH")
      )
    }
  }
})

testthat::test_that("test_usage takes pairwise and explicit contrasts", {
  testthat::skip_if_not_installed("edgeR")
  x <- usage_experiment(n_groups = 3L)
  pairwise <- test_usage(x$iso, group = "condition", contrasts = "pairwise")
  testthat::expect_setequal(
    unique(pairwise$contrast),
    c("b vs a", "c vs a", "c vs b")
  )
  c_vs_b <- pairwise[pairwise$contrast == "c vs b", ]
  genes <- unique(c_vs_b[c_vs_b$gene_qvalue < 0.05, "gene_id"])$gene_id
  # Both sets of genes differ between b and c
  testthat::expect_gte(sum(c(x$switch_genes, x$c_genes) %in% genes), 12L)

  one <- test_usage(x$iso, group = "condition", contrasts = c("c", "b"))
  testthat::expect_identical(unique(one$contrast), "c vs b")
  testthat::expect_equal(
    one$usage_log2FC[order(one$transcript_id)],
    c_vs_b$usage_log2FC[order(c_vs_b$transcript_id)]
  )
  two <- test_usage(
    x$iso,
    group = "condition",
    contrasts = list(c("b", "a"), c("c", "b"))
  )
  testthat::expect_identical(unique(two$contrast), c("b vs a", "c vs b"))

  # Only some of the levels
  sub <- test_usage(x$iso, group = "condition", levels = c("a", "c"))
  testthat::expect_identical(unique(sub$contrast), "c vs a")

  for (bad in list("all", c("a", "z"), c("a", "a"), list(c("a", "b", "c")))) {
    testthat::expect_error(
      test_usage(x$iso, group = "condition", contrasts = bad),
      class = "isoformic_invalid_contrasts"
    )
  }
  testthat::expect_error(
    test_usage(x$iso, group = "condition", levels = "a"),
    class = "isoformic_group_levels"
  )
})

testthat::test_that("the omnibus test finds a difference among any of the groups", {
  x <- usage_experiment(n_groups = 3L)
  for (method in c("diffsplice", "satuRn", "DEXSeq")) {
    testthat::skip_if_not_installed(switch(
      method,
      diffsplice = "edgeR",
      satuRn = "satuRn",
      DEXSeq = "DEXSeq"
    ))
    res <- suppressWarnings(test_usage(
      x$iso,
      group = "condition",
      contrasts = "omnibus",
      method = method
    ))
    testthat::expect_identical(unique(res$contrast), "omnibus")
    testthat::expect_true(all(is.na(res$usage_log2FC)))
    testthat::expect_false(anyDuplicated(res$transcript_id) > 0L)
    genes <- unique(res[res$gene_qvalue < 0.05, "gene_id"])$gene_id
    # Both switch sets differ among the groups
    testthat::expect_gte(sum(c(x$switch_genes, x$c_genes) %in% genes), 12L)
    testthat::expect_lte(sum(!genes %in% c(x$switch_genes, x$c_genes)), 4L)
  }
})

testthat::test_that("with two groups the omnibus test is the single contrast", {
  x <- usage_experiment()
  for (method in c("diffsplice", "DEXSeq")) {
    testthat::skip_if_not_installed(
      if (method == "DEXSeq") "DEXSeq" else "edgeR"
    )
    one <- suppressWarnings(test_usage(
      x$iso,
      group = "condition",
      method = method
    ))
    omni <- suppressWarnings(test_usage(
      x$iso,
      group = "condition",
      contrasts = "omnibus",
      method = method
    ))
    one <- one[order(one$transcript_id), ]
    omni <- omni[order(omni$transcript_id), ]
    testthat::expect_identical(omni$transcript_id, one$transcript_id)
    testthat::expect_equal(omni$pvalue, one$pvalue, tolerance = 1e-8)
  }
})

testthat::test_that("omnibus_usage combines the contrasts with Simes", {
  res <- tibble::tibble(
    transcript_id = c("t1", "t1", "t2", "t2"),
    contrast = c("b vs a", "c vs a", "b vs a", "c vs a"),
    usage_log2FC = c(1, -1, 0, 0),
    pvalue = c(0.01, 0.5, 0.2, 0.3)
  )
  omni <- omnibus_usage(res)
  testthat::expect_equal(omni$pvalue, c(0.02, 0.3))
  testthat::expect_true(all(is.na(omni$usage_log2FC)))
})
