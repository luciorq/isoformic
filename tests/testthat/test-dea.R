testthat::test_that("de_tx<- and de_gene<- standardise column names and validate IDs", {
  iso <- mini_experiment(with_dea = FALSE)
  de_tx(iso) <- data.frame(
    transcript_id = c("tx1", "tx2"),
    log2FoldChange = c(1, -1),
    PValue = c(0.01, 0.2),
    FDR = c(0.02, 0.3)
  )
  det <- de_tx(iso)
  testthat::expect_equal(
    colnames(det),
    c("transcript_id", "log2FC", "pvalue", "qvalue", "det_sig")
  )
  testthat::expect_equal(det$det_sig[det$transcript_id == "tx1"], "yes")

  # row names as IDs, GENCODE suffixes, qvalue computed when missing
  deg <- data.frame(logFC = c(0.5, 2), pvalue = c(0.5, 0.001))
  rownames(deg) <- c("gene1", "gene2")
  de_gene(iso) <- deg
  res <- de_gene(iso)
  testthat::expect_true("qvalue" %in% colnames(res))
  testthat::expect_equal(
    res$qvalue[res$gene_id == "gene2"],
    stats::p.adjust(c(0.5, 0.001), "BH")[2]
  )

  testthat::expect_error(
    {
      de_tx(iso) <- data.frame(transcript_id = "nope", log2FC = 1, pvalue = 0.1)
    },
    class = "isoformic_invalid_assay_transcript_ids"
  )
  testthat::expect_error(
    {
      de_tx(iso) <- data.frame(transcript_id = "tx1", pvalue = 0.1)
    },
    class = "isoformic_dea_invalid_det"
  )
  testthat::expect_error(
    {
      de_gene(iso) <- data.frame(log2FC = 1, pvalue = 0.1)
    },
    class = "isoformic_invalid_dea"
  )
  testthat::expect_error(
    de_tx(mini_experiment(with_dea = FALSE)),
    class = "isoformic_missing_dea_results"
  )
})

testthat::test_that("cutoffs drive significance flags and can be overridden per call", {
  iso <- mini_experiment()
  testthat::expect_equal(
    de_cutoffs(iso),
    list(log2fc = 1, pvalue = 0.05, use_fdr = TRUE)
  )
  testthat::expect_equal(de_tx(iso)$det_sig, c("yes", "yes", "no"))
  testthat::expect_equal(
    de_tx(iso, log2fc_cutoff = 1.8)$det_sig,
    c("yes", "no", "no")
  )
  testthat::expect_equal(nrow(de_tx(iso, sig_only = TRUE)), 2L)
  de_cutoffs(iso) <- list(log2fc = 1.8, use_fdr = FALSE)
  testthat::expect_equal(de_cutoffs(iso)$log2fc, 1.8)
  testthat::expect_equal(de_tx(iso)$det_sig, c("yes", "no", "no"))
  testthat::expect_error(
    {
      de_cutoffs(iso) <- list(nope = 1)
    },
    class = "isoformic_invalid_cutoffs"
  )
  testthat::expect_error(
    {
      de_cutoffs(iso) <- list(pvalue = 2)
    },
    class = "isoformic_invalid_cutoffs"
  )
  testthat::expect_equal(resolve_cutoffs(list(log2FC = 2))$log2fc, 2)
})

testthat::test_that("de_table combines gene and transcript results in long format", {
  iso <- mini_experiment()
  tbl <- de_table(iso)
  testthat::expect_equal(nrow(tbl), 5L)
  testthat::expect_equal(sum(tbl$feature_type == "gene"), 2L)
  testthat::expect_true(all(
    c(
      "feature_id",
      "feature_name",
      "feature_type",
      "gene_id",
      "gene_name",
      "log2FC",
      "is_de"
    ) %in%
      colnames(tbl)
  ))
  testthat::expect_equal(tbl$is_de[tbl$feature_id == "tx3"], "no")
  testthat::expect_equal(tbl$feature_name[tbl$feature_id == "gene1"], "GENEA")
})

testthat::test_that("join_DEG_DET and is_deg_sig keep working on loose tables", {
  deg_table <- data.frame(
    gene_id = c("gene1", "gene2"),
    gene_name = c("GeneA", "GeneB"),
    log2FC = c(1.5, -2.0),
    pvalue = c(0.01, 0.04)
  )
  det_table <- data.frame(
    transcript_id = c("tx1", "tx2", "tx3", "tx4"),
    transcript_name = c("T1", "T2", "T3", "T4"),
    transcript_type = c(
      "protein_coding",
      "lncRNA",
      "processed_transcript",
      "TEC"
    ),
    log2FC = c(1.2, -1.8, 0.5, 3),
    pvalue = c(0.02, 0.03, 0.2, 0.001)
  )
  res <- join_DEG_DET(deg_table, det_table, logfc_cut = 1, pval_cut = 0.05)
  testthat::expect_equal(nrow(res), 5L) # TEC transcript filtered out
  testthat::expect_equal(res$is_de, c("yes", "yes", "yes", "yes", "no"))
  testthat::expect_false("tx4" %in% res$feature_id)
  all_types <- join_DEG_DET(deg_table, det_table, transcript_types = NULL)
  testthat::expect_equal(nrow(all_types), 6L)

  flagged <- is_deg_sig(c("T1", "T3"), det_table)
  testthat::expect_equal(flagged$DEG_sig, c("YES", "NO", "YES", "NO"))
  testthat::expect_equal(flagged$transcript_name, det_table$transcript_name)
})

testthat::test_that("de_table keeps the annotation names when the results carry their own", {
  iso <- mini_experiment(with_dea = FALSE)
  de_tx(iso) <- data.frame(
    transcript_id = c("tx1", "tx2"),
    transcript_name = c("mine1", "mine2"),
    gene_id = c("gene1", "gene1"),
    gene_name = c("other", "other"),
    log2FC = c(1, -1),
    pvalue = c(0.01, 0.2),
    qvalue = c(0.02, 0.3)
  )
  de_gene(iso) <- data.frame(
    gene_id = c("gene1", "gene2"),
    gene_name = c("mine", "mine"),
    log2FC = c(0.5, 2),
    pvalue = c(0.5, 0.001),
    qvalue = c(0.6, 0.01)
  )
  res <- de_table(iso)
  testthat::expect_false(any(grepl("\\.(x|y)$", colnames(res))))
  genes <- iso@annot_data_genes
  testthat::expect_identical(
    res$gene_name[res$feature_id == "gene2"],
    genes$gene_name[genes$gene_id == "gene2"]
  )
  testthat::expect_no_error(plot_log2fc(iso, feature = genes$gene_name[1L]))
})

testthat::test_that("join_DEG_DET and the stored cutoffs agree at the boundary", {
  boundary <- data.frame(
    gene_id = "gene1",
    gene_name = "GeneA",
    log2FC = 1,
    pvalue = 0.05,
    qvalue = 0.5
  )
  det_table <- data.frame(
    transcript_id = "tx1",
    transcript_name = "T1",
    transcript_type = "protein_coding",
    log2FC = 1,
    pvalue = 0.05,
    qvalue = 0.5
  )
  res <- join_DEG_DET(boundary, det_table, logfc_cut = 1, pval_cut = 0.05)
  testthat::expect_identical(res$is_de, c("yes", "yes"))
  fdr <- join_DEG_DET(
    boundary,
    det_table,
    logfc_cut = 1,
    pval_cut = 0.05,
    use_fdr = TRUE
  )
  testthat::expect_identical(fdr$is_de, c("no", "no"))

  iso <- mini_experiment(with_dea = FALSE)
  de_tx(iso) <- data.frame(
    transcript_id = "tx1",
    log2FC = 1,
    pvalue = 0.05,
    qvalue = 0.05
  )
  testthat::expect_identical(de_tx(iso)$det_sig, "yes")
})
