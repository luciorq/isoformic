testthat::test_that("isoform_fraction divides by the gene total", {
  iso <- mini_experiment(with_dea = FALSE)
  frac <- isoform_fraction(iso)
  testthat::expect_equal(dim(frac), dim(iso))
  testthat::expect_equal(frac["tx3", ], c(s1 = 1, s2 = 1, s3 = 1, s4 = 1))
  testthat::expect_equal(unname(frac["tx1", "s1"]), 10 / 30)
  testthat::expect_equal(unname(colSums(frac[c("tx1", "tx2"), ])), rep(1, 4))
  zero <- mini_assay()
  zero["tx3", ] <- 0
  assay(iso, "abundance") <- zero
  testthat::expect_true(all(is.na(isoform_fraction(iso)["tx3", ])))
  testthat::expect_error(
    isoform_fraction(iso, assay = "nope"),
    class = "isoformic_assay_not_found"
  )
})

testthat::test_that("tx_diversity summarises isoform usage per gene", {
  iso <- mini_experiment(with_dea = FALSE)
  res <- tx_diversity(iso, min_expr = 5, min_samples = 2)
  testthat::expect_equal(nrow(res), 2L)
  a <- res[res$gene_name == "GENEA", ]
  testthat::expect_equal(a$n_tx_annotated, 2L)
  testthat::expect_equal(a$n_tx_assayed, 2L)
  testthat::expect_equal(a$n_tx_expressed, 2L)
  testthat::expect_equal(a$n_tx_types, 2L)
  testthat::expect_equal(a$tx_types, "protein_coding;retained_intron")
  testthat::expect_equal(a$dominant_transcript_id, "tx2")
  testthat::expect_true(a$effective_isoforms > 1 && a$effective_isoforms <= 2)
  b <- res[res$gene_name == "GENEB", ]
  testthat::expect_equal(b$effective_isoforms, 1)
  testthat::expect_equal(b$dominant_fraction, 1)
  testthat::expect_equal(nrow(tx_diversity(iso, gene_name = "GENEB")), 1L)
})

testthat::test_that("vectorised diversity helpers agree with a naive implementation", {
  iso <- example_isoformic()
  frac <- isoform_fraction(iso)
  rd <- row_data(iso)
  mf <- rowMeans(frac, na.rm = TRUE)
  mf[is.nan(mf)] <- NA
  naive_eff <- function(p) {
    p <- p[!is.na(p) & p > 0]
    if (!length(p)) {
      return(NA_real_)
    }
    p <- p / sum(p)
    exp(-sum(p * log(p)))
  }
  naive <- tibble::tibble(gene_id = rd$gene_id, tx = rownames(frac), mf = mf) |>
    dplyr::group_by(gene_id) |>
    dplyr::summarise(
      dom = tx[if (all(is.na(mf))) 1L else which.max(mf)],
      eff = naive_eff(mf),
      .groups = "drop"
    )
  res <- tx_diversity(iso)
  cmp <- dplyr::inner_join(naive, res, by = "gene_id")
  testthat::expect_equal(nrow(cmp), nrow(res))
  testthat::expect_equal(cmp$dom, cmp$dominant_transcript_id)
  testthat::expect_equal(cmp$eff, cmp$effective_isoforms)
})

testthat::test_that("dominant_isoform works overall and per group", {
  iso <- mini_experiment(with_dea = FALSE)
  overall <- dominant_isoform(iso)
  testthat::expect_equal(nrow(overall), 2L)
  testthat::expect_false("dominant_switch" %in% colnames(overall))
  testthat::expect_equal(
    overall$dominant_transcript_id[overall$gene_id == "gene1"],
    "tx2"
  )

  m <- mini_assay()
  m["tx1", c("s3", "s4")] <- c(100, 100) # tx1 dominates in condition b
  assay(iso, "abundance") <- m
  by_group <- dominant_isoform(iso, group = "condition")
  testthat::expect_equal(nrow(by_group), 4L)
  testthat::expect_true("condition" %in% colnames(by_group))
  g1 <- by_group[by_group$gene_id == "gene1", ]
  testthat::expect_equal(g1$dominant_transcript_id[g1$condition == "a"], "tx2")
  testthat::expect_equal(g1$dominant_transcript_id[g1$condition == "b"], "tx1")
  testthat::expect_true(all(g1$dominant_switch))
  testthat::expect_false(any(by_group$dominant_switch[
    by_group$gene_id == "gene2"
  ]))
  testthat::expect_error(
    dominant_isoform(iso, group = "nope"),
    class = "isoformic_group_not_found"
  )
})

testthat::test_that("isoform_switches classifies genes from DE results", {
  iso <- mini_experiment()
  res <- isoform_switches(iso)
  testthat::expect_equal(nrow(res), 1L)
  testthat::expect_equal(res$gene_name, "GENEA")
  testthat::expect_equal(res$n_det_sig, 2L)
  testthat::expect_equal(res$n_det_up, 1L)
  testthat::expect_equal(res$n_det_down, 1L)
  testthat::expect_equal(res$switch_type, "opposite_direction")
  testthat::expect_equal(res$deg_sig, "yes")
  testthat::expect_equal(res$det_sig_types, "protein_coding;retained_intron")

  iso@dea$deg <- NULL
  res2 <- isoform_switches(iso)
  testthat::expect_true(is.na(res2$deg_sig))
  de_tx(iso) <- data.frame(
    transcript_id = c("tx1", "tx2"),
    log2FC = c(2, 1.5),
    pvalue = c(0.001, 0.01)
  )
  testthat::expect_equal(isoform_switches(iso)$switch_type, "det_not_deg")
})

testthat::test_that("exon_diversity counts constitutive and alternative exons", {
  iso <- mini_experiment(with_dea = FALSE)
  res <- exon_diversity(iso)
  a <- res[res$gene_name == "GENEA", ]
  testthat::expect_equal(a$n_tx, 2L)
  testthat::expect_equal(a$n_exons, 3L)
  testthat::expect_equal(a$n_constitutive, 0L)
  testthat::expect_equal(a$n_alternative, 3L)
  testthat::expect_equal(a$max_exons_per_tx, 2L)
  b <- res[res$gene_name == "GENEB", ]
  testthat::expect_equal(b$n_constitutive, 1L)
  testthat::expect_equal(nrow(exon_diversity(iso, gene_name = "GENEB")), 1L)
})

testthat::test_that("dominant_isoform reports whether the dominant isoform is canonical", {
  annot_dir <- mini_annotation("mini_tags.gtf")
  m <- matrix(c(9, 1, 8, 2), 2, dimnames = list(c("tx1", "tx2"), c("a", "b")))
  iso <- IsoformicExperiment(
    annot_path = annot_dir,
    assays = list(abundance = m)
  )
  res <- dominant_isoform(iso)
  testthat::expect_equal(res$dominant_transcript_id, "tx1")
  testthat::expect_true(res$dominant_is_canonical)
  m2 <- m[c("tx2", "tx1"), ]
  rownames(m2) <- c("tx1", "tx2")
  assay(iso, "abundance") <- m2
  testthat::expect_false(dominant_isoform(iso)$dominant_is_canonical)
})

testthat::test_that("haplotype_summary compares haplotypes of the same gene", {
  gff <- testthat::test_path("testdata", "mini_liftoff.gff3")
  merged <- suppressMessages(merge_annotations(
    c(mat = gff, pat = gff),
    output_path = withr::local_tempdir(),
    annot_name = "diploid"
  ))
  ids <- c(
    "XM_047436352.1_mat",
    "XM_047436352.1_pat",
    "CHM13_T0000001_mat",
    "CHM13_T0000001_pat"
  )
  m <- matrix(
    c(
      30,
      10,
      4,
      4, # sample a
      20,
      20,
      8,
      0 # sample b
    ),
    nrow = 4,
    dimnames = list(ids, c("a", "b"))
  )
  iso <- IsoformicExperiment(
    annot_path = merged,
    assays = list(abundance = m),
    col_data = data.frame(sample_id = c("a", "b"), condition = c("x", "y"))
  )
  res <- haplotype_summary(iso)
  testthat::expect_equal(nrow(res), 4L)
  loc <- res[res$gene_key == "LOC112268260", ]
  testthat::expect_equal(loc$haplotype, c("mat", "pat"))
  testthat::expect_equal(loc$haplotype_fraction, c(25 / 40, 15 / 40))
  testthat::expect_false(any(res$haplotype_isoform_switch))
  cat_gene <- res[res$gene_key == "ENSG00000235146.2", ]
  testthat::expect_equal(cat_gene$gene_name, c("AC114498.1", "AC114498.1"))

  by_group <- haplotype_summary(iso, group = "condition")
  testthat::expect_true("condition" %in% colnames(by_group))
  testthat::expect_equal(nrow(by_group), 8L)
  y <- by_group[
    by_group$condition == "y" & by_group$gene_key == "ENSG00000235146.2",
  ]
  testthat::expect_equal(y$haplotype_fraction, c(1, 0))
  testthat::expect_error(
    haplotype_summary(mini_experiment()),
    class = "isoformic_not_haplotype_resolved"
  )
})
