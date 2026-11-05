expect_ggplot <- function(p) {
  testthat::expect_s3_class(p, "ggplot")
  built <- ggplot2::ggplot_build(p)
  testthat::expect_true(length(built$data) > 0)
  invisible(built)
}

testthat::test_that("plot_log2fc works on objects and tables, with facets for several genes", {
  iso <- mini_experiment()
  p <- plot_log2fc(iso, "GENEA")
  built <- expect_ggplot(p)
  testthat::expect_equal(nrow(built$data[[2]]), 3L) # gene + 2 transcripts
  p2 <- plot_log2fc(iso, c("GENEA", "GENEB"))
  expect_ggplot(p2)
  testthat::expect_s3_class(p2$facet, "FacetWrap")
  p3 <- plot_log2fc(de_table(iso), "gene2", feature_column = "gene_id")
  expect_ggplot(p3)
  expect_ggplot(plot_log2FC(de_table(iso), "GENEB"))
  testthat::expect_error(
    plot_log2fc(iso, "NOPE"),
    class = "isoformic_empty_plot_data"
  )
  # A requested feature without data does not turn the plot into a facet
  p4 <- plot_log2fc(iso, c("GENEA", "NOPE"))
  testthat::expect_s3_class(p4$facet, "FacetNull")
  testthat::expect_identical(ggplot2::get_labs(p4)$title, "GENEA")
  testthat::expect_error(
    plot_log2fc(iso, "GENEA", feature_column = "nope"),
    class = "isoformic_plot_log2fc_invalid_feature_column"
  )
  testthat::expect_error(
    plot_log2fc(data.frame(x = 1), "a"),
    class = "isoformic_invalid_de_table"
  )
})

testthat::test_that("profile_data and plot_tx_expr summarise expression per group", {
  iso <- mini_experiment()
  prof <- profile_data(iso, gene_name = "GENEA", group = "condition")
  testthat::expect_equal(attr(prof, "group_var"), "condition")
  testthat::expect_equal(levels(prof$group), c("a", "b"))
  gene_a <- prof[prof$feature_type == "gene" & prof$group == "a", ]
  testthat::expect_equal(
    gene_a$mean_expr,
    mean(colSums(mini_assay()[c("tx1", "tx2"), c("s1", "s2")]))
  )
  testthat::expect_equal(prof$is_de[prof$feature_id == "tx1"], c("yes", "yes"))
  testthat::expect_equal(prof$is_de[prof$feature_id == "tx2"], c("yes", "yes"))
  testthat::expect_equal(unique(prof$n_samples), 2L)
  prof_lv <- profile_data(iso, "GENEA", "condition", levels = c("b", "a"))
  testthat::expect_equal(levels(prof_lv$group), c("b", "a"))
  testthat::expect_error(
    profile_data(iso, "NOPE", "condition"),
    class = "isoformic_gene_not_found"
  )

  expect_ggplot(plot_tx_expr(iso, c("GENEA", "GENEB"), group = "condition"))
  expect_ggplot(plot_tx_expr(prof))
  expect_ggplot(plot_tx_expr(
    iso,
    "GENEA",
    group = "condition",
    log_scale = FALSE,
    label = FALSE
  ))
  testthat::expect_error(
    plot_tx_expr(data.frame(x = 1)),
    class = "isoformic_invalid_profile_data"
  )
})

testthat::test_that("prepare_profile_data matches profile_data on loose tables", {
  iso <- mini_experiment()
  prof <- prepare_profile_data(
    txi_transcript = assay(iso, "abundance"),
    sample_metadata = col_data(iso),
    tx_to_gene = tx_annot(iso),
    de_result_gene = de_gene(iso),
    de_result_transcript = de_tx(iso),
    var = "condition",
    gene_name = "GENEA"
  )
  ref <- profile_data(iso, "GENEA", "condition")
  testthat::expect_equal(prof$mean_expr, ref$mean_expr)
  testthat::expect_equal(prof$is_de, ref$is_de)
  testthat::expect_equal(prof$feature_name, ref$feature_name)

  # tximport-style list and data frame input, no DE tables
  df <- as.data.frame(assay(iso, "abundance"))
  df$transcript_id <- rownames(df)
  prof2 <- prepare_profile_data(
    txi_transcript = list(abundance = assay(iso), counts = assay(iso)),
    sample_metadata = col_data(iso),
    tx_to_gene = tx_annot(iso),
    var = "condition"
  )
  testthat::expect_true(all(prof2$is_de == "no"))
  prof3 <- prepare_profile_data(
    df,
    col_data(iso),
    tx_annot(iso),
    var = "condition",
    var_levels = c("b", "a")
  )
  testthat::expect_equal(levels(prof3$group), c("b", "a"))
  testthat::expect_error(
    prepare_profile_data(
      df,
      col_data(iso),
      tx_annot(iso)[, 1:2],
      var = "condition"
    ),
    class = "isoformic_invalid_tx_to_gene"
  )
})

testthat::test_that("plot_genomic_context draws transcript structures", {
  iso <- mini_experiment()
  p <- plot_genomic_context(iso, "GENEA")
  expect_ggplot(p)
  layer_classes <- vapply(p$layers, function(l) class(l$geom)[1], character(1))
  testthat::expect_true("GeomRect" %in% layer_classes)
  testthat::expect_true("GeomSegment" %in% layer_classes)
  p2 <- plot_genomic_context(
    iso,
    c("GENEA", "GENEB"),
    assay = NULL,
    highlight_de = FALSE
  )
  testthat::expect_s3_class(p2$facet, "FacetWrap")
  expect_ggplot(p2)

  # table-based and legacy alias
  ex <- exon_table(iso, "GENEA")
  expect_ggplot(plot_genomic_context(ex))
  legacy <- data.frame(
    tx_id = c("t1", "t1", "t2"),
    exon_left = c(1L, 50L, 10L),
    exon_right = c(20L, 80L, 70L),
    transcript_type = "protein_coding",
    transcript_name = c("T1", "T1", "T2")
  )
  expect_ggplot(plot_tx_context(legacy))
  testthat::expect_error(
    plot_genomic_context(data.frame(a = 1)),
    class = "isoformic_invalid_exon_table"
  )
  testthat::expect_warning(
    plot_genomic_context(iso, c("GENEA", "NOPE")),
    "not found"
  )
})

testthat::test_that("plot_enrichment draws a dot plot", {
  df <- data.frame(
    pathway = rep(c("P1", "P2", "P3"), 2),
    NES = c(2, -1.5, 0.5, 1, -2, 0.2),
    padj = c(0.01, 0.02, 0.5, 0.03, 0.01, 0.6),
    size = c(10, 20, 30, 10, 20, 30),
    experiment = rep(c("protein_coding", "unproductive"), each = 3)
  )
  expect_ggplot(plot_enrichment(df, n = 2))
  expect_ggplot(plot_enrichment(
    df,
    experiments = "unproductive",
    nes_cutoff = 1
  ))
  testthat::expect_error(
    plot_enrichment(df, nes_cutoff = 10),
    class = "isoformic_empty_plot_data"
  )
  testthat::expect_error(
    plot_enrichment(df[, 1:2]),
    class = "isoformic_invalid_enrichment"
  )
})

testthat::test_that("tx_type_palette is a complete named colour vector", {
  pal <- tx_type_palette()
  testthat::expect_type(pal, "character")
  testthat::expect_true(all(grepl("^#[0-9a-fA-F]{6}$", pal)))
  testthat::expect_true(all(
    c(
      "gene",
      "protein_coding",
      "retained_intron",
      "nonsense_mediated_decay"
    ) %in%
      names(pal)
  ))
  testthat::expect_equal(
    complete_palette(pal, c("gene", "new_type"))[["new_type"]],
    "#9e9e9e"
  )
})
