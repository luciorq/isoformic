cds_experiment <- function(envir = parent.frame()) {
  annot <- mini_annotation("mini_cds.gff3", envir = envir)
  ids <- c("tA1", "tA2", "tA3", "tA4", "tB1")
  m <- matrix(
    c(
      50,
      20,
      10,
      5,
      8, # sample a (x)
      40,
      25,
      20,
      6,
      9, # sample b (x)
      10,
      10,
      60,
      30,
      7, # sample c (y)
      12,
      8,
      70,
      25,
      6 # sample d (y)
    ),
    nrow = 5,
    dimnames = list(ids, c("a", "b", "c", "d"))
  )
  IsoformicExperiment(
    annot_path = annot,
    assays = list(abundance = m),
    col_data = data.frame(
      sample_id = c("a", "b", "c", "d"),
      condition = c("x", "x", "y", "y")
    )
  )
}

testthat::test_that("CDS features are aggregated by the parser", {
  tx <- prepare_annotation(mini_annotation("mini_cds.gff3"))$transcripts
  tx_a1 <- tx[tx$transcript_id == "tA1", ]
  testthat::expect_equal(tx_a1$cds_start, 1100L)
  testthat::expect_equal(tx_a1$cds_end, 8300L)
  testthat::expect_equal(tx_a1$cds_length, 101L + 201L + 201L + 301L)
  testthat::expect_true(is.na(tx$cds_length[tx$transcript_id == "tA4"]))
  testthat::expect_true(tx_a1$is_canonical)
})

testthat::test_that("transcript_features computes ORF, UTR, NMD, ends and events", {
  iso <- cds_experiment()
  f <- transcript_features(iso)
  f <- f[match(c("tA1", "tA2", "tA3", "tA4", "tB1"), f$transcript_id), ]
  testthat::expect_equal(f$n_exons, c(4L, 3L, 5L, 3L, 2L))
  testthat::expect_equal(f$tx_length[1], 201L + 201L + 201L + 1001L)
  testthat::expect_equal(f$has_cds, c(TRUE, TRUE, TRUE, FALSE, FALSE))
  testthat::expect_equal(f$orf_length[1:3], c(804L, 603L, 353L))
  # tA1: 5' UTR = 100 nt (1000..1099), stop at the last exon -> no NMD
  testthat::expect_equal(f$utr5_length[1], 100)
  testthat::expect_equal(f$utr3_length[1], 9000 - 8300)
  testthat::expect_false(f$nmd_predicted[1])
  # tA3: ORF ends in exon 3 of 5, far upstream of the last junction -> NMD
  testthat::expect_true(f$nmd_predicted[3])
  testthat::expect_true(f$stop_to_last_junction[3] > 50)
  testthat::expect_equal(f$canonical_transcript_id[1:4], rep("tA1", 4))
  testthat::expect_equal(f$orf_fraction_of_canonical[2], 603 / 804)
  testthat::expect_true(f$cds_differs_from_canonical[2])
  testthat::expect_false(f$cds_differs_from_canonical[1])
  # ends: minus-strand gene has its TSS at the highest coordinate
  testthat::expect_equal(f$tss_pos[5], 21000L)
  testthat::expect_equal(f$tes_pos[5], 20000L)
  testthat::expect_equal(f$tss_group[1:4], rep(1L, 4))
  # splicing events relative to tA1
  testthat::expect_equal(f$n_skipped_exons[2], 1L) # exon 5000-5200 skipped
  testthat::expect_equal(f$n_novel_exons[3], 1L) # 4000-4100 is a novel (poison) exon
  testthat::expect_equal(f$n_retained_introns[4], 1L) # one exon spanning the intron between canonical exons 2 and 3
  testthat::expect_false(f$alt_first_exon[2])
  testthat::expect_false(f$alt_last_exon[2])
  ev <- splicing_events(iso, gene_name = "GENEA")
  testthat::expect_equal(nrow(ev), 4L)
  # annotation-directory method gives the same result
  f2 <- transcript_features(iso@annot_path)
  testthat::expect_equal(sort(f2$transcript_id), sort(f$transcript_id))
})

testthat::test_that("transcript_productivity grades isoforms and productive_fraction summarises genes", {
  iso <- cds_experiment()
  p <- transcript_productivity(iso)
  p <- p[match(c("tA1", "tA2", "tA3", "tA4", "tB1"), p$transcript_id), ]
  testthat::expect_equal(
    as.character(p$productivity),
    c(
      "productive",
      "altered_productive",
      "unproductive",
      "unproductive",
      "noncoding"
    )
  )
  testthat::expect_match(p$reason[3], "stop codon|biotype")
  testthat::expect_true(all(p$source == "annotation"))
  # biotype can be switched off: tA4 has no ORF in a coding gene -> still unproductive
  p2 <- transcript_productivity(iso, use_biotype = FALSE)
  testthat::expect_equal(
    as.character(p2$productivity[p2$transcript_id == "tA4"]),
    "unproductive"
  )

  pf <- productive_fraction(iso, group = "condition")
  a <- pf[pf$gene_id == "gA", ]
  testthat::expect_equal(
    a$productive_fraction_x,
    mean(c((50 + 20) / 85, (40 + 25) / 91))
  )
  testthat::expect_equal(a$productive_fraction_y, mean(c(20 / 110, 20 / 115)))
  testthat::expect_equal(
    a$delta_productive_fraction,
    a$productive_fraction_y - a$productive_fraction_x
  )
  testthat::expect_false("gB" %in% pf$gene_id) # non-coding genes are excluded
  mat <- productive_fraction(iso, as_matrix = TRUE)
  testthat::expect_equal(dim(mat), c(1L, 4L))
  pf_strict <- productive_fraction(iso, include_altered = FALSE)
  testthat::expect_equal(
    pf_strict$productive_fraction,
    mean(c(50 / 85, 40 / 91, 10 / 110, 12 / 115))
  )

  pe <- poison_exons(iso)
  testthat::expect_equal(nrow(pe), 1L)
  testthat::expect_equal(c(pe$start_pos, pe$end_pos), c(4000L, 4100L))
  testthat::expect_equal(pe$transcript_ids, "tA3")
  testthat::expect_s3_class(plot_genomic_context(iso, "GENEA"), "ggplot")
  testthat::expect_s3_class(
    plot_genomic_context(iso, "GENEA", highlight_poison_exons = FALSE),
    "ggplot"
  )
  testthat::expect_named(
    productivity_palette(),
    c(
      "productive",
      "altered_productive",
      "unproductive",
      "noncoding",
      "uncertain"
    )
  )
})

testthat::test_that("run_productivity_enrichment ranks genes by change in productive fraction", {
  testthat::skip_if_not_installed("fgsea")
  iso <- example_isoformic()
  sets <- read_gmt(system.file(
    "extdata",
    "reactome_subset.gmt",
    package = "isoformic"
  ))
  res <- suppressWarnings(run_productivity_enrichment(
    iso,
    sets,
    group = "condition",
    min_size = 3
  ))
  testthat::expect_true(all(
    c("pathway", "NES", "padj", "experiment") %in% colnames(res)
  ))
  testthat::expect_equal(unique(res$experiment), "productive_fraction")
  testthat::expect_error(
    run_productivity_enrichment(iso[, 1:3], sets, group = "condition"),
    class = "isoformic_group_levels"
  )
})

testthat::test_that("add_annotation_columns attaches evidence and overrides productivity", {
  iso <- cds_experiment()
  evidence <- data.frame(
    transcript_id = c("tA2", "tA3", "zzz"),
    riboseq_translated = c(TRUE, FALSE, NA),
    productivity_evidence = c("productive", "productive", NA)
  )
  testthat::expect_message(
    add_annotation_columns(iso, evidence),
    "Added 2 columns"
  )
  rd <- row_data(iso)
  testthat::expect_true(rd$riboseq_translated[rd$transcript_id == "tA2"])
  p <- transcript_productivity(iso)
  testthat::expect_equal(
    as.character(p$productivity[p$transcript_id == "tA3"]),
    "productive"
  )
  testthat::expect_equal(p$source[p$transcript_id == "tA3"], "evidence")
  testthat::expect_message(
    add_annotation_columns(
      iso,
      data.frame(gene_id = "gA", gtex_dominant = "tA1"),
      level = "gene",
      prefix = "ref_"
    ),
    "ref_gtex_dominant"
  )
  testthat::expect_equal(gene_data(iso)$ref_gtex_dominant, c("tA1", NA))
  testthat::expect_error(
    add_annotation_columns(iso, data.frame(x = 1), by = "nope"),
    class = "isoformic_invalid_evidence"
  )
  testthat::expect_error(
    add_annotation_columns(iso, data.frame(transcript_id = "none", v = 1)),
    class = "isoformic_invalid_evidence"
  )
})

testthat::test_that("FANTOM CAT GTF attributes are parsed and companion tables attached", {
  annot_dir <- mini_annotation("mini_fantom.gtf")
  annot <- prepare_annotation(annot_dir)
  g <- annot$genes[
    match(c("CATG00000000004.1", "ENSG00000999999.1"), annot$genes$gene_id),
  ]
  testthat::expect_equal(
    g$gene_type[g$gene_id == "CATG00000000004.1"],
    "lncRNA_divergent"
  )
  testthat::expect_equal(g$gene_class, c("lncRNA_divergent", "coding_mRNA"))
  testthat::expect_equal(g$gene_category, c("p_lncRNA_divergent", NA))
  testthat::expect_equal(g$dhs_type, c("DHS_promoter", "DHS_promoter"))
  testthat::expect_equal(g$coding_status, c("nonCoding", "coding"))
  tx <- annot$transcripts[
    match(
      c("ENCT00000005542.1", "ENST00000999999.1"),
      annot$transcripts$transcript_id
    ),
  ]
  testthat::expect_equal(
    tx$transcript_type,
    c("lncRNA_divergent", "protein_coding")
  )
  testthat::expect_equal(
    tx$transcript_name,
    c("ENCT00000005542.1", "GENEX-201")
  )
  testthat::expect_equal(tx$tie_score, c(60.31233, 80.5))
  testthat::expect_equal(tx$coding_status, c("nonCoding", "coding"))
  testthat::expect_true(all(tx$transcript_type %in% names(tx_type_palette())))

  mapping <- data.frame(
    geneID = c("CATG00000000004.1", "ENSG00000999999.1"),
    transcriptID = c("ENCT00000005542.1", "ENST00000999999.1"),
    CAGEClusterID = c("chr1:47082699..47082703,+", "chr1:109990..110010,-")
  )
  info <- data.frame(
    geneID = "CATG00000000004.1",
    geneClass = "lncRNA_divergent",
    TIEScore_tier = "tier3_robust",
    splicingType = "type1:unspliced",
    directionality = -0.96
  )
  suppressMessages(add_fantom_cat_tables(annot_dir, mapping, info))
  tx2 <- read_annotation_table(annot_dir, "transcripts")
  tx2 <- tx2[
    match(c("ENCT00000005542.1", "ENST00000999999.1"), tx2$transcript_id),
  ]
  testthat::expect_equal(tx2$cage_cluster_id[1], "chr1:47082699..47082703,+")
  testthat::expect_equal(tx2$cage_cluster_start[1], 47082699L)
  testthat::expect_equal(tx2$cage_cluster_seqid[2], "chr1")
  g2 <- read_annotation_table(annot_dir, "genes")
  g2 <- g2[match(c("CATG00000000004.1", "ENSG00000999999.1"), g2$gene_id), ]
  testthat::expect_equal(g2$tie_score_tier, c("tier3_robust", NA))
  testthat::expect_equal(g2$splicing_type[1], "type1:unspliced")
  testthat::expect_false("gene_class.x" %in% colnames(g2))
  testthat::expect_error(
    add_fantom_cat_tables(annot_dir, data.frame(x = 1)),
    class = "isoformic_invalid_evidence"
  )
})
