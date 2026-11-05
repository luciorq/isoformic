seq_experiment <- function(envir = parent.frame(), file = "mini_cds.gff3") {
  annot <- mini_annotation(file, envir = envir)
  extract_transcript_sequences(
    annot,
    genome_fasta = mini_genome_fasta(envir = envir)
  )
  tx <- read_parquet(fs::path(
    annot,
    "annot_data_transcripts.parquet"
  ))$transcript_id
  mat <- matrix(
    c(10, 5, 2, 1, 3, 10, 4, 3, 2, 3, 2, 1, 8, 6, 3, 1, 2, 9, 7, 3),
    nrow = 5,
    dimnames = list(c("tA1", "tA2", "tA3", "tA4", "tB1"), paste0("s", 1:4))
  )
  if (length(tx) > 5L) {
    extra <- matrix(
      c(4, 2, 1, 3, 1, 2, 2, 3, 1, 2, 4, 5),
      nrow = 3,
      dimnames = list(c("tA5", "tA6", "tB2"), paste0("s", 1:4))
    )
    mat <- rbind(mat, extra)
  }
  IsoformicExperiment(
    experiment_name = "seq",
    annot_path = annot,
    assays = list(abundance = mat, counts = round(mat * 10)),
    col_data = data.frame(
      sample_id = paste0("s", 1:4),
      condition = c("a", "a", "b", "b"),
      stringsAsFactors = FALSE
    )
  )
}

testthat::test_that("sequence-aware productivity uses projected ORFs and NMD classes", {
  iso <- seq_experiment()
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
  testthat::expect_identical(p$source[1:2], c("sequence", "sequence"))
  testthat::expect_match(p$reason[2], "in-frame internal change")
  testthat::expect_identical(
    p$orf_source,
    c("annotation", "annotation", "annotation", "projected", "longest")
  )
  testthat::expect_identical(p$nmd_class[3:4], c("ptc", "ptc"))
  # Without biotypes the sequence rules alone recover the NMD and retained-intron isoforms
  rules <- transcript_productivity(iso, use_biotype = FALSE)
  rules <- rules[match(c("tA3", "tA4"), rules$transcript_id), ]
  testthat::expect_identical(
    as.character(rules$productivity),
    c("unproductive", "unproductive")
  )
  testthat::expect_match(rules$reason[1], "premature stop codon")
  testthat::expect_identical(rules$source, c("sequence", "sequence"))
  # use_sequence = FALSE falls back to the annotation rules
  p_annot <- transcript_productivity(iso, use_sequence = FALSE)
  testthat::expect_true(all(p_annot$source == "annotation"))
  testthat::expect_false("orf_source" %in% colnames(p_annot))
  # A start-proximal stop becomes uncertain
  iso2 <- seq_experiment()
  p3 <- transcript_productivity(iso2, use_biotype = FALSE, nmd_rule = 50)
  testthat::expect_identical(
    as.character(p3$productivity[p3$transcript_id == "tA3"]),
    "unproductive"
  )
})

testthat::test_that("productivity_fractions splits expression by class and mechanism", {
  iso <- seq_experiment()
  pf <- productivity_fractions(iso, group = "condition")
  a <- pf[pf$gene_id == "gA", ]
  testthat::expect_setequal(
    as.character(a$productivity),
    c("productive", "altered_productive", "nmd", "retained_intron")
  )
  testthat::expect_equal(sum(a$fraction), 1)
  testthat::expect_equal(sum(a$fraction_a), 1)
  testthat::expect_equal(sum(a$fraction_b), 1)
  nmd <- a[a$productivity == "nmd", ]
  testthat::expect_equal(nmd$fraction_a, mean(c(2 / 18, 3 / 19)))
  testthat::expect_equal(nmd$delta_fraction, nmd$fraction_b - nmd$fraction_a)
  b <- pf[pf$gene_id == "gB", ]
  testthat::expect_identical(as.character(b$productivity), "noncoding")
  testthat::expect_equal(b$fraction, 1)
})

testthat::test_that("test_productive_fraction fits the two-group quasi-binomial model", {
  iso <- seq_experiment()
  res <- test_productive_fraction(iso, group = "condition")
  testthat::expect_identical(res$gene_id, "gA")
  testthat::expect_equal(res$productive_fraction_a, (150 + 140) / (180 + 190))
  testthat::expect_equal(res$productive_fraction_b, (30 + 30) / (170 + 190))
  testthat::expect_lt(res$pvalue, 0.05)
  testthat::expect_equal(res$df, 2)
  testthat::expect_true(res$log_odds_ratio < 0)
  testthat::expect_equal(res$padj, res$pvalue)
  testthat::expect_error(
    test_productive_fraction(iso, group = "sample_id"),
    class = "isoformic_group_levels"
  )
  # The closed form matches glm(quasibinomial) on the pooled counts
  counts <- assay(iso, "counts")
  prod <- c(TRUE, TRUE, FALSE, FALSE)
  y <- colSums(counts[1:4, ][prod, ])
  n <- colSums(counts[1:4, ])
  fit <- stats::glm(
    cbind(y, n - y) ~ factor(c("a", "a", "b", "b")),
    family = stats::quasibinomial()
  )
  testthat::expect_equal(
    unname(stats::coef(fit)[2]),
    res$log_odds_ratio,
    tolerance = 0.05
  )
})

testthat::test_that("two_group_proportion_test handles degenerate input", {
  num <- matrix(c(0, 0, 5, 5), nrow = 1, dimnames = list("g", paste0("s", 1:4)))
  den <- matrix(
    c(10, 10, 10, 10),
    nrow = 1,
    dimnames = list("g", paste0("s", 1:4))
  )
  groups <- factor(c("a", "a", "b", "b"))
  res <- two_group_proportion_test(num, den, groups)
  testthat::expect_equal(res$proportion_a, 0)
  testthat::expect_equal(res$proportion_b, 0.5)
  testthat::expect_true(is.finite(res$log_odds_ratio))
  testthat::expect_equal(res$dispersion, 1)
  zero <- two_group_proportion_test(num * 0, den * 0, groups)
  testthat::expect_true(is.na(zero$pvalue))
  testthat::expect_error(
    two_group_proportion_test(
      num,
      den,
      factor(c("a", "b", "c", "c"))
    ),
    class = "isoformic_group_levels"
  )
})

testthat::test_that("productive_fraction honours min_expr, inferential replicates and a reference", {
  iso <- seq_experiment()
  base <- productive_fraction(iso, group = "condition")
  floored <- productive_fraction(iso, group = "condition", min_expr = 3)
  # In sample s1 only tA1 (10) and tA2 (5) pass the floor: fully productive
  mat <- productive_fraction(iso, as_matrix = TRUE, min_expr = 3)
  testthat::expect_equal(unname(mat["gA", "s1"]), 1)
  testthat::expect_true(
    floored$productive_fraction_a > base$productive_fraction_a
  )
  set.seed(1)
  reps <- lapply(1:20, function(i) {
    m <- assay(iso, "abundance")
    m * matrix(stats::runif(length(m), 0.8, 1.2), nrow = nrow(m))
  })
  with_reps <- productive_fraction(iso, group = "condition", inf_reps = reps)
  testthat::expect_true(all(
    c("productive_fraction_sd", "delta_productive_fraction_sd", "delta_z") %in%
      colnames(with_reps)
  ))
  testthat::expect_true(with_reps$productive_fraction_sd > 0)
  testthat::expect_true(abs(with_reps$delta_z) > 1)
  iso@metadata$inf_reps <- reps
  from_meta <- productive_fraction(iso, group = "condition")
  testthat::expect_equal(
    from_meta$productive_fraction_sd,
    with_reps$productive_fraction_sd
  )
  testthat::expect_error(
    productive_fraction(iso, inf_reps = list(matrix(1, 2, 2))),
    class = "isoformic_invalid_inf_reps"
  )
  testthat::expect_error(
    productive_fraction(iso, inf_reps = list()),
    class = "isoformic_invalid_inf_reps"
  )
  reference <- tibble::tibble(
    gene_id = "gA",
    reference_tissue = "Synthetic",
    reference_productive_fraction = 0.9,
    reference_dominant_productivity = "productive"
  )
  with_ref <- productive_fraction(
    iso,
    group = "condition",
    reference = reference
  )
  testthat::expect_equal(with_ref$reference_productive_fraction, 0.9)
  testthat::expect_equal(
    with_ref$delta_vs_reference,
    with_ref$productive_fraction_b - 0.9
  )
  no_group <- productive_fraction(iso, reference = reference)
  testthat::expect_equal(
    no_group$delta_vs_reference,
    no_group$productive_fraction - 0.9
  )
})

testthat::test_that("dominant_isoform honours min_expr", {
  iso <- seq_experiment()
  d <- dominant_isoform(iso, group = "condition", min_expr = 100)
  testthat::expect_true(all(is.na(d$dominant_fraction[d$gene_id == "gA"])))
  d2 <- dominant_isoform(iso, group = "condition", min_expr = 0)
  testthat::expect_false(anyNA(d2$dominant_fraction[d2$gene_id == "gA"]))
})

testthat::test_that("benchmark_productivity measures the rules against biotypes", {
  iso <- seq_experiment()
  bm <- benchmark_productivity(iso)
  testthat::expect_named(bm, c("summary", "confusion"))
  s <- bm$summary
  testthat::expect_equal(
    s$value[s$metric == "nmd_biotype_called_unproductive"],
    1
  )
  testthat::expect_equal(
    s$value[s$metric == "retained_intron_called_unproductive"],
    1
  )
  testthat::expect_equal(s$value[s$metric == "canonical_called_productive"], 1)
  testthat::expect_equal(s$n[s$metric == "canonical_called_productive"], 1L)
  testthat::expect_equal(s$value[s$metric == "lncrna_called_noncoding"], 1)
  testthat::expect_identical(unique(s$source), "sequence")
  testthat::expect_true(all(
    c("transcript_type", "productivity", "n") %in% colnames(bm$confusion)
  ))
  bm2 <- benchmark_productivity(iso, use_sequence = FALSE)
  testthat::expect_identical(unique(bm2$summary$source), "annotation")
})

testthat::test_that("evidence columns override computed calls in priority order", {
  iso <- seq_experiment()
  annot <- iso@annot_path
  add_annotation_columns(
    annot,
    data.frame(
      transcript_id = c("tA1", "tA3", "tA4"),
      appris_principal = c(TRUE, TRUE, FALSE)
    )
  )
  clear_cache(iso)
  p <- transcript_productivity(iso)
  testthat::expect_identical(p$source[p$transcript_id == "tA1"], "evidence")
  testthat::expect_identical(
    as.character(p$productivity[p$transcript_id == "tA3"]),
    "productive"
  )
  add_annotation_columns(
    annot,
    data.frame(transcript_id = "tA3", nmd_validated_exon = TRUE)
  )
  clear_cache(iso)
  p <- transcript_productivity(iso)
  testthat::expect_identical(
    as.character(p$productivity[p$transcript_id == "tA3"]),
    "unproductive"
  )
  testthat::expect_match(p$reason[p$transcript_id == "tA3"], "validated")
  add_annotation_columns(
    annot,
    data.frame(transcript_id = "tA3", productivity_evidence = "uncertain")
  )
  clear_cache(iso)
  p <- transcript_productivity(iso)
  testthat::expect_identical(
    as.character(p$productivity[p$transcript_id == "tA3"]),
    "uncertain"
  )
  p_off <- transcript_productivity(iso, use_evidence = FALSE)
  testthat::expect_false(any(p_off$source == "evidence"))
})

testthat::test_that("test_productive_fraction contrasts two of three groups", {
  iso <- seq_experiment()
  sample_table <- col_data(iso)
  sample_table$condition <- c("a", "a", "b", "other")
  col_data(iso) <- sample_table
  testthat::expect_error(
    test_productive_fraction(iso, group = "condition"),
    class = "isoformic_group_levels"
  )
  res <- test_productive_fraction(
    iso,
    group = "condition",
    levels = c("a", "b")
  )
  testthat::expect_identical(res$gene_id, "gA")
  testthat::expect_equal(res$productive_fraction_a, (150 + 140) / (180 + 190))
  testthat::expect_equal(res$productive_fraction_b, 30 / 170)
})
