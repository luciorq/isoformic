seq_experiment <- function(envir = parent.frame(), file = "mini_cds.gff3") {
  annot <- mini_annotation(file, envir = envir)
  suppressMessages(extract_transcript_sequences(
    annot,
    genome_fasta = mini_genome_fasta(envir = envir)
  ))
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

testthat::test_that("protein_coding_LoF is an unproductive biotype", {
  testthat::expect_true("protein_coding_LoF" %in% unproductive_biotypes())
  iso <- seq_experiment()
  tx_file <- fs::path(iso@annot_path, "annot_data_transcripts.parquet")
  tx <- read_parquet(tx_file)
  tx$transcript_type[tx$transcript_id == "tA1"] <- "protein_coding_LoF"
  write_parquet(tx, tx_file)
  iso <- IsoformicExperiment(
    experiment_name = "lof",
    annot_path = iso@annot_path,
    assays = list(abundance = tx_assay(iso, "abundance")),
    col_data = as.data.frame(col_data(iso))
  )
  p <- transcript_productivity(iso)
  at <- p$transcript_id == "tA1"
  testthat::expect_identical(as.character(p$productivity[at]), "unproductive")
  testthat::expect_identical(p$reason[at], "unproductive biotype")
  rules <- transcript_productivity(iso, use_biotype = FALSE)
  testthat::expect_identical(
    as.character(rules$productivity[rules$transcript_id == "tA1"]),
    "productive"
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
  counts <- tx_assay(iso, "counts")
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
    m <- tx_assay(iso, "abundance")
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
  suppressMessages(add_annotation_columns(
    annot,
    data.frame(
      transcript_id = c("tA1", "tA3", "tA4"),
      appris_principal = c(TRUE, TRUE, FALSE)
    )
  ))
  clear_cache(iso)
  p <- transcript_productivity(iso)
  testthat::expect_identical(p$source[p$transcript_id == "tA1"], "evidence")
  testthat::expect_identical(
    as.character(p$productivity[p$transcript_id == "tA3"]),
    "productive"
  )
  suppressMessages(add_annotation_columns(
    annot,
    data.frame(transcript_id = "tA3", nmd_validated_exon = TRUE)
  ))
  clear_cache(iso)
  p <- transcript_productivity(iso)
  testthat::expect_identical(
    as.character(p$productivity[p$transcript_id == "tA3"]),
    "unproductive"
  )
  testthat::expect_match(p$reason[p$transcript_id == "tA3"], "validated")
  suppressMessages(add_annotation_columns(
    annot,
    data.frame(transcript_id = "tA3", productivity_evidence = "uncertain")
  ))
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

testthat::test_that("min_gene_count drops genes with few reads in any group", {
  iso <- seq_experiment()
  # Mean count of gA: 185 in group a, 180 in group b
  kept <- productive_fraction(iso, group = "condition", min_gene_count = 180)
  testthat::expect_identical(kept$gene_id, "gA")
  dropped <- productive_fraction(
    iso,
    group = "condition",
    min_gene_count = 181
  )
  testthat::expect_identical(nrow(dropped), 0L)
  # Without group, the mean over all the samples (182.5)
  testthat::expect_identical(
    nrow(productive_fraction(iso, min_gene_count = 182)),
    1L
  )
  testthat::expect_identical(
    nrow(productive_fraction(iso, as_matrix = TRUE, min_gene_count = 183)),
    0L
  )
  testthat::expect_identical(
    nrow(test_productive_fraction(
      iso,
      group = "condition",
      min_gene_count = 181
    )),
    0L
  )
  empty <- test_productive_fraction(
    iso,
    group = "condition",
    min_gene_count = 181
  )
  testthat::expect_type(empty$pvalue, "double")
  testthat::expect_type(empty$productive_fraction_a, "double")
  testthat::expect_error(
    productive_fraction(iso, min_gene_count = -1),
    class = "isoformic_invalid_argument"
  )
  testthat::expect_error(
    productive_fraction(iso, min_gene_count = NA),
    class = "isoformic_invalid_argument"
  )
  # Off by default; without a counts assay a filter is an error
  no_counts <- IsoformicExperiment(
    experiment_name = "no_counts",
    annot_path = iso@annot_path,
    assays = list(abundance = tx_assay(iso, "abundance")),
    col_data = col_data(iso)
  )
  testthat::expect_identical(
    productive_fraction(no_counts, group = "condition")$gene_id,
    "gA"
  )
  testthat::expect_identical(
    nrow(test_productive_fraction(iso, group = "condition")),
    1L
  )
  testthat::expect_error(
    productive_fraction(no_counts, min_gene_count = 10),
    class = "isoformic_assay_not_found"
  )
})

testthat::test_that("two_group_proportion_test does not test genes without usage information", {
  groups <- factor(c("a", "a", "b", "b"))
  den <- matrix(
    c(10, 10, 10, 10, 10, 10, 40, 40, 10, 10, 40, 40, 10, 10, 0, 0),
    nrow = 4,
    byrow = TRUE
  )
  # Row 1: never productive; row 2: always productive with a change in the
  # total; row 3: a real change; row 4: no reads in group b
  num <- rbind(c(0, 0, 0, 0), c(10, 10, 40, 40), c(9, 8, 4, 6), c(5, 5, 0, 0))
  res <- two_group_proportion_test(num, den, groups, min_class_count = 1)
  testthat::expect_equal(res$log_odds_ratio[1:2], c(0, 0))
  testthat::expect_true(all(is.na(res$pvalue[c(1, 2, 4)])))
  testthat::expect_true(all(is.na(res$statistic[c(1, 2, 4)])))
  testthat::expect_true(is.na(res$log_odds_ratio[4]))
  testthat::expect_false(is.na(res$pvalue[3]))
  # Only the testable row enters the BH adjustment
  testthat::expect_equal(res$padj[3], res$pvalue[3])
})

testthat::test_that("run_productivity_enrichment ranks by the test statistic", {
  testthat::skip_if_not_installed("fgsea")
  sets <- list(
    up = paste0("G", 1:10),
    down = paste0("G", 31:40)
  )
  tested <- tibble::tibble(
    gene_id = paste0("g", 1:40),
    gene_name = paste0("G", 1:40),
    delta_productive_fraction = seq(0.4, -0.38, by = -0.02),
    statistic = seq(8, -7.6, by = -0.4),
    pvalue = 2 * stats::pt(-abs(seq(8, -7.6, by = -0.4)), df = 4)
  )
  set.seed(1)
  # fgsea warns that some p-values are overestimated on so few genes
  res <- suppressWarnings(run_productivity_enrichment(
    tested,
    sets,
    rank_by = "statistic",
    min_size = 5
  ))
  testthat::expect_identical(unique(res$rank_by), "statistic")
  testthat::expect_true(res$NES[res$pathway == "up"] > 0)
  testthat::expect_true(res$NES[res$pathway == "down"] < 0)
  testthat::expect_error(
    run_productivity_enrichment(
      tested[, c("gene_id", "gene_name", "delta_productive_fraction")],
      sets,
      rank_by = "statistic"
    ),
    class = "isoformic_invalid_productive_fraction"
  )
})

testthat::test_that("a trace of reads in one class does not make a usage change", {
  groups <- factor(rep(c("a", "b"), each = 3))
  # 8-fold more reads in group b, all productive but for a trace
  den <- matrix(
    c(100, 110, 120, 800, 850, 900),
    nrow = 3,
    ncol = 6,
    byrow = TRUE
  )
  trace <- c(0.02, 2, 6)
  num <- den
  num[, 6] <- num[, 6] - trace
  res <- two_group_proportion_test(num, den, groups)
  # 0.02 and 2 reads (below 5) are not tested; 6 reads are
  testthat::expect_true(all(is.na(res$statistic[1:2])))
  testthat::expect_equal(res$log_odds_ratio[1:2], c(0, 0))
  testthat::expect_false(is.na(res$pvalue[3]))
  # Through test_productive_fraction(): the unproductive isoforms of gA hold
  # 0.02 reads in one sample only, while the gene grows 8-fold
  iso <- seq_experiment()
  counts <- tx_assay(iso, "counts")
  counts[] <- 0
  counts[c("tA1", "tA2"), ] <- cbind(
    c(60, 40),
    c(60, 50),
    c(480, 400),
    c(500, 400)
  )
  counts["tA3", "s4"] <- 0.02
  counts["tB1", ] <- 30
  traced <- IsoformicExperiment(
    experiment_name = "traced",
    annot_path = iso@annot_path,
    assays = list(abundance = tx_assay(iso, "abundance"), counts = counts),
    col_data = col_data(iso)
  )
  res_iso <- test_productive_fraction(traced, group = "condition")
  testthat::expect_true(is.na(res_iso$pvalue[res_iso$gene_id == "gA"]))
})

testthat::test_that("the statistic ranking compares genes with different degrees of freedom", {
  z <- signed_normal_quantile(c(-4.32, 3.47, 0, 2), c(0.012, 0.074, 1, 0))
  testthat::expect_equal(z[1], -stats::qnorm(0.006, lower.tail = FALSE))
  testthat::expect_equal(z[2], stats::qnorm(0.037, lower.tail = FALSE))
  testthat::expect_equal(z[3], 0)
  # A p-value of zero takes the largest finite score
  testthat::expect_equal(z[4], abs(z[1]))
  testthat::expect_true(abs(z[1]) > abs(z[2]))
})

testthat::test_that("poison exon candidates match a grouped summary", {
  withr::local_seed(2)
  n <- 400L
  exons <- tibble::tibble(
    transcript_id = sample(paste0("t", 1:60), n, replace = TRUE),
    seqid = "chr1",
    start_pos = sample(c(100L, 200L, 300L, 400L, 500L), n, replace = TRUE),
    end_pos = 0L,
    strand = sample(c("+", "-"), n, replace = TRUE)
  )
  exons$end_pos <- exons$start_pos + 50L
  exons$gene_id <- paste0(
    "g",
    as.integer(sub("t", "", exons$transcript_id, fixed = TRUE)) %% 4L
  )
  exons$productivity <- sample(
    c("unproductive", "productive", "noncoding", NA),
    n,
    replace = TRUE,
    prob = c(0.6, 0.2, 0.15, 0.05)
  )
  ref <- exons |>
    dplyr::summarise(
      n_unproductive = sum(.data$productivity == "unproductive"),
      n_other = sum(.data$productivity != "unproductive"),
      transcript_ids = paste(sort(unique(.data$transcript_id)), collapse = ";"),
      .by = c("gene_id", "seqid", "start_pos", "end_pos", "strand")
    ) |>
    dplyr::filter(.data$n_unproductive > 0L, .data$n_other == 0L)
  testthat::expect_gt(nrow(ref), 0L)
  testthat::expect_identical(poison_exon_candidates(exons), ref)
  testthat::expect_identical(
    nrow(poison_exon_candidates(exons[0L, ])),
    0L
  )
})

testthat::test_that("poison_exons reads only the exons it needs", {
  iso <- seq_experiment()
  invisible(transcript_productivity(iso))
  # Drop the exon table that the productivity calls loaded
  rm(
    list = paste0("exons::", attr(iso, "annot_path", exact = TRUE)),
    envir = annotation_cache
  )
  testthat::expect_false(annot_table_cached(iso, "exons"))
  cold <- poison_exons(iso)
  testthat::expect_false(annot_table_cached(iso, "exons"))
  invisible(iso@annot_data_exons)
  testthat::expect_true(annot_table_cached(iso, "exons"))
  testthat::expect_identical(poison_exons(iso), cold)
  testthat::expect_identical(
    poison_exons(iso, gene_name = "GENEA"),
    cold[cold$gene_name %in% "GENEA", ]
  )
})
