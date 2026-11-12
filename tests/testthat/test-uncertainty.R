testthat::test_that("read_inferential_replicates reads Salmon draws exactly", {
  fx <- salmon_uncertainty_fixture()
  reps <- read_inferential_replicates(fx$dirs)
  testthat::expect_length(reps, 30L)
  testthat::expect_identical(attr(reps, "samp_type"), "gibbs")
  testthat::expect_identical(dimnames(reps[[1L]]), dimnames(fx$draws)[1:2])
  testthat::expect_equal(reps[[7L]], fx$draws[seq_len(8L), seq_len(6L), 7L])
  few <- read_inferential_replicates(
    fx$dirs,
    transcripts = c("tA3", "tB1"),
    n_draws = 5L
  )
  testthat::expect_length(few, 5L)
  testthat::expect_identical(rownames(few[[1L]]), c("tA3", "tB1"))
  testthat::expect_equal(few[[2L]]["tB1", "s4"], fx$draws["tB1", "s4", 2L])
  testthat::expect_message(
    read_inferential_replicates(fx$dirs, n_draws = 99L),
    "Only 30"
  )
  testthat::expect_error(
    read_inferential_replicates(fx$dirs, transcripts = "nope"),
    class = "isoformic_no_annotated_transcripts"
  )
  # Quantifications without replicates
  quants <- list.files(
    testthat::test_path("..", "..", "inst", "extdata", "mini_quants"),
    full.names = TRUE
  )
  if (length(quants)) {
    testthat::expect_error(
      read_inferential_replicates(quants[1L]),
      class = "isoformic_no_inferential_replicates"
    )
  }
  # A truncated file is detected
  bad <- fs::path(fx$dirs[["s1"]], "aux_info", "bootstrap", "bootstraps.gz")
  con <- gzfile(bad, "wb")
  writeBin(c(1, 2, 3), con)
  close(con)
  testthat::expect_error(
    read_inferential_replicates(fx$dirs["s1"]),
    class = "isoformic_invalid_inferential_replicates"
  )
})

testthat::test_that("equivalence_support reports the unique read support", {
  fx <- salmon_uncertainty_fixture()
  eq <- equivalence_support(fx$dirs)
  testthat::expect_equal(nrow(eq), 8L)
  one <- function(id, col) eq[[col]][eq$transcript_id == id]
  testthat::expect_equal(one("tA1", "unique_fraction"), 40 / 160)
  testthat::expect_equal(one("tA3", "unique_fraction"), 20 / 30)
  testthat::expect_equal(one("tA6", "unique_fraction"), 1)
  testthat::expect_equal(one("tA5", "n_classes"), 0)
  testthat::expect_true(is.na(one("tA5", "unique_fraction")))
  testthat::expect_equal(one("tA1", "n_classes"), 2)
  testthat::expect_equal(one("tA1", "n_fragments"), 6 * 160)
  testthat::expect_equal(one("tA1", "n_samples"), 6L)
  testthat::expect_identical(
    equivalence_support(fx$dirs, transcripts = "tA3")$transcript_id,
    "tA3"
  )
  fs::file_delete(fs::path(fx$dirs, "aux_info", "eq_classes.txt.gz"))
  testthat::expect_error(
    equivalence_support(fx$dirs),
    class = "isoformic_no_equivalence_classes"
  )
})

testthat::test_that("salmon_run_qc summarises meta_info.json", {
  fx <- salmon_uncertainty_fixture()
  qc <- salmon_run_qc(fx$dirs)
  testthat::expect_identical(qc$sample_id, fx$samples)
  testthat::expect_identical(unique(qc$salmon_version), "2.8.0")
  testthat::expect_identical(unique(qc$library_type), "IU")
  testthat::expect_equal(qc$percent_mapped[3L], 30)
  testthat::expect_identical(qc$qc_flag[3L], "low_mapping_rate")
  testthat::expect_true(all(is.na(qc$qc_flag[-3L])))
  testthat::expect_true(all(qc$em_converged))
  testthat::expect_identical(unique(qc$replicate_type), "gibbs")
  testthat::expect_equal(unique(qc$num_replicates), 30)
  # No replicates is flagged
  quants <- list.files(
    testthat::test_path("..", "..", "inst", "extdata", "mini_quants"),
    full.names = TRUE
  )
  if (length(quants)) {
    testthat::expect_match(salmon_run_qc(quants[1L])$qc_flag, "no_replicates")
  }
})

testthat::test_that("add_uncertainty measures quantification uncertainty and ambiguity", {
  x <- uncertainty_experiment()
  iso <- add_uncertainty(x$iso, paths = x$fixture$dirs)
  unc <- tx_uncertainty(iso)
  get <- function(id, col) unc[[col]][unc$transcript_id == id]
  # tA5 is not expressed: its draws are not read
  testthat::expect_false("tA5" %in% unc$transcript_id)
  testthat::expect_false("tA5" %in% rownames(iso@metadata$inf_reps[[1L]]))
  testthat::expect_true(all(unc$assessed))
  # A higher expression floor leaves low transcripts unassessed
  floor <- tx_uncertainty(add_uncertainty(iso, min_count = 75))
  testthat::expect_setequal(
    floor$transcript_id[floor$assessed],
    c("tA1", "tA2")
  )
  testthat::expect_true(is.na(floor$mean_inf_rv[floor$transcript_id == "tA3"]))
  # Inferential relative variance orders the isoforms by how well they are determined
  testthat::expect_lt(get("tA3", "mean_inf_rv"), 0.5)
  testthat::expect_lt(get("tA6", "mean_inf_rv"), 0.5)
  testthat::expect_gt(get("tA4", "mean_inf_rv"), 3)
  testthat::expect_gt(get("tA1", "mean_inf_rv"), 1)
  testthat::expect_true(get("tA4", "uncertain"))
  testthat::expect_false(get("tA3", "uncertain"))
  # tA1 and tA2 share their reads: draws are perfectly anti-correlated
  pairs <- isoform_ambiguity(iso)
  amb <- pairs[pairs$ambiguous, ]
  testthat::expect_equal(nrow(amb), 1L)
  testthat::expect_setequal(
    c(amb$transcript_id_1, amb$transcript_id_2),
    c("tA1", "tA2")
  )
  testthat::expect_lt(amb$mean_cor, -0.99)
  testthat::expect_equal(amb$n_samples, 6L)
  testthat::expect_identical(
    isoform_ambiguity(iso, ambiguous_only = TRUE)$gene_id,
    "gA"
  )
  testthat::expect_identical(get("tA1", "ambiguity_group"), "gA_amb1")
  testthat::expect_identical(get("tA2", "ambiguity_group"), "gA_amb1")
  testthat::expect_true(is.na(get("tA3", "ambiguity_group")))
  testthat::expect_equal(get("tA1", "n_ambiguous_partners"), 1L)
  testthat::expect_equal(get("tA3", "n_ambiguous_partners"), 0L)
  testthat::expect_lt(get("tA1", "min_sibling_cor"), -0.99)
  # Independent isoforms are not correlated
  tb <- pairs[pairs$gene_id == "gB", ]
  testthat::expect_lt(abs(tb$mean_cor), 0.3)
  # Equivalence classes and run metadata are stored
  testthat::expect_equal(get("tA1", "unique_fraction"), 0.25)
  testthat::expect_equal(get("tA3", "unique_fraction"), 20 / 30)
  testthat::expect_true("salmon_percent_mapped" %in% colnames(col_data(iso)))
  testthat::expect_identical(
    col_data(iso)$salmon_qc_flag[3L],
    "low_mapping_rate"
  )
  testthat::expect_equal(iso@metadata$uncertainty_n_draws, 30L)
  testthat::expect_identical(iso@metadata$uncertainty_replicate_type, "gibbs")
  # Sample IDs must match
  testthat::expect_error(
    add_uncertainty(
      x$iso,
      paths = stats::setNames(x$fixture$dirs, paste0("z", 1:6))
    ),
    class = "isoformic_sample_mismatch"
  )
  # Stored tables are used without the draws
  lean <- add_uncertainty(x$iso, paths = x$fixture$dirs, keep_draws = FALSE)
  testthat::expect_null(lean@metadata$inf_reps)
  testthat::expect_equal(tx_uncertainty(lean)$mean_inf_rv, unc$mean_inf_rv)
  testthat::expect_error(
    de_confidence(lean, group = "condition"),
    class = "isoformic_no_inferential_replicates"
  )
  # Recomputing from the stored draws with another cutoff
  again <- add_uncertainty(iso, inf_rv_cutoff = 100)
  testthat::expect_false(any(tx_uncertainty(again)$uncertain, na.rm = TRUE))
  testthat::expect_equal(again@metadata$uncertainty_inf_rv_cutoff, 100)
  testthat::expect_equal(
    tx_uncertainty(again)$unique_fraction,
    unc$unique_fraction
  )
  testthat::expect_identical(
    tx_uncertainty(iso, gene_name = "GENEB")$transcript_id,
    c("tB1", "tB2")
  )
})

testthat::test_that("without uncertainty tables the accessors explain what to do", {
  x <- uncertainty_experiment()
  testthat::expect_error(
    tx_uncertainty(x$iso),
    class = "isoformic_no_uncertainty"
  )
  testthat::expect_error(
    isoform_ambiguity(x$iso),
    class = "isoformic_no_uncertainty"
  )
  testthat::expect_error(
    de_confidence(x$iso),
    class = "isoformic_no_uncertainty"
  )
  testthat::expect_error(
    add_uncertainty(x$iso),
    class = "isoformic_no_inferential_replicates"
  )
  testthat::expect_error(
    isoform_switches(x$iso, confident_only = TRUE),
    class = "isoformic_no_uncertainty"
  )
})

testthat::test_that("read_salmon(inferential = TRUE) attaches the uncertainty tables", {
  x <- uncertainty_experiment(inferential = TRUE)
  testthat::expect_false(is.null(x$iso@dea$uncertainty))
  testthat::expect_false(is.null(x$iso@dea$ambiguity))
  testthat::expect_length(x$iso@metadata$inf_reps, 30L)
  testthat::expect_true("salmon_percent_mapped" %in% colnames(col_data(x$iso)))
})

testthat::test_that("de_confidence repeats the contrast on every draw", {
  x <- uncertainty_experiment()
  iso <- add_uncertainty(x$iso, paths = x$fixture$dirs, group = "condition")
  conf <- de_confidence(iso)
  get <- function(id, col) conf[[col]][conf$feature_id == id]
  testthat::expect_true(all(conf$feature_type == "transcript"))
  # Well determined, strongly changed isoform: high confidence
  testthat::expect_identical(as.character(get("tA3", "confidence")), "high")
  testthat::expect_gt(get("tA3", "lfc_lo"), 1)
  testthat::expect_equal(get("tA3", "sign_agreement"), 1)
  testthat::expect_gt(get("tA3", "frac_draws_sig"), 0.5)
  testthat::expect_true(is.na(get("tA3", "reason")))
  testthat::expect_true(get("tA3", "confident_de"))
  testthat::expect_equal(get("tA3", "det_log2FC"), 1.8)
  testthat::expect_equal(get("tA3", "n_det_sig"), 1L)
  # Ambiguous and uncertain isoform: low confidence even though its call is stable
  testthat::expect_identical(as.character(get("tA1", "confidence")), "low")
  testthat::expect_match(get("tA1", "reason"), "ambiguous")
  testthat::expect_match(get("tA1", "reason"), "uncertain quantification")
  testthat::expect_false(get("tA1", "confident_de"))
  testthat::expect_gt(get("tA1", "lfc_median"), 0)
  testthat::expect_lt(get("tA2", "lfc_median"), 0)
  # Overdispersed isoform without a real change: unstable fold-change
  testthat::expect_lt(get("tA4", "sign_agreement"), 0.95)
  testthat::expect_match(get("tA4", "reason"), "unstable fold-change")
  testthat::expect_identical(as.character(get("tA4", "confidence")), "low")
  testthat::expect_false(get("tA4", "confident_de"))
  # Not significant in the primary test
  testthat::expect_equal(get("tA6", "n_det_sig"), 0L)
  testthat::expect_false(get("tA6", "confident_de"))
  testthat::expect_s3_class(conf$confidence, "ordered")
  testthat::expect_identical(iso@metadata$uncertainty_group, "condition")
  testthat::expect_identical(iso@metadata$uncertainty_levels, "a|b")
  # Computed on demand, with the reference level swapped
  rev <- de_confidence(iso, group = "condition", levels = c("b", "a"))
  testthat::expect_equal(
    rev$lfc_median[rev$feature_id == "tA3"],
    -get("tA3", "lfc_median")
  )
  gene_b <- de_confidence(iso, gene_name = "GENEB")
  testthat::expect_identical(unique(gene_b$gene_id), "gB")
  testthat::expect_error(
    de_confidence(iso, group = "sample_id"),
    class = "isoformic_group_levels"
  )
  testthat::expect_error(
    de_confidence(iso, group = "condition", levels = c("a", "zzz")),
    class = "isoformic_group_levels"
  )
})

testthat::test_that("collapse_groups tests ambiguous isoforms as one feature", {
  x <- uncertainty_experiment()
  iso <- add_uncertainty(x$iso, paths = x$fixture$dirs)
  grp <- de_confidence(iso, group = "condition", collapse_groups = TRUE)
  testthat::expect_true("gA_amb1" %in% grp$feature_id)
  testthat::expect_false(any(c("tA1", "tA2") %in% grp$feature_id))
  g <- grp[grp$feature_id == "gA_amb1", ]
  testthat::expect_identical(g$feature_type, "group")
  testthat::expect_setequal(
    strsplit(g$transcript_ids, ";", fixed = TRUE)[[1L]],
    c("tA1", "tA2")
  )
  testthat::expect_false(g$ambiguous)
  # Their sum does not change between the groups
  testthat::expect_lt(abs(g$lfc_median), 0.5)
  testthat::expect_equal(g$n_det_sig, 2L)
  testthat::expect_lt(g$mean_inf_rv, 0.5)
  testthat::expect_true(is.na(g$det_log2FC))
})

testthat::test_that("library_size in col_data sets the scaling of the fold-changes", {
  x <- uncertainty_experiment()
  iso <- add_uncertainty(x$iso, paths = x$fixture$dirs)
  base <- de_confidence(iso, group = "condition")
  sample_table <- col_data(iso)
  sample_table$library_size <- ifelse(
    sample_table$condition == "b",
    2e6,
    1e6
  )
  col_data(iso) <- sample_table
  scaled <- de_confidence(iso, group = "condition")
  # A larger library size for group b scales its counts down, so the
  # fold-change of the well-determined tA3 falls by about one log2 unit
  shift <- scaled$lfc_median[scaled$feature_id == "tA3"] -
    base$lfc_median[base$feature_id == "tA3"]
  testthat::expect_lt(shift, -0.5)
  testthat::expect_gt(shift, -1.5)
})

testthat::test_that("uncertainty tables survive write_isoformic without the draws", {
  x <- uncertainty_experiment()
  iso <- add_uncertainty(x$iso, paths = x$fixture$dirs, group = "condition")
  out <- withr::local_tempdir()
  testthat::expect_no_warning(
    written <- write_isoformic(iso, fs::path(out, "exp"), overwrite = TRUE)
  )
  testthat::expect_true(fs::file_exists(fs::path(
    out,
    "exp",
    "dea",
    "uncertainty.parquet"
  )))
  back <- read_isoformic(fs::path(out, "exp"))
  testthat::expect_null(back@metadata$inf_reps)
  testthat::expect_equal(
    tx_uncertainty(back)$mean_inf_rv,
    tx_uncertainty(iso)$mean_inf_rv
  )
  testthat::expect_equal(
    as.data.frame(isoform_ambiguity(back))$mean_cor,
    as.data.frame(isoform_ambiguity(iso))$mean_cor
  )
  testthat::expect_equal(
    de_confidence(back)$sign_agreement,
    de_confidence(iso)$sign_agreement
  )
  testthat::expect_equal(back@metadata$uncertainty_n_draws, "30")
  testthat::expect_equal(uncertainty_settings(back)$inf_rv_cutoff, 0.5)
})

testthat::test_that("isoform_switches and dominant_isoform use the confidence", {
  x <- uncertainty_experiment()
  iso <- add_uncertainty(x$iso, paths = x$fixture$dirs, group = "condition")
  sw <- isoform_switches(iso)
  ga <- sw[sw$gene_id == "gA", ]
  testthat::expect_identical(ga$switch_type, "opposite_direction")
  testthat::expect_equal(ga$n_det_sig, 3L)
  testthat::expect_equal(ga$n_det_confident, 1L)
  conf_sw <- isoform_switches(iso, confident_only = TRUE)
  ga <- conf_sw[conf_sw$gene_id == "gA", ]
  testthat::expect_equal(ga$n_det_sig, 1L)
  testthat::expect_identical(ga$switch_type, "det_not_deg")
  dom <- dominant_isoform(iso, group = "condition")
  testthat::expect_true(all(
    c("dominant_inf_rv", "dominant_ambiguous") %in% colnames(dom)
  ))
  # Without tables the columns are absent
  testthat::expect_false(
    "dominant_inf_rv" %in% colnames(dominant_isoform(x$iso))
  )
  testthat::expect_false(
    "n_det_confident" %in% colnames(isoform_switches(x$iso))
  )
})

testthat::test_that("uncertainty_report and uncertainty_methods summarise the settings", {
  x <- uncertainty_experiment()
  iso <- add_uncertainty(x$iso, paths = x$fixture$dirs, group = "condition")
  rep <- uncertainty_report(iso)
  testthat::expect_named(
    rep,
    c("settings", "replicates", "run", "transcripts", "genes", "de", "switches")
  )
  testthat::expect_equal(rep$replicates$n_draws, 30L)
  testthat::expect_equal(rep$replicates$resampling_df, 6L * 29L)
  testthat::expect_identical(rep$replicates$replicate_type, "gibbs")
  testthat::expect_gte(rep$replicates$median_overdispersion, 1)
  testthat::expect_equal(sum(rep$transcripts$n_assessed), 7L)
  testthat::expect_equal(
    rep$genes$n_ambiguous_groups[rep$genes$gene_id == "gA"],
    1L
  )
  testthat::expect_equal(nrow(rep$run), 6L)
  de <- rep$de
  testthat::expect_identical(
    as.character(de$confidence),
    c("high", "moderate", "low")
  )
  testthat::expect_equal(sum(de$n_significant), 3L)
  testthat::expect_equal(de$n_significant[de$confidence == "high"], 1L)
  testthat::expect_equal(de$n_ambiguous[de$confidence == "low"], 2L)
  testthat::expect_equal(rep$switches$n_genes, c(1L, 1L))
  txt <- uncertainty_methods(iso)
  testthat::expect_match(txt, "Salmon 2.8.0", fixed = TRUE)
  testthat::expect_match(txt, "30 Gibbs samples", fixed = TRUE)
  testthat::expect_match(txt, "174 degrees of freedom", fixed = TRUE)
  testthat::expect_match(txt, "Turro", fixed = TRUE)
  testthat::expect_match(txt, "b versus a", fixed = TRUE)
  testthat::expect_match(txt, "InfRV of at least 0.5", fixed = TRUE)
  testthat::expect_match(txt, "95%", fixed = TRUE)
})

testthat::test_that("the uncertainty plots build", {
  x <- uncertainty_experiment()
  iso <- add_uncertainty(x$iso, paths = x$fixture$dirs, group = "condition")
  p1 <- plot_de_confidence(iso, gene_name = "GENEA")
  testthat::expect_s3_class(p1, "ggplot")
  testthat::expect_no_error(ggplot2::ggplot_build(p1))
  testthat::expect_setequal(as.character(p1$data$confidence), c("high", "low"))
  p2 <- plot_isoform_correlation(iso, gene_name = c("GENEA", "GENEB"))
  testthat::expect_s3_class(p2, "ggplot")
  testthat::expect_no_error(ggplot2::ggplot_build(p2))
  testthat::expect_equal(nrow(p2$data[p2$data$mean_cor == 1, ]), 7L)
  p3 <- plot_inf_rv(iso)
  testthat::expect_s3_class(p3, "ggplot")
  testthat::expect_no_error(ggplot2::ggplot_build(p3))
  testthat::expect_error(
    plot_de_confidence(iso, gene_name = "nope"),
    class = "isoformic_no_uncertainty"
  )
  testthat::expect_named(confidence_palette(), c("high", "moderate", "low"))
})

testthat::test_that("productive_fraction accepts replicates of expressed transcripts only", {
  x <- uncertainty_experiment()
  iso <- add_uncertainty(x$iso, paths = x$fixture$dirs)
  pf <- productive_fraction(iso, group = "condition")
  testthat::expect_true("productive_fraction_sd" %in% colnames(pf))
})

testthat::test_that("salmon_quant draws 30 Gibbs samples by default and can skip them", {
  testthat::skip_if_not_installed("condathis")
  captured <- new.env()
  mock_run <- function(...) {
    assign("args", unlist(list(...)), envir = captured)
    invisible(NULL)
  }
  index_dir <- withr::local_tempdir()
  testthat::local_mocked_bindings(
    env_exists = function(...) TRUE,
    run = mock_run,
    .package = "condathis"
  )
  salmon_quant(input_r1 = "reads.fq", index_path = index_dir)
  args <- captured$args
  testthat::expect_identical(
    unname(args[which(args == "--numGibbsSamples") + 1L]),
    "30"
  )
  salmon_quant(input_r1 = "reads.fq", index_path = index_dir, num_gibbs = 0)
  testthat::expect_false("--numGibbsSamples" %in% captured$args)
})

testthat::test_that("replicates convert to a long table and back without loss", {
  fx <- salmon_uncertainty_fixture()
  reps <- read_inferential_replicates(fx$dirs, n_draws = 7L)
  tab <- inferential_replicates_table(reps)
  testthat::expect_identical(
    colnames(tab),
    c("transcript_id", "sample_id", "draw", "count")
  )
  testthat::expect_equal(
    nrow(tab),
    nrow(reps[[1L]]) * ncol(reps[[1L]]) * 7L
  )
  back <- as_inferential_replicates(tab, samp_type = "gibbs")
  testthat::expect_identical(attr(back, "samp_type"), "gibbs")
  attr(reps, "samp_type") <- NULL
  attr(back, "samp_type") <- NULL
  testthat::expect_equal(back, reps)
  # Row order of the table only changes the order of the dimnames
  shuffled <- as_inferential_replicates(tab[rev(seq_len(nrow(tab))), ])
  aligned <- lapply(shuffled, function(m) {
    m[rownames(reps[[1L]]), colnames(reps[[1L]])]
  })
  testthat::expect_equal(aligned, unclass(reps))
  testthat::expect_error(
    as_inferential_replicates(tab[-1L, ]),
    class = "isoformic_invalid_inferential_replicates"
  )
  testthat::expect_error(
    as_inferential_replicates(tab[, c("transcript_id", "count")]),
    class = "isoformic_invalid_inferential_replicates"
  )
  testthat::expect_error(
    inferential_replicates_table(list()),
    class = "isoformic_invalid_inferential_replicates"
  )
})

testthat::test_that("inferential_overdispersion follows the quasi-Poisson estimator of catchSalmon", {
  x <- uncertainty_experiment()
  reps <- read_inferential_replicates(x$fixture$dirs)
  iso <- x$iso
  iso@metadata$inf_reps <- reps
  od <- inferential_overdispersion(iso)
  testthat::expect_identical(
    colnames(od),
    c("transcript_id", "overdispersion", "df")
  )
  testthat::expect_true(all(od$overdispersion >= 1))

  # The estimator by hand: pooled Pearson statistic, then moderated
  arr <- simplify2array(unname(lapply(reps, unclass)))
  m <- apply(arr, c(1L, 2L), mean)
  ss <- apply(arr, c(1L, 2L), function(v) sum((v - mean(v))^2))
  ok <- m > 0
  num <- rowSums(ifelse(ok, ss / m, 0))
  df <- rowSums(ok) * (dim(arr)[3L] - 1L)
  phi <- num / df
  has <- df > 0
  prior <- max(
    1,
    stats::median(phi[has]) / stats::qf(0.5, stats::median(df[has]), 3)
  )
  expected <- pmax((3 * prior + df * phi) / (3 + df), 1)
  expected[!has] <- prior
  testthat::expect_equal(od$overdispersion, unname(expected))
  testthat::expect_equal(attr(od, "prior"), prior)

  # A noisy, ambiguous isoform has more overdispersion than a clean one
  by_id <- stats::setNames(od$overdispersion, od$transcript_id)
  testthat::expect_gt(by_id[["tA1"]], by_id[["tA3"]])

  testthat::skip_if_not_installed("edgeR")
  testthat::skip_if_not_installed("jsonlite")
  testthat::skip_if_not_installed("readr")
  cs <- suppressMessages(
    edgeR::catchSalmon(unname(x$fixture$dirs), verbose = FALSE)
  )
  testthat::expect_equal(
    unname(by_id[rownames(cs$counts)]),
    unname(cs$annotation$Overdispersion)
  )
})

testthat::test_that("scale_counts divides the counts by the overdispersion", {
  x <- uncertainty_experiment()
  iso <- add_uncertainty(x$iso, paths = x$fixture$dirs)
  scaled <- scale_counts(iso)
  testthat::expect_true("counts_scaled" %in% assay_names(scaled))
  od <- inferential_overdispersion(iso)
  counts <- assay(iso, "counts")
  ratio <- counts / assay(scaled, "counts_scaled")
  testthat::expect_equal(
    unname(ratio[od$transcript_id, 1L]),
    od$overdispersion
  )
  testthat::expect_error(
    scale_counts(iso, assay = "nope"),
    class = "isoformic_assay_not_found"
  )
})

testthat::test_that("the inferential relative variance is the one of fishpond", {
  moments <- list(
    mean = matrix(c(100, 100, 0, 2), ncol = 1L),
    var = matrix(c(100, 310, 0, 1), ncol = 1L)
  )
  # no excess over counting noise gives the shift, never less
  testthat::expect_equal(
    inferential_rv(moments)[, 1L],
    c(0.01, 210 / 105 + 0.01, 0.01, 0.01)
  )
  testthat::skip_if_not_installed("fishpond")
  testthat::skip_if_not_installed("SummarizedExperiment")
  withr::local_seed(1L)
  arr <- array(
    stats::rpois(4L * 3L * 30L, 50) * rep(c(1, 1, 3, 0.2), 90L),
    dim = c(4L, 3L, 30L),
    dimnames = list(paste0("t", 1:4), paste0("s", 1:3), NULL)
  )
  mom <- draw_moments(arr)
  reps <- stats::setNames(
    lapply(seq_len(30L), function(k) draw_slice(arr, k)),
    paste0("infRep", seq_len(30L))
  )
  se <- SummarizedExperiment::SummarizedExperiment(
    c(list(counts = mom$mean), reps)
  )
  expected <- SummarizedExperiment::mcols(fishpond::computeInfRV(se))$meanInfRV
  testthat::expect_equal(
    unname(rowMeans(inferential_rv(mom))),
    unname(expected)
  )
})

testthat::test_that("too few replicates over the study are reported", {
  testthat::expect_equal(resampling_df(6L, 30L), 174L)
  testthat::expect_message(
    inform_resampling_df(resampling_df(2L, 10L), "bootstrap"),
    class = "isoformic_few_inferential_replicates"
  )
  testthat::expect_message(
    inform_resampling_df(resampling_df(2L, 10L), "bootstrap"),
    "Gibbs samples"
  )
  testthat::expect_no_message(inform_resampling_df(
    resampling_df(3L, 30L),
    "gibbs"
  ))
  testthat::expect_no_message(inform_resampling_df(NULL))
})

testthat::test_that("collapse_gain is the Terminus score of a pair", {
  x <- uncertainty_experiment()
  iso <- add_uncertainty(x$iso, paths = x$fixture$dirs)
  pairs <- isoform_ambiguity(iso)
  unc <- tx_uncertainty(iso)
  rv <- function(id) unc$mean_inf_rv[unc$transcript_id == id]
  # tA1 + tA2 is constant over the draws, so the sum has the floor value
  shared <- pairs[pairs$ambiguous, ]
  testthat::expect_equal(
    shared$collapse_gain,
    0.01 - (rv("tA1") + rv("tA2")) / 2
  )
  testthat::expect_lt(shared$collapse_gain, -1)
  # Independent isoforms gain nothing from being summed
  independent <- pairs[pairs$gene_id == "gB", ]
  testthat::expect_lt(abs(independent$collapse_gain), 0.2)
  # Same value as computing the three InfRV from the draws
  draws <- x$fixture$draws
  one <- function(v) max(stats::var(v) - mean(v), 0) / (mean(v) + 5) + 0.01
  by_hand <- mean(vapply(
    seq_len(6L),
    function(s) {
      one(draws["tA3", s, ] + draws["tA4", s, ]) -
        (one(draws["tA3", s, ]) + one(draws["tA4", s, ])) / 2
    },
    numeric(1L)
  ))
  pair <- pairs[
    pairs$transcript_id_1 == "tA3" & pairs$transcript_id_2 == "tA4",
  ]
  testthat::expect_equal(pair$collapse_gain, by_hand)
})

testthat::test_that("low_information marks uncertain isoforms without unique fragments", {
  x <- uncertainty_experiment()
  iso <- add_uncertainty(x$iso, paths = x$fixture$dirs)
  unc <- tx_uncertainty(iso)
  testthat::expect_identical(
    unc$low_information,
    unc$unique_fraction == 0 & unc$uncertain
  )
  no_classes <- add_uncertainty(
    x$iso,
    paths = x$fixture$dirs,
    equivalence = FALSE
  )
  testthat::expect_false(
    "low_information" %in% colnames(tx_uncertainty(no_classes))
  )
})

testthat::test_that("cor_cutoff = 'adaptive' derives the cutoff from the pairs", {
  pairs <- tibble::tibble(
    transcript_id_1 = c("a", "a", "b", "d"),
    transcript_id_2 = c("b", "c", "c", "e"),
    mean_cor = c(-0.9, 0.2, 0.1, 0.4)
  )
  # Largest correlation per isoform: a 0.2, b 0.1, c 0.2, d 0.4, e 0.4
  testthat::expect_equal(
    resolve_cor_cutoff("adaptive", pairs),
    -stats::quantile(c(0.2, 0.1, 0.2, 0.4, 0.4), 0.975, names = FALSE)
  )
  testthat::expect_identical(resolve_cor_cutoff(-0.3, pairs), -0.3)
  testthat::expect_error(
    resolve_cor_cutoff("auto", pairs),
    class = "isoformic_invalid_cor_cutoff"
  )
  # Nothing to derive it from: the default, with a message
  testthat::expect_message(
    fallback <- resolve_cor_cutoff("adaptive", pairs[0L, ]),
    class = "isoformic_no_adaptive_cutoff"
  )
  testthat::expect_identical(fallback, -0.5)
  x <- uncertainty_experiment()
  iso <- suppressMessages(add_uncertainty(
    x$iso,
    paths = x$fixture$dirs,
    cor_cutoff = "adaptive"
  ))
  testthat::expect_identical(
    iso@metadata$uncertainty_cor_cutoff_rule,
    "adaptive"
  )
  testthat::expect_lt(iso@metadata$uncertainty_cor_cutoff, 0)
  pairs <- isoform_ambiguity(iso)
  testthat::expect_identical(
    pairs$ambiguous,
    pairs$mean_cor <= iso@metadata$uncertainty_cor_cutoff
  )
  fixed <- add_uncertainty(x$iso, paths = x$fixture$dirs)
  testthat::expect_identical(
    fixed@metadata$uncertainty_cor_cutoff_rule,
    "fixed"
  )
  testthat::expect_match(uncertainty_methods(fixed), "with a fixed cutoff")
})
