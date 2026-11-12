# `group` + `levels` choose the samples in the same way in every function:
# the result is the one of the object subset to the samples of those levels.

three_group_example <- function() {
  iso <- example_isoformic()
  sample_table <- col_data(iso)
  sample_table$g3 <- rep(c("a", "b", "c"), length.out = nrow(sample_table))
  col_data(iso) <- sample_table
  # Warm the cache: one example transcript warns about its sequence length
  suppressWarnings(transcript_productivity(iso))
  iso
}

example_draws <- function(iso, n = 8L) {
  set.seed(42)
  draws <- lapply(seq_len(n), function(k) {
    m <- assay(iso, "abundance")
    m * matrix(stats::runif(length(m), 0.7, 1.3), nrow = nrow(m))
  })
  attr(draws, "samp_type") <- "gibbs"
  draws
}

testthat::test_that("group_factor drops, orders and checks the levels", {
  values <- c("a", "b", "c", "a", NA)
  res <- group_factor(values, group = "g", levels = c("c", "a"))
  testthat::expect_identical(res$keep, c(TRUE, FALSE, TRUE, TRUE, FALSE))
  testthat::expect_identical(levels(res$groups), c("c", "a"))
  testthat::expect_identical(as.character(res$groups), c("a", "c", "a"))
  # Without levels: order of appearance, missing values dropped with a message
  testthat::expect_message(
    res <- group_factor(values, group = "g"),
    class = "isoformic_group_missing_values"
  )
  testthat::expect_identical(levels(res$groups), c("a", "b", "c"))
  testthat::expect_identical(sum(res$keep), 4L)
  # Levels of a factor without samples are dropped
  res <- group_factor(factor(c("b", "a"), levels = c("a", "zz", "b")), "g")
  testthat::expect_identical(levels(res$groups), c("a", "b"))
  testthat::expect_error(
    group_factor(values, group = "g", levels = c("a", "zzz")),
    class = "isoformic_group_levels"
  )
  testthat::expect_error(
    group_factor(values, group = "g", levels = c("a", "a")),
    class = "isoformic_group_levels"
  )
  testthat::expect_error(
    group_factor(values, group = "g", levels = "a", min_levels = 2L),
    class = "isoformic_group_levels"
  )
  testthat::expect_error(
    suppressMessages(
      group_factor(values, group = "g", min_levels = 2L, max_levels = 2L)
    ),
    class = "isoformic_group_levels"
  )
})

testthat::test_that("levels give the result of the object subset to those samples", {
  iso <- three_group_example()
  keep <- col_data(iso)$sample_id[col_data(iso)$g3 %in% c("a", "c")]
  sub <- iso[, keep]
  lv <- c("a", "c")

  testthat::expect_equal(
    productive_fraction(iso, group = "g3", levels = lv),
    productive_fraction(sub, group = "g3")
  )
  testthat::expect_identical(
    colnames(productive_fraction(
      iso,
      group = "g3",
      levels = lv,
      as_matrix = TRUE
    )),
    keep
  )
  testthat::expect_equal(
    productivity_fractions(iso, group = "g3", levels = lv),
    productivity_fractions(sub, group = "g3")
  )
  testthat::expect_equal(
    test_productive_fraction(iso, group = "g3", levels = lv),
    test_productive_fraction(sub, group = "g3")
  )
  testthat::expect_equal(
    apa_usage(iso, group = "g3", levels = lv),
    apa_usage(sub, group = "g3")
  )
  testthat::expect_identical(
    colnames(apa_usage(iso, group = "g3", levels = lv, as_matrix = TRUE)),
    keep
  )
  testthat::expect_equal(
    test_apa(iso, group = "g3", levels = lv),
    test_apa(sub, group = "g3")
  )
  testthat::expect_equal(
    dominant_isoform(iso, group = "g3", levels = lv),
    dominant_isoform(sub, group = "g3")
  )
  gene <- gene_names(iso)[[1L]]
  prof <- profile_data(iso, gene, group = "g3", levels = lv)
  testthat::expect_equal(prof, profile_data(sub, gene, group = "g3"))
  testthat::expect_identical(levels(prof$group), lv)
  testthat::expect_identical(unique(prof$n_samples), c(9L, 8L))

  # The order of `levels` is the direction of the change
  rev_lv <- productive_fraction(iso, group = "g3", levels = rev(lv))
  fwd_lv <- productive_fraction(iso, group = "g3", levels = lv)
  testthat::expect_equal(
    rev_lv$delta_productive_fraction[match(fwd_lv$gene_id, rev_lv$gene_id)],
    -fwd_lv$delta_productive_fraction
  )
})

testthat::test_that("every function rejects the same invalid levels", {
  iso <- three_group_example()
  gene <- gene_names(iso)[[1L]]
  calls <- list(
    function(...) productive_fraction(iso, ...),
    function(...) productivity_fractions(iso, ...),
    function(...) test_productive_fraction(iso, ...),
    function(...) apa_usage(iso, ...),
    function(...) test_apa(iso, ...),
    function(...) dominant_isoform(iso, ...),
    function(...) profile_data(iso, gene, ...)
  )
  # test_usage() checks its backend (Suggests) before the groups
  if (rlang::is_installed(c("limma", "edgeR"))) {
    calls <- c(calls, function(...) test_usage(iso, ...))
  }
  for (call in calls) {
    testthat::expect_error(
      call(group = "g3", levels = c("a", "zzz")),
      class = "isoformic_group_levels"
    )
    testthat::expect_error(
      call(group = "nope", levels = c("a", "b")),
      class = "isoformic_group_not_found"
    )
  }
  # `levels` without `group` is an error, not silently ignored
  for (call in calls[c(1L, 2L, 4L, 6L)]) {
    testthat::expect_error(
      call(levels = c("a", "b")),
      class = "isoformic_group_levels"
    )
  }
})

testthat::test_that("samples without a group value and unused factor levels are dropped", {
  iso <- three_group_example()
  sample_table <- col_data(iso)
  sample_table$g3[sample_table$g3 == "c"] <- NA
  col_data(iso) <- sample_table
  keep <- sample_table$sample_id[!is.na(sample_table$g3)]
  testthat::expect_message(
    res <- productive_fraction(iso, group = "g3"),
    class = "isoformic_group_missing_values"
  )
  testthat::expect_equal(res, productive_fraction(iso[, keep], group = "g3"))
  testthat::expect_no_message(
    productive_fraction(iso, group = "g3", levels = c("a", "b"))
  )

  # A subset keeps the factor levels of the sample table; they are not groups
  iso <- three_group_example()
  sample_table <- col_data(iso)
  sample_table$g3 <- factor(sample_table$g3, levels = c("a", "b", "c"))
  col_data(iso) <- sample_table
  sub <- iso[, sample_table$sample_id[sample_table$g3 != "b"]]
  res <- productive_fraction(sub, group = "g3")
  testthat::expect_false("productive_fraction_b" %in% colnames(res))
  testthat::expect_false(all(is.na(res$delta_productive_fraction)))
  testthat::expect_identical(
    unique(dominant_isoform(sub, group = "g3")$g3),
    c("a", "c")
  )
})

testthat::test_that("inferential summaries follow the levels", {
  iso <- three_group_example()
  iso@metadata$inf_reps <- example_draws(iso)
  keep <- col_data(iso)$sample_id[col_data(iso)$g3 %in% c("a", "b")]
  with_levels <- productive_fraction(iso, group = "g3", levels = c("a", "b"))
  on_subset <- productive_fraction(iso[, keep], group = "g3")
  testthat::expect_true(
    all(
      c("productive_fraction_sd", "delta_z") %in% colnames(with_levels)
    )
  )
  testthat::expect_equal(with_levels, on_subset)
  # A change of exactly zero in every draw is a z-score of zero
  flat <- with_levels$delta_productive_fraction %in%
    0 &
    with_levels$delta_productive_fraction_sd %in% 0
  testthat::expect_true(all(with_levels$delta_z[flat] == 0))
})

testthat::test_that("test_usage scales counts with the samples of the test", {
  testthat::skip_if_not_installed("limma")
  testthat::skip_if_not_installed("edgeR")
  x <- usage_experiment(n_groups = 3L)
  set.seed(7)
  draws <- lapply(1:10, function(k) {
    m <- assay(x$iso, "counts")
    m[] <- stats::rpois(length(m), as.vector(m) + 1)
    m
  })
  attr(draws, "samp_type") <- "gibbs"
  x$iso@metadata$inf_reps <- draws
  sample_table <- col_data(x$iso)
  lv <- unique(sample_table$condition)[1:2]
  keep <- sample_table$sample_id[sample_table$condition %in% lv]
  with_levels <- test_usage(x$iso, group = "condition", levels = lv)
  on_subset <- test_usage(x$iso[, keep], group = "condition")
  testthat::expect_true(unique(with_levels$scaled))
  testthat::expect_equal(with_levels, on_subset)
})

testthat::test_that("run_productivity_enrichment ranks by delta or delta_z", {
  testthat::skip_if_not_installed("fgsea")
  iso <- three_group_example()
  sets <- read_gmt(system.file(
    "extdata",
    "reactome_subset.gmt",
    package = "isoformic"
  ))
  lv <- c("a", "b")
  # Three groups in the object: the contrast is chosen with `levels`
  set.seed(1)
  res <- suppressWarnings(run_productivity_enrichment(
    iso,
    sets,
    group = "g3",
    levels = lv,
    min_size = 3
  ))
  testthat::expect_identical(unique(res$rank_by), "delta")
  testthat::expect_error(
    run_productivity_enrichment(
      iso,
      sets,
      group = "g3",
      levels = lv,
      rank_by = "delta_z"
    ),
    class = "isoformic_no_inferential_replicates"
  )
  testthat::expect_error(
    run_productivity_enrichment(iso, list(1, 2), group = "g3", levels = lv),
    class = "isoformic_invalid_genesets"
  )
  testthat::expect_error(
    run_productivity_enrichment(iso, sets, group = "g3", rank_by = "nope")
  )
  draws <- example_draws(iso)
  set.seed(1)
  res_z <- suppressWarnings(run_productivity_enrichment(
    iso,
    sets,
    group = "g3",
    levels = lv,
    min_size = 3,
    rank_by = "delta_z",
    inf_reps = draws
  ))
  testthat::expect_identical(unique(res_z$rank_by), "delta_z")
  testthat::expect_true(all(
    c("pathway", "NES", "padj", "experiment") %in% colnames(res_z)
  ))
  # The arguments of productive_fraction() reach it
  set.seed(1)
  strict <- suppressWarnings(run_productivity_enrichment(
    iso,
    sets,
    group = "g3",
    levels = lv,
    min_size = 3,
    include_altered = FALSE,
    min_expr = 1
  ))
  testthat::expect_s3_class(strict, "tbl_df")
  testthat::expect_false(isTRUE(all.equal(strict, res)))
  # As on the object subset to the samples of the levels
  keep <- col_data(iso)$sample_id[col_data(iso)$g3 %in% lv]
  set.seed(1)
  on_subset <- suppressWarnings(run_productivity_enrichment(
    iso[, keep],
    sets,
    group = "g3",
    min_size = 3
  ))
  testthat::expect_equal(res, on_subset)
  # One replicate gives no standard deviation, so no delta_z
  testthat::expect_error(
    run_productivity_enrichment(
      iso,
      sets,
      group = "g3",
      levels = lv,
      rank_by = "delta_z",
      inf_reps = draws[1L]
    ),
    class = "isoformic_no_inferential_replicates"
  )
})

testthat::test_that("delta_z ranks the genes with a certain change at the extreme", {
  testthat::skip_if_not_installed("fgsea")
  iso <- three_group_example()
  lv <- c("a", "b")
  # The same draw twice: every change has a standard deviation of zero
  draws <- rep(list(assay(iso, "abundance")), 2L)
  pf <- productive_fraction(iso, group = "g3", levels = lv, inf_reps = draws)
  changed <- !pf$delta_productive_fraction %in% c(0, NA)
  testthat::expect_true(any(changed))
  testthat::expect_true(all(is.na(pf$delta_z[changed])))
  ranked <- NULL
  testthat::local_mocked_bindings(
    fgsea_multilevel = function(pathways, stats, ...) {
      ranked <<- stats
      data.frame(pathway = character(), NES = numeric(), padj = numeric())
    }
  )
  run_productivity_enrichment(
    iso,
    list(set = "A"),
    group = "g3",
    levels = lv,
    rank_by = "delta_z",
    inf_reps = draws
  )
  genes <- pf$gene_name[changed & !is.na(pf$gene_name)]
  testthat::expect_true(all(genes %in% names(ranked)))
  testthat::expect_identical(
    unname(sign(ranked[genes[[1L]]])),
    sign(pf$delta_productive_fraction[match(genes[[1L]], pf$gene_name)])
  )
  testthat::expect_true(all(abs(ranked[genes]) == max(abs(ranked))))
})

testthat::test_that("a single level, numeric levels and profile tables are handled", {
  iso <- three_group_example()
  keep <- col_data(iso)$sample_id[col_data(iso)$g3 == "a"]
  res <- dominant_isoform(iso, group = "g3", levels = "a")
  testthat::expect_equal(res, dominant_isoform(iso[, keep], group = "g3"))
  testthat::expect_identical(unique(res$g3), "a")

  # A level given as a number is the value of the column, whatever its text
  res <- group_factor(c(100000L, 5L, 100000L), group = "dose", levels = 1e5)
  testthat::expect_identical(res$keep, c(TRUE, FALSE, TRUE))
  testthat::expect_identical(levels(res$groups), "100000")
  testthat::expect_error(
    group_factor(c(NA, NA), group = "g"),
    class = "isoformic_group_levels"
  )

  # `levels` chooses the groups of a profile table as well
  gene <- gene_names(iso)[[1L]]
  prof <- profile_data(iso, gene, group = "g3")
  plot <- plot_tx_expr(prof, levels = c("c", "a"))
  testthat::expect_identical(levels(plot$data$group), c("c", "a"))
  testthat::expect_error(
    plot_tx_expr(prof, levels = c("a", "zzz")),
    class = "isoformic_group_levels"
  )
  tables <- list(
    de_tx = de_tx(iso),
    tx_to_gene = tx_annot(iso)
  )
  testthat::expect_error(
    prepare_profile_data(
      txi_transcript = assay(iso, "abundance"),
      tx_to_gene = tables$tx_to_gene,
      sample_metadata = col_data(iso),
      de_result_gene = NULL,
      de_result_transcript = NULL,
      var = "g3",
      var_levels = c("a", "zzz")
    ),
    class = "isoformic_group_levels"
  )
})

testthat::test_that("the contrast of de_confidence follows the levels, its tables do not", {
  x <- uncertainty_experiment(inferential = TRUE)
  iso <- x$iso
  sample_table <- col_data(iso)
  sample_table$g3 <- c("a", "a", "b", "b", "c", "c")
  col_data(iso) <- sample_table
  testthat::expect_error(
    de_confidence(iso, group = "g3"),
    class = "isoformic_group_levels"
  )
  conf <- de_confidence(iso, group = "g3", levels = c("a", "b"))
  testthat::expect_true(all(c("feature_id", "confidence") %in% colnames(conf)))
  # `levels` without `group` is an error here as well
  testthat::expect_error(
    de_confidence(iso, levels = c("a", "b")),
    class = "isoformic_group_levels"
  )
  testthat::expect_error(
    add_uncertainty(iso, levels = c("a", "b")),
    class = "isoformic_group_levels"
  )
  with_conf <- add_uncertainty(iso, group = "g3", levels = c("a", "b"))
  testthat::expect_identical(
    with_conf@dea[["uncertainty"]]$mean_inf_rv,
    add_uncertainty(iso)@dea[["uncertainty"]]$mean_inf_rv
  )
})

testthat::test_that("collapsed ambiguity groups follow a transcript subset", {
  x <- uncertainty_experiment(inferential = TRUE)
  iso <- add_uncertainty(x$iso, group = "condition", collapse_groups = TRUE)
  conf <- de_confidence(iso)
  groups <- conf[conf$feature_type == "group", ]
  testthat::skip_if(nrow(groups) == 0L, "The fixture has no ambiguity group")
  members <- strsplit(groups$transcript_ids[[1L]], ";", fixed = TRUE)[[1L]]
  others <- setdiff(rownames(iso), members)
  group_id <- groups$feature_id[[1L]]
  kept <- function(ids) de_confidence(suppressMessages(iso[ids, ]))$feature_id
  testthat::expect_false(group_id %in% kept(others))
  testthat::expect_false(group_id %in% kept(c(members[[1L]], others)))
  testthat::expect_true(group_id %in% kept(c(members, others[[1L]])))
  # A member left alone is a transcript again, flagged as ambiguous
  alone <- suppressMessages(iso[c(members[[1L]], others), ])
  res <- de_confidence(alone, group = "condition", collapse_groups = TRUE)
  row <- res[res$feature_id == members[[1L]], ]
  testthat::expect_identical(row$feature_type, "transcript")
  testthat::expect_true(row$ambiguous)
  testthat::expect_false(group_id %in% res$feature_id)
})

testthat::test_that("subsetting carries the replicates and the uncertainty tables", {
  x <- uncertainty_experiment(inferential = TRUE)
  iso <- x$iso
  draws <- iso@metadata$inf_reps
  testthat::expect_message(
    sub <- iso[, c("s5", "s1", "s2")],
    class = "isoformic_stale_uncertainty"
  )
  sub_draws <- sub@metadata$inf_reps
  testthat::expect_length(sub_draws, length(draws))
  testthat::expect_identical(colnames(sub_draws[[1L]]), c("s5", "s1", "s2"))
  testthat::expect_identical(
    sub_draws[[3L]],
    draws[[3L]][, c("s5", "s1", "s2"), drop = FALSE]
  )
  testthat::expect_identical(
    attr(sub_draws, "samp_type"),
    attr(draws, "samp_type")
  )
  # A reorder keeps every sample: no message
  testthat::expect_no_message(iso[, rev(colnames(iso))])

  rows <- suppressMessages(iso[c("tA2", "tA1"), ])
  row_draws <- rows@metadata$inf_reps
  testthat::expect_setequal(rownames(row_draws[[1L]]), c("tA1", "tA2"))
  testthat::expect_identical(colnames(row_draws[[1L]]), colnames(iso))
  testthat::expect_setequal(
    tx_uncertainty(rows)$transcript_id,
    c("tA1", "tA2")
  )
  pairs <- rows@dea[["ambiguity"]]
  testthat::expect_true(all(
    c(pairs$transcript_id_1, pairs$transcript_id_2) %in% c("tA1", "tA2")
  ))
  testthat::expect_identical(dim(draws_array(rows))[1:2], c(2L, 6L))
  testthat::expect_setequal(
    inferential_overdispersion(rows)$transcript_id,
    c("tA1", "tA2")
  )
  # Draws hold the expressed transcripts only: unexpressed rows are no error
  testthat::expect_false("tA5" %in% rownames(draws[[1L]]))
  with_unexpressed <- suppressMessages(iso[c("tA5", "tA1"), ])
  testthat::expect_identical(
    rownames(with_unexpressed@metadata$inf_reps[[1L]]),
    "tA1"
  )
  # No kept transcript has draws: the object has no replicates left
  none <- suppressMessages(iso["tA5", ])
  testthat::expect_null(none@metadata$inf_reps)
  testthat::expect_error(
    inferential_overdispersion(none),
    class = "isoformic_no_inferential_replicates"
  )
  # The tables are recomputed from the subset draws on request
  redone <- add_uncertainty(sub)
  testthat::expect_false(identical(
    redone@dea[["uncertainty"]]$mean_inf_rv,
    iso@dea[["uncertainty"]]$mean_inf_rv
  ))
})

testthat::test_that("subscripts are read as IDs and checked", {
  iso <- mini_experiment(with_dea = FALSE)
  testthat::expect_identical(
    colnames(iso[, factor("s2", levels = c("s2", "s1"))]),
    "s2"
  )
  testthat::expect_identical(rownames(iso[-1, ]), c("tx2", "tx3"))
  testthat::expect_identical(
    colnames(iso[, c(FALSE, TRUE, TRUE, FALSE)]),
    c("s2", "s3")
  )
  testthat::expect_error(iso[, "nope"], class = "isoformic_invalid_subset")
  testthat::expect_error(iso[9, ], class = "isoformic_invalid_subset")
  testthat::expect_error(iso[, c(1, 1)], class = "isoformic_invalid_subset")
  testthat::expect_error(iso[, NA], class = "isoformic_invalid_subset")
  testthat::expect_error(iso[0, ], class = "isoformic_invalid_subset")
  testthat::expect_error(
    iso[rownames(iso) %in% "absent", ],
    class = "isoformic_invalid_subset"
  )
  testthat::expect_error(
    iso[, c(TRUE, TRUE, TRUE, TRUE, FALSE)],
    class = "isoformic_invalid_subset"
  )
})

testthat::test_that("subsetting keeps stored gene assays in step", {
  iso <- mini_experiment(with_dea = FALSE)
  stored <- gene_assay(iso, "abundance") * 3
  # Stored in another column order than the transcript assays
  gene_assays(iso) <- list(abundance = stored[, rev(colnames(stored))])
  sub <- iso[, c("s1", "s3")]
  testthat::expect_identical(
    gene_assay(sub, "abundance"),
    stored[, c("s1", "s3")]
  )
  one_gene <- iso[gene_ids_of(iso, rownames(iso)) == "gene1", ]
  testthat::expect_identical(
    rownames(gene_assay(one_gene, "abundance")),
    "gene1"
  )
  # Stored gene assays of some samples only, or of none of the kept genes
  part <- iso
  gene_assays(part) <- list(abundance = stored[, c("s1", "s2")])
  testthat::expect_identical(
    colnames(attr(part[1, ], "gene_assays")$abundance),
    c("s1", "s2")
  )
  testthat::expect_length(attr(part[, c("s3", "s4")], "gene_assays"), 0L)
  # The cached fractions follow a change of the stored gene assays
  derived <- mini_experiment(with_dea = FALSE)
  before <- isoform_fraction(derived)
  gene_assays(derived) <- list(abundance = stored)
  testthat::expect_equal(isoform_fraction(derived), before / 3)
  # Gene and transcript matrices are matched by sample, not by position
  frac <- isoform_fraction(iso)
  testthat::expect_equal(
    frac["tx1", "s1"],
    assay(iso, "abundance")["tx1", "s1"] / stored["gene1", "s1"]
  )
})

testthat::test_that("running moments match the standard deviation over the draws", {
  set.seed(3)
  x <- matrix(stats::rnorm(60, mean = 1e3), nrow = 6)
  x[1L, ] <- NA
  x[2L, -1L] <- NA
  x[3L, c(2L, 5L)] <- NA
  x[4L, ] <- 7
  stats_run <- running_moments(nrow(x))
  for (k in seq_len(ncol(x))) {
    stats_run <- running_moments_add(stats_run, x[, k])
  }
  testthat::expect_equal(
    running_moments_sd(stats_run),
    apply(x, 1L, stats::sd, na.rm = TRUE)
  )
  testthat::expect_identical(running_moments_sd(stats_run)[[4L]], 0)
})

testthat::test_that("productive_fraction reads the draws one at a time with the same result", {
  iso <- three_group_example()
  lv <- c("c", "a")
  draws <- example_draws(iso)
  # Draws with fewer transcripts, in another order, with one unknown row
  draws <- lapply(draws, function(m) {
    m <- m[rev(seq_len(nrow(m)))[-(1:5)], rev(seq_len(ncol(m)))]
    rbind(m, unknown_tx = 1)
  })
  pf <- productive_fraction(
    iso,
    group = "g3",
    levels = lv,
    inf_reps = draws,
    min_expr = 2
  )
  prod <- transcript_productivity(iso)
  samples <- col_data(iso)$sample_id[col_data(iso)$g3 %in% lv]
  groups <- col_data(iso)$g3[match(samples, col_data(iso)$sample_id)]
  tx_ids <- rownames(assay(iso, "abundance"))
  gene_ids <- prod$gene_id[match(tx_ids, prod$transcript_id)]
  is_prod <- prod$productivity[match(tx_ids, prod$transcript_id)] %in%
    c("productive", "altered_productive")
  fraction_of <- function(m) {
    full <- matrix(0, length(tx_ids), length(samples))
    dimnames(full) <- list(tx_ids, samples)
    shared <- intersect(tx_ids, rownames(m))
    full[shared, ] <- m[shared, samples]
    full[full < 2] <- 0
    frac <- rowsum(full * is_prod, gene_ids) / rowsum(full, gene_ids)
    frac[!is.finite(frac)] <- NA_real_
    frac[pf$gene_id, , drop = FALSE]
  }
  delta_of <- function(frac) {
    last <- rowMeans(frac[, groups == "a", drop = FALSE], na.rm = TRUE)
    first <- rowMeans(frac[, groups == "c", drop = FALSE], na.rm = TRUE)
    unname(last - first)
  }
  point <- fraction_of(assay(iso, "abundance"))
  testthat::expect_equal(pf$delta_productive_fraction, delta_of(point))
  draw_fracs <- lapply(draws, fraction_of)
  draw_delta <- vapply(draw_fracs, delta_of, numeric(nrow(pf)))
  expected_sd <- apply(draw_delta, 1L, stats::sd, na.rm = TRUE)
  testthat::expect_equal(pf$delta_productive_fraction_sd, expected_sd)
  sq_dev <- Reduce(`+`, lapply(draw_fracs, function(f) (f - point)^2))
  expected_sample_sd <- rowMeans(sqrt(sq_dev / length(draws)), na.rm = TRUE)
  expected_sample_sd[is.nan(expected_sample_sd)] <- NA_real_
  testthat::expect_equal(pf$productive_fraction_sd, unname(expected_sample_sd))
  testthat::expect_true(any(pf$productive_fraction_sd > 0, na.rm = TRUE))
})

testthat::test_that("inf_reps = FALSE skips the stored replicates", {
  iso <- three_group_example()
  lv <- c("a", "b")
  plain <- productive_fraction(iso, group = "g3", levels = lv)
  with_draws <- iso
  with_draws@metadata$inf_reps <- example_draws(iso)
  testthat::expect_true(
    "delta_z" %in%
      colnames(productive_fraction(with_draws, group = "g3", levels = lv))
  )
  testthat::expect_identical(
    productive_fraction(
      with_draws,
      group = "g3",
      levels = lv,
      inf_reps = FALSE
    ),
    plain
  )
  testthat::expect_error(
    productive_fraction(iso, inf_reps = TRUE),
    class = "isoformic_invalid_inf_reps"
  )
  # The default ranking does not read the replicates: unusable ones are not
  # even looked at
  broken <- iso
  broken@metadata$inf_reps <- list(matrix(1, 2, 2))
  testthat::expect_error(
    productive_fraction(broken, group = "g3", levels = lv),
    class = "isoformic_invalid_inf_reps"
  )
  ranked <- list()
  testthat::local_mocked_bindings(
    fgsea_multilevel = function(pathways, stats, ...) {
      ranked[[length(ranked) + 1L]] <<- stats
      data.frame(pathway = character(), NES = numeric(), padj = numeric())
    }
  )
  testthat::skip_if_not_installed("fgsea")
  run_productivity_enrichment(
    broken,
    list(set = "A"),
    group = "g3",
    levels = lv
  )
  run_productivity_enrichment(iso, list(set = "A"), group = "g3", levels = lv)
  testthat::expect_identical(ranked[[1L]], ranked[[2L]])
  testthat::expect_error(
    run_productivity_enrichment(
      with_draws,
      list(set = "A"),
      group = "g3",
      levels = lv,
      rank_by = "delta_z",
      inf_reps = FALSE
    ),
    class = "isoformic_no_inferential_replicates"
  )
  testthat::expect_error(
    run_productivity_enrichment(iso, list(set = "A")),
    class = "isoformic_group_levels"
  )
})

testthat::test_that("run_productivity_enrichment takes a productive_fraction table", {
  testthat::skip_if_not_installed("fgsea")
  iso <- three_group_example()
  lv <- c("a", "b")
  draws <- example_draws(iso)
  pf <- productive_fraction(iso, group = "g3", levels = lv, inf_reps = draws)
  ranked <- list()
  testthat::local_mocked_bindings(
    fgsea_multilevel = function(pathways, stats, ...) {
      ranked[[length(ranked) + 1L]] <<- stats
      data.frame(pathway = names(pathways), NES = 1, padj = 1)
    }
  )
  for (rank_by in c("delta", "delta_z")) {
    from_object <- run_productivity_enrichment(
      iso,
      list(set = "A"),
      group = "g3",
      levels = lv,
      rank_by = rank_by,
      inf_reps = draws
    )
    from_table <- run_productivity_enrichment(
      pf,
      list(set = "A"),
      rank_by = rank_by
    )
    testthat::expect_identical(from_table, from_object)
    testthat::expect_identical(from_table$rank_by, rank_by)
    testthat::expect_identical(from_table$experiment, "productive_fraction")
  }
  testthat::expect_identical(ranked[[1L]], ranked[[2L]])
  testthat::expect_identical(ranked[[3L]], ranked[[4L]])
  testthat::expect_false(identical(ranked[[1L]], ranked[[3L]]))
  err <- testthat::expect_error(
    run_productivity_enrichment(
      productive_fraction(iso, group = "g3", levels = lv),
      list(set = "A"),
      rank_by = "delta_z"
    ),
    class = "isoformic_no_inferential_replicates"
  )
  testthat::expect_no_match(conditionMessage(err), "inf_reps", fixed = TRUE)
  testthat::expect_error(
    run_productivity_enrichment(productive_fraction(iso), list(set = "A")),
    class = "isoformic_invalid_productive_fraction"
  )
})

testthat::test_that("draws without rows, integer draws and infinite values are handled", {
  sums <- gene_sums(
    matrix(1, 1, 2, dimnames = list("t", c("a", "b"))),
    2L,
    c("g1", "g2")
  )
  testthat::expect_identical(unname(sums[, "a"]), c(0, 1))
  testthat::expect_true(all(
    gene_sums(sums[0, , drop = FALSE], integer(), "g") == 0
  ))
  tx <- c("t1", "t2", "t3")
  fraction <- function(draw) {
    draw_productive_fraction(
      draw,
      tx_ids = tx,
      sample_ids = c("a", "b"),
      gene_idx = c(1L, 1L, 2L),
      gene_levels = c("g1", "g2"),
      is_prod = c(TRUE, FALSE, TRUE)
    )
  }
  big <- matrix(
    c(2000000000L, 2000000000L, 5L, 15L),
    nrow = 2,
    dimnames = list(c("t1", "t2"), c("a", "b"))
  )
  # No integer overflow; a gene without any row in the draw has no fraction
  testthat::expect_equal(fraction(big)["g1", ], c(a = 0.5, b = 0.25))
  testthat::expect_true(all(is.na(fraction(big)["g2", ])))
  # Only the unproductive transcript is in the draw: fraction 0, not missing
  testthat::expect_equal(
    unname(fraction(big["t2", , drop = FALSE])["g1", ]),
    c(0, 0)
  )
  # An infinite unproductive value gives no fraction
  inf_draw <- big * 1
  inf_draw["t2", "a"] <- Inf
  testthat::expect_identical(unname(fraction(inf_draw)["g1", "a"]), NA_real_)
  testthat::expect_true(all(is.na(fraction(big[0, , drop = FALSE]))))
  testthat::expect_error(
    fraction(unname(big)),
    class = "isoformic_invalid_inf_reps"
  )
  iso <- three_group_example()
  draws <- example_draws(iso, n = 3L)
  draws[[1L]] <- draws[[1L]][0, , drop = FALSE]
  pf <- productive_fraction(iso, group = "g3", inf_reps = draws)
  testthat::expect_true(any(!is.na(pf$delta_productive_fraction_sd)))
})
