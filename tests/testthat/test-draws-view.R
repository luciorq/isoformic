# x[i, j] keeps the inferential replicates as a view of the draws of x. The
# reference below is the copy that x[i, j] stored before: every draw
# restricted to the kept transcripts (in their order) and samples, and the
# moments restricted to the kept samples.
copied_draws <- function(reps, tx_ids, sample_ids) {
  out <- lapply(reps, function(m) {
    rows <- match(tx_ids, rownames(m))
    m[rows[!is.na(rows)], sample_ids, drop = FALSE]
  })
  attr(out, "samp_type") <- attr(reps, "samp_type")
  mom <- attr(reps, "moments")
  if (!is.null(mom)) {
    mom$mean <- mom$mean[, sample_ids, drop = FALSE]
    mom$var <- mom$var[, sample_ids, drop = FALSE]
    attr(out, "moments") <- mom
  }
  out
}

# The subset with the copied draws in place of the view
with_copied_draws <- function(sub, reps) {
  sub@metadata$inf_reps <- copied_draws(
    reps,
    rownames(tx_assay(sub)),
    colnames(tx_assay(sub))
  )
  sub
}

testthat::test_that("a subset shares the draws and reads them through its selection", {
  x <- uncertainty_experiment(inferential = TRUE)
  iso <- x$iso
  draws <- iso@metadata$inf_reps
  sub <- suppressMessages(iso[c("tA2", "tB1", "tA1"), c("s5", "s1", "s2")])
  reps <- sub@metadata$inf_reps
  testthat::expect_identical(
    attr(reps, "selection"),
    list(
      transcripts = c("tA2", "tB1", "tA1"),
      samples = c("s5", "s1", "s2")
    )
  )
  # The matrices are those of the parent, not copies of the kept part
  for (k in seq_along(draws)) {
    testthat::expect_identical(reps[[k]], draws[[k]])
  }
  if (isTRUE(capabilities("profmem"))) {
    address <- function(m) {
      out <- tracemem(m)
      untracemem(m)
      out
    }
    testthat::expect_identical(address(reps[[1L]]), address(draws[[1L]]))
  }
  # The memory of the subset does not grow with the draws: it holds the
  # parent matrices (all their rows) and the IDs of its selection
  testthat::expect_identical(
    vapply(reps, nrow, integer(1L)),
    vapply(draws, nrow, integer(1L))
  )
  # The accessor gives the copy that x[i, j] stored before
  testthat::expect_identical(
    inf_reps_materialize(reps),
    copied_draws(draws, c("tA2", "tB1", "tA1"), c("s5", "s1", "s2"))
  )
  read_draw <- inf_reps_reader(reps)
  testthat::expect_identical(
    read_draw(2L),
    draws[[2L]][c("tA2", "tB1", "tA1"), c("s5", "s1", "s2")]
  )
  # Without a selection the draws are read as they are
  testthat::expect_identical(inf_reps_materialize(draws), draws)
  testthat::expect_identical(inf_reps_reader(draws)(1L), draws[[1L]])
  testthat::expect_identical(inf_reps_moments(draws), attr(draws, "moments"))
  # The long tables follow the selection
  table <- inferential_replicates_table(reps)
  testthat::expect_identical(
    unique(table$transcript_id),
    c("tA2", "tB1", "tA1")
  )
  testthat::expect_identical(unique(table$sample_id), c("s5", "s1", "s2"))
  testthat::expect_identical(
    table,
    inferential_replicates_table(inf_reps_materialize(reps))
  )
  moments <- inferential_moments_table(reps)
  testthat::expect_identical(unique(moments$sample_id), c("s5", "s1", "s2"))
})

testthat::test_that("nested subsets compose their selections", {
  x <- uncertainty_experiment(inferential = TRUE)
  iso <- x$iso
  draws <- iso@metadata$inf_reps
  outer <- suppressMessages(iso[c("tA3", "tA2", "tB1", "tA1"), -3L])
  inner <- suppressMessages(outer[c("tA1", "tA3"), c("s6", "s1")])
  reps <- inner@metadata$inf_reps
  # Still the matrices of the first object
  testthat::expect_identical(reps[[1L]], draws[[1L]])
  testthat::expect_identical(
    attr(reps, "selection"),
    list(transcripts = c("tA1", "tA3"), samples = c("s6", "s1"))
  )
  direct <- suppressMessages(iso[c("tA1", "tA3"), c("s6", "s1")])
  testthat::expect_identical(
    inf_reps_materialize(reps),
    inf_reps_materialize(direct@metadata$inf_reps)
  )
  # The same as subsetting the copy of the first subset
  copy_outer <- with_copied_draws(outer, draws)
  copy_inner <- suppressMessages(copy_outer[c("tA1", "tA3"), c("s6", "s1")])
  testthat::expect_identical(
    inf_reps_materialize(reps),
    inf_reps_materialize(copy_inner@metadata$inf_reps)
  )
  testthat::expect_identical(
    inferential_overdispersion(inner),
    inferential_overdispersion(copy_inner)
  )
})

testthat::test_that("a subset gives the results of a copy of its draws", {
  x <- uncertainty_experiment(inferential = TRUE)
  iso <- x$iso
  draws <- iso@metadata$inf_reps
  subsets <- list(
    samples = suppressMessages(iso[, c("s6", "s1", "s2", "s4")]),
    transcripts = suppressMessages(iso[c("tA3", "tA1", "tA2", "tB1"), ]),
    both = suppressMessages(iso[c("tA2", "tA1", "tA3"), c("s4", "s1", "s5")])
  )
  for (nm in names(subsets)) {
    sub <- subsets[[nm]]
    copy <- with_copied_draws(sub, draws)
    testthat::expect_identical(
      inf_reps_materialize(sub@metadata$inf_reps),
      copy@metadata$inf_reps,
      label = nm
    )
    testthat::expect_identical(
      tx_assay(scale_counts(sub), "counts_scaled"),
      tx_assay(scale_counts(copy), "counts_scaled"),
      label = nm
    )
    testthat::expect_identical(
      inferential_overdispersion(sub),
      inferential_overdispersion(copy),
      label = nm
    )
    testthat::expect_identical(
      draws_array(sub),
      draws_array(copy),
      label = nm
    )
    testthat::expect_identical(
      de_confidence(sub, group = "condition"),
      de_confidence(copy, group = "condition"),
      label = nm
    )
    testthat::expect_identical(
      productive_fraction(sub, group = "condition"),
      productive_fraction(copy, group = "condition"),
      label = nm
    )
    # Without the moments the overdispersion comes from the draws themselves
    no_mom <- sub
    attr(no_mom@metadata$inf_reps, "moments") <- NULL
    no_mom_copy <- copy
    attr(no_mom_copy@metadata$inf_reps, "moments") <- NULL
    testthat::expect_identical(
      inferential_overdispersion(no_mom),
      inferential_overdispersion(no_mom_copy),
      label = nm
    )
    unc <- suppressMessages(add_uncertainty(sub))
    unc_copy <- suppressMessages(add_uncertainty(copy))
    testthat::expect_identical(unc@dea, unc_copy@dea, label = nm)
    # add_uncertainty() keeps the view
    testthat::expect_identical(
      attr(unc@metadata$inf_reps, "selection"),
      attr(sub@metadata$inf_reps, "selection")
    )
  }
})

testthat::test_that("test_usage of a subset with draws gives the results of a copy", {
  testthat::skip_if_not_installed("edgeR")
  testthat::skip_if_not_installed("limma")
  x <- usage_experiment()
  iso <- x$iso
  draws <- withr::with_seed(1L, {
    lapply(1:10, function(k) {
      m <- tx_assay(iso, "counts")
      m[] <- stats::rpois(length(m), as.vector(m) + 1)
      m
    })
  })
  attr(draws, "samp_type") <- "gibbs"
  iso@metadata$inf_reps <- draws
  keep_tx <- rownames(tx_assay(iso))[-seq(3L, nrow(tx_assay(iso)), by = 7L)]
  sub <- iso[keep_tx, c(6L, 1L, 2L, 4L, 5L)]
  copy <- with_copied_draws(sub, draws)
  res <- test_usage(sub, group = "condition")
  testthat::expect_true(unique(res$scaled))
  testthat::expect_identical(res, test_usage(copy, group = "condition"))
})

testthat::test_that("a subset without kept draws keeps only the moments", {
  x <- uncertainty_experiment(inferential = TRUE)
  iso <- x$iso
  # Draws of the expressed transcripts only (as add_uncertainty() keeps),
  # with the moments of every transcript
  draws <- iso@metadata$inf_reps
  testthat::expect_false(is.null(attr(draws, "moments")))
  cut <- lapply(draws, function(m) m[c("tA1", "tA2"), , drop = FALSE])
  attributes(cut) <- attributes(draws)
  iso@metadata$inf_reps <- cut
  tx <- c("tA3", "tA4", "tA5", "tA6", "tB1", "tB2")
  sub <- suppressMessages(iso[tx, c("s1", "s2", "s4", "s5")])
  reps <- sub@metadata$inf_reps
  # Plain zero-row copies with the moments of the kept samples: the subset
  # does not hold the draws of the parent
  testthat::expect_null(attr(reps, "selection"))
  testthat::expect_identical(
    unique(vapply(reps, dim, integer(2L), USE.NAMES = FALSE)[1L, ]),
    0L
  )
  testthat::expect_identical(
    colnames(attr(reps, "moments")$mean),
    c("s1", "s2", "s4", "s5")
  )
  # The functions that need the draws see an object without them
  none <- sub
  none@metadata$inf_reps <- NULL
  testthat::expect_identical(
    productive_fraction(sub, group = "condition"),
    productive_fraction(none, group = "condition")
  )
  testthat::expect_error(
    add_uncertainty(sub),
    class = "isoformic_no_inferential_replicates"
  )
  testthat::expect_error(
    draws_view(sub),
    class = "isoformic_no_inferential_replicates"
  )
  # The overdispersion comes from the moments
  od <- inferential_overdispersion(sub)
  testthat::expect_setequal(od$transcript_id, tx)
  testthat::expect_identical(
    tx_assay(scale_counts(sub), "counts_scaled"),
    tx_assay(scale_counts(with_copied_draws(sub, cut)), "counts_scaled")
  )
  out <- fs::path(withr::local_tempdir(), "sub")
  write_isoformic(sub, out)
  testthat::expect_identical(read_isoformic(out)@metadata$inf_reps, reps)
})

testthat::test_that("list operations on the draws of a subset keep its selection", {
  x <- uncertainty_experiment(inferential = TRUE)
  iso <- x$iso
  draws <- iso@metadata$inf_reps
  tx <- c("tA1", "tA2", "tA3", "tB1")
  samples <- c("s1", "s2", "s4", "s5")
  sub <- suppressMessages(iso[tx, samples])
  copy <- with_copied_draws(sub, draws)
  edits <- list(
    rounded = function(reps) lapply(reps, round),
    thinned = function(reps) reps[1:10]
  )
  for (nm in names(edits)) {
    view <- sub
    view@metadata$inf_reps <- edits[[nm]](view@metadata$inf_reps)
    # The attribute is gone, the matrices are still those of the parent
    testthat::expect_null(attr(view@metadata$inf_reps, "selection"))
    testthat::expect_gt(nrow(view@metadata$inf_reps[[1L]]), length(tx))
    old <- copy
    old@metadata$inf_reps <- edits[[nm]](old@metadata$inf_reps)
    testthat::expect_identical(
      draws_view(view)$ids,
      draws_view(old)$ids,
      label = nm
    )
    testthat::expect_identical(draws_array(view), draws_array(old), label = nm)
    testthat::expect_identical(
      inferential_overdispersion(view),
      inferential_overdispersion(old),
      label = nm
    )
    testthat::expect_identical(
      productive_fraction(view, group = "condition"),
      productive_fraction(old, group = "condition"),
      label = nm
    )
    unc <- suppressMessages(add_uncertainty(view))
    unc_old <- suppressMessages(add_uncertainty(old))
    testthat::expect_identical(unc@dea, unc_old@dea, label = nm)
    testthat::expect_setequal(unc@dea[["uncertainty"]]$transcript_id, tx)
    out <- fs::path(withr::local_tempdir(), nm)
    write_isoformic(view, out)
    back <- read_isoformic(out)@metadata$inf_reps
    testthat::expect_identical(dim(back[[1L]]), c(4L, 4L), label = nm)
    testthat::expect_equal(
      lapply(back, identity),
      old@metadata$inf_reps,
      label = nm
    )
  }
  # Draws that only cover the object are read as they are
  testthat::expect_identical(object_inf_reps(iso), draws)
  testthat::expect_identical(object_inf_reps(copy), copy@metadata$inf_reps)
})

testthat::test_that("the table helpers read the draws of x after a list operation", {
  iso <- uncertainty_experiment(inferential = TRUE)$iso
  sub <- suppressMessages(iso[c("tA2", "tA1"), c("s1", "s2")])
  reps <- sub@metadata$inf_reps
  n_draws <- length(reps)
  expected <- inferential_replicates_table(reps)
  testthat::expect_identical(nrow(expected), 2L * 2L * n_draws)
  testthat::expect_identical(inferential_replicates_table(x = sub), expected)
  # lapply() and `[` drop the selection: the list alone gives the parent
  plain <- lapply(reps, identity)
  testthat::expect_null(attr(plain, "selection"))
  testthat::expect_identical(
    nrow(inferential_replicates_table(plain)),
    length(iso@metadata$inf_reps[[1L]]) * n_draws
  )
  testthat::expect_identical(
    inferential_replicates_table(plain, x = sub),
    expected
  )
  first <- reps[1:5]
  testthat::expect_identical(
    inferential_replicates_table(first, x = sub),
    expected[expected$draw <= 5L, ]
  )
  # The moments: the kept samples of the unmodified list or of x
  moments <- inferential_moments_table(reps)
  testthat::expect_identical(unique(moments$sample_id), c("s1", "s2"))
  testthat::expect_identical(inferential_moments_table(x = sub), moments)
  # A list operation drops the moments too: an error, not other samples
  testthat::expect_error(
    inferential_moments_table(plain, x = sub),
    class = "isoformic_no_replicate_moments"
  )
  testthat::expect_error(
    inferential_replicates_table(reps, x = data.frame()),
    class = "isoformic_invalid_inferential_replicates"
  )
})
