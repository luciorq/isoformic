# Reference: a transcript is identifiable when its unit vector is in the row
# space of the class by transcript matrix
brute_identifiable <- function(tid, cls, n_tx) {
  a <- matrix(0, max(cls), n_tx)
  a[cbind(cls, tid)] <- 1
  dec <- svd(a, nu = 0L)
  rank <- sum(dec$d > 1e-10 * dec$d[[1L]])
  res <- rowSums(dec$v[, seq_len(rank), drop = FALSE]^2) > 1 - 1e-8
  res[tabulate(tid, n_tx) == 0L] <- NA
  res
}

random_classes <- function(n_tx, n_cls, max_size) {
  size <- sample.int(max_size, n_cls, replace = TRUE)
  tid <- unlist(lapply(size, function(s) sort(sample.int(n_tx, s))))
  list(tid = tid, cls = rep.int(seq_len(n_cls), size))
}

write_eq_classes <- function(dir, tx, classes) {
  fs::dir_create(fs::path(dir, "aux_info"))
  writeLines(
    c(length(tx), length(classes), tx, classes),
    fs::path(dir, "aux_info", "eq_classes.txt")
  )
  as.character(dir)
}

testthat::test_that("hand-made classes give the expected identifiability", {
  # 1 and 2 always together; 3 alone; 4 and 5 separated by a class of 4
  tid <- c(1L, 2L, 3L, 4L, 5L, 4L)
  cls <- c(1L, 1L, 2L, 3L, 3L, 4L)
  testthat::expect_identical(
    identifiable_from_classes(tid, cls, n_tx = 6L),
    c(FALSE, FALSE, TRUE, TRUE, TRUE, NA)
  )
  # {1,2}, {2,3}, {1,3}: full rank although no class has one member
  tid <- c(1L, 2L, 2L, 3L, 1L, 3L)
  cls <- c(1L, 1L, 2L, 2L, 3L, 3L)
  testthat::expect_identical(
    identifiable_from_classes(tid, cls, n_tx = 3L),
    c(TRUE, TRUE, TRUE)
  )
  # {1,2}, {3,4}, {1,3}, {2,4}: a cycle of four with rank three
  tid <- c(1L, 2L, 3L, 4L, 1L, 3L, 2L, 4L)
  cls <- c(1L, 1L, 2L, 2L, 3L, 3L, 4L, 4L)
  testthat::expect_identical(
    identifiable_from_classes(tid, cls, n_tx = 4L),
    c(FALSE, FALSE, FALSE, FALSE)
  )
  testthat::expect_identical(
    identifiable_from_classes(integer(), integer(), n_tx = 2L),
    c(NA, NA)
  )
})

testthat::test_that("the reduction agrees with the full decomposition", {
  withr::local_seed(11L)
  for (i in seq_len(30L)) {
    n_tx <- sample(5:40, 1L)
    x <- random_classes(n_tx, n_cls = sample(3:30, 1L), max_size = 4L)
    testthat::expect_identical(
      identifiable_from_classes(x$tid, x$cls, n_tx),
      brute_identifiable(x$tid, x$cls, n_tx)
    )
  }
})

testthat::test_that("the iterative solver agrees with the dense one", {
  testthat::skip_if_not_installed("Matrix")
  withr::local_seed(12L)
  for (i in seq_len(10L)) {
    n_tx <- sample(20:60, 1L)
    x <- random_classes(n_tx, n_cls = sample(15:50, 1L), max_size = 4L)
    testthat::expect_identical(
      identifiable_from_classes(x$tid, x$cls, n_tx, max_dense = 1L),
      brute_identifiable(x$tid, x$cls, n_tx)
    )
  }
  # A block that does not converge is reported and left undetermined
  x <- random_classes(40L, n_cls = 30L, max_size = 4L)
  testthat::expect_message(
    res <- identifiable_from_classes(
      x$tid,
      x$cls,
      40L,
      max_dense = 1L,
      max_iterations = 1L
    ),
    class = "isoformic_identifiability_not_converged"
  )
  testthat::expect_true(anyNA(res[tabulate(x$tid, 40L) > 0L]))
})

testthat::test_that("equivalence_identifiability() summarises samples", {
  tmp <- withr::local_tempdir()
  tx <- c("tx1|g1", "tx2|g1", "tx3|g2", "tx4|g2")
  dirs <- c(
    s1 = write_eq_classes(
      fs::path(tmp, "s1"),
      tx,
      c("2\t0\t1\t10", "1\t2\t5")
    ),
    # weights are accepted and ignored; tx1 gets a class of its own
    s2 = write_eq_classes(
      fs::path(tmp, "s2"),
      tx,
      c("2\t0\t1\t0.4\t0.6\t10", "1\t0\t1\t3", "1\t2\t5")
    )
  )
  res <- equivalence_identifiability(dirs)
  testthat::expect_identical(res$transcript_id, c("tx1", "tx2", "tx3"))
  testthat::expect_identical(res$n_samples, c(2L, 2L, 2L))
  testthat::expect_identical(res$n_identifiable, c(1L, 1L, 2L))
  testthat::expect_identical(res$identifiable, c(FALSE, FALSE, TRUE))
  testthat::expect_equal(res$identifiable_fraction, c(0.5, 0.5, 1))
  kept <- equivalence_identifiability(
    dirs,
    transcripts = "tx3",
    ignore_after_bar = TRUE
  )
  testthat::expect_identical(kept$transcript_id, "tx3")
  testthat::expect_identical(
    equivalence_identifiability(dirs, ignore_after_bar = FALSE)$transcript_id,
    tx[1:3]
  )
  empty <- fs::dir_create(fs::path(tmp, "empty"))
  testthat::expect_error(
    equivalence_identifiability(c(s = as.character(empty))),
    class = "isoformic_no_equivalence_classes"
  )
})
