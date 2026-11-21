testthat::test_that("write_isoformic can overwrite the directory it was read from", {
  iso <- example_isoformic()
  dir <- withr::local_tempdir()
  out <- file.path(dir, "exp")
  iso <- write_isoformic(iso, out, copy_annotation = TRUE)
  again <- read_isoformic(out)
  # Same directory, with and without copying the annotation
  testthat::expect_no_error(
    again <- write_isoformic(again, overwrite = TRUE, copy_annotation = TRUE)
  )
  testthat::expect_no_error(
    again <- write_isoformic(again, overwrite = TRUE)
  )
  back <- read_isoformic(out)
  testthat::expect_identical(dim(back), dim(iso))
  testthat::expect_true(fs::dir_exists(fs::path(out, "annotation")))
  testthat::expect_length(fs::dir_ls(dir, glob = "*.writing-*"), 0L)
})

testthat::test_that("write_isoformic refuses to delete other directories", {
  iso <- example_isoformic()
  dir <- withr::local_tempdir()
  writeLines("keep me", file.path(dir, "notes.txt"))
  testthat::expect_error(
    write_isoformic(iso, dir, overwrite = TRUE),
    class = "isoformic_not_experiment_dir"
  )
  testthat::expect_true(file.exists(file.path(dir, "notes.txt")))
})

# Inferential replicates of the example object: draws for the first half of
# the transcripts, moments for all of them, as read_inferential_replicates()
# records them
example_inf_reps <- function(iso, n_draws = 3L) {
  counts <- tx_assay(iso, "counts")
  kept <- counts[seq_len(nrow(counts) %/% 2L), , drop = FALSE]
  withr::with_seed(1L, {
    reps <- lapply(seq_len(n_draws), function(k) {
      kept * stats::runif(length(kept), 0.5, 1.5)
    })
    mom_mean <- counts * stats::runif(length(counts), 0.9, 1.1)
  })
  attr(reps, "samp_type") <- "gibbs"
  attr(reps, "moments") <- list(
    mean = mom_mean,
    var = mom_mean * 0.25,
    n = as.integer(n_draws)
  )
  reps
}

testthat::test_that("write_isoformic and read_isoformic round-trip the inferential replicates", {
  iso <- example_isoformic()
  iso@metadata$inf_reps <- example_inf_reps(iso)
  out <- fs::path(withr::local_tempdir(), "exp")
  iso <- write_isoformic(iso, out, copy_annotation = TRUE)
  testthat::expect_true(all(fs::file_exists(fs::path(
    out,
    "inf_reps",
    c("draw_0001.parquet", "draw_0002.parquet", "draw_0003.parquet")
  ))))
  testthat::expect_true(all(fs::file_exists(fs::path(
    out,
    "inf_reps",
    "moments",
    c("mean.parquet", "var.parquet")
  ))))
  # One wide table per draw: transcript_id and one column per sample
  draw <- read_parquet(fs::path(out, "inf_reps", "draw_0002.parquet"))
  testthat::expect_identical(
    colnames(draw),
    c("transcript_id", colnames(iso@metadata$inf_reps[[2L]]))
  )
  back <- read_isoformic(out)
  testthat::expect_identical(back@metadata$inf_reps, iso@metadata$inf_reps)
  # The keys that describe the replicates do not stay in the metadata
  testthat::expect_false(any(
    c("inf_reps_n_draws", "inf_reps_samp_type", "inf_reps_moments_n") %in%
      names(back@metadata)
  ))
  testthat::expect_identical(
    setdiff(names(back@metadata), "inf_reps"),
    setdiff(names(iso@metadata), "inf_reps")
  )
  # Writing the object read back gives the same replicates again
  again <- fs::path(withr::local_tempdir(), "exp")
  write_isoformic(back, again)
  testthat::expect_identical(
    read_isoformic(again)@metadata$inf_reps,
    iso@metadata$inf_reps
  )
})

testthat::test_that("subsets keep their inferential replicates on disk", {
  iso <- example_isoformic()
  iso@metadata$inf_reps <- example_inf_reps(iso)
  tx <- rownames(iso@metadata$inf_reps[[1L]])
  samples <- col_data(iso)$sample_id
  sub <- iso[tx[seq(1L, length(tx), by = 3L)], samples[c(5L, 1L, 9L, 2L)]]
  # The subset shares the draws: only the selected rows and columns are
  # written, and they read back as plain matrices of those
  kept <- inf_reps_materialize(sub@metadata$inf_reps)
  testthat::expect_identical(colnames(kept[[1L]]), samples[c(5L, 1L, 9L, 2L)])
  testthat::expect_identical(rownames(kept[[1L]]), tx[seq(1L, length(tx), 3L)])
  out <- fs::path(withr::local_tempdir(), "sub")
  write_isoformic(sub, out)
  draw <- read_parquet(fs::path(out, "inf_reps", "draw_0001.parquet"))
  testthat::expect_identical(draw$transcript_id, rownames(kept[[1L]]))
  testthat::expect_identical(
    colnames(draw),
    c("transcript_id", samples[c(5L, 1L, 9L, 2L)])
  )
  back <- read_isoformic(out)@metadata$inf_reps
  testthat::expect_identical(back, kept)
  testthat::expect_null(attr(back, "selection"))
  testthat::expect_identical(
    colnames(attr(back, "moments")$mean),
    samples[c(5L, 1L, 9L, 2L)]
  )
  # A subset of the read subset gives what a subset of the parent gives
  inner <- tx[seq(1L, length(tx), by = 6L)]
  testthat::expect_identical(
    inf_reps_materialize(
      read_isoformic(out)[inner, samples[c(9L, 5L)]]@metadata$inf_reps
    ),
    inf_reps_materialize(sub[inner, samples[c(9L, 5L)]]@metadata$inf_reps)
  )
})

testthat::test_that("replicate moments without draws round-trip", {
  iso <- example_isoformic()
  iso@metadata$inf_reps <- example_inf_reps(iso)
  # Zero-row draws that keep their moments
  draws <- iso@metadata$inf_reps
  reps <- lapply(draws, function(m) m[0L, , drop = FALSE])
  attributes(reps) <- attributes(draws)
  testthat::expect_false(is.null(attr(reps, "moments")))
  iso@metadata$inf_reps <- reps
  out <- fs::path(withr::local_tempdir(), "moments_only")
  write_isoformic(iso, out)
  testthat::expect_identical(read_isoformic(out)@metadata$inf_reps, reps)
  # A subset of transcripts without draws keeps zero-row draws and moments
  iso@metadata$inf_reps <- draws
  no_draws <- setdiff(rownames(tx_assay(iso)), rownames(draws[[1L]]))
  sub <- iso[no_draws, ]
  sub_reps <- sub@metadata$inf_reps
  testthat::expect_identical(nrow(sub_reps[[1L]]), 0L)
  testthat::expect_false(is.null(attr(sub_reps, "moments")))
  out_sub <- fs::path(withr::local_tempdir(), "no_draws")
  write_isoformic(sub, out_sub)
  testthat::expect_identical(
    read_isoformic(out_sub)@metadata$inf_reps,
    sub_reps
  )
})

testthat::test_that("replicates without attributes and of integer counts round-trip", {
  iso <- example_isoformic()
  counts <- tx_assay(iso, "counts")
  reps <- lapply(1:2, function(k) {
    m <- round(counts) + k
    storage.mode(m) <- "integer"
    m
  })
  iso@metadata$inf_reps <- reps
  out <- fs::path(withr::local_tempdir(), "plain")
  write_isoformic(iso, out)
  back <- read_isoformic(out)@metadata$inf_reps
  testthat::expect_identical(back, reps)
  testthat::expect_null(attr(back, "samp_type"))
  testthat::expect_null(attr(back, "moments"))
})

testthat::test_that("objects and directories without replicates read without them", {
  iso <- example_isoformic()
  iso@metadata$note <- "kept"
  out <- fs::path(withr::local_tempdir(), "none")
  write_isoformic(iso, out)
  testthat::expect_false(fs::dir_exists(fs::path(out, "inf_reps")))
  back <- read_isoformic(out)
  testthat::expect_null(back@metadata$inf_reps)
  testthat::expect_identical(back@metadata$note, "kept")
  # A directory of an earlier version (the bundled example) has no
  # inf_reps/ and keeps its metadata as it is
  bundled <- system.file("extdata", "example_experiment", package = "isoformic")
  testthat::expect_false(fs::dir_exists(fs::path(bundled, "inf_reps")))
  old <- read_isoformic(bundled, annot_path = iso@annot_path)
  testthat::expect_null(old@metadata$inf_reps)
})

testthat::test_that("replicates that cannot be written are dropped with a warning", {
  iso <- example_isoformic()
  iso@metadata$inf_reps <- list(as.data.frame(tx_assay(iso, "counts")))
  out <- fs::path(withr::local_tempdir(), "bad")
  testthat::expect_warning(
    write_isoformic(iso, out),
    class = "isoformic_inf_reps_not_written"
  )
  testthat::expect_false(fs::dir_exists(fs::path(out, "inf_reps")))
  testthat::expect_null(read_isoformic(out)@metadata$inf_reps)
})

testthat::test_that("read_isoformic reports a missing draw file", {
  iso <- example_isoformic()
  iso@metadata$inf_reps <- example_inf_reps(iso, n_draws = 2L)
  out <- fs::path(withr::local_tempdir(), "broken")
  write_isoformic(iso, out)
  fs::file_delete(fs::path(out, "inf_reps", "draw_0002.parquet"))
  testthat::expect_error(
    read_isoformic(out),
    class = "isoformic_invalid_experiment_dir"
  )
})

testthat::test_that("metadata entries named like the replicate keys stay without replicates", {
  iso <- example_isoformic()
  iso@metadata$inf_reps_n_draws <- "mine"
  out <- fs::path(withr::local_tempdir(), "collide_none")
  write_isoformic(iso, out)
  testthat::expect_identical(
    read_isoformic(out)@metadata$inf_reps_n_draws,
    "mine"
  )
  # With replicates the key describes them, and the user entry is replaced
  iso@metadata$inf_reps <- example_inf_reps(iso, n_draws = 2L)
  out <- fs::path(withr::local_tempdir(), "collide_reps")
  testthat::expect_warning(
    write_isoformic(iso, out),
    class = "isoformic_metadata_replaced"
  )
  back <- read_isoformic(out)@metadata
  testthat::expect_null(back$inf_reps_n_draws)
  testthat::expect_identical(back$inf_reps, iso@metadata$inf_reps)
})

testthat::test_that("the storage mode of the moments round-trips", {
  iso <- example_isoformic()
  reps <- example_inf_reps(iso, n_draws = 2L)
  mom <- attr(reps, "moments")
  storage.mode(mom$mean) <- "integer"
  mom$n <- 2
  attr(reps, "moments") <- mom
  iso@metadata$inf_reps <- reps
  out <- fs::path(withr::local_tempdir(), "mom_types")
  write_isoformic(iso, out)
  back <- read_isoformic(out)@metadata$inf_reps
  testthat::expect_identical(back, reps)
  testthat::expect_type(attr(back, "moments")$mean, "integer")
  testthat::expect_type(attr(back, "moments")$n, "double")
  testthat::expect_identical(encode_moments_n(30L), "30")
  testthat::expect_identical(decode_moments_n(encode_moments_n(30L)), 30L)
  for (n in c(30, 2.5, 1e20, NaN, NA_real_)) {
    testthat::expect_identical(decode_moments_n(encode_moments_n(n)), n)
  }
  testthat::expect_identical(
    decode_moments_n(encode_moments_n(NA_integer_)),
    NA_integer_
  )
})

testthat::test_that("properties of the replicates that are not written give a warning", {
  iso <- example_isoformic()
  base <- example_inf_reps(iso, n_draws = 2L)
  named <- base
  names(named) <- c("a", "b")
  extra <- base
  attr(extra, "foo") <- "bar"
  dn_names <- base
  names(dimnames(dn_names[[1L]])) <- c("transcript", "sample")
  for (reps in list(named, extra, dn_names)) {
    iso@metadata$inf_reps <- reps
    out <- fs::path(withr::local_tempdir(), "changed")
    testthat::expect_warning(
      write_isoformic(iso, out),
      class = "isoformic_inf_reps_changed"
    )
    back <- read_isoformic(out)@metadata$inf_reps
    testthat::expect_identical(back, base)
  }
  # The replicates of the package's readers give no warning
  iso@metadata$inf_reps <- base
  out <- fs::path(withr::local_tempdir(), "unchanged")
  testthat::expect_no_warning(write_isoformic(iso, out))
})

testthat::test_that("small draws do not pay for a full collection each", {
  iso <- example_isoformic()
  m <- tx_assay(iso, "counts")[1:2, ]
  reps <- rep(list(m), 60L)
  iso@metadata$inf_reps <- reps
  out <- fs::path(withr::local_tempdir(), "many")
  gc_calls <- new.env(parent = emptyenv())
  gc_calls$n <- 0L
  testthat::local_mocked_bindings(
    full_gc = function() {
      gc_calls$n <- gc_calls$n + 1L
    }
  )
  write_isoformic(iso, out)
  testthat::expect_identical(read_isoformic(out)@metadata$inf_reps, reps)
  testthat::expect_identical(gc_calls$n, 0L)
})

testthat::test_that("gc_pacer collects once the cells reach the threshold", {
  gc_calls <- new.env(parent = emptyenv())
  gc_calls$n <- 0L
  testthat::local_mocked_bindings(
    full_gc = function() {
      gc_calls$n <- gc_calls$n + 1L
    }
  )
  collect <- gc_pacer(threshold = 100)
  for (k in 1:9) {
    collect(10)
  }
  testthat::expect_identical(gc_calls$n, 0L)
  collect(10)
  testthat::expect_identical(gc_calls$n, 1L)
  collect(1000)
  testthat::expect_identical(gc_calls$n, 2L)
  testthat::expect_identical(matrix_cells(matrix(0, 3, 4)), 15)
})

testthat::test_that("print shows data_path relative to the working directory", {
  iso <- example_isoformic()
  dir <- withr::local_tempdir()
  withr::local_dir(dir)
  iso <- write_isoformic(iso, fs::path(dir, "data", "exp"))
  testthat::expect_true(fs::is_absolute_path(iso@data_path))
  testthat::expect_output(print(iso), "on disk    : data/exp", fixed = TRUE)
  outside <- withr::local_tempdir()
  withr::local_dir(outside)
  testthat::expect_output(
    print(iso),
    paste0("on disk    : ", iso@data_path),
    fixed = TRUE
  )
})
