testthat::test_that("read_terminus_groups reads the clusters.txt that Terminus writes", {
  # Five lines of the output of terminus collapse 0.1.0 on Salmon 2.8.0 output
  file <- testthat::test_path("testdata", "terminus_clusters.txt")
  groups <- read_terminus_groups(file)
  testthat::expect_identical(
    colnames(groups),
    c("group_id", "transcript_id", "group_size")
  )
  testthat::expect_identical(nrow(groups), 10L)
  testthat::expect_identical(unique(groups$group_size), 2L)
  testthat::expect_identical(attr(groups, "source"), "terminus")
  testthat::expect_identical(
    groups$transcript_id[groups$group_id == "NewTr459365"],
    c("ENST00000831217.1", "ENST00000831165.1")
  )
  # The folder of a sample, and the folder that holds the sample folder
  out <- withr::local_tempdir()
  sample_dir <- fs::path(out, "sample1")
  fs::dir_create(sample_dir)
  fs::file_copy(file, fs::path(sample_dir, "clusters.txt"))
  testthat::expect_equal(read_terminus_groups(sample_dir), groups)
  testthat::expect_equal(read_terminus_groups(out), groups)
})

testthat::test_that("read_terminus_groups strips bars and rejects other files", {
  clusters <- withr::local_tempfile(fileext = ".txt")
  writeLines(c("G1,tA1|gA|x,tA2|gA|x,tA3|gA|x"), clusters)
  groups <- read_terminus_groups(clusters)
  testthat::expect_identical(groups$transcript_id, c("tA1", "tA2", "tA3"))
  testthat::expect_identical(unique(groups$group_size), 3L)
  raw <- read_terminus_groups(clusters, ignore_after_bar = FALSE)
  testthat::expect_identical(raw$transcript_id[1L], "tA1|gA|x")

  writeLines(c("G1,tA1"), clusters)
  testthat::expect_error(
    read_terminus_groups(clusters),
    class = "isoformic_invalid_terminus_groups"
  )
  writeLines(character(), clusters)
  testthat::expect_error(
    read_terminus_groups(clusters),
    class = "isoformic_invalid_terminus_groups"
  )
  empty <- withr::local_tempdir()
  testthat::expect_error(
    read_terminus_groups(empty),
    class = "isoformic_invalid_terminus_groups"
  )
})

testthat::test_that("external groups replace the correlation groups in add_uncertainty", {
  x <- uncertainty_experiment()
  base <- add_uncertainty(x$iso, paths = x$fixture$dirs)
  unc <- tx_uncertainty(base)
  testthat::expect_false(is.na(unc$ambiguity_group[unc$transcript_id == "tA1"]))
  testthat::expect_identical(
    base@metadata$uncertainty_ambiguity_source,
    "correlation"
  )

  clusters <- withr::local_tempfile(fileext = ".txt")
  writeLines("NewTr1,tA3,tA4", clusters)
  iso <- add_uncertainty(x$iso, paths = x$fixture$dirs, groups = clusters)
  unc <- tx_uncertainty(iso)
  testthat::expect_identical(
    unc$ambiguity_group[unc$transcript_id == "tA3"],
    "NewTr1"
  )
  testthat::expect_identical(
    unc$ambiguity_group[unc$transcript_id == "tA4"],
    "NewTr1"
  )
  testthat::expect_true(is.na(unc$ambiguity_group[unc$transcript_id == "tA1"]))
  testthat::expect_identical(
    iso@metadata$uncertainty_ambiguity_source,
    "terminus"
  )
  testthat::expect_match(uncertainty_methods(iso), "Terminus", fixed = TRUE)

  # The confidence classes and collapse_groups follow the new groups
  conf <- de_confidence(iso, group = "condition", collapse_groups = TRUE)
  testthat::expect_true("NewTr1" %in% conf$feature_id)
  testthat::expect_false(any(c("tA3", "tA4") %in% conf$feature_id))
  plain <- de_confidence(iso, group = "condition")
  testthat::expect_true(all(plain$ambiguous[
    plain$feature_id %in% c("tA3", "tA4")
  ]))
  testthat::expect_false(any(plain$ambiguous[
    plain$feature_id %in% c("tA1", "tA2")
  ]))

  # A table works too, and the source is then external
  table_groups <- data.frame(transcript_id = c("tA3", "tA4"), group_id = "G")
  iso2 <- add_uncertainty(x$iso, paths = x$fixture$dirs, groups = table_groups)
  testthat::expect_identical(
    iso2@metadata$uncertainty_ambiguity_source,
    "external"
  )
  testthat::expect_match(
    uncertainty_methods(iso2),
    "external grouping",
    fixed = TRUE
  )
  testthat::expect_error(
    add_uncertainty(
      x$iso,
      paths = x$fixture$dirs,
      groups = data.frame(
        transcript_id = c("tA3", "tA3"),
        group_id = c("G", "H")
      )
    ),
    class = "isoformic_invalid_groups"
  )
  testthat::expect_error(
    add_uncertainty(x$iso, paths = x$fixture$dirs, groups = data.frame(a = 1)),
    class = "isoformic_invalid_groups"
  )
})

# Salmon-like directories with replicates and equivalence classes, with or
# without the class weights
terminus_input <- function(root, names = c("s1", "s2"), weights = TRUE) {
  vapply(
    names,
    function(nm) {
      dir <- fs::path(root, nm)
      fs::dir_create(fs::path(dir, "aux_info", "bootstrap"))
      writeLines("Name\tLength", fs::path(dir, "quant.sf"))
      con <- gzfile(
        fs::path(dir, "aux_info", "bootstrap", "bootstraps.gz"),
        "wb"
      )
      writeBin(c(1, 2, 3), con)
      close(con)
      class_line <- if (weights) "2\t0\t1\t0.4\t0.6\t10" else "2\t0\t1\t10"
      con <- gzfile(fs::path(dir, "aux_info", "eq_classes.txt.gz"), "w")
      writeLines(c("2", "1", "tA1", "tA2", class_line), con)
      close(con)
      as.character(dir)
    },
    character(1L)
  )
}

testthat::test_that("run_terminus groups every sample, then builds the consensus", {
  testthat::skip_if_not_installed("condathis")
  root <- withr::local_tempdir()
  dirs <- terminus_input(root)
  out <- fs::path(root, "terminus")
  calls <- new.env(parent = emptyenv())
  calls$value <- list()
  testthat::local_mocked_bindings(
    env_exists = function(...) TRUE,
    run = function(cmd, ...) {
      args <- c(...)
      calls$value[[length(calls$value) + 1L]] <- list(
        cmd = cmd,
        args = args,
        env = list(...)$env_name
      )
      out_dir <- args[[which(args == "--out") + 1L]]
      if (identical(args[[1L]], "group")) {
        sample <- basename(args[[which(args == "--dir") + 1L]])
        fs::dir_create(fs::path(out_dir, sample))
        writeLines("0,1", fs::path(out_dir, sample, "groups.txt"))
        writeLines("big", fs::path(out_dir, sample, "delta.log"))
        writeLines("log", fs::path(out_dir, sample, "infrv.log"))
      } else {
        for (sample in c("s1", "s2")) {
          writeLines(
            "NewTr1,tA1,tA2",
            fs::path(out_dir, sample, "clusters.txt")
          )
        }
      }
      invisible(NULL)
    },
    .package = "condathis"
  )
  groups <- suppressMessages(run_terminus(
    dirs,
    output_dir = out,
    consensus_thresh = 0.7,
    env_name = "my-terminus"
  ))
  testthat::expect_identical(groups$transcript_id, c("tA1", "tA2"))
  testthat::expect_identical(attr(groups, "source"), "terminus")
  testthat::expect_identical(attr(groups, "output_dir"), as.character(out))
  # Two group calls, one per sample, then one collapse on both
  testthat::expect_length(calls$value, 3L)
  testthat::expect_identical(
    vapply(calls$value, function(x) x$args[[1L]], ""),
    c("group", "group", "collapse")
  )
  testthat::expect_true(all(
    vapply(calls$value, function(x) x$cmd, "") == "terminus"
  ))
  testthat::expect_true(all(
    vapply(calls$value, function(x) x$env, "") == "my-terminus"
  ))
  collapse <- calls$value[[3L]]$args
  testthat::expect_true(all(unname(dirs) %in% collapse))
  testthat::expect_identical(
    collapse[[which(collapse == "--consensus-thresh") + 1L]],
    "0.7"
  )
  # The large logs are gone, the groups and the consensus remain
  testthat::expect_false(any(fs::file_exists(fs::path(
    out,
    c("s1", "s2"),
    "delta.log"
  ))))
  testthat::expect_false(any(fs::file_exists(fs::path(
    out,
    c("s1", "s2"),
    "infrv.log"
  ))))
  testthat::expect_true(all(fs::file_exists(fs::path(
    out,
    c("s1", "s2"),
    "clusters.txt"
  ))))
})

testthat::test_that("run_terminus checks its input before running anything", {
  testthat::skip_if_not_installed("condathis")
  root <- withr::local_tempdir()
  ran <- new.env(parent = emptyenv())
  ran$value <- FALSE
  testthat::local_mocked_bindings(
    env_exists = function(...) TRUE,
    run = function(...) {
      ran$value <- TRUE
      invisible(NULL)
    },
    .package = "condathis"
  )
  no_weights <- terminus_input(fs::path(root, "a"), weights = FALSE)
  testthat::expect_error(
    run_terminus(no_weights, output_dir = fs::path(root, "out")),
    class = "isoformic_no_equivalence_weights"
  )
  with_weights <- terminus_input(fs::path(root, "b"))
  fs::file_delete(fs::path(
    with_weights[[1L]],
    "aux_info",
    "bootstrap",
    "bootstraps.gz"
  ))
  testthat::expect_error(
    run_terminus(with_weights, output_dir = fs::path(root, "out")),
    class = "isoformic_no_inferential_replicates"
  )
  fs::file_delete(fs::path(with_weights[[2L]], "aux_info", "eq_classes.txt.gz"))
  testthat::expect_error(
    run_terminus(with_weights[2L], output_dir = fs::path(root, "out")),
    class = "isoformic_no_equivalence_classes"
  )
  same_name <- c(
    terminus_input(fs::path(root, "c"), names = "s1"),
    terminus_input(fs::path(root, "d"), names = "s1")
  )
  testthat::expect_error(
    run_terminus(
      stats::setNames(same_name, c("x", "y")),
      output_dir = fs::path(root, "out")
    ),
    class = "isoformic_terminus_input"
  )
  oarfish <- fs::path(
    fs::path_package("isoformic", "extdata", "oarfish"),
    "control_1"
  )
  testthat::expect_error(
    run_terminus(oarfish, output_dir = fs::path(root, "out")),
    class = "isoformic_terminus_input"
  )
  testthat::expect_false(ran$value)
})

testthat::test_that("run_terminus groups samples in parallel and reads plain-text classes", {
  testthat::skip_if_not_installed("condathis")
  testthat::skip_on_os("windows")
  root <- withr::local_tempdir()
  dirs <- terminus_input(root, names = c("s1", "s2", "s3"))
  # Salmon can write the classes uncompressed
  eq <- fs::path(dirs[["s3"]], "aux_info", "eq_classes.txt.gz")
  writeLines(
    readLines(gzfile(eq)),
    fs::path(dirs[["s3"]], "aux_info", "eq_classes.txt")
  )
  fs::file_delete(eq)
  out <- fs::path(root, "terminus")
  testthat::local_mocked_bindings(
    env_exists = function(...) TRUE,
    run = function(cmd, ...) {
      args <- c(...)
      out_dir <- args[[which(args == "--out") + 1L]]
      if (identical(args[[1L]], "group")) {
        sample <- basename(args[[which(args == "--dir") + 1L]])
        fs::dir_create(fs::path(out_dir, sample))
        writeLines("big", fs::path(out_dir, sample, "delta.log"))
      } else {
        for (sample in c("s1", "s2", "s3")) {
          writeLines(
            "NewTr1,tA1,tA2",
            fs::path(out_dir, sample, "clusters.txt")
          )
        }
      }
      invisible(NULL)
    },
    .package = "condathis"
  )
  # R CMD check caps forked workers at 2 (`_R_CHECK_LIMIT_CORES_`).
  groups <- suppressMessages(run_terminus(dirs, output_dir = out, workers = 2L))
  testthat::expect_identical(groups$group_id, c("NewTr1", "NewTr1"))
  testthat::expect_false(any(fs::file_exists(fs::path(
    out,
    c("s1", "s2", "s3"),
    "delta.log"
  ))))
  kept <- suppressMessages(
    run_terminus(dirs, output_dir = fs::path(root, "kept"), keep_logs = TRUE)
  )
  testthat::expect_true(all(fs::file_exists(fs::path(
    root,
    "kept",
    c("s1", "s2", "s3"),
    "delta.log"
  ))))
  testthat::expect_identical(nrow(kept), 2L)
})

testthat::test_that("run_terminus reports failed samples and environment errors", {
  testthat::skip_if_not_installed("condathis")
  testthat::skip_on_os("windows")
  root <- withr::local_tempdir()
  dirs <- terminus_input(root, names = c("s1", "s2"))
  testthat::local_mocked_bindings(
    env_exists = function(...) TRUE,
    run = function(cmd, ...) {
      args <- c(...)
      if (
        identical(args[[1L]], "group") &&
          grepl("s2$", args[[which(args == "--dir") + 1L]])
      ) {
        stop("terminus crashed", call. = FALSE)
      }
      invisible(NULL)
    },
    .package = "condathis"
  )
  testthat::expect_error(
    run_terminus(dirs, output_dir = fs::path(root, "out"), workers = 2L),
    class = "isoformic_terminus_failed"
  )
})

testthat::test_that("check_terminus_env creates the environment and reports failures", {
  testthat::skip_if_not_installed("condathis")
  created <- new.env(parent = emptyenv())
  created$args <- NULL
  testthat::local_mocked_bindings(
    env_exists = function(...) FALSE,
    create_env = function(...) {
      created$args <- list(...)
      invisible(NULL)
    },
    .package = "condathis"
  )
  check_terminus_env("my-env")
  testthat::expect_identical(created$args$env_name, "my-env")
  testthat::expect_identical(created$args$package, "bioconda::terminus==0.1.0")
  testthat::local_mocked_bindings(
    create_env = function(...) stop("no network", call. = FALSE),
    .package = "condathis"
  )
  testthat::expect_error(
    check_terminus_env("my-env"),
    class = "isoformic_terminus_env_create_error"
  )
})
