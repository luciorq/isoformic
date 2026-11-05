oarfish_example <- function() {
  dir <- fs::path_package("isoformic", "extdata", "oarfish")
  prefixes <- as.character(fs::path(
    dir,
    c("control_1", "control_2", "treated_1")
  ))
  names(prefixes) <- basename(prefixes)
  list(
    prefixes = prefixes,
    annot = fs::path_package(
      "isoformic",
      "extdata",
      "example_experiment",
      "annotation"
    )
  )
}

testthat::test_that("read_oarfish builds the assays from the real Oarfish output", {
  x <- oarfish_example()
  iso <- read_oarfish(x$prefixes, annot_path = x$annot)
  testthat::expect_identical(dim(iso), c(72L, 3L))
  quant <- read.delim(paste0(x$prefixes[["control_1"]], ".quant"))
  counts <- assay(iso, "counts")
  testthat::expect_equal(
    unname(counts[, "control_1"]),
    quant$num_reads
  )
  testthat::expect_equal(unname(assay(iso, "length")[, "control_1"]), quant$len)
  # TPM is computed from reads per base and sums to a million per sample
  testthat::expect_equal(unname(colSums(assay(iso, "abundance"))), rep(1e6, 3L))
  testthat::expect_identical(
    rownames(counts)[1L],
    sub("\\|.*$", "", quant$tname[1L])
  )
  # Directories with one .quant file and the .quant path itself are accepted
  one <- withr::local_tempdir()
  fs::file_copy(paste0(x$prefixes[["control_1"]], ".quant"), one)
  testthat::expect_identical(
    dim(read_oarfish(one, annot_path = x$annot)),
    c(72L, 1L)
  )
  testthat::expect_identical(
    dim(read_oarfish(
      paste0(x$prefixes[["control_2"]], ".quant"),
      annot_path = x$annot
    )),
    c(72L, 1L)
  )
})

testthat::test_that("Oarfish replicates, support and run metadata use the shared readers", {
  x <- oarfish_example()
  reps <- read_inferential_replicates(x$prefixes)
  testthat::expect_length(reps, 30L)
  testthat::expect_identical(attr(reps, "samp_type"), "bootstrap")
  testthat::expect_identical(dim(reps[[1L]]), c(72L, 3L))
  testthat::expect_identical(colnames(reps[[1L]]), names(x$prefixes))
  raw <- read_parquet(paste0(x$prefixes[["treated_1"]], ".infreps.pq"))
  testthat::expect_equal(
    unname(reps[[5L]][, "treated_1"]),
    raw[["bootstrap.4"]]
  )
  few <- read_inferential_replicates(x$prefixes, n_draws = 5L)
  testthat::expect_length(few, 5L)

  support <- equivalence_support(x$prefixes)
  ambig <- read.delim(paste0(x$prefixes[["control_1"]], ".ambig_info.tsv"))
  one <- equivalence_support(x$prefixes[["control_1"]])
  testthat::expect_equal(
    one$unique_fraction,
    ifelse(
      ambig$total_reads > 0,
      ambig$unique_reads / ambig$total_reads,
      NA_real_
    )
  )
  testthat::expect_true(all(is.na(support$n_classes)))
  testthat::expect_identical(support$n_samples[1L], 3L)
  testthat::expect_identical(
    equivalence_support(
      x$prefixes,
      transcripts = support$transcript_id[1:3]
    )$transcript_id,
    support$transcript_id[1:3]
  )

  qc <- salmon_run_qc(x$prefixes)
  testthat::expect_identical(attr(qc, "tool"), "oarfish")
  testthat::expect_identical(qc$num_replicates, rep(30, 3L))
  testthat::expect_equal(qc$percent_mapped, rep(100, 3L))
  testthat::expect_true(all(is.na(qc$qc_flag)))
})

testthat::test_that("add_uncertainty works on Oarfish output and names the tool", {
  x <- oarfish_example()
  iso <- read_oarfish(x$prefixes, annot_path = x$annot, inferential = TRUE)
  unc <- tx_uncertainty(iso)
  testthat::expect_gt(nrow(unc), 20L)
  testthat::expect_true(all(
    c("mean_inf_rv", "unique_fraction") %in% colnames(unc)
  ))
  testthat::expect_true(all(
    unc$unique_fraction[!is.na(unc$unique_fraction)] <= 1
  ))
  testthat::expect_true("oarfish_percent_mapped" %in% colnames(col_data(iso)))
  testthat::expect_false(any(grepl("^salmon_", colnames(col_data(iso)))))
  testthat::expect_match(uncertainty_methods(iso), "Oarfish", fixed = TRUE)
  testthat::expect_match(uncertainty_methods(iso), "bootstrap", fixed = TRUE)
})

testthat::test_that("Oarfish readers fail clearly on missing or mixed input", {
  x <- oarfish_example()
  copy <- withr::local_tempdir()
  fs::file_copy(paste0(x$prefixes[["control_1"]], ".quant"), copy)
  prefix <- file.path(copy, "control_1")
  testthat::expect_error(
    read_inferential_replicates(prefix),
    class = "isoformic_no_inferential_replicates"
  )
  testthat::expect_error(
    equivalence_support(prefix),
    class = "isoformic_no_equivalence_classes"
  )
  testthat::expect_error(
    salmon_run_qc(prefix),
    class = "isoformic_no_meta_info"
  )
  salmon_dir <- fs::path_package(
    "isoformic",
    "extdata",
    "mini_quants",
    "SRR11498039_quant"
  )
  testthat::expect_error(
    read_inferential_replicates(c(
      a = x$prefixes[["control_1"]],
      b = salmon_dir
    )),
    class = "isoformic_mixed_quantifiers"
  )
})

oarfish_sc_example <- function() {
  dir <- fs::path_package("isoformic", "extdata", "oarfish_sc")
  list(
    prefix = as.character(fs::path(dir, "out")),
    dir = as.character(dir),
    cells = utils::read.delim(fs::path(dir, "cells.tsv")),
    annot = fs::path_package(
      "isoformic",
      "extdata",
      "example_experiment",
      "annotation"
    )
  )
}

testthat::test_that("read_oarfish_sc_counts reads the real single-cell output of Oarfish", {
  x <- oarfish_sc_example()
  counts <- read_oarfish_sc_counts(x$prefix)
  testthat::expect_identical(
    colnames(counts),
    c("barcode", "transcript_id", "count")
  )
  mtx <- readLines(paste0(x$prefix, ".count.mtx"))
  dims <- as.integer(strsplit(mtx[3L], " ")[[1L]])
  testthat::expect_identical(nrow(counts), dims[[3L]])
  testthat::expect_setequal(unique(counts$barcode), x$cells$barcode)
  testthat::expect_false(any(grepl("|", counts$transcript_id, fixed = TRUE)))
  # 60 cells of 150 reads, all assigned
  testthat::expect_equal(sum(counts$count), 9000, tolerance = 1e-3)
  # The same table from the .count.mtx path and from the folder
  testthat::expect_equal(
    read_oarfish_sc_counts(paste0(x$prefix, ".count.mtx")),
    counts
  )
  testthat::expect_equal(read_oarfish_sc_counts(x$dir), counts)
  raw <- read_oarfish_sc_counts(x$prefix, ignore_after_bar = FALSE)
  testthat::expect_true(all(grepl("|", raw$transcript_id, fixed = TRUE)))
})

testthat::test_that("read_oarfish_sc sums the cells of each sample", {
  x <- oarfish_sc_example()
  iso <- read_oarfish_sc(x$prefix, cells = x$cells, annot_path = x$annot)
  testthat::expect_identical(dim(iso), c(72L, 3L))
  cd <- col_data(iso)
  testthat::expect_identical(
    cd$sample_id,
    c("control_1", "control_2", "treated_1")
  )
  testthat::expect_identical(cd$n_cells, rep(20L, 3L))
  testthat::expect_equal(cd$n_reads, rep(3000, 3L), tolerance = 1e-3)
  counts <- read_oarfish_sc_counts(x$prefix)
  counts$sample_id <- x$cells$sample_id[match(counts$barcode, x$cells$barcode)]
  by_hand <- tapply(
    counts$count,
    list(counts$transcript_id, counts$sample_id),
    sum
  )
  pseudo <- assay(iso, "counts")
  ids <- rownames(by_hand)
  testthat::expect_equal(
    unname(pseudo[ids, "treated_1"]),
    unname(by_hand[ids, "treated_1"])
  )
  testthat::expect_equal(unname(colSums(assay(iso, "abundance"))), rep(1e6, 3L))
  # The simulation raises three isoforms four-fold in the treated cells; the
  # one with enough reads to tell (RBPJ-203, about 80 in each control) shows it
  fasta <- readLines(gzfile(fs::path_package(
    "isoformic",
    "extdata",
    "example_transcripts.fa.gz"
  )))
  headers <- sub("^>([^|]*)\\|.*", "\\1", fasta[startsWith(fasta, ">")])
  boosted <- headers[c(2L, 5L, 9L)]
  control <- rowMeans(pseudo[boosted, c("control_1", "control_2")])
  enough <- boosted[control >= 20]
  testthat::expect_length(enough, 1L)
  testthat::expect_gt(pseudo[enough, "treated_1"] / control[[enough]], 3)

  # Metadata, subsets of cells and minimum numbers of cells
  meta <- data.frame(
    sample_id = c("control_1", "control_2", "treated_1"),
    condition = c("control", "control", "treated")
  )
  iso2 <- read_oarfish_sc(
    x$prefix,
    cells = x$cells,
    annot_path = x$annot,
    col_data = meta
  )
  testthat::expect_identical(col_data(iso2)$condition, meta$condition)
  some <- x$cells[c(1:5, 41:60), ]
  iso3 <- suppressMessages(
    read_oarfish_sc(
      x$prefix,
      cells = some,
      annot_path = x$annot,
      min_cells = 6L
    )
  )
  testthat::expect_identical(col_data(iso3)$sample_id, "treated_1")
  testthat::expect_message(
    read_oarfish_sc(x$prefix, cells = some, annot_path = x$annot),
    "Dropping 35 barcodes"
  )
})

testthat::test_that("read_oarfish_sc rejects bad input", {
  x <- oarfish_sc_example()
  testthat::expect_error(
    read_oarfish_sc(x$prefix, cells = data.frame(a = 1), annot_path = x$annot),
    class = "isoformic_invalid_oarfish_sc"
  )
  testthat::expect_error(
    read_oarfish_sc(
      x$prefix,
      cells = rbind(x$cells, x$cells[1L, ]),
      annot_path = x$annot
    ),
    class = "isoformic_invalid_oarfish_sc"
  )
  testthat::expect_error(
    read_oarfish_sc(
      x$prefix,
      cells = x$cells,
      annot_path = x$annot,
      min_cells = 100L
    ),
    class = "isoformic_invalid_oarfish_sc"
  )
  broken <- withr::local_tempdir()
  for (ext in c(".barcodes.txt", ".features.txt")) {
    fs::file_copy(paste0(x$prefix, ext), fs::path(broken, paste0("out", ext)))
  }
  mtx <- readLines(paste0(x$prefix, ".count.mtx"))
  writeLines(mtx[-length(mtx)], fs::path(broken, "out.count.mtx"))
  testthat::expect_error(
    read_oarfish_sc_counts(fs::path(broken, "out")),
    class = "isoformic_invalid_oarfish_sc"
  )
  writeLines(
    c("%%MatrixMarket matrix coordinate real general", "2 2 0"),
    fs::path(broken, "out.count.mtx")
  )
  testthat::expect_error(
    read_oarfish_sc_counts(fs::path(broken, "out")),
    class = "isoformic_invalid_oarfish_sc"
  )
  testthat::expect_error(
    read_oarfish_sc_counts(withr::local_tempdir()),
    class = "isoformic_invalid_oarfish_sc"
  )
})

testthat::test_that("read_oarfish_sc handles transcripts missing from the annotation", {
  x <- oarfish_sc_example()
  copy <- withr::local_tempdir()
  for (ext in c(".count.mtx", ".barcodes.txt", ".features.txt")) {
    fs::file_copy(paste0(x$prefix, ext), fs::path(copy, paste0("out", ext)))
  }
  features <- readLines(fs::path(copy, "out.features.txt"))
  features[[1L]] <- "NOT_IN_ANNOTATION|x|"
  writeLines(features, fs::path(copy, "out.features.txt"))
  testthat::expect_error(
    read_oarfish_sc(
      fs::path(copy, "out"),
      cells = x$cells,
      annot_path = x$annot
    ),
    class = "isoformic_invalid_assay_transcript_ids"
  )
  iso <- suppressMessages(read_oarfish_sc(
    fs::path(copy, "out"),
    cells = x$cells,
    annot_path = x$annot,
    drop_unannotated = TRUE
  ))
  testthat::expect_identical(nrow(iso), 71L)
  testthat::expect_true(all(is.finite(assay(iso, "abundance"))))
})

testthat::test_that("the shared readers accept every form of Oarfish path", {
  x <- oarfish_example()
  by_prefix <- read_inferential_replicates(x$prefixes[1:2])
  quant_files <- stats::setNames(
    paste0(x$prefixes[1:2], ".quant"),
    names(x$prefixes)[1:2]
  )
  testthat::expect_equal(read_inferential_replicates(quant_files), by_prefix)
  # A directory with one .quant file per sample
  root <- withr::local_tempdir()
  dirs <- vapply(
    names(x$prefixes)[1:2],
    function(s) {
      d <- fs::path(root, s)
      fs::dir_create(d)
      for (ext in c(
        ".quant",
        ".infreps.pq",
        ".ambig_info.tsv",
        ".meta_info.json"
      )) {
        fs::file_copy(paste0(x$prefixes[[s]], ext), fs::path(d, paste0(s, ext)))
      }
      as.character(d)
    },
    character(1L)
  )
  testthat::expect_equal(read_inferential_replicates(dirs), by_prefix)
  testthat::expect_equal(
    equivalence_support(dirs),
    equivalence_support(x$prefixes[1:2])
  )
})

testthat::test_that("Oarfish readers report inconsistent or incomplete input", {
  x <- oarfish_example()
  root <- withr::local_tempdir()
  # A sample quantified against fewer transcripts
  short <- fs::path(root, "short")
  quant <- readLines(paste0(x$prefixes[["control_2"]], ".quant"))
  writeLines(quant[-2L], paste0(short, ".quant"))
  testthat::expect_error(
    read_oarfish(
      c(a = x$prefixes[["control_1"]], b = short),
      annot_path = x$annot
    ),
    class = "isoformic_inconsistent_transcripts"
  )
  testthat::expect_error(
    read_oarfish(
      c(a = x$prefixes[["control_1"]], a = x$prefixes[["control_2"]]),
      annot_path = x$annot
    ),
    class = "isoformic_duplicated_sample_ids"
  )
  # Replicates with another number of rows than the .quant file
  fs::file_copy(
    paste0(x$prefixes[["control_2"]], ".infreps.pq"),
    paste0(short, ".infreps.pq")
  )
  testthat::expect_error(
    read_inferential_replicates(short),
    class = "isoformic_invalid_inferential_replicates"
  )
  # A folder with two .quant files is ambiguous
  two <- fs::path(root, "two")
  fs::dir_create(two)
  fs::file_copy(
    paste0(x$prefixes[["control_1"]], ".quant"),
    fs::path(two, "a.quant")
  )
  fs::file_copy(
    paste0(x$prefixes[["control_2"]], ".quant"),
    fs::path(two, "b.quant")
  )
  testthat::expect_error(
    read_oarfish(two, annot_path = x$annot),
    class = "isoformic_no_oarfish_quant"
  )
  # Unannotated transcripts can be dropped
  odd <- fs::path(root, "odd")
  quant[[2L]] <- sub("^[^\t]*", "NOT_IN_ANNOTATION", quant[[2L]])
  writeLines(quant, paste0(odd, ".quant"))
  iso <- suppressMessages(read_oarfish(
    odd,
    annot_path = x$annot,
    drop_unannotated = TRUE
  ))
  testthat::expect_identical(nrow(iso), 71L)
  # Run flags
  qc <- salmon_run_qc(x$prefixes, min_percent_mapped = 101)
  testthat::expect_identical(qc$qc_flag, rep("low_mapping_rate", 3L))
})

testthat::test_that("read_matrix_market rejects files without a size line", {
  bad <- withr::local_tempfile(fileext = ".mtx")
  writeLines(
    c("%%MatrixMarket matrix coordinate real general", "% only comments"),
    bad
  )
  testthat::expect_error(
    read_matrix_market(bad),
    class = "isoformic_invalid_oarfish_sc"
  )
  writeLines(
    c("%%MatrixMarket matrix coordinate real general", "not a size line"),
    bad
  )
  testthat::expect_error(
    read_matrix_market(bad),
    class = "isoformic_invalid_oarfish_sc"
  )
})
