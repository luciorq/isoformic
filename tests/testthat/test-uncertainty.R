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

testthat::test_that("replicate_moments sums the first n replicates exactly", {
  withr::local_seed(3L)
  counts <- matrix(
    stats::rgamma(7L * 9L, shape = 2, rate = 0.01),
    nrow = 7L,
    dimnames = list(paste0("t", 1:7), NULL)
  )
  # The reference: rowSums() of a copy of the first n columns
  reference <- function(counts, n) {
    kept <- counts[, seq_len(n), drop = FALSE]
    s1 <- rowSums(kept)
    s2 <- rowSums(kept * kept)
    mean_draw <- s1 / n
    var_draw <- if (n > 1L) {
      pmax(s2 - n * mean_draw^2, 0) / (n - 1L)
    } else {
      mean_draw * 0
    }
    list(mean = unname(mean_draw), var = unname(var_draw))
  }
  for (n in c(1L, 4L, 9L)) {
    testthat::expect_identical(
      replicate_moments(counts, n),
      reference(counts, n)
    )
  }
  testthat::expect_identical(
    replicate_moments(counts[0L, , drop = FALSE], 4L),
    list(mean = numeric(), var = numeric())
  )
})

testthat::test_that("read_inferential_replicates matches transcripts listed in another order", {
  fx <- salmon_uncertainty_fixture()
  same <- read_inferential_replicates(fx$dirs, transcripts = c("tB1", "tA3"))
  # Rewrite s2 with its transcripts in another order
  perm <- c(8L, 3L, 1L, 6L, 2L, 7L, 5L, 4L)
  tx <- dimnames(fx$draws)[[1L]][perm]
  boot <- fs::path(fx$dirs[["s2"]], "aux_info", "bootstrap")
  names_con <- gzfile(fs::path(boot, "names.tsv.gz"), "wb")
  writeBin(charToRaw(paste0(paste(tx, collapse = "\t"), "\n")), names_con)
  close(names_con)
  boot_con <- gzfile(fs::path(boot, "bootstraps.gz"), "wb")
  writeBin(
    as.vector(fx$draws[perm, "s2", ]),
    boot_con,
    size = 8L,
    endian = "little"
  )
  close(boot_con)
  shuffled <- read_inferential_replicates(
    fx$dirs,
    transcripts = c("tB1", "tA3")
  )
  testthat::expect_identical(shuffled, same)
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

# The line-by-line parser that the readers replaced
eq_classes_by_line <- function(file) {
  lines <- readLines(file)
  n_tx <- as.integer(lines[[1L]])
  n_classes <- as.integer(lines[[2L]])
  tx <- lines[2L + seq_len(n_tx)]
  parts <- strsplit(lines[2L + n_tx + seq_len(n_classes)], "\t", fixed = TRUE)
  size <- vapply(parts, function(p) as.integer(p[[1L]]), integer(1L))
  count <- vapply(parts, function(p) as.numeric(p[[length(p)]]), numeric(1L))
  tid <- unlist(lapply(parts, function(p) {
    as.integer(p[1L + seq_len(as.integer(p[[1L]]))]) + 1L
  }))
  cls <- rep.int(seq_along(size), size)
  list(
    transcripts = tx,
    members = data.frame(
      tid = as.integer(tid),
      size = size[cls],
      count = count[cls],
      class = cls
    )
  )
}

support_by_line <- function(dirs) {
  ec <- lapply(dirs, function(d) {
    eq_classes_by_line(eq_classes_file(d))
  })
  ids <- sub("\\|.*$", "", ec[[1L]]$transcripts)
  total <- unique_total <- n_classes <- numeric(length(ids))
  for (e in ec) {
    m <- e$members
    tid <- match(sub("\\|.*$", "", e$transcripts), ids)[m$tid]
    for (i in seq_along(tid)) {
      total[tid[i]] <- total[tid[i]] + m$count[i]
      unique_total[tid[i]] <- unique_total[tid[i]] +
        m$count[i] * (m$size[i] == 1L)
      n_classes[tid[i]] <- n_classes[tid[i]] + 1
    }
  }
  tibble::tibble(
    transcript_id = ids,
    n_classes = n_classes / length(ec),
    n_fragments = total,
    n_unique_fragments = unique_total,
    unique_fraction = ifelse(total > 0, unique_total / total, NA_real_),
    n_samples = length(ec)
  )
}

testthat::test_that("the equivalence class readers match a line-by-line parser", {
  fx <- salmon_uncertainty_fixture()
  tmp <- withr::local_tempdir()
  tx <- c("tx1|g1", "tx2|g1", "tx3|g2", "tx4|g2", "tx5|g3")
  write_eq <- function(name, classes, gz = FALSE) {
    dir <- fs::path(tmp, name)
    fs::dir_create(fs::path(dir, "aux_info"))
    file <- fs::path(dir, "aux_info", "eq_classes.txt")
    if (gz) {
      file <- paste0(file, ".gz")
    }
    con <- if (gz) gzfile(file, "w") else file(file, "w")
    writeLines(c(length(tx), length(classes), tx, classes), con)
    close(con)
    as.character(dir)
  }
  dirs <- c(
    # Plain, with weights (--dumpEqWeights), with weights compressed
    plain = write_eq("plain", c("2\t0\t1\t10", "1\t2\t5", "3\t4\t2\t0\t7")),
    weights = write_eq(
      "weights",
      c("2\t0\t1\t0.4\t0.6\t10", "1\t2\t1\t5", "3\t4\t2\t0\t0.2\t0.3\t0.5\t7")
    ),
    weights_gz = write_eq(
      "weights_gz",
      c("1\t1\t1\t2", "2\t3\t4\t0.5\t0.5\t1e3"),
      gz = TRUE
    ),
    empty = write_eq("empty", character()),
    # Hand-written or post-processed files may end lines with a tab
    trailing_tab = write_eq(
      "trailing_tab",
      c("2\t0\t1\t10\t", "1\t2\t5\t", "2\t3\t4\t0.5\t0.5\t7")
    )
  )
  # The fixture mixes lines with and without weights
  for (dir in c(fx$dirs, dirs)) {
    expected <- eq_classes_by_line(eq_classes_file(dir))
    testthat::expect_identical(read_eq_classes(dir), expected)
    members <- read_eq_class_members(dir)
    testthat::expect_identical(
      members$members,
      expected$members[, c("tid", "class")]
    )
  }
  testthat::expect_null(read_eq_classes(fs::path(tmp, "none")))
  for (set in list(fx$dirs, dirs[1:3], dirs)) {
    testthat::expect_equal(
      equivalence_support(set),
      support_by_line(set)
    )
  }
  # A class line cut short is an error, not a shifted parse
  bad <- write_eq("bad", c("2\t0\t1\t10", "3\t1"))
  testthat::expect_error(
    read_eq_class_members(bad),
    class = "isoformic_invalid_equivalence_classes"
  )
  testthat::expect_error(
    equivalence_support(bad),
    class = "isoformic_invalid_equivalence_classes"
  )
  not_numeric <- write_eq("not_numeric", c("2\t0\t1\tx"))
  testthat::expect_error(
    equivalence_support(not_numeric),
    class = "isoformic_invalid_equivalence_classes"
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
  # No replicates is flagged, alone or with another flag
  meta_file <- fs::path(fx$dirs[[1L]], "aux_info", "meta_info.json")
  meta <- readLines(meta_file)
  meta <- sub('"num_bootstraps": [0-9]+', '"num_bootstraps": 0', meta)
  writeLines(meta, meta_file)
  qc0 <- salmon_run_qc(fx$dirs)
  testthat::expect_equal(qc0$num_replicates[1L], 0)
  testthat::expect_identical(qc0$qc_flag[1L], "no_replicates")
  testthat::expect_true(is.na(qc0$qc_flag[2L]))
  meta_file <- fs::path(fx$dirs[[3L]], "aux_info", "meta_info.json")
  meta <- readLines(meta_file)
  meta <- sub('"num_bootstraps": [0-9]+', '"num_bootstraps": 0', meta)
  writeLines(meta, meta_file)
  testthat::expect_identical(
    salmon_run_qc(fx$dirs)$qc_flag[3L],
    "low_mapping_rate;no_replicates"
  )
  # The bundled mini quantifications map to 15 genes only, with 20 draws
  quants <- list.files(
    fs::path_package("isoformic", "extdata", "mini_quants"),
    full.names = TRUE
  )
  mini <- salmon_run_qc(quants[1L])
  testthat::expect_equal(mini$num_replicates, 20)
  testthat::expect_identical(mini$qc_flag, "low_mapping_rate")
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

testthat::test_that("add_uncertainty keeps the cached features and productivity", {
  x <- uncertainty_experiment()
  calls <- new.env(parent = emptyenv())
  calls$features <- 0L
  calls$productivity <- 0L
  real_features <- compute_transcript_features
  real_classify <- classify_productivity
  testthat::local_mocked_bindings(
    compute_transcript_features = function(...) {
      calls$features <- calls$features + 1L
      real_features(...)
    },
    classify_productivity = function(...) {
      calls$productivity <- calls$productivity + 1L
      real_classify(...)
    }
  )
  feats <- transcript_features(x$iso)
  prod <- transcript_productivity(x$iso)
  testthat::expect_identical(
    c(features = calls$features, productivity = calls$productivity),
    c(features = 1L, productivity = 1L)
  )
  iso <- suppressMessages(
    add_uncertainty(x$iso, paths = x$fixture$dirs, group = "condition")
  )
  # The cached tables are reused, by the new object and by the original one
  testthat::expect_identical(transcript_features(iso), feats)
  testthat::expect_identical(transcript_productivity(iso), prod)
  testthat::expect_identical(transcript_productivity(x$iso), prod)
  testthat::expect_identical(
    c(features = calls$features, productivity = calls$productivity),
    c(features = 1L, productivity = 1L)
  )
  # And they are those computed from scratch on the new object
  clear_cache(iso)
  testthat::expect_identical(transcript_features(iso), feats)
  testthat::expect_identical(transcript_productivity(iso), prod)
  testthat::expect_identical(
    c(features = calls$features, productivity = calls$productivity),
    c(features = 2L, productivity = 2L)
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
    isoform_switches(x$iso, drop_unstable = TRUE),
    class = "isoformic_no_uncertainty"
  )
  testthat::expect_error(
    isoform_switches(x$iso, confident_only = TRUE),
    class = "isoformic_defunct_argument"
  )
  testthat::expect_error(
    isoform_switches(x$iso, drop_unstable = NA),
    class = "isoformic_invalid_argument"
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
  testthat::expect_false(any(
    c("confidence", "confident_de", "reason") %in% colnames(conf)
  ))
  testthat::expect_identical(
    colnames(conf),
    c(
      "feature_id",
      "feature_type",
      "gene_id",
      "gene_name",
      "transcript_ids",
      "transcript_name",
      "transcript_type",
      "mean_count",
      "lfc_median",
      "lfc_lo",
      "lfc_hi",
      "lfc_sd",
      "sign_agreement",
      "unstable",
      "frac_draws_effect",
      "frac_draws_sig",
      "mean_inf_rv",
      "ambiguous",
      "flags",
      "det_log2FC",
      "det_qvalue",
      "n_det_sig",
      "unstable_de"
    )
  )
  testthat::expect_identical(conf$unstable, conf$sign_agreement < 0.95)
  # Well determined, strongly changed isoform: stable direction, no flag
  testthat::expect_false(get("tA3", "unstable"))
  testthat::expect_gt(get("tA3", "lfc_lo"), 1)
  testthat::expect_equal(get("tA3", "sign_agreement"), 1)
  testthat::expect_gt(get("tA3", "frac_draws_sig"), 0.5)
  testthat::expect_true(is.na(get("tA3", "flags")))
  testthat::expect_false(get("tA3", "unstable_de"))
  testthat::expect_equal(get("tA3", "det_log2FC"), 1.8)
  testthat::expect_equal(get("tA3", "n_det_sig"), 1L)
  # Ambiguous and uncertain isoform: the descriptive flags are listed, but
  # they do not make its stable call unstable
  testthat::expect_false(get("tA1", "unstable"))
  testthat::expect_match(get("tA1", "flags"), "ambiguous")
  testthat::expect_match(get("tA1", "flags"), "uncertain quantification")
  testthat::expect_false(get("tA1", "unstable_de"))
  testthat::expect_gt(get("tA1", "lfc_median"), 0)
  testthat::expect_lt(get("tA2", "lfc_median"), 0)
  # Overdispersed isoform without a real change: unstable direction
  testthat::expect_lt(get("tA4", "sign_agreement"), 0.95)
  testthat::expect_true(get("tA4", "unstable"))
  testthat::expect_match(get("tA4", "flags"), "^unstable direction")
  testthat::expect_false(get("tA4", "unstable_de"))
  # Not significant in the primary test
  testthat::expect_equal(get("tA6", "n_det_sig"), 0L)
  testthat::expect_false(get("tA6", "unstable_de"))
  # The threshold is the one given
  strict <- de_confidence(iso, group = "condition", min_agreement = 1)
  testthat::expect_identical(strict$unstable, strict$sign_agreement < 1)
  # Below min_count the direction is not assessed: unstable is NA and the
  # flags do not mention it; such a call is never counted as unstable
  high_min <- add_uncertainty(iso, group = "condition", min_count = 40)
  na_conf <- de_confidence(high_min)
  low <- na_conf$mean_count < 40
  testthat::expect_true(any(low))
  testthat::expect_true(all(is.na(na_conf$unstable[low])))
  testthat::expect_false(any(
    grepl("unstable direction", na_conf$flags[low], fixed = TRUE)
  ))
  testthat::expect_true(all(!is.na(na_conf$unstable[!low])))
  testthat::expect_false(any(na_conf$unstable_de[low]))
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

testthat::test_that("uncertainty tables and draws survive write_isoformic", {
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
  testthat::expect_identical(back@metadata$inf_reps, iso@metadata$inf_reps)
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
  testthat::expect_identical(
    de_confidence(back)$unstable,
    de_confidence(iso)$unstable
  )
  testthat::expect_equal(back@metadata$uncertainty_n_draws, "30")
  testthat::expect_equal(uncertainty_settings(back)$inf_rv_cutoff, 0.5)
  # Without the draws only the tables are written
  no_draws <- add_uncertainty(
    x$iso,
    paths = x$fixture$dirs,
    keep_draws = FALSE
  )
  write_isoformic(no_draws, fs::path(out, "tables"))
  testthat::expect_false(fs::dir_exists(fs::path(out, "tables", "inf_reps")))
  testthat::expect_null(
    read_isoformic(fs::path(out, "tables"))@metadata$inf_reps
  )
})

testthat::test_that("fold-change tables of earlier development versions are upgraded", {
  x <- uncertainty_experiment()
  iso <- add_uncertainty(x$iso, paths = x$fixture$dirs, group = "condition")
  current <- de_confidence(iso)
  # The layout written before the confidence class was removed
  legacy <- current[,
    setdiff(colnames(current), c("unstable", "unstable_de", "flags"))
  ]
  legacy$confidence <- factor(
    ifelse(current$unstable, "low", "high"),
    levels = c("high", "moderate", "low"),
    ordered = TRUE
  )
  legacy$reason <- sub(
    "unstable direction",
    "unstable fold-change",
    current$flags,
    fixed = TRUE
  )
  legacy$confident_de <- current$n_det_sig > 0L & !current$unstable
  old <- iso
  dea <- old@dea
  dea[["confidence"]] <- legacy
  old@dea <- dea
  # In memory, silently
  testthat::expect_no_message(upgraded <- de_confidence(old))
  testthat::expect_identical(upgraded$unstable, upgraded$sign_agreement < 0.95)
  testthat::expect_false(any(
    c("confidence", "confident_de", "reason") %in% colnames(upgraded)
  ))
  testthat::expect_equal(upgraded, current)
  sw <- isoform_switches(old)
  testthat::expect_true("n_det_unstable" %in% colnames(sw))
  # From disk, with one message; the class may be stored as text
  out <- withr::local_tempdir()
  for (as_text in c(FALSE, TRUE)) {
    if (as_text) {
      dea[["confidence"]]$confidence <- as.character(legacy$confidence)
      old@dea <- dea
    }
    dir <- fs::path(out, paste0("legacy_", as_text))
    write_isoformic(old, dir)
    testthat::expect_message(
      back <- read_isoformic(dir),
      class = "isoformic_legacy_confidence"
    )
    testthat::expect_false(any(
      c("confidence", "confident_de", "reason") %in%
        colnames(back@dea[["confidence"]])
    ))
    testthat::expect_equal(
      as.data.frame(de_confidence(back)),
      as.data.frame(current),
      ignore_attr = TRUE
    )
  }
})

testthat::test_that("legacy tables lose the uncertain flag of features that were not assessed", {
  x <- uncertainty_experiment()
  iso <- add_uncertainty(x$iso, paths = x$fixture$dirs, group = "condition")
  # From the stored draws, so that the features below `min_count` are kept
  iso <- add_uncertainty(iso, group = "condition", min_count = 45)
  current <- de_confidence(iso)
  low <- is.na(current$sign_agreement)
  testthat::expect_true(any(low))
  testthat::expect_true(any(!low))
  # Older versions kept the InfRV, and could set the uncertain flag, below
  # `min_count`
  legacy <- current[,
    setdiff(colnames(current), c("unstable", "unstable_de", "flags"))
  ]
  legacy$mean_inf_rv[low] <- 0.9
  legacy$reason <- sub(
    "unstable direction",
    "unstable fold-change",
    current$flags,
    fixed = TRUE
  )
  legacy$reason[low] <- "uncertain quantification"
  legacy$confidence <- factor(
    ifelse(low | current$unstable %in% TRUE, "low", "high"),
    levels = c("high", "moderate", "low"),
    ordered = TRUE
  )
  legacy$confident_de <- current$n_det_sig > 0L & !current$unstable %in% TRUE
  old <- iso
  dea <- old@dea
  dea[["confidence"]] <- legacy
  old@dea <- dea
  stored <- uncertainty_settings(iso)$min_agreement
  upgraded <- de_confidence(old)
  testthat::expect_identical(
    upgraded$flags,
    de_confidence(old, min_agreement = stored)$flags
  )
  testthat::expect_true(all(is.na(upgraded$flags[low])))
  testthat::expect_true(all(is.na(upgraded$mean_inf_rv[low])))
  testthat::expect_equal(upgraded, current)
  testthat::expect_equal(uncertainty_report(old)$de, uncertainty_report(iso)$de)
  # The same from disk
  out <- withr::local_tempdir()
  dir <- fs::path(out, "legacy")
  write_isoformic(old, dir)
  testthat::expect_message(
    back <- read_isoformic(dir),
    class = "isoformic_legacy_confidence"
  )
  testthat::expect_equal(
    as.data.frame(de_confidence(back)),
    as.data.frame(current),
    ignore_attr = TRUE
  )
})

testthat::test_that("draws read for an object without counts cover only its rows", {
  x <- uncertainty_experiment()
  sub <- x$iso[c("tA1", "tA2", "tB1"), ]
  tx_assay(sub, "counts") <- NULL
  sub <- add_uncertainty(sub, paths = x$fixture$dirs)
  testthat::expect_setequal(
    rownames(sub@metadata$inf_reps[[1L]]),
    c("tA1", "tA2", "tB1")
  )
  unc <- tx_uncertainty(sub)
  testthat::expect_setequal(unc$transcript_id, c("tA1", "tA2", "tB1"))
  # The tables recomputed from the stored draws match those from the paths
  again <- add_uncertainty(sub)
  testthat::expect_equal(
    tx_uncertainty(again)$mean_inf_rv,
    unc$mean_inf_rv
  )
  # The draws written are the draws in memory
  out <- withr::local_tempdir()
  write_isoformic(sub, fs::path(out, "exp"))
  back <- read_isoformic(fs::path(out, "exp"))
  testthat::expect_identical(back@metadata$inf_reps, sub@metadata$inf_reps)
  testthat::expect_equal(
    tx_uncertainty(back)$mean_inf_rv,
    unc$mean_inf_rv
  )
})

testthat::test_that("isoform_switches and dominant_isoform use the stability", {
  x <- uncertainty_experiment()
  iso <- add_uncertainty(x$iso, paths = x$fixture$dirs, group = "condition")
  sw <- isoform_switches(iso)
  ga <- sw[sw$gene_id == "gA", ]
  testthat::expect_identical(ga$switch_type, "opposite_direction")
  testthat::expect_equal(ga$n_det_sig, 3L)
  testthat::expect_equal(ga$n_det_unstable, 0L)
  testthat::expect_false("n_det_confident" %in% colnames(sw))
  # A significant call with an unstable direction (tA4) is counted, and
  # drop_unstable = TRUE leaves out exactly that call
  de_tx(iso) <- data.frame(
    transcript_id = c("tA1", "tA2", "tA3", "tA4", "tB1"),
    log2FC = c(1.2, -1.2, 1.8, 1.5, 0.1),
    pvalue = c(0.001, 0.001, 0.0001, 0.001, 0.6),
    qvalue = c(0.01, 0.01, 0.001, 0.01, 0.7)
  )
  sw <- isoform_switches(iso)
  ga <- sw[sw$gene_id == "gA", ]
  testthat::expect_equal(ga$n_det_sig, 4L)
  testthat::expect_equal(ga$n_det_unstable, 1L)
  # The stored table follows the new results as well
  conf <- de_confidence(iso)
  testthat::expect_equal(conf$n_det_sig[conf$feature_id == "tA4"], 1L)
  testthat::expect_true(conf$unstable_de[conf$feature_id == "tA4"])
  testthat::expect_equal(conf$det_log2FC[conf$feature_id == "tA4"], 1.5)
  de <- uncertainty_report(iso)$de
  testthat::expect_equal(de$n_significant[de$flag == "unstable"], 4L)
  testthat::expect_equal(de$n_flagged[de$flag == "unstable"], 1L)
  testthat::expect_identical(
    colnames(sw)[match("n_det_sig", colnames(sw)) + 1L],
    "n_det_unstable"
  )
  kept <- isoform_switches(iso, drop_unstable = TRUE)
  ga_kept <- kept[kept$gene_id == "gA", ]
  testthat::expect_equal(ga_kept$n_det_sig, ga$n_det_sig - ga$n_det_unstable)
  testthat::expect_identical(ga_kept$switch_type, "opposite_direction")
  testthat::expect_error(
    isoform_switches(iso, confident_only = FALSE),
    class = "isoformic_defunct_argument"
  )
  # An old positional call does not silently mean drop_unstable
  testthat::expect_error(
    isoform_switches(iso, NULL, TRUE),
    class = "isoformic_defunct_argument"
  )
  dom <- dominant_isoform(iso, group = "condition")
  testthat::expect_true(all(
    c("dominant_inf_rv", "dominant_ambiguous") %in% colnames(dom)
  ))
  # Without tables the columns are absent
  testthat::expect_false(
    "dominant_inf_rv" %in% colnames(dominant_isoform(x$iso))
  )
  testthat::expect_false(any(
    c("n_det_unstable", "n_det_confident") %in%
      colnames(isoform_switches(x$iso))
  ))
})

testthat::test_that("members of a collapsed group take the direction of the group", {
  x <- uncertainty_experiment()
  iso <- add_uncertainty(
    x$iso,
    paths = x$fixture$dirs,
    group = "condition",
    collapse_groups = TRUE
  )
  # gA_amb1 sums tA1 and tA2, both significant; tA3 is significant alone
  sw <- isoform_switches(iso)
  testthat::expect_equal(sw$n_det_sig[sw$gene_id == "gA"], 3L)
  testthat::expect_equal(sw$n_det_unstable[sw$gene_id == "gA"], 0L)
  dea <- iso@dea
  g <- dea$confidence$feature_type == "group"
  testthat::expect_identical(dea$confidence$feature_id[g], "gA_amb1")
  dea$confidence$sign_agreement[g] <- 0.5
  dea$confidence$unstable[g] <- TRUE
  iso@dea <- dea
  conf <- de_confidence(iso)
  testthat::expect_true(conf$unstable_de[conf$feature_id == "gA_amb1"])
  sw <- isoform_switches(iso)
  testthat::expect_equal(sw$n_det_sig[sw$gene_id == "gA"], 3L)
  testthat::expect_equal(sw$n_det_unstable[sw$gene_id == "gA"], 2L)
  kept <- isoform_switches(iso, drop_unstable = TRUE)
  testthat::expect_equal(kept$n_det_sig[kept$gene_id == "gA"], 1L)
  # The report counts the group once, as one significant feature
  de <- uncertainty_report(iso)$de
  testthat::expect_equal(de$n_significant, rep(2L, 3L))
  testthat::expect_equal(de$n_flagged[de$flag == "unstable"], 1L)
  testthat::expect_identical(
    uncertainty_report(iso)$switches$n_opposite_direction,
    c(1L, 0L)
  )
})

testthat::test_that("uncertainty_report counts switches without a gene-level table", {
  x <- uncertainty_experiment()
  iso <- add_uncertainty(x$iso, paths = x$fixture$dirs, group = "condition")
  with_genes <- uncertainty_report(iso)$switches
  dea <- iso@dea
  dea[["deg"]] <- NULL
  iso@dea <- dea
  sw <- uncertainty_report(iso)$switches
  testthat::expect_s3_class(sw, "tbl_df")
  testthat::expect_identical(nrow(sw), 2L)
  testthat::expect_identical(sw$n_genes, with_genes$n_genes)
  testthat::expect_identical(
    sw$n_opposite_direction,
    with_genes$n_opposite_direction
  )
  testthat::expect_identical(sw$n_det_not_deg, c(NA_integer_, NA_integer_))
})

testthat::test_that("de_confidence applies new thresholds to the stored table", {
  x <- uncertainty_experiment()
  iso <- add_uncertainty(x$iso, paths = x$fixture$dirs, group = "condition")
  stored <- de_confidence(iso)
  testthat::expect_identical(
    de_confidence(iso, min_agreement = 0.95, inf_rv_cutoff = 0.5),
    stored
  )
  strict <- de_confidence(iso, min_agreement = 1)
  testthat::expect_identical(strict$unstable, strict$sign_agreement < 1)
  testthat::expect_identical(
    strict$unstable_de,
    strict$n_det_sig > 0L & strict$sign_agreement < 1
  )
  testthat::expect_true(any(strict$unstable & !stored$unstable))
  # tA6 agrees in 29 of 30 replicates
  testthat::expect_match(
    strict$flags[strict$feature_id == "tA6"],
    "unstable direction",
    fixed = TRUE
  )
  loose <- de_confidence(iso, inf_rv_cutoff = Inf)
  testthat::expect_false(any(grepl(
    "uncertain quantification",
    loose$flags,
    fixed = TRUE
  )))
  testthat::expect_identical(loose$unstable, stored$unstable)
})

testthat::test_that("invalid thresholds are rejected", {
  x <- uncertainty_experiment()
  iso <- add_uncertainty(x$iso, paths = x$fixture$dirs, group = "condition")
  bad_agreement <- list("a", c(0.5, 0.9), NA_real_, 1.5, -0.1)
  for (bad in bad_agreement) {
    testthat::expect_error(
      de_confidence(iso, min_agreement = bad),
      class = "isoformic_invalid_argument"
    )
    testthat::expect_error(
      de_confidence(iso, group = "condition", min_agreement = bad),
      class = "isoformic_invalid_argument"
    )
    testthat::expect_error(
      add_uncertainty(iso, min_agreement = bad),
      class = "isoformic_invalid_argument"
    )
  }
  for (bad in list("z", -1, NA_real_, c(1, 2))) {
    testthat::expect_error(
      de_confidence(iso, inf_rv_cutoff = bad),
      class = "isoformic_invalid_argument"
    )
    testthat::expect_error(
      add_uncertainty(iso, inf_rv_cutoff = bad),
      class = "isoformic_invalid_argument"
    )
  }
  # A valid value marks exactly the features below it
  res <- de_confidence(iso, min_agreement = 0.85)
  testthat::expect_identical(res$unstable, res$sign_agreement < 0.85)
  testthat::expect_identical(
    de_confidence(iso, min_agreement = 1L)$unstable,
    res$sign_agreement < 1
  )
})

testthat::test_that("features below min_count are not flagged as uncertain", {
  x <- uncertainty_experiment()
  iso <- add_uncertainty(x$iso, paths = x$fixture$dirs, group = "condition")
  high <- add_uncertainty(iso, group = "condition", min_count = 300)
  unc <- tx_uncertainty(high)
  testthat::expect_false(any(unc$assessed))
  conf <- de_confidence(high)
  testthat::expect_true(all(is.na(conf$sign_agreement)))
  testthat::expect_true(all(is.na(conf$mean_inf_rv)))
  no_uncertain <- function(flags) {
    !any(grepl("uncertain quantification", flags, fixed = TRUE))
  }
  testthat::expect_true(no_uncertain(conf$flags))
  testthat::expect_true(no_uncertain(
    de_confidence(high, inf_rv_cutoff = 0.1)$flags
  ))
  testthat::expect_true(no_uncertain(
    de_confidence(high, group = "condition")$flags
  ))
  # A table stored before this change keeps the InfRV of such features
  old <- de_confidence(iso)
  old$sign_agreement[] <- NA_real_
  testthat::expect_true(no_uncertain(
    rethreshold_confidence(
      old,
      list(min_agreement = 0.95, inf_rv_cutoff = 0)
    )$flags
  ))
})

testthat::test_that("uncertainty_report counts switches without the stability", {
  x <- uncertainty_experiment()
  i0 <- add_uncertainty(x$iso, paths = x$fixture$dirs)
  sw <- uncertainty_report(i0)$switches
  testthat::expect_s3_class(sw, "tbl_df")
  testthat::expect_identical(nrow(sw), 2L)
  all_sw <- isoform_switches(i0)
  testthat::expect_identical(sw$n_genes[[1L]], nrow(all_sw))
  testthat::expect_identical(
    sw$n_opposite_direction[[1L]],
    sum(all_sw$switch_type == "opposite_direction")
  )
  testthat::expect_identical(
    sw$n_det_not_deg[[1L]],
    sum(all_sw$switch_type == "det_not_deg")
  )
  testthat::expect_true(is.na(sw$n_genes[[2L]]))
  testthat::expect_true(is.na(sw$n_opposite_direction[[2L]]))
  testthat::expect_true(is.na(sw$n_det_not_deg[[2L]]))
  # With the stability both rows are filled
  ig <- add_uncertainty(i0, group = "condition")
  testthat::expect_false(anyNA(uncertainty_report(ig)$switches$n_genes))
  # Without a transcript-level table there is no switch summary
  dea <- i0@dea
  dea[["det"]] <- NULL
  i0@dea <- dea
  testthat::expect_null(uncertainty_report(i0)$switches)
})

testthat::test_that("a missing value in the replicates is reported", {
  x <- uncertainty_experiment()
  iso <- add_uncertainty(x$iso, paths = x$fixture$dirs, keep_draws = TRUE)
  meta <- iso@metadata
  meta$inf_reps[[1L]]["tA4", 2L] <- NA
  iso@metadata <- meta
  testthat::expect_error(
    de_confidence(iso, group = "condition"),
    class = "isoformic_missing_draws"
  )
  testthat::expect_error(
    add_uncertainty(iso, group = "condition"),
    class = "isoformic_missing_draws"
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
  testthat::expect_true(
    all(c("q25_inf_rv", "q75_inf_rv") %in% colnames(rep$transcripts))
  )
  testthat::expect_true(all(
    rep$transcripts$q25_inf_rv <= rep$transcripts$median_inf_rv &
      rep$transcripts$median_inf_rv <= rep$transcripts$q75_inf_rv
  ))
  de <- rep$de
  testthat::expect_identical(de$flag, c("unstable", "uncertain", "ambiguous"))
  testthat::expect_named(
    de,
    c("flag", "rule", "n_significant", "n_assessed", "n_flagged", "pct_flagged")
  )
  testthat::expect_identical(de$rule[1L], "sign_agreement < 0.95")
  testthat::expect_identical(de$rule[2L], "mean_inf_rv >= 0.5")
  testthat::expect_equal(de$n_significant, rep(3L, 3L))
  testthat::expect_equal(de$n_assessed, rep(3L, 3L))
  testthat::expect_equal(de$n_flagged, c(0L, 2L, 2L))
  testthat::expect_equal(de$pct_flagged, 100 * c(0, 2, 2) / 3)
  testthat::expect_identical(
    rep$switches$set,
    c("all significant isoforms", "without unstable directions")
  )
  testthat::expect_equal(rep$switches$n_genes, c(1L, 1L))
  txt <- uncertainty_methods(iso)
  testthat::expect_match(txt, "Salmon 2.8.0", fixed = TRUE)
  testthat::expect_match(txt, "30 Gibbs samples", fixed = TRUE)
  testthat::expect_match(txt, "174 degrees of freedom", fixed = TRUE)
  testthat::expect_match(txt, "Turro", fixed = TRUE)
  testthat::expect_match(txt, "b versus a", fixed = TRUE)
  testthat::expect_match(txt, "InfRV of at least 0.5", fixed = TRUE)
  testthat::expect_match(txt, "95%", fixed = TRUE)
  testthat::expect_match(txt, "sign agreement", fixed = TRUE)
  testthat::expect_match(txt, "unstable direction", fixed = TRUE)
  testthat::expect_no_match(txt, "high confidence", fixed = TRUE)
  testthat::expect_no_match(txt, "moderate", fixed = TRUE)
  testthat::expect_no_match(txt, "identified each", fixed = TRUE)
})

testthat::test_that("the uncertainty plots build", {
  x <- uncertainty_experiment()
  iso <- add_uncertainty(x$iso, paths = x$fixture$dirs, group = "condition")
  p1 <- plot_de_confidence(iso, gene_name = "GENEA")
  testthat::expect_s3_class(p1, "ggplot")
  testthat::expect_no_error(ggplot2::ggplot_build(p1))
  testthat::expect_identical(
    levels(p1$data$direction),
    c("stable", "unstable")
  )
  conf_a <- de_confidence(iso, gene_name = "GENEA")
  testthat::expect_identical(
    as.character(p1$data$direction[match(
      conf_a$feature_id,
      p1$data$feature_id
    )]),
    ifelse(conf_a$unstable, "unstable", "stable")
  )
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
  testthat::expect_false(
    exists("confidence_palette", envir = asNamespace("isoformic"))
  )
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
  testthat::expect_null(attr(back, "moments"))
  # The moments of the replicates go to a table of their own
  mom_tab <- inferential_moments_table(reps)
  testthat::expect_identical(
    colnames(mom_tab),
    c("transcript_id", "sample_id", "n_draws", "mean", "var")
  )
  with_mom <- as_inferential_replicates(tab, moments = mom_tab)
  testthat::expect_equal(attr(with_mom, "moments"), attr(reps, "moments"))
  # Row order of the table only changes the order of the dimnames
  rev_mom <- attr(
    as_inferential_replicates(
      tab,
      moments = mom_tab[rev(seq_len(nrow(mom_tab))), ]
    ),
    "moments"
  )
  testthat::expect_equal(
    rev_mom$var[rownames(reps[[1L]]), colnames(reps[[1L]])],
    attr(reps, "moments")$var
  )
  testthat::expect_error(
    as_inferential_replicates(tab, moments = mom_tab[-1L, ]),
    class = "isoformic_invalid_inferential_replicates"
  )
  testthat::expect_error(
    as_inferential_replicates(
      tab[tab$draw <= 6L, ],
      moments = mom_tab
    ),
    class = "isoformic_invalid_inferential_replicates"
  )
  testthat::expect_error(
    inferential_moments_table(back),
    class = "isoformic_no_replicate_moments"
  )
  attr(reps, "samp_type") <- NULL
  attr(back, "samp_type") <- NULL
  attr(reps, "moments") <- NULL
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

testthat::test_that("canonical and shuffled replicate tables give the same draws", {
  fx <- salmon_uncertainty_fixture()
  reps <- read_inferential_replicates(fx$dirs, n_draws = 5L)
  tab <- inferential_replicates_table(reps)
  ids <- unique(tab$transcript_id)
  samples <- unique(tab$sample_id)
  testthat::expect_identical(
    canonical_replicates_layout(tab),
    list(ids = ids, samples = samples, draws = 1:5)
  )
  # Two rows of the second draw swapped: same cells, general path
  swapped <- tab
  rows <- length(ids) * length(samples) + 1:2
  swapped[rows, ] <- tab[rev(rows), ]
  testthat::expect_null(canonical_replicates_layout(swapped))
  # Draws out of order, or a draw block of another size, are not canonical
  testthat::expect_null(
    canonical_replicates_layout(tab[order(-tab$draw, seq_len(nrow(tab))), ])
  )
  testthat::expect_null(canonical_replicates_layout(tab[-1L, ]))
  fast <- as_inferential_replicates(tab)
  general <- as_inferential_replicates(swapped)
  testthat::expect_identical(general, fast)
  # Integer counts and factor columns give the same double matrices
  int_tab <- tab
  int_tab$count <- as.integer(round(int_tab$count))
  int_tab$transcript_id <- factor(int_tab$transcript_id, levels = ids)
  as_int <- as_inferential_replicates(int_tab)
  testthat::expect_type(as_int[[1L]], "double")
  testthat::expect_identical(
    as_int[[2L]],
    round(unclass(fast)[[2L]])
  )
})

testthat::test_that("a repeated cell in a replicate or moments table is an error", {
  tab <- tibble::tibble(
    transcript_id = rep(c("t1", "t2"), 4L),
    sample_id = rep(rep(c("s1", "s2"), each = 2L), 2L),
    draw = rep(1:2, each = 4L),
    count = as.numeric(1:8)
  )
  ok <- as_inferential_replicates(tab)
  testthat::expect_identical(ok[[2L]]["t2", "s1"], 6)
  # Row 2 (t2, s1, draw 1) replaced by a copy of row 1: as many rows, one
  # cell twice and one missing
  dup <- tab
  dup[2L, ] <- tab[1L, ]
  testthat::expect_error(
    as_inferential_replicates(dup),
    class = "isoformic_invalid_inferential_replicates",
    regexp = "more than one row"
  )
  mom <- tibble::tibble(
    transcript_id = c("t1", "t2", "t1", "t2"),
    sample_id = c("s1", "s1", "s2", "s2"),
    n_draws = 2L,
    mean = 1:4,
    var = 1:4
  )
  testthat::expect_equal(moments_from_table(mom)$mean[["t2", "s2"]], 4)
  mom_dup <- mom
  mom_dup[2L, ] <- mom[1L, ]
  testthat::expect_error(
    moments_from_table(mom_dup),
    class = "isoformic_invalid_inferential_replicates",
    regexp = "more than one row"
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
  # edgeR 4.99 (Bioconductor devel) takes the directories as `sample.dirs`
  # (its first argument is now `parent.dir`) and returns a DGEList by
  # default; earlier versions take them as `paths` and return a list
  dirs <- unname(x$fixture$dirs)
  cs_args <- if ("sample.dirs" %in% names(formals(edgeR::catchSalmon))) {
    list(sample.dirs = dirs, DGEList = FALSE, divide = FALSE)
  } else {
    list(paths = dirs)
  }
  cs <- suppressMessages(
    do.call(edgeR::catchSalmon, c(cs_args, list(verbose = FALSE)))
  )
  testthat::expect_equal(
    unname(by_id[rownames(cs$counts)]),
    unname(cs$annotation$Overdispersion)
  )
})

testthat::test_that("the overdispersion does not depend on the transcripts whose replicates are kept", {
  x <- uncertainty_experiment()
  all_reps <- read_inferential_replicates(x$fixture$dirs)
  full <- x$iso
  full@metadata$inf_reps <- all_reps
  attr(full@metadata$inf_reps, "moments") <- NULL
  od_full <- inferential_overdispersion(full)

  # add_uncertainty() keeps the transcripts with a mean count of at least
  # min_count: tA5 (mean 0.2) has no stored replicates
  kept <- add_uncertainty(x$iso, paths = x$fixture$dirs, min_count = 1)
  testthat::expect_false("tA5" %in% rownames(kept@metadata$inf_reps[[1L]]))
  od_kept <- inferential_overdispersion(kept)
  testthat::expect_equal(attr(od_kept, "prior"), attr(od_full, "prior"))
  testthat::expect_equal(
    od_kept$overdispersion[match(od_full$transcript_id, od_kept$transcript_id)],
    od_full$overdispersion
  )
  # A higher min_count keeps fewer transcripts and changes nothing
  fewer <- add_uncertainty(x$iso, paths = x$fixture$dirs, min_count = 60)
  testthat::expect_lt(
    nrow(fewer@metadata$inf_reps[[1L]]),
    nrow(kept@metadata$inf_reps[[1L]])
  )
  testthat::expect_equal(
    inferential_overdispersion(fewer)$overdispersion,
    od_kept$overdispersion
  )
  scaled <- scale_counts(fewer)
  counts <- tx_assay(fewer, "counts")
  ratio <- counts / tx_assay(scaled, "counts_scaled")
  ratio <- ratio[counts[, 1L] > 0, 1L]
  by_id <- stats::setNames(od_full$overdispersion, od_full$transcript_id)
  testthat::expect_true("tA5" %in% names(ratio))
  testthat::expect_equal(unname(ratio), unname(by_id[names(ratio)]))

  # Without the moments the estimates cover the stored transcripts only
  testthat::expect_true("tA5" %in% od_kept$transcript_id)
  no_mom <- kept
  attr(no_mom@metadata$inf_reps, "moments") <- NULL
  testthat::expect_false(
    "tA5" %in% inferential_overdispersion(no_mom)$transcript_id
  )

  # Subsetting the samples keeps the moments of those samples only
  two <- suppressMessages(kept[, c("s1", "s4")])
  mom <- inf_reps_moments(two@metadata$inf_reps)
  testthat::expect_identical(colnames(mom$mean), c("s1", "s4"))
  testthat::expect_true("tA5" %in% rownames(mom$mean))

  # A subset made only of transcripts without stored replicates keeps the
  # moments, so their overdispersion is still available
  only_a5 <- suppressMessages(kept["tA5", ])
  testthat::expect_false(is.null(only_a5@metadata$inf_reps))
  testthat::expect_false(is.null(attr(only_a5@metadata$inf_reps, "moments")))
  testthat::expect_equal(
    tx_assay(scale_counts(only_a5), "counts_scaled")["tA5", ],
    tx_assay(scale_counts(kept), "counts_scaled")["tA5", ]
  )
  od_a5 <- inferential_overdispersion(only_a5)
  testthat::expect_identical(od_a5$transcript_id, "tA5")
  testthat::expect_equal(
    od_a5$overdispersion,
    od_kept$overdispersion[od_kept$transcript_id == "tA5"]
  )
  # Without the moments such a subset has no replicates left
  no_mom_a5 <- suppressMessages(no_mom["tA5", ])
  testthat::expect_null(no_mom_a5@metadata$inf_reps)
})

testthat::test_that("the moments use the replicates kept in every sample", {
  many <- salmon_uncertainty_fixture(n_draws = 30L)
  few <- salmon_uncertainty_fixture(n_draws = 10L)
  dirs <- c(many$dirs[c("s1", "s2")], few$dirs[c("s3", "s4")])
  reps <- read_inferential_replicates(dirs, transcripts = c("tA1", "tA3"))
  testthat::expect_length(reps, 10L)
  mom <- attr(reps, "moments")
  testthat::expect_identical(mom$n, 10L)
  testthat::expect_identical(rownames(mom$mean), rownames(many$draws))
  all_reps <- read_inferential_replicates(dirs)
  expected <- draw_moments(simplify2array(unname(lapply(all_reps, unclass))))
  testthat::expect_equal(unname(mom$mean), unname(expected$mean))
  testthat::expect_equal(unname(mom$var), unname(expected$var))

  # The same with the samples with fewer draws first
  dirs_rev <- c(few$dirs[c("s3", "s4")], many$dirs[c("s1", "s2")])
  reps_rev <- read_inferential_replicates(
    dirs_rev,
    transcripts = c("tA1", "tA3")
  )
  testthat::expect_length(reps_rev, 10L)
  mom_rev <- attr(reps_rev, "moments")
  testthat::expect_identical(mom_rev$n, 10L)
  all_rev <- read_inferential_replicates(dirs_rev)
  expected_rev <- draw_moments(simplify2array(unname(lapply(all_rev, unclass))))
  testthat::expect_equal(unname(mom_rev$mean), unname(expected_rev$mean))
  testthat::expect_equal(unname(mom_rev$var), unname(expected_rev$var))

  # n_draws: the moments use the same first draws as the replicates
  reps_7 <- read_inferential_replicates(
    dirs,
    transcripts = c("tA1", "tA3"),
    n_draws = 7L
  )
  testthat::expect_length(reps_7, 7L)
  mom_7 <- attr(reps_7, "moments")
  testthat::expect_identical(mom_7$n, 7L)
  all_7 <- read_inferential_replicates(dirs, n_draws = 7L)
  expected_7 <- draw_moments(simplify2array(unname(lapply(all_7, unclass))))
  testthat::expect_equal(unname(mom_7$mean), unname(expected_7$mean))
  testthat::expect_equal(unname(mom_7$var), unname(expected_7$var))
})

testthat::test_that("scale_counts divides the counts by the overdispersion", {
  x <- uncertainty_experiment()
  iso <- add_uncertainty(x$iso, paths = x$fixture$dirs)
  scaled <- scale_counts(iso)
  testthat::expect_true("counts_scaled" %in% assay_names(scaled))
  od <- inferential_overdispersion(iso)
  counts <- tx_assay(iso, "counts")
  ratio <- counts / tx_assay(scaled, "counts_scaled")
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

testthat::test_that("identifiability is reported and does not change the direction flag", {
  x <- uncertainty_experiment()
  plain <- add_uncertainty(x$iso, paths = x$fixture$dirs, group = "condition")
  testthat::expect_false("identifiable" %in% colnames(de_confidence(plain)))
  ident <- equivalence_identifiability(x$fixture$dirs)
  iso <- add_uncertainty(
    x$iso,
    paths = x$fixture$dirs,
    group = "condition",
    identifiability = TRUE
  )
  unc <- tx_uncertainty(iso)
  testthat::expect_true(
    all(c("identifiable_fraction", "identifiable") %in% colnames(unc))
  )
  testthat::expect_identical(
    unc$identifiable,
    ident$identifiable[match(unc$transcript_id, ident$transcript_id)]
  )
  # The table of equivalence_identifiability() gives the same result
  from_table <- add_uncertainty(
    x$iso,
    paths = x$fixture$dirs,
    group = "condition",
    identifiability = ident
  )
  testthat::expect_identical(de_confidence(from_table), de_confidence(iso))
  # A transcript marked non-identifiable is flagged as such; its direction
  # flag does not change
  ident$identifiable[ident$transcript_id == "tA3"] <- FALSE
  flagged <- add_uncertainty(
    plain,
    group = "condition",
    identifiability = ident
  )
  conf <- de_confidence(flagged)
  at <- conf$feature_id == "tA3"
  testthat::expect_false(conf$identifiable[at])
  testthat::expect_match(conf$flags[at], "not identifiable")
  plain_conf <- de_confidence(plain)
  testthat::expect_identical(
    conf$unstable,
    plain_conf$unstable[match(conf$feature_id, plain_conf$feature_id)]
  )
  testthat::expect_identical(
    conf$unstable_de,
    plain_conf$unstable_de[match(conf$feature_id, plain_conf$feature_id)]
  )
  testthat::expect_identical(
    colnames(conf)[match("ambiguous", colnames(conf)) + 1:2],
    c("identifiable", "flags")
  )
  de <- uncertainty_report(flagged)$de
  testthat::expect_identical(
    de$flag,
    c("unstable", "uncertain", "ambiguous", "not_identifiable")
  )
  testthat::expect_equal(de$n_flagged[de$flag == "not_identifiable"], 1L)
  testthat::expect_match(
    paste(uncertainty_methods(flagged), collapse = " "),
    "structurally identifiable from the equivalence classes of every sample"
  )
  testthat::expect_match(
    paste(uncertainty_methods(flagged), collapse = " "),
    "not the same as, the range-of-optima non-identifiability"
  )
  # The flag is kept when the tables are recomputed from the stored draws
  again <- add_uncertainty(flagged, group = "condition")
  testthat::expect_identical(
    de_confidence(again)$identifiable,
    conf$identifiable
  )
  # Collapsed groups are not assessed
  grouped <- de_confidence(flagged, group = "condition", collapse_groups = TRUE)
  testthat::expect_true(
    all(is.na(grouped$identifiable[grouped$feature_type == "group"]))
  )
  testthat::expect_error(
    add_uncertainty(plain, identifiability = TRUE),
    class = "isoformic_identifiability_needs_paths"
  )
  testthat::expect_error(
    add_uncertainty(plain, identifiability = data.frame(transcript_id = "tA1")),
    class = "isoformic_invalid_identifiability"
  )
})

testthat::test_that("row_stats matches median, quantile and sd row by row", {
  set.seed(11)
  for (d in c(1L, 2L, 7L, 30L)) {
    m <- matrix(round(stats::rnorm(40L * d), 1L), nrow = 40L)
    m[1L, ] <- 2.5
    m[2L, ] <- c(0, rep(1, d - 1L))
    m[3L, 1L] <- NA
    if (d > 1L) {
      m[4L, ] <- NA
      m[5L, 2L] <- NaN
    }
    m[6L, 1L] <- -Inf
    rownames(m) <- paste0("f", seq_len(nrow(m)))
    probs <- c(0.025, 0.5, 0.975)
    res <- row_stats(m, probs)
    testthat::expect_identical(res$median, unname(apply(m, 1L, stats::median)))
    testthat::expect_identical(res$quantiles, row_quantiles(m, probs))
    testthat::expect_equal(
      res$sd,
      unname(apply(m, 1L, stats::sd)),
      tolerance = 1e-12
    )
  }
  empty <- row_stats(matrix(numeric(), 0L, 5L), 0.5)
  testthat::expect_identical(empty$median, numeric())
  testthat::expect_identical(dim(empty$quantiles), c(1L, 0L))
})

testthat::test_that("the fold-change of de_confidence is a difference of mean log counts", {
  # one feature, group a = (0, 0, 100), group b = (30, 30, 30): the log
  # ratio of the arithmetic means is negative, the documented statistic (the
  # difference of the group means of log2(count + 5)) is positive
  counts <- matrix(c(0, 0, 100, 30, 30, 30), nrow = 1L)
  in_b <- c(FALSE, FALSE, FALSE, TRUE, TRUE, TRUE)
  w <- welch_matrix(log2(counts + 5), in_b)
  expected <- log2(35) - mean(log2(c(5, 5, 105)))
  testthat::expect_equal(w$lfc, expected, tolerance = 1e-12)
  testthat::expect_gt(w$lfc, 1.3)
  testthat::expect_lt(log2(35) - log2(mean(c(0, 0, 100)) + 5), 0)
})
