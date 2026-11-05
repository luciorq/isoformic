testthat::test_that("read_treeterminus reads the real Newick output of TreeTerminus", {
  # Three trees of `treeterminus group --mean_inf true` on six Salmon 2.8.0
  # samples; leaves are positions in quant.sf
  file <- testthat::test_path("testdata", "treeterminus_trees.txt")
  transcripts <- paste0("ENST", sprintf("%08d", 0:300000), ".1|gene|x|")
  trees <- read_treeterminus(file, transcripts = transcripts)
  testthat::expect_identical(
    colnames(trees),
    c("tree_id", "node_id", "transcript_id", "parent_id")
  )
  testthat::expect_identical(attr(trees, "source"), "treeterminus")
  leaves <- trees[!is.na(trees$transcript_id), ]
  nodes <- trees[!is.na(trees$node_id), ]
  # 3 + 15 + 7 leaves; a binary tree of n leaves has n - 1 internal nodes
  testthat::expect_identical(as.vector(table(leaves$tree_id)), c(3L, 15L, 7L))
  testthat::expect_identical(as.vector(table(nodes$tree_id)), c(2L, 14L, 6L))
  testthat::expect_identical(sum(is.na(nodes$parent_id)), 3L)
  testthat::expect_true(all(nodes$parent_id < nodes$node_id, na.rm = TRUE))
  # Positions are 0-based, and names are cut at the first bar
  first <- leaves[leaves$tree_id == 1L, ]
  testthat::expect_setequal(
    first$transcript_id,
    sprintf("ENST%08d.1", c(285370, 285371, 285365))
  )
  inner <- first$parent_id[first$transcript_id == "ENST00285370.1"]
  outer <- first$parent_id[first$transcript_id == "ENST00285365.1"]
  testthat::expect_identical(nodes$parent_id[nodes$node_id == inner], outer)
  # The folder form finds the file
  dir <- withr::local_tempdir()
  fs::file_copy(file, fs::path(dir, "cluster_nwk.txt"))
  testthat::expect_equal(
    read_treeterminus(dir, transcripts = transcripts),
    trees
  )
})

testthat::test_that("read_treeterminus names leaves from a Salmon directory and checks its input", {
  x <- salmon_uncertainty_fixture()
  quant <- read.delim(fs::path(x$dirs[[1L]], "quant.sf"))
  nwk <- withr::local_tempfile(fileext = ".txt")
  writeLines("((0,1),2);", nwk)
  trees <- read_treeterminus(nwk, salmon_dir = x$dirs[[1L]])
  testthat::expect_setequal(
    trees$transcript_id[!is.na(trees$transcript_id)],
    sub("\\|.*$", "", quant$Name[1:3])
  )
  testthat::expect_error(
    read_treeterminus(nwk),
    class = "isoformic_invalid_treeterminus"
  )
  testthat::expect_error(
    read_treeterminus(nwk, transcripts = "only_one"),
    class = "isoformic_invalid_treeterminus"
  )
  writeLines("((0,1),2;", nwk)
  testthat::expect_error(
    read_treeterminus(nwk, transcripts = letters),
    class = "isoformic_invalid_treeterminus"
  )
  writeLines(character(), nwk)
  testthat::expect_error(
    read_treeterminus(nwk, transcripts = letters),
    class = "isoformic_invalid_treeterminus"
  )
  testthat::expect_error(
    read_treeterminus(withr::local_tempdir()),
    class = "isoformic_invalid_treeterminus"
  )
  # Named leaves and branch lengths
  writeLines("((tA1:0.1,tA2:0.2):0.5,tA3);", nwk)
  named <- read_treeterminus(nwk)
  testthat::expect_setequal(
    named$transcript_id[!is.na(named$transcript_id)],
    c("tA1", "tA2", "tA3")
  )
})

testthat::test_that("cut_transcript_trees groups uncertain isoforms at the node the reads determine", {
  x <- uncertainty_experiment()
  iso <- add_uncertainty(x$iso, paths = x$fixture$dirs)
  nwk <- withr::local_tempfile(fileext = ".txt")
  # tA1 and tA2 are anti-correlated: each is uncertain, their sum is not;
  # tA4 is overdispersed and its tree never gets below the cutoff
  writeLines(c("((tA1,tA2),tA3);", "(tA4,tB1);", "(tB2,tA6);"), nwk)
  trees <- read_treeterminus(nwk)
  cut <- cut_transcript_trees(trees, iso)
  testthat::expect_identical(attr(cut, "source"), "treeterminus")
  testthat::expect_identical(attr(cut, "cutoff"), 0.5)
  testthat::expect_setequal(cut$transcript_id, c("tA1", "tA2"))
  testthat::expect_identical(unique(cut$group_size), 2L)
  testthat::expect_lt(unique(cut$group_inf_rv), 0.5)
  # The roots of unresolved trees, on request
  rooted <- cut_transcript_trees(trees, iso, unresolved = "root")
  testthat::expect_setequal(rooted$transcript_id, c("tA1", "tA2", "tA4", "tB1"))
  testthat::expect_length(
    unique(rooted$group_id[rooted$transcript_id %in% c("tA4", "tB1")]),
    1L
  )
  # A well-determined leaf inside a chosen node joins it: the groups never overlap
  writeLines("(((tA1,tA2),tA3),tB2);", nwk)
  nested <- cut_transcript_trees(read_treeterminus(nwk), iso, cutoff = 0.5)
  testthat::expect_false(anyDuplicated(nested$transcript_id) > 0L)

  # add_uncertainty uses the cut, and the methods text says so
  grouped <- add_uncertainty(iso, groups = cut)
  testthat::expect_identical(
    grouped@metadata$uncertainty_ambiguity_source,
    "treeterminus"
  )
  unc <- tx_uncertainty(grouped)
  testthat::expect_identical(
    unique(unc$ambiguity_group[unc$transcript_id %in% c("tA1", "tA2")]),
    unique(cut$group_id)
  )
  testthat::expect_match(
    uncertainty_methods(grouped),
    "TreeTerminus",
    fixed = TRUE
  )
  conf <- de_confidence(grouped, group = "condition", collapse_groups = TRUE)
  testthat::expect_true(unique(cut$group_id) %in% conf$feature_id)

  testthat::expect_error(
    cut_transcript_trees(data.frame(a = 1), iso),
    class = "isoformic_invalid_treeterminus"
  )
  writeLines("(zz1,zz2);", nwk)
  testthat::expect_error(
    cut_transcript_trees(read_treeterminus(nwk), iso),
    class = "isoformic_invalid_treeterminus"
  )
})
