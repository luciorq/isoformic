testthat::test_that("IsoformicExperiment is an S7 class with a working empty constructor", {
  testthat::expect_true(S7::S7_inherits(IsoformicExperiment))
  iso <- IsoformicExperiment()
  testthat::expect_s7_class(iso, IsoformicExperiment)
  testthat::expect_match(iso@experiment_name, "^isoformic_")
  testthat::expect_null(iso@annot_path)
  testthat::expect_equal(dim(iso), c(0L, 0L))
  testthat::expect_equal(dimnames(iso), list(NULL, NULL))
  testthat::expect_equal(nrow(annot_data_transcripts(iso)), 0L)
  testthat::expect_named(tx_type_palette(), names(iso@tx_type_palette))
  testthat::expect_output(print(iso), "<IsoformicExperiment>")
})

testthat::test_that("construction validates names, assays, col_data and annotation", {
  annot <- mini_annotation()
  testthat::expect_error(
    IsoformicExperiment(experiment_name = c("a", "b")),
    class = "isoformic_not_string"
  )
  testthat::expect_error(
    IsoformicExperiment(assays = list(abundance = mini_assay())),
    class = "isoformic_missing_annotation"
  )
  testthat::expect_error(
    IsoformicExperiment(annot_path = annot, assays = list(mini_assay())),
    class = "isoformic_invalid_assays"
  )
  bad_rows <- mini_assay()
  rownames(bad_rows)[1] <- "unknown_tx"
  testthat::expect_error(
    IsoformicExperiment(
      annot_path = annot,
      assays = list(abundance = bad_rows)
    ),
    class = "isoformic_invalid_assay_transcript_ids"
  )
  testthat::expect_error(
    IsoformicExperiment(
      annot_path = annot,
      assays = list(abundance = mini_assay()),
      col_data = data.frame(sample_id = c("s1", "s2"))
    ),
    class = "isoformic_invalid_assay_sample_ids"
  )
  testthat::expect_error(
    IsoformicExperiment(
      annot_path = annot,
      assays = list(abundance = mini_assay(), counts = mini_assay()[, 1:2])
    ),
    class = "isoformic_invalid_assays"
  )
  testthat::expect_error(
    IsoformicExperiment(annot_path = withr::local_tempdir()),
    class = "isoformic_annot_path_dont_exist"
  )
})

testthat::test_that("assays and col_data are standardised", {
  annot <- mini_annotation()
  df <- as.data.frame(mini_assay())
  df$transcript_id <- rownames(df)
  cd <- data.frame(condition = c("a", "a", "b", "b"))
  rownames(cd) <- c("s1", "s2", "s3", "s4")
  iso <- IsoformicExperiment(
    annot_path = annot,
    assays = list(abundance = df),
    col_data = cd
  )
  testthat::expect_true(is.matrix(assay(iso)))
  testthat::expect_equal(rownames(iso), c("tx1", "tx2", "tx3"))
  testthat::expect_equal(colnames(col_data(iso))[1], "sample_id")
  testthat::expect_equal(col_data(iso)$sample_id, c("s1", "s2", "s3", "s4"))

  # a single matrix is accepted as the abundance assay
  iso2 <- IsoformicExperiment(annot_path = annot, assays = mini_assay())
  testthat::expect_equal(assay_names(iso2), "abundance")
  testthat::expect_equal(col_names(iso2), c("s1", "s2", "s3", "s4"))
  testthat::expect_equal(colnames(col_data(iso2)), "sample_id")
})

testthat::test_that("accessors return the expected tables and matrices", {
  iso <- mini_experiment(with_dea = FALSE)
  testthat::expect_equal(assay_names(iso), c("abundance", "counts"))
  testthat::expect_equal(assay(iso, "counts"), mini_assay() * 2)
  testthat::expect_equal(assay(iso, 2), mini_assay() * 2)
  testthat::expect_error(
    assay(iso, "nope"),
    class = "isoformic_assay_not_found"
  )
  testthat::expect_equal(row_names(iso), c("tx1", "tx2", "tx3"))
  testthat::expect_equal(col_names(iso), c("s1", "s2", "s3", "s4"))
  testthat::expect_equal(dim(iso), c(3L, 4L))
  testthat::expect_equal(nrow(iso), 3L)

  rd <- row_data(iso)
  testthat::expect_equal(rd$transcript_id, c("tx1", "tx2", "tx3"))
  testthat::expect_equal(rd$gene_name, c("GENEA", "GENEA", "GENEB"))
  testthat::expect_false(any(
    c("seqid", "start_pos", "end_pos", "strand") %in% colnames(rd)
  ))
  testthat::expect_equal(nrow(tx_annot(iso)), 3L)
  testthat::expect_equal(
    colnames(tx_to_gene(iso)),
    c("transcript_id", "gene_id")
  )
  testthat::expect_equal(nrow(annot_data_genes(iso)), 2L)
  testthat::expect_equal(nrow(annot_data_exons(iso)), 4L)
  testthat::expect_equal(annot_metadata(iso)$annot_name, "mini")
  testthat::expect_true(all(
    c("gene_name", "gene_start", "gene_end") %in% colnames(annot_data(iso))
  ))
  testthat::expect_true(
    "exon_start" %in% colnames(annot_data(iso, with_exons = TRUE))
  )
  testthat::expect_equal(nrow(exon_table(iso, "GENEA")), 3L)
})

testthat::test_that("gene assays are derived by summing isoforms, cached and overridable", {
  iso <- mini_experiment(with_dea = FALSE)
  ga <- gene_assays(iso)
  testthat::expect_named(ga, c("abundance", "counts"))
  testthat::expect_equal(
    ga$abundance["gene1", ],
    colSums(mini_assay()[c("tx1", "tx2"), ])
  )
  testthat::expect_equal(ga$abundance["gene2", ], mini_assay()["tx3", ])
  testthat::expect_equal(gene_names(iso), c("gene1", "gene2"))
  testthat::expect_equal(gene_data(iso)$gene_name, c("GENEA", "GENEB"))
  testthat::expect_true("gene_assays" %in% ls(attr(iso, ".cache")))

  # modifying the transcript assays invalidates the derived gene assays
  assay(iso, "abundance") <- mini_assay() * 10
  testthat::expect_equal(
    gene_assay(iso, "abundance")["gene2", ],
    mini_assay()["tx3", ] * 10
  )

  # explicit gene assays win
  custom <- matrix(
    1,
    2,
    4,
    dimnames = list(c("gene1", "gene2"), c("s1", "s2", "s3", "s4"))
  )
  gene_assays(iso) <- list(abundance = custom)
  testthat::expect_equal(gene_assay(iso, "abundance"), custom)
  testthat::expect_error(
    {
      gene_assays(iso) <- list(
        abundance = matrix(
          1,
          1,
          4,
          dimnames = list("nope", c("s1", "s2", "s3", "s4"))
        )
      )
    },
    class = "isoformic_invalid_assay_gene_ids"
  )
})

testthat::test_that("length assays are summarised as abundance-weighted means", {
  tx2gene <- data.frame(
    transcript_id = c("tx1", "tx2", "tx3"),
    gene_id = c("g", "g", "h")
  )
  abundance <- matrix(
    c(1, 3, 5, 0, 0, 2),
    3,
    2,
    dimnames = list(c("tx1", "tx2", "tx3"), c("a", "b"))
  )
  length_mat <- matrix(
    c(100, 300, 50, 100, 300, 50),
    3,
    2,
    dimnames = dimnames(abundance)
  )
  res <- summarize_to_gene(
    list(abundance = abundance, length = length_mat),
    tx2gene
  )
  testthat::expect_equal(res$abundance["g", ], c(a = 4, b = 0))
  testthat::expect_equal(res$length["g", "a"], (1 * 100 + 3 * 300) / 4)
  testthat::expect_equal(res$length["g", "b"], 200) # plain mean when abundance is zero
  single <- summarize_to_gene(abundance, tx2gene)
  testthat::expect_equal(single["h", ], c(a = 5, b = 2))
})

testthat::test_that("subsetting behaves like SummarizedExperiment", {
  iso <- mini_experiment(with_dea = FALSE)
  sub <- iso[c("tx1", "tx3"), c("s1", "s4")]
  testthat::expect_equal(dim(sub), c(2L, 2L))
  testthat::expect_equal(rownames(sub), c("tx1", "tx3"))
  testthat::expect_equal(col_data(sub)$sample_id, c("s1", "s4"))
  testthat::expect_equal(
    assay(sub, "counts"),
    (mini_assay() * 2)[c("tx1", "tx3"), c("s1", "s4")]
  )
  testthat::expect_equal(dim(iso[, 1:2]), c(3L, 2L))
  testthat::expect_equal(dim(iso[1:2, ]), c(2L, 4L))
  testthat::expect_equal(dim(iso[c(TRUE, FALSE, TRUE), ]), c(2L, 4L))
  testthat::expect_equal(gene_names(sub), c("gene1", "gene2"))
  testthat::expect_error(
    IsoformicExperiment()[1, ],
    class = "isoformic_no_assays"
  )
})

testthat::test_that("print summarises all levels", {
  iso <- mini_experiment()
  out <- paste(capture.output(print(iso)), collapse = "\n")
  testthat::expect_match(out, "2 genes, 3 transcripts, 4 exons")
  testthat::expect_match(out, "samples    : 4 \\[condition\\]")
  testthat::expect_match(out, "abundance, counts")
  testthat::expect_match(out, "det \\(3 transcripts\\), deg \\(2 genes\\)")
})

testthat::test_that("annotation tables are cached by file signature and refreshed on change", {
  annot <- mini_annotation()
  iso <- IsoformicExperiment(
    annot_path = annot,
    assays = list(abundance = mini_assay())
  )
  clear_cache(iso, annotation = TRUE)
  testthat::expect_equal(nrow(annot_data_genes(iso)), 2L)
  key <- paste0("genes::", annot)
  testthat::expect_true(key %in% ls(annotation_cache))
  # Rewrite the genes table with an extra gene: the cache must notice
  genes <- read_parquet(fs::path(annot, "annot_data_genes.parquet"))
  genes <- rbind(genes, genes[1, ])
  genes$gene_id[3] <- "gene3"
  Sys.sleep(0.05)
  write_parquet(genes, fs::path(annot, "annot_data_genes.parquet"))
  testthat::expect_equal(nrow(annot_data_genes(iso)), 3L)
  clear_cache(iso, annotation = TRUE)
  testthat::expect_length(ls(annotation_cache), 0L)
})

testthat::test_that("per-object cache is keyed by inputs so copies never go stale", {
  iso <- mini_experiment(with_dea = FALSE)
  frac1 <- isoform_fraction(iso)
  iso2 <- iso
  assay(iso2, "abundance") <- mini_assay() * 2
  frac2 <- isoform_fraction(iso2)
  testthat::expect_equal(frac1, frac2) # fractions are scale invariant
  testthat::expect_equal(
    gene_assay(iso2, "abundance")["gene2", ],
    mini_assay()["tx3", ] * 2
  )
  testthat::expect_equal(
    gene_assay(iso, "abundance")["gene2", ],
    mini_assay()["tx3", ]
  )
})

testthat::test_that("tx_type_palette property validates and defaults", {
  iso <- IsoformicExperiment()
  testthat::expect_error(
    {
      iso@tx_type_palette <- c("#ffffff", "#000000")
    },
    class = "isoformic_invalid_tx_type_palette"
  )
  iso@tx_type_palette <- c(gene = "#000000")
  testthat::expect_equal(iso@tx_type_palette, c(gene = "#000000"))
  iso@tx_type_palette <- character()
  testthat::expect_equal(iso@tx_type_palette, tx_type_palette())
})

testthat::test_that("write_isoformic and read_isoformic round-trip the object", {
  iso <- mini_experiment()
  iso@metadata <- list(source = "test", n = 3L)
  iso@gsea <- list(
    fake = data.frame(
      pathway = "p",
      NES = 1.2,
      padj = 0.01,
      size = 3L,
      experiment = "x"
    )
  )
  dir <- fs::path(withr::local_tempdir(), "exp")
  iso <- write_isoformic(iso, dir, copy_annotation = TRUE)
  testthat::expect_equal(iso@data_path, as.character(fs::path_abs(dir)))
  testthat::expect_true(fs::file_exists(fs::path(dir, "experiment.parquet")))
  testthat::expect_true(fs::file_exists(fs::path(
    dir,
    "assays",
    "abundance.parquet"
  )))
  testthat::expect_true(fs::file_exists(fs::path(dir, "dea", "det.parquet")))
  testthat::expect_true(fs::file_exists(fs::path(dir, "gsea", "fake.parquet")))
  testthat::expect_true(is_annotation_dir(fs::path(dir, "annotation")))

  iso2 <- read_isoformic(dir)
  testthat::expect_equal(iso2@experiment_name, "mini")
  testthat::expect_equal(assay(iso2, "abundance"), assay(iso, "abundance"))
  testthat::expect_equal(assay(iso2, "counts"), assay(iso, "counts"))
  testthat::expect_equal(col_data(iso2), col_data(iso))
  testthat::expect_equal(iso2@metadata$source, "test")
  testthat::expect_equal(de_tx(iso2), de_tx(iso))
  testthat::expect_equal(de_gene(iso2), de_gene(iso))
  testthat::expect_named(iso2@gsea, "fake")
  testthat::expect_equal(
    iso2@annot_path,
    as.character(fs::path(dir, "annotation"))
  )

  testthat::expect_error(
    write_isoformic(iso, dir),
    class = "isoformic_dir_exists"
  )
  testthat::expect_error(
    read_isoformic(withr::local_tempdir()),
    class = "isoformic_invalid_experiment_dir"
  )

  # without copying, the original annotation path is stored and can be overridden
  dir2 <- fs::path(withr::local_tempdir(), "exp2")
  write_isoformic(iso, dir2)
  testthat::expect_false(fs::dir_exists(fs::path(dir2, "annotation")))
  iso3 <- read_isoformic(dir2, annot_path = iso@annot_path)
  testthat::expect_equal(nrow(row_data(iso3)), 3L)
})

testthat::test_that("as_isoformic converts tximport-style lists and strips GENCODE suffixes", {
  annot <- mini_annotation()
  ab <- mini_assay()
  rownames(ab) <- paste0(rownames(ab), "|ENSG|extra")
  txi <- list(abundance = ab, counts = ab * 2, countsFromAbundance = "no")
  iso <- as_isoformic(txi, annot_path = annot, col_data = mini_col_data())
  testthat::expect_equal(rownames(iso), c("tx1", "tx2", "tx3"))
  testthat::expect_equal(assay_names(iso), c("abundance", "counts"))
  testthat::expect_identical(as_isoformic(iso, annot_path = annot), iso)
  testthat::expect_error(
    as_isoformic("x", annot_path = annot),
    class = "isoformic_convert_error"
  )
  testthat::expect_error(
    as_isoformic(list(a = 1), annot_path = annot),
    class = "isoformic_invalid_assays"
  )

  extra <- rbind(mini_assay(), unknown = c(1, 1, 1, 1))
  testthat::expect_message(
    iso2 <- as_isoformic(
      list(abundance = extra),
      annot_path = annot,
      drop_unannotated = TRUE
    ),
    "Dropping 1 transcript"
  )
  testthat::expect_equal(nrow(iso2), 3L)
})

testthat::test_that("as_isoformic and as_summarized_experiment round-trip through SummarizedExperiment", {
  testthat::skip_if_not_installed("SummarizedExperiment")
  iso <- mini_experiment(with_dea = FALSE)
  se <- as_summarized_experiment(iso)
  testthat::expect_s4_class(se, "SummarizedExperiment")
  testthat::expect_equal(dim(se), c(3L, 4L))
  testthat::expect_equal(
    as.data.frame(SummarizedExperiment::rowData(se))$gene_name,
    c("GENEA", "GENEA", "GENEB")
  )
  se_gene <- as_summarized_experiment(iso, level = "gene")
  testthat::expect_equal(dim(se_gene), c(2L, 4L))
  testthat::expect_equal(SummarizedExperiment::assayNames(se_gene)[1], "counts")
  testthat::expect_false(
    "strand" %in% colnames(SummarizedExperiment::rowData(se_gene))
  )
  testthat::expect_s4_class(se_gene[1, ], "SummarizedExperiment")
  back <- as_isoformic(se, annot_path = iso@annot_path)
  testthat::expect_equal(assay(back, "abundance"), assay(iso, "abundance"))
  testthat::expect_equal(col_data(back)$condition, col_data(iso)$condition)
})

testthat::test_that("example_isoformic loads the bundled experiment", {
  iso <- example_isoformic()
  testthat::expect_s7_class(iso, IsoformicExperiment)
  testthat::expect_true(all(
    c("RBPJ", "EGFR", "PNCK", "FLT1") %in% row_data(iso)$gene_name
  ))
  testthat::expect_equal(ncol(iso), 25L)
  testthat::expect_equal(assay_names(iso), c("abundance", "counts"))
  testthat::expect_false(is.null(iso@dea$det))
  testthat::expect_null(iso@data_path)
})

testthat::test_that("read_salmon reads quant.sf files into an experiment", {
  quant_dirs <- list.files(
    system.file("extdata", "mini_quants", package = "isoformic"),
    full.names = TRUE
  )
  testthat::skip_if(length(quant_dirs) == 0L)
  quant <- read_tabular(fs::path(quant_dirs[1], "quant.sf"))
  ids <- sub("\\|.*$", "", quant$Name[1:5])
  gff <- withr::local_tempfile(fileext = ".gff3")
  lines <- c(
    "##gff-version 3",
    "chr1\tX\tgene\t1\t1000\t.\t+\t.\tID=g1;gene_id=g1;gene_name=G1;gene_type=protein_coding;"
  )
  for (i in seq_along(ids)) {
    lines <- c(
      lines,
      sprintf(
        paste0(
          "chr1\tX\ttranscript\t%d\t%d\t.\t+\t.\tID=%s;Parent=g1;gene_id=g1;",
          "transcript_id=%s;transcript_name=G1-%d;transcript_type=protein_coding;"
        ),
        i * 10,
        i * 10 + 100,
        ids[i],
        ids[i],
        i
      ),
      sprintf(
        "chr1\tX\texon\t%d\t%d\t.\t+\t.\tID=e%d;Parent=%s;transcript_id=%s;exon_number=1;exon_id=e%d;",
        i * 10,
        i * 10 + 50,
        i,
        ids[i],
        ids[i],
        i
      )
    )
  }
  writeLines(lines, gff)
  annot <- suppressMessages(prepare_isoformic_annotation(
    gff,
    output_path = withr::local_tempdir()
  ))
  names(quant_dirs) <- sub("_quant$", "", basename(quant_dirs))
  testthat::expect_message(
    iso <- read_salmon(quant_dirs, annot_path = annot, drop_unannotated = TRUE),
    "Dropping"
  )
  testthat::expect_equal(nrow(iso), 5L)
  testthat::expect_length(quant_dirs, ncol(iso))
  testthat::expect_equal(assay_names(iso), c("abundance", "counts", "length"))
  testthat::expect_equal(col_names(iso), names(quant_dirs))
  testthat::expect_equal(
    unname(assay(iso, "abundance")[ids[1], names(quant_dirs)[1]]),
    quant$TPM[1]
  )
  testthat::expect_error(
    read_salmon(quant_dirs, annot_path = annot),
    class = "isoformic_invalid_assay_transcript_ids"
  )
  dup <- quant_dirs[c(1, 1)]
  names(dup) <- c("a", "a")
  testthat::expect_error(
    read_salmon(dup, annot_path = annot),
    class = "isoformic_duplicated_sample_ids"
  )
})

testthat::test_that("the assay signature is computed when the assays are set", {
  iso <- mini_experiment(with_dea = FALSE)
  signature <- attr(iso, "assays_hash", exact = TRUE)
  testthat::expect_type(signature, "character")
  testthat::expect_identical(
    signature,
    rlang::hash(attr(iso, "assays", exact = TRUE))
  )
  # A change of the assays changes the signature, and so the cached values
  changed <- iso
  assays <- attr(changed, "assays", exact = TRUE)
  assays[[1L]][1L, 1L] <- assays[[1L]][1L, 1L] + 1
  changed@assays <- assays
  testthat::expect_false(
    identical(attr(changed, "assays_hash", exact = TRUE), signature)
  )
  testthat::expect_false(identical(
    gene_assay(changed, names(assays)[[1L]]),
    gene_assay(iso, names(assays)[[1L]])
  ))
  # Subsetting rebuilds the object and its signature
  sub <- iso[1:2, ]
  testthat::expect_identical(
    attr(sub, "assays_hash", exact = TRUE),
    rlang::hash(attr(sub, "assays", exact = TRUE))
  )
})
