testthat::test_that("prepare_isoformic_annotation parses GFF3 into parquet tables", {
  annot_dir <- mini_annotation("mini.gff")
  testthat::expect_true(is_annotation_dir(annot_dir))
  testthat::expect_setequal(
    list.files(annot_dir),
    c(
      "annot_data_exons.parquet",
      "annot_data_genes.parquet",
      "annot_data_transcripts.parquet",
      "annot_metadata.parquet",
      "annot_seq_regions.parquet"
    )
  )
  genes <- read_parquet(fs::path(annot_dir, "annot_data_genes.parquet"))
  tx <- read_parquet(fs::path(annot_dir, "annot_data_transcripts.parquet"))
  exons <- read_parquet(fs::path(annot_dir, "annot_data_exons.parquet"))
  testthat::expect_equal(nrow(genes), 2L)
  testthat::expect_equal(nrow(tx), 3L)
  testthat::expect_equal(nrow(exons), 4L)
  testthat::expect_equal(genes$gene_name, c("GENEA", "GENEB"))
  testthat::expect_equal(tx$n_exons[tx$transcript_id == "tx1"], 2L)
  testthat::expect_equal(tx$tx_length[tx$transcript_id == "tx1"], 202L)
  testthat::expect_equal(
    tx$transcript_type[tx$transcript_id == "tx2"],
    "retained_intron"
  )
  testthat::expect_type(exons$exon_number, "integer")
  meta <- kv_to_list(read_parquet(fs::path(
    annot_dir,
    "annot_metadata.parquet"
  )))
  testthat::expect_equal(meta$annot_name, "mini")
  testthat::expect_equal(meta$file_type, "gff")
  testthat::expect_equal(meta[["gff-version"]], "3")
  seq_regions <- read_parquet(fs::path(annot_dir, "annot_seq_regions.parquet"))
  testthat::expect_equal(seq_regions$seqid, c("chr1", "chr2"))
})

testthat::test_that("GTF input yields the same tables as GFF3", {
  gff <- prepare_annotation(mini_annotation("mini.gff"))
  gtf <- prepare_annotation(mini_annotation("mini.gtf"))
  testthat::expect_equal(gff$genes, gtf$genes)
  testthat::expect_equal(gff$transcripts, gtf$transcripts)
  testthat::expect_equal(gff$exons, gtf$exons)
  testthat::expect_equal(gtf$metadata$file_type, "gtf")
})

testthat::test_that("prepare_isoformic_annotation caches and can overwrite", {
  out <- withr::local_tempdir()
  gff <- testthat::test_path("testdata", "mini.gff")
  dir1 <- suppressMessages(prepare_isoformic_annotation(gff, output_path = out))
  testthat::expect_message(
    dir2 <- prepare_isoformic_annotation(gff, output_path = out),
    "cached"
  )
  testthat::expect_identical(dir1, dir2)
  testthat::expect_message(
    prepare_isoformic_annotation(gff, output_path = out, overwrite = TRUE),
    "Parsing"
  )
  testthat::expect_error(
    prepare_isoformic_annotation("/no/such.gff"),
    class = "isoformic_annot_file_dont_exist"
  )
  bad <- withr::local_tempfile(fileext = ".txt")
  writeLines("x", bad)
  testthat::expect_error(
    prepare_isoformic_annotation(bad),
    class = "isoformic_unsupported_file_type"
  )
})

testthat::test_that("read_annotation_header parses pragmas, comments and sequence regions", {
  tmp <- withr::local_tempfile(fileext = ".gff3")
  writeLines(
    c(
      "##gff-version 3",
      "#description: test annotation",
      "#provider: TEST",
      "##sequence-region chr1 1 1000",
      "##sequence-region chr2 1 500",
      "chr1\tsrc\tgene\t1\t10\t.\t+\t.\tID=g"
    ),
    tmp
  )
  header <- read_annotation_header(tmp)
  testthat::expect_equal(header$metadata$description, "test annotation")
  testthat::expect_equal(header$metadata$provider, "TEST")
  testthat::expect_equal(header$metadata[["gff-version"]], "3")
  testthat::expect_equal(header$seq_regions$seqid, c("chr1", "chr2"))
  testthat::expect_equal(header$seq_regions$end_pos, c(1000L, 500L))
})

testthat::test_that("prepare_exon_annotation and build_exon_table select genes", {
  annot_dir <- mini_annotation()
  res <- prepare_exon_annotation(c("GENEA", "GENEB"), annot_dir)
  testthat::expect_equal(nrow(res), 4L)
  testthat::expect_setequal(unique(res$transcript_id), c("tx1", "tx2", "tx3"))
  testthat::expect_true(all(
    c(
      "transcript_name",
      "transcript_type",
      "gene_name",
      "start_pos",
      "end_pos"
    ) %in%
      colnames(res)
  ))
  one <- prepare_exon_annotation("GENEB", annot_dir)
  testthat::expect_equal(nrow(one), 1L)
  testthat::expect_warning(
    prepare_exon_annotation(c("GENEB", "NOPE"), annot_dir),
    "not found"
  )
})

testthat::test_that("read_annotation_table selects rows in DuckDB as %in% does", {
  annot_dir <- mini_annotation("mini_cds.gff3")
  tx <- read_annotation_table(annot_dir, "transcripts")
  wanted <- c(rev(unique(tx$gene_id)), tx$gene_id[1L], "NOPE")
  testthat::expect_identical(
    read_annotation_table(
      annot_dir,
      "transcripts",
      filter = list(gene_id = wanted)
    ),
    tx[tx$gene_id %in% wanted, ]
  )
  one <- tx$gene_id[nrow(tx)]
  testthat::expect_identical(
    read_annotation_table(
      annot_dir,
      "transcripts",
      columns = c("transcript_id", "gene_id"),
      filter = list(gene_id = one)
    ),
    tx[tx$gene_id %in% one, c("transcript_id", "gene_id")]
  )
  testthat::expect_identical(
    read_annotation_table(
      annot_dir,
      "transcripts",
      filter = list(gene_id = character())
    ),
    tx[0L, ]
  )
  testthat::expect_error(
    read_annotation_table(annot_dir, "transcripts", filter = list("x")),
    class = "isoformic_invalid_argument"
  )
})

testthat::test_that("prepare_exon_annotation and exon_table read only the selected genes", {
  annot_dir <- mini_annotation("mini_cds.gff3")
  full <- prepare_annotation(annot_dir)
  for (genes in list("GENEA", c("GENEB", "GENEA"), full$genes$gene_id[2L])) {
    testthat::expect_identical(
      prepare_exon_annotation(genes, annot_dir),
      build_exon_table(full$genes, full$transcripts, full$exons, genes)
    )
  }
  testthat::expect_warning(
    none <- prepare_exon_annotation("NOPE", annot_dir),
    "not found"
  )
  testthat::expect_identical(
    none,
    suppressWarnings(
      build_exon_table(full$genes, full$transcripts, full$exons, "NOPE")
    )
  )
  iso <- example_isoformic()
  clear_cache(iso, annotation = TRUE)
  from_dir <- exon_table(iso, c("RBPJ", "EGFR"))
  testthat::expect_false(annot_table_cached(iso, "exons"))
  invisible(iso@annot_data_exons)
  testthat::expect_true(annot_table_cached(iso, "exons"))
  testthat::expect_identical(exon_table(iso, c("RBPJ", "EGFR")), from_dir)
})

testthat::test_that("make_tx_to_gene parses FASTA headers and annotation directories", {
  lines <- c(
    ">ENST00000001.1|ENSG00000001.1|OTTHUMG00000001.1|OTTHUMT00000001.1|GENEA-201|GENEA|1000|protein_coding|",
    "ACGTACGTACGT",
    ">ENST00000002.1|ENSG00000001.1|OTTHUMG00000001.1|OTTHUMT00000002.1|GENEA-202|GENEA|500|retained_intron|",
    "ACGTACGT"
  )
  fasta <- withr::local_tempfile(lines = lines, fileext = ".fa")
  res <- make_tx_to_gene(fasta)
  testthat::expect_s3_class(res, "tbl_df")
  testthat::expect_equal(colnames(res), tx_to_gene_cols())
  testthat::expect_equal(
    res$transcript_id,
    c("ENST00000001.1", "ENST00000002.1")
  )
  testthat::expect_equal(res$tx_length, c(1000L, 500L))
  testthat::expect_equal(
    res$transcript_type,
    c("protein_coding", "retained_intron")
  )

  from_annot <- make_tx_to_gene(mini_annotation())
  testthat::expect_equal(colnames(from_annot), tx_to_gene_cols())
  testthat::expect_equal(from_annot$gene_name, c("GENEA", "GENEA", "GENEB"))
  testthat::expect_true(all(is.na(from_annot$havana_gene_id)))

  from_gff <- make_tx_to_gene(
    testthat::test_path("testdata", "mini.gff"),
    output_path = withr::local_tempdir()
  )
  testthat::expect_equal(from_gff$transcript_id, from_annot$transcript_id)
  testthat::expect_error(
    make_tx_to_gene("/nonexistent.fa"),
    class = "isoformic_annot_file_dont_exist"
  )
})

testthat::test_that("resolve_annot_path accepts directories and files, rejects others", {
  annot_dir <- mini_annotation()
  testthat::expect_identical(
    resolve_annot_path(annot_dir),
    as.character(fs::path_abs(annot_dir))
  )
  out <- withr::local_tempdir()
  testthat::expect_true(is_annotation_dir(suppressMessages(
    resolve_annot_path(
      testthat::test_path("testdata", "mini.gff"),
      output_path = out
    )
  )))
  testthat::expect_error(
    resolve_annot_path("/no/such/path"),
    class = "isoformic_annot_path_dont_exist"
  )
  testthat::expect_error(
    resolve_annot_path(withr::local_tempdir()),
    class = "isoformic_annot_path_dont_exist"
  )
})

testthat::test_that("assert helpers raise typed errors", {
  testthat::expect_true(assert_file_exists(testthat::test_path(
    "testdata",
    "mini.gff"
  )))
  testthat::expect_error(
    assert_file_exists("/nope"),
    class = "isoformic_annot_file_dont_exist"
  )
  testthat::expect_error(
    assert_dir_exists("/nope"),
    class = "isoformic_dir_dont_exist"
  )
  testthat::expect_error(check_string(1), class = "isoformic_not_string")
  testthat::expect_null(check_string(NULL, allow_null = TRUE))
})

testthat::test_that("check_installed reports availability", {
  rlang::local_interactive(FALSE)
  testthat::expect_true(check_installed("testthat"))
  testthat::expect_false(check_installed("someNonExistentPackage12345"))
})

testthat::test_that("NCBI RefSeq style GFF3 (T2T-CHM13) is parsed with prefixes stripped", {
  annot <- prepare_annotation(mini_annotation("mini_ncbi.gff3"))
  testthat::expect_setequal(
    annot$genes$gene_id,
    c("OR4F16", "LOC127239154", "LOC999")
  )
  testthat::expect_equal(
    annot$genes$gene_type[annot$genes$gene_id == "OR4F16"],
    "protein_coding"
  )
  tx <- annot$transcripts
  testthat::expect_setequal(
    tx$transcript_id,
    c("NM_001005277.1", "NR_182074.1", "LOC999")
  )
  testthat::expect_equal(
    tx$transcript_type[tx$transcript_id == "NM_001005277.1"],
    "protein_coding"
  )
  testthat::expect_equal(
    tx$transcript_type[tx$transcript_id == "NR_182074.1"],
    "lncRNA"
  )
  # a pseudogene carrying exons directly becomes its own transcript
  testthat::expect_equal(tx$gene_id[tx$transcript_id == "LOC999"], "LOC999")
  testthat::expect_equal(
    tx$transcript_type[tx$transcript_id == "LOC999"],
    "pseudogene"
  )
  testthat::expect_equal(tx$n_exons[tx$transcript_id == "NM_001005277.1"], 2L)
  testthat::expect_equal(nrow(annot$exons), 5L)
  testthat::expect_true(all(annot$exons$transcript_id %in% tx$transcript_id))
  testthat::expect_equal(annot$metadata[["genome-build"]], "T2T-CHM13v2.0")
  testthat::expect_equal(annot$seq_regions$end_pos, 248387328L)
})

testthat::test_that("Ensembl Rapid Release style GFF3 is parsed", {
  annot <- prepare_annotation(mini_annotation("mini_ensembl.gff3"))
  testthat::expect_setequal(
    annot$genes$gene_id,
    c("ENSG05220060769", "ENSG05220062570")
  )
  testthat::expect_equal(
    annot$genes$gene_name[annot$genes$gene_id == "ENSG05220060769"],
    "OR4F29"
  )
  testthat::expect_equal(
    annot$genes$gene_name[annot$genes$gene_id == "ENSG05220062570"],
    "ENSG05220062570"
  )
  testthat::expect_equal(annot$genes$gene_type, c("protein_coding", "lncRNA"))
  tx <- annot$transcripts
  testthat::expect_equal(
    tx$gene_id[tx$transcript_id == "ENST05220237111"],
    "ENSG05220062570"
  )
  testthat::expect_equal(tx$transcript_type, c("protein_coding", "lncRNA"))
  exons <- annot$exons[annot$exons$transcript_id == "ENST05220237111", ]
  testthat::expect_equal(exons$exon_number, c(1L, 2L)) # from `rank`
  testthat::expect_equal(exons$exon_id, c("ENSE05220757331", "ENSE05220757330"))
})

testthat::test_that("T2T Liftoff/CAT GFF3 keeps source identifiers and infers types", {
  annot <- prepare_annotation(mini_annotation("mini_liftoff.gff3"))
  genes <- annot$genes
  testthat::expect_equal(
    genes$source_gene_id[genes$gene_id == "CHM13_G0000001"],
    "ENSG00000235146.2"
  )
  testthat::expect_true(is.na(genes$source_gene_id[
    genes$gene_id == "LOC112268260"
  ]))
  tx <- annot$transcripts
  testthat::expect_equal(
    tx$transcript_id[tx$gene_id == "LOC112268260"],
    "XM_047436352.1"
  )
  testthat::expect_equal(
    tx$transcript_type[tx$transcript_id == "XM_047436352.1"],
    "protein_coding"
  )
  testthat::expect_equal(
    tx$source_transcript_id[tx$transcript_id == "CHM13_T0000001"],
    "ENST00000423796.1"
  )
  testthat::expect_equal(
    tx$transcript_name[tx$transcript_id == "CHM13_T0000001"],
    "AC114498.1-201"
  )
  testthat::expect_equal(tx$n_exons, c(2L, 2L))
  # the annotation directory works as the backbone of an experiment
  m <- matrix(
    c(5, 10, 1, 2),
    2,
    dimnames = list(c("XM_047436352.1", "CHM13_T0000001"), c("a", "b"))
  )
  iso <- IsoformicExperiment(
    annot_path = mini_annotation("mini_liftoff.gff3"),
    assays = list(abundance = m)
  )
  testthat::expect_equal(
    row_data(iso)$source_transcript_id,
    c(NA, "ENST00000423796.1")
  )
})

testthat::test_that("transcript tags and the canonical flag are parsed from GFF3 and GTF", {
  gff <- withr::local_tempfile(fileext = ".gff3")
  writeLines(
    c(
      "##gff-version 3",
      "chr1\tX\tgene\t100\t5000\t.\t+\t.\tID=gene1;gene_id=gene1;gene_name=GENEA;gene_type=protein_coding;",
      paste0(
        "chr1\tX\ttranscript\t100\t5000\t.\t+\t.\tID=tx1;Parent=gene1;gene_id=gene1;",
        "transcript_id=tx1;transcript_type=protein_coding;tag=basic,MANE_Select,CCDS;"
      ),
      "chr1\tX\texon\t100\t200\t.\t+\t.\tID=e1;Parent=tx1;transcript_id=tx1;exon_number=1;",
      paste0(
        "chr1\tX\ttranscript\t600\t5000\t.\t+\t.\tID=tx2;Parent=gene1;",
        "gene_id=gene1;transcript_id=tx2;transcript_type=retained_intron;tag=basic;"
      ),
      "chr1\tX\texon\t600\t700\t.\t+\t.\tID=e2;Parent=tx2;transcript_id=tx2;exon_number=1;",
      paste0(
        "chr1\tX\ttranscript\t700\t5000\t.\t+\t.\tID=tx3;",
        "Parent=gene1;gene_id=gene1;transcript_id=tx3;transcript_type=lncRNA;"
      ),
      "chr1\tX\texon\t700\t800\t.\t+\t.\tID=e3;Parent=tx3;transcript_id=tx3;exon_number=1;"
    ),
    gff
  )
  out <- withr::local_tempdir()
  tx <- prepare_annotation(suppressMessages(prepare_isoformic_annotation(
    gff,
    output_path = out
  )))$transcripts
  testthat::expect_equal(tx$tags, c("basic,MANE_Select,CCDS", "basic", NA))
  testthat::expect_equal(tx$is_canonical, c(TRUE, FALSE, FALSE))
  testthat::expect_true(all(is.na(tx$haplotype)))

  gtf <- prepare_annotation(mini_annotation("mini_tags.gtf"))$transcripts
  testthat::expect_equal(gtf$tags, c("basic,MANE_Select,CCDS", "basic"))
  testthat::expect_equal(gtf$is_canonical, c(TRUE, FALSE))
})

testthat::test_that("StringTie GTF keeps reference identifiers and leaves novel types missing", {
  annot <- prepare_annotation(mini_annotation("mini_stringtie.gtf"))
  tx <- annot$transcripts
  testthat::expect_equal(tx$transcript_id, c("STRG.1.1", "STRG.1.2"))
  testthat::expect_equal(tx$source_transcript_id, c("ENST00000001.1", NA))
  testthat::expect_equal(tx$source_gene_id, c("ENSG00000001.1", NA))
  testthat::expect_true(all(is.na(tx$transcript_type)))
  testthat::expect_equal(tx$n_exons, c(2L, 2L))
  # the gene row is synthesised from the transcripts' parent
  testthat::expect_equal(annot$genes$gene_id, "STRG.1")
  testthat::expect_equal(annot$genes$gene_name, "GENEA")
})

testthat::test_that("add_sqanti_classification assigns structural categories", {
  annot_dir <- mini_annotation("mini_stringtie.gtf")
  classification <- data.frame(
    isoform = c("STRG.1.1", "STRG.1.2", "unrelated"),
    structural_category = c(
      "full-splice_match",
      "novel_not_in_catalog",
      "genic"
    ),
    associated_gene = c("GENEA", "GENEA", "X"),
    associated_transcript = c("ENST00000001.1", "novel", "novel")
  )
  testthat::expect_message(
    add_sqanti_classification(annot_dir, classification),
    "2 of 2"
  )
  tx <- read_annotation_table(annot_dir, "transcripts")
  testthat::expect_equal(
    tx$transcript_type,
    c("full-splice_match", "novel_not_in_catalog")
  )
  testthat::expect_equal(tx$source_transcript_id, c("ENST00000001.1", NA))
  testthat::expect_true(all(tx$transcript_type %in% names(tx_type_palette())))

  m <- matrix(
    c(5, 1, 2, 8),
    2,
    dimnames = list(c("STRG.1.1", "STRG.1.2"), c("a", "b"))
  )
  iso <- IsoformicExperiment(
    annot_path = annot_dir,
    assays = list(abundance = m)
  )
  testthat::expect_equal(
    row_data(iso)$transcript_type,
    c("full-splice_match", "novel_not_in_catalog")
  )
  testthat::expect_s3_class(plot_genomic_context(iso, "GENEA"), "ggplot")
  testthat::expect_error(
    add_sqanti_classification(annot_dir, data.frame(isoform = "x")),
    class = "isoformic_invalid_sqanti"
  )
  testthat::expect_error(
    add_sqanti_classification(
      annot_dir,
      data.frame(isoform = "nope", structural_category = "genic")
    ),
    class = "isoformic_invalid_sqanti"
  )
})

testthat::test_that("merge_annotations builds a haplotype-resolved annotation", {
  gff <- testthat::test_path("testdata", "mini_liftoff.gff3")
  out <- withr::local_tempdir()
  merged <- suppressMessages(merge_annotations(
    c(hap1 = gff, hap2 = gff),
    output_path = out,
    annot_name = "diploid"
  ))
  testthat::expect_true(is_annotation_dir(merged))
  annot <- prepare_annotation(merged)
  testthat::expect_equal(
    sort(unique(annot$transcripts$haplotype)),
    c("hap1", "hap2")
  )
  testthat::expect_true(all(grepl(
    "_hap[12]$",
    annot$transcripts$transcript_id
  )))
  testthat::expect_true(all(grepl("_hap[12]$", annot$genes$gene_id)))
  testthat::expect_true(all(
    annot$exons$transcript_id %in% annot$transcripts$transcript_id
  ))
  testthat::expect_equal(annot$metadata$haplotypes, "hap1;hap2")
  testthat::expect_equal(
    unique(annot$transcripts$source_transcript_id),
    c(NA, "ENST00000423796.1")
  )
  testthat::expect_message(
    merge_annotations(
      c(hap1 = gff, hap2 = gff),
      output_path = out,
      annot_name = "diploid"
    ),
    "cached"
  )
  testthat::expect_error(
    merge_annotations(c(gff, gff)),
    class = "isoformic_invalid_haplotypes"
  )
  testthat::expect_error(
    merge_annotations(c(a = gff)),
    class = "isoformic_invalid_haplotypes"
  )
})

testthat::test_that("exon and CDS attributes are parsed with placeholders and prefixes", {
  gff <- withr::local_tempfile(fileext = ".gff3")
  writeLines(
    c(
      "##gff-version 3",
      "chr1\tX\tgene\t100\t900\t.\t+\t.\tID=gene:g1;gene_id=g1;gene_type=protein_coding",
      "chr1\tX\tmRNA\t100\t900\t.\t+\t.\tID=transcript:t1;Parent=gene:g1",
      "chr1\tX\texon\t100\t200\t.\t+\t.\tParent=transcript:t1;exon_id=NA;rank=1",
      "chr1\tX\texon\t300\t400\t.\t+\t.\tParent=transcript:t1;exon_id=e2;exon_number=.;rank=2",
      "chr1\tX\texon\t800\t900\t.\t+\t.\tParent=transcript:t1;exon_id=__na;exon_number=3",
      "chr1\tX\tCDS\t150\t200\t.\t+\t0\tParent=transcript:t1",
      "chr1\tX\tCDS\t300\t400\t.\t+\t1\tParent=transcript:t1;exon_id=e2",
      "chr1\tX\tfive_prime_UTR\t100\t149\t.\t+\t.\tParent=transcript:t1",
      "chr2\tX\tgene\t10\t50\t.\t-\t.\tID=g2;Name=GENE2"
    ),
    gff
  )
  annot <- prepare_annotation(suppressMessages(prepare_isoformic_annotation(
    gff,
    output_path = withr::local_tempdir()
  )))
  testthat::expect_equal(annot$exons$transcript_id, rep("t1", 3L))
  testthat::expect_equal(annot$exons$exon_id, c(NA, "e2", NA))
  testthat::expect_equal(annot$exons$exon_number, 1:3)
  tx <- annot$transcripts
  testthat::expect_equal(tx$transcript_id, "t1")
  testthat::expect_equal(tx$gene_id, "g1")
  testthat::expect_equal(tx$transcript_type, "protein_coding")
  testthat::expect_equal(c(tx$cds_start, tx$cds_end), c(150L, 400L))
  testthat::expect_equal(tx$cds_length, 51L + 101L)
  testthat::expect_equal(c(tx$n_exons, tx$tx_length), c(3L, 101L * 3L))
  testthat::expect_equal(annot$genes$gene_id, c("g1", "g2"))
  testthat::expect_equal(annot$genes$gene_name, c("g1", "GENE2"))
  testthat::expect_equal(annot$seq_regions$end_pos, c(900L, 50L))
})

testthat::test_that("GTF UTR rows give gene attributes and sequence regions", {
  gtf <- withr::local_tempfile(fileext = ".gtf")
  writeLines(
    c(
      paste0(
        "chrA\tST\ttranscript\t100\t500\t.\t+\t.\t",
        "gene_id \"G1\"; transcript_id \"T1\";"
      ),
      paste0(
        "chrA\tST\texon\t100\t500\t.\t+\t.\t",
        "gene_id \"G1\"; transcript_id \"T1\"; exon_number \"1\";"
      ),
      paste0(
        "chrA\tST\tUTR\t100\t150\t.\t+\t.\t",
        "gene_id \"G1\"; transcript_id \"T1\"; gene_name \"GN1\"; ",
        "gene_type \"lncRNA\";"
      ),
      paste0(
        "chrZ\tST\tUTR\t10\t20\t.\t+\t.\t",
        "gene_id \"G9\"; transcript_id \"T9\";"
      )
    ),
    gtf
  )
  annot <- prepare_annotation(suppressMessages(prepare_isoformic_annotation(
    gtf,
    output_path = withr::local_tempdir()
  )))
  testthat::expect_equal(annot$genes$gene_id, "G1")
  testthat::expect_equal(annot$genes$gene_name, "GN1")
  testthat::expect_equal(annot$genes$gene_type, "lncRNA")
  testthat::expect_equal(annot$transcripts$transcript_type, "lncRNA")
  testthat::expect_equal(annot$seq_regions$seqid, c("chrA", "chrZ"))
})

testthat::test_that("merge_annotations keeps every row and column, in order", {
  out <- withr::local_tempdir()
  gff <- system.file("extdata", "example.gff3.gz", package = "isoformic")
  fasta <- system.file(
    "extdata",
    "example_transcripts.fa.gz",
    package = "isoformic"
  )
  dir1 <- suppressMessages(prepare_isoformic_annotation(
    gff,
    output_path = out,
    annot_name = "hap_a"
  ))
  suppressMessages(extract_transcript_sequences(dir1, transcript_fasta = fasta))
  dir2 <- suppressMessages(prepare_isoformic_annotation(
    testthat::test_path("testdata", "mini_liftoff.gff3"),
    output_path = out,
    annot_name = "hap_b",
    haplotype = "old_label"
  ))
  # A table written by an older version, without some columns
  tx2 <- read_parquet(annot_file_path(dir2, "transcripts"))
  write_parquet(
    tx2[, setdiff(colnames(tx2), c("tie_score", "haplotype"))],
    annot_file_path(dir2, "transcripts")
  )
  dir3 <- fs::path(out, "hap_c")
  fs::dir_copy(dir1, dir3)

  expected <- function(dirs, suffix) {
    tables <- c("genes", "transcripts", "exons", "seq_regions", "sequences")
    parts <- lapply(names(dirs), function(hap) {
      sfx <- if (suffix) paste0("_", hap) else ""
      lapply(stats::setNames(nm = tables), function(type) {
        file <- annot_file_path(dirs[[hap]], type)
        if (!fs::file_exists(file)) {
          return(NULL)
        }
        x <- read_parquet(file)
        if (!identical(type, "sequences")) {
          x$haplotype <- hap
        }
        for (col in intersect(c("gene_id", "transcript_id"), colnames(x))) {
          x[[col]] <- paste0(x[[col]], sfx)
        }
        if ("exon_id" %in% colnames(x)) {
          x$exon_id <- ifelse(is.na(x$exon_id), NA, paste0(x$exon_id, sfx))
        }
        x
      })
    })
    lapply(stats::setNames(nm = tables), function(type) {
      dplyr::bind_rows(lapply(parts, function(p) p[[type]]))
    })
  }

  dirs <- c(hap1 = dir1, hap2 = dir2)
  merged <- suppressMessages(merge_annotations(
    dirs,
    output_path = out,
    annot_name = "merged"
  ))
  want <- expected(dirs, suffix = TRUE)
  for (type in c("genes", "transcripts", "exons", "seq_regions")) {
    testthat::expect_identical(
      read_parquet(annot_file_path(merged, type)),
      want[[type]]
    )
  }
  testthat::expect_true("tie_score" %in% colnames(want$transcripts))
  # Sequences are merged only when every haplotype has them
  testthat::expect_false(has_sequences(merged))

  dirs <- c(hap1 = dir1, hap2 = as.character(dir3))
  merged <- suppressMessages(merge_annotations(
    dirs,
    output_path = out,
    annot_name = "merged_seq",
    id_suffix = FALSE
  ))
  want <- expected(dirs, suffix = FALSE)
  for (type in c("genes", "transcripts", "exons", "seq_regions", "sequences")) {
    testthat::expect_identical(
      read_parquet(annot_file_path(merged, type)),
      want[[type]]
    )
  }
  meta <- kv_to_list(read_annotation_table(merged, "metadata"))
  testthat::expect_equal(meta$id_suffix, "FALSE")
  testthat::expect_equal(meta$hap2.annot_name, "hap_a")

  # The sequences are copied with one thread; the connection gets its
  # number of threads back
  con <- duckdb_connect()
  DBI::dbExecute(con, "SET threads = 3")
  merge_annotation_table(
    con,
    files = as.character(annot_file_path(dirs, "sequences")),
    haplotypes = names(dirs),
    suffixes = NULL,
    type = "sequences",
    output_file = file.path(out, "sequences.parquet")
  )
  testthat::expect_identical(
    read_parquet(file.path(out, "sequences.parquet")),
    want$sequences
  )
  testthat::expect_equal(
    DBI::dbGetQuery(con, "SELECT current_setting('threads') AS n")$n,
    3
  )
})

testthat::test_that("sequence regions declared between features are all kept", {
  gff <- withr::local_tempfile(fileext = ".gff3")
  attr <- function(id, parent = NULL, extra = "") {
    paste0(
      "ID=",
      id,
      if (!is.null(parent)) paste0(";Parent=", parent),
      extra
    )
  }
  row <- function(seqid, type, start, end, attributes) {
    paste(
      seqid,
      "HAVANA",
      type,
      start,
      end,
      ".",
      "+",
      ".",
      attributes,
      sep = "\t"
    )
  }
  writeLines(
    c(
      "##gff-version 3",
      "##sequence-region chr1 1 1000",
      row(
        "chr1",
        "gene",
        10,
        200,
        attr("g1", extra = ";gene_id=g1;gene_name=A;gene_type=lncRNA")
      ),
      row(
        "chr1",
        "transcript",
        10,
        200,
        attr("t1", "g1", ";gene_id=g1;transcript_id=t1;transcript_type=lncRNA")
      ),
      row(
        "chr1",
        "exon",
        10,
        200,
        attr("exon:t1:1", "t1", ";gene_id=g1;transcript_id=t1;exon_number=1")
      ),
      "##sequence-region chr2 1 500",
      row(
        "chr2",
        "gene",
        5,
        90,
        attr("g2", extra = ";gene_id=g2;gene_name=B;gene_type=lncRNA")
      ),
      row(
        "chr2",
        "transcript",
        5,
        90,
        attr("t2", "g2", ";gene_id=g2;transcript_id=t2;transcript_type=lncRNA")
      ),
      row(
        "chr2",
        "exon",
        5,
        90,
        attr("exon:t2:1", "t2", ";gene_id=g2;transcript_id=t2;exon_number=1")
      ),
      # a sequence without a pragma keeps the largest feature end
      row(
        "chr3",
        "gene",
        5,
        70,
        attr("g3", extra = ";gene_id=g3;gene_name=C;gene_type=lncRNA")
      ),
      row(
        "chr3",
        "transcript",
        5,
        70,
        attr("t3", "g3", ";gene_id=g3;transcript_id=t3;transcript_type=lncRNA")
      ),
      row(
        "chr3",
        "exon",
        5,
        70,
        attr("exon:t3:1", "t3", ";gene_id=g3;transcript_id=t3;exon_number=1")
      )
    ),
    gff
  )
  regions <- read_sequence_regions(gff)
  testthat::expect_equal(regions$seqid, c("chr1", "chr2"))
  testthat::expect_equal(regions$end_pos, c(1000L, 500L))
  annot_dir <- suppressMessages(prepare_isoformic_annotation(
    gff,
    output_path = withr::local_tempdir()
  ))
  seq_regions <- read_parquet(fs::path(annot_dir, "annot_seq_regions.parquet"))
  testthat::expect_equal(seq_regions$seqid, c("chr1", "chr2", "chr3"))
  testthat::expect_equal(seq_regions$end_pos, c(1000L, 500L, 70L))
})
