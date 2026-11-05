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
