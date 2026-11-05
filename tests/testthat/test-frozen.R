testthat::test_that("the frozen dataset registry is consistent", {
  reg <- frozen_datasets()
  testthat::expect_setequal(
    reg$name,
    c(
      "appris",
      "appris_domains",
      "riboseq_orfs",
      "nmd_exons",
      "polyasite",
      "gtex_transcript_tpm"
    )
  )
  testthat::expect_true(all(grepl("\\.parquet$", reg$file)))
  testthat::expect_true(all(nzchar(reg$version)))
  for (b in reg$builder) {
    testthat::expect_true(
      is.function(get(b, envir = asNamespace("isoformic"))),
      info = b
    )
  }
  testthat::expect_error(
    fetch_frozen_data("nope"),
    class = "isoformic_unknown_frozen_dataset"
  )
  testthat::expect_error(
    read_frozen_data("nope"),
    class = "isoformic_unknown_frozen_dataset"
  )
})

testthat::test_that("fetch_frozen_data returns cached files and verifies their checksum", {
  dir <- withr::local_tempdir()
  entry <- frozen_dataset_entry("nmd_exons")
  path <- fs::path(dir, entry$file)
  write_parquet(data.frame(seqid = "chr1", start_pos = 1L), path)
  testthat::expect_identical(
    fetch_frozen_data("nmd_exons", output_path = dir, verify = FALSE),
    as.character(path)
  )
  if (!is.na(entry$md5)) {
    testthat::expect_warning(
      fetch_frozen_data("nmd_exons", output_path = dir),
      "MD5"
    )
  }
  res <- suppressWarnings(
    read_frozen_data("nmd_exons", output_path = dir, build = "never")
  )
  testthat::expect_equal(nrow(res), 1L)
})

testthat::test_that("build_frozen_appris converts the APPRIS score table", {
  dir <- withr::local_tempdir()
  src <- fs::path(dir, "sources")
  fs::dir_create(src)
  header <- paste(
    "Ensembl Gene ID",
    "Gene name (HGNC)",
    "Transcript ID",
    "Translation ID",
    "Is translated?",
    "Transcript type",
    "Not found tag",
    "CCDS ID",
    "Transcript support level",
    "Protein length",
    "Functional residues (firestar)",
    "Structure score (Matador)",
    "Conservation (CORSAIR)",
    "Domain Score (SPADE)",
    "Trans-membrane helices (THUMP)",
    "Signal sequence (CRASH)",
    "Trifid Score",
    "Peptides",
    "APPRIS score",
    "APPRIS Annotation",
    "MANE",
    sep = "\t"
  )
  rows <- c(
    paste(
      "ENSG1",
      "GENEA",
      "ENST00000000001.3",
      "ENSP1",
      "TRANSLATION",
      "protein_coding",
      "-",
      "CCDS1",
      "1",
      "100",
      "3",
      "10.5",
      "2",
      "120.0",
      "0",
      "-,-",
      "0.5",
      "-",
      "8.1",
      "PRINCIPAL:1",
      "MANE_Select",
      sep = "\t"
    ),
    paste(
      "ENSG1",
      "GENEA",
      "ENST00000000002.1",
      "-",
      "NO_TRANSLATION",
      "retained_intron",
      "",
      "",
      "",
      "",
      "",
      "",
      "",
      "",
      "",
      "",
      "",
      "",
      "",
      "",
      "",
      sep = "\t"
    )
  )
  writeLines(c(header, rows), fs::path(src, "appris_data.appris.txt"))
  path <- build_frozen_appris(output_path = dir, sources_dir = src)
  res <- read_parquet(path)
  testthat::expect_equal(nrow(res), 2L)
  testthat::expect_identical(
    res$transcript_id,
    c("ENST00000000001", "ENST00000000002")
  )
  testthat::expect_identical(res$appris_principal, c(TRUE, FALSE))
  testthat::expect_equal(res$spade_domain_score[1], 120)
  testthat::expect_identical(res$mane[1], "MANE_Select")
})

testthat::test_that("build_frozen_appris_domains parses the SPADE GTF", {
  dir <- withr::local_tempdir()
  src <- fs::path(dir, "sources")
  fs::dir_create(src)
  gtf <- fs::path(src, "appris_method.spade.gtf.gz")
  con <- gzfile(gtf, "w")
  writeLines(
    c(
      paste(
        "1",
        "SPADE",
        "domain",
        "2558396",
        "2559875",
        "41.6",
        "+",
        ".",
        paste0(
          'gene_id "ENSG00000000001.2"; ',
          'note "hmm_acc:PF00020.17,hmm_name:TNFR_c6,evalue:4.3e-11,pep_start:78,pep_end:119"; ',
          'transcript_id "ENST00000000001.3"'
        ),
        sep = "\t"
      ),
      paste(
        "chrX",
        "SPADE",
        "domain",
        "100",
        "200",
        "5",
        "-",
        ".",
        'note "hmm_acc:PF2,hmm_name:Dom2,evalue:1e-3,pep_start:1,pep_end:20"; gene_id "ENSG2"; transcript_id "ENST2"',
        sep = "\t"
      )
    ),
    con
  )
  close(con)
  path <- build_frozen_appris_domains(output_path = dir, sources_dir = src)
  res <- read_parquet(path)
  testthat::expect_equal(nrow(res), 2L)
  testthat::expect_identical(res$seqid, c("chr1", "chrX"))
  testthat::expect_identical(res$transcript_id, c("ENST00000000001", "ENST2"))
  testthat::expect_identical(res$pfam_name, c("TNFR_c6", "Dom2"))
  testthat::expect_equal(res$pep_end, c(119L, 20L))
  testthat::expect_equal(res$evalue[1], 4.3e-11)
})

testthat::test_that("build_frozen_riboseq_orfs merges phase 1 and phase 2 BED files", {
  dir <- withr::local_tempdir()
  src <- fs::path(dir, "sources")
  fs::dir_create(src)
  phase1 <- paste(
    "chr1",
    "826870",
    "829027",
    "c1norep1",
    "0",
    "+",
    "826870",
    "829027",
    "0",
    "2",
    "53,25,",
    "0,2132,",
    "c1norep1",
    "none",
    "none",
    "0,2,",
    "lncRNA",
    "ENSG00000228794.1",
    "LINC01128",
    "lncRNA",
    "lncRNA",
    "MKKPLP",
    "ENST00000623808.2,ENST00000416570.1",
    "ENSG00000228794",
    "no",
    "PMID:32139545",
    sep = "\t"
  )
  phase2 <- paste(
    "chr1",
    "826699",
    "826726",
    "c1riboseqorf662:LINC00115:lncRNA",
    "0",
    "-",
    "826699",
    "826726",
    "128,128,128",
    "1",
    "27",
    "0",
    sep = "\t"
  )
  writeLines(phase1, fs::path(src, "Ribo-seq_ORFs.bed"))
  writeLines(phase2, fs::path(src, "Ribo-seq_ORFs.comprehensive.bed"))
  path <- build_frozen_riboseq_orfs(output_path = dir, sources_dir = src)
  res <- read_parquet(path)
  testthat::expect_equal(nrow(res), 2L)
  testthat::expect_equal(res$phase, c(1L, 2L))
  testthat::expect_identical(res$block_starts[1], "826871;829003")
  testthat::expect_identical(res$block_ends[1], "826923;829027")
  testthat::expect_identical(
    res$transcript_ids[1],
    "ENST00000623808;ENST00000416570"
  )
  testthat::expect_identical(res$gene_id[1], "ENSG00000228794")
  testthat::expect_identical(res$orf_id[2], "c1riboseqorf662")
  testthat::expect_identical(res$gene_name[2], "LINC00115")
  testthat::expect_equal(res$start_pos[2], 826700L)
  testthat::expect_identical(res$replicated, c(FALSE, NA))
})

testthat::test_that("build_frozen_polyasite parses the atlas BED", {
  dir <- withr::local_tempdir()
  src <- fs::path(dir, "sources")
  fs::dir_create(src)
  bed <- fs::path(src, "atlas.clusters.2.0.GRCh38.96.bed.gz")
  con <- gzfile(bed, "w")
  writeLines(
    c(
      paste(
        "1",
        "16442",
        "16452",
        "1:16450:-",
        "0.2690",
        "-",
        "0.15",
        "4",
        "0.2690",
        "IN",
        "CATAAA@-18@16467",
        sep = "\t"
      ),
      paste(
        "1",
        "29535",
        "29539",
        "1:29536:-",
        "0.0042",
        "-",
        "0.00",
        "1",
        "0.0042",
        "AI",
        "NA",
        sep = "\t"
      )
    ),
    con
  )
  close(con)
  path <- build_frozen_polyasite(output_path = dir, sources_dir = src)
  res <- read_parquet(path)
  testthat::expect_identical(res$seqid, c("chr1", "chr1"))
  testthat::expect_equal(res$cleavage_pos, c(16450L, 29536L))
  testthat::expect_equal(res$start_pos, c(16443L, 29536L))
  testthat::expect_identical(res$site_class, c("IN", "AI"))
  testthat::expect_identical(res$signals, c("CATAAA@-18@16467", NA))
})

testthat::test_that("build_frozen_gtex computes median TPM per tissue", {
  dir <- withr::local_tempdir()
  src <- fs::path(dir, "sources")
  fs::dir_create(src)
  gct <- fs::path(
    src,
    "GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct.gz"
  )
  con <- gzfile(gct, "w")
  writeLines(
    c(
      "#1.2",
      "2\t4",
      paste("transcript_id", "gene_id", "S1", "S2", "S3", "S4", sep = "\t"),
      paste(
        "ENST00000000001.1",
        "ENSG00000000001.2",
        "1",
        "3",
        "10",
        "20",
        sep = "\t"
      ),
      paste(
        "ENST00000000002.1",
        "ENSG00000000001.2",
        "0",
        "0",
        "5",
        "0",
        sep = "\t"
      )
    ),
    con
  )
  close(con)
  writeLines(
    c(
      paste("SAMPID", "SMTSD", "SMAFRZE", sep = "\t"),
      paste("S1", "Liver", "RNASEQ", sep = "\t"),
      paste("S2", "Liver", "RNASEQ", sep = "\t"),
      paste("S3", "Lung", "RNASEQ", sep = "\t"),
      paste("S4", "Lung", "EXCLUDE", sep = "\t")
    ),
    fs::path(src, "GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt")
  )
  path <- build_frozen_gtex(output_path = dir, sources_dir = src)
  res <- read_parquet(path)
  testthat::expect_equal(nrow(res), 4L)
  liver <- res[res$tissue == "Liver" & res$transcript_id == "ENST00000000001", ]
  testthat::expect_equal(liver$median_tpm, 2)
  testthat::expect_equal(liver$n_samples, 2L)
  lung <- res[res$tissue == "Lung" & res$transcript_id == "ENST00000000002", ]
  testthat::expect_equal(lung$median_tpm, 5)
  testthat::expect_equal(lung$n_samples, 1L)
  testthat::expect_identical(unique(res$gene_id), "ENSG00000000001")
})

testthat::test_that("build_frozen_nmd_exons reads the supplementary table", {
  testthat::skip_if_not_installed("readxl")
  testthat::skip_if_not_installed("writexl")
  dir <- withr::local_tempdir()
  src <- fs::path(dir, "sources")
  fs::dir_create(src)
  tbl <- data.frame(
    GeneID = "MSTRG.1",
    chr = "chr11",
    strand = "-",
    exonStart_0base = 6473083,
    exonEnd = 6473248,
    upstreamES = 6465564,
    upstreamEE = 6465732,
    downstreamES = 6473790,
    downstreamEE = 6473923,
    NMD_isoforms = "ENST00000359518.5",
    non_NMD_isoforms = "ENST00000345851.2",
    gene_id = "ENSG00000110171.1",
    gene_name = "TRIM3"
  )
  writexl::write_xlsx(tbl, fs::path(src, "13059_2021_2439_MOESM2_ESM.xlsx"))
  path <- build_frozen_nmd_exons(output_path = dir, sources_dir = src)
  res <- read_parquet(path)
  testthat::expect_equal(res$start_pos, 6473084L)
  testthat::expect_identical(res$nmd_transcript_ids, "ENST00000359518")
  testthat::expect_identical(res$gene_id, "ENSG00000110171")
})

testthat::test_that("frozen_mirrors lists user, GitHub and Zenodo bases in order", {
  m <- frozen_mirrors(zenodo_record = NA)
  testthat::expect_match(
    m[length(m)],
    "github.com/luciorq/isoformic/releases/download/frozen-data/"
  )
  withr::local_options(isoformic.frozen_mirrors = "https://example.org/mirror")
  m <- frozen_mirrors(
    release_tag = "tag-x",
    zenodo_record = "https://zenodo.org/records/123456"
  )
  testthat::expect_identical(m[1], "https://example.org/mirror")
  testthat::expect_match(m[2], "tag-x/$")
  testthat::expect_identical(m[3], "https://zenodo.org/records/123456/files/")
  testthat::expect_identical(
    frozen_mirror_url(m[3], "a.parquet"),
    "https://zenodo.org/records/123456/files/a.parquet?download=1"
  )
  testthat::expect_identical(
    frozen_mirror_url("https://example.org/mirror", "a.parquet"),
    "https://example.org/mirror/a.parquet"
  )
  testthat::expect_length(frozen_mirrors(zenodo_record = "7"), 3L)
})

testthat::test_that("bundle_frozen_data assembles the snapshot with checksums and a README", {
  cache <- withr::local_tempdir()
  for (name in frozen_datasets()$name) {
    entry <- frozen_dataset_entry(name)
    write_parquet(data.frame(id = name), fs::path(cache, entry$file))
  }
  out <- fs::path(withr::local_tempdir(), "bundle")
  testthat::expect_message(
    bundle_frozen_data(
      out,
      output_path = cache,
      build = "never",
      verify = FALSE
    ),
    "Bundled"
  )
  files <- fs::path_file(fs::dir_ls(out))
  testthat::expect_true(all(frozen_datasets()$file %in% files))
  testthat::expect_true(all(
    c(
      "MD5SUMS",
      "README.md",
      "hexamer_logodds.parquet",
      "coding_model_coefficients.parquet"
    ) %in%
      files
  ))
  md5 <- read.table(fs::path(out, "MD5SUMS"), col.names = c("md5", "file"))
  testthat::expect_equal(nrow(md5), nrow(frozen_datasets()) + 2L)
  testthat::expect_identical(
    unname(tools::md5sum(as.character(fs::path(out, md5$file[1])))),
    md5$md5[1]
  )
  readme <- readLines(fs::path(out, "README.md"))
  testthat::expect_true(any(grepl("e113v49", readme, fixed = TRUE)))
})

testthat::test_that("a download that fails its checksum is discarded and the next mirror is used", {
  root <- withr::local_tempdir()
  bad <- fs::path(root, "bad")
  good <- fs::path(root, "good")
  fs::dir_create(c(bad, good))
  write_parquet(data.frame(a = 1L), fs::path(good, "toy.parquet"))
  write_parquet(data.frame(a = 2L), fs::path(bad, "toy.parquet"))
  good_md5 <- unname(tools::md5sum(as.character(fs::path(good, "toy.parquet"))))
  testthat::local_mocked_bindings(
    frozen_dataset_entry = function(name) {
      list(file = "toy.parquet", md5 = good_md5, version = "test")
    }
  )
  out <- fs::path(root, "out")
  mirrors <- paste0("file://", c(bad, good))
  res <- suppressMessages(
    fetch_frozen_data(
      "toy",
      output_path = out,
      build = "never",
      mirrors = mirrors
    )
  )
  testthat::expect_identical(
    unname(tools::md5sum(res)),
    good_md5
  )
  # Only a bad mirror: nothing is cached and the call fails
  fs::file_delete(res)
  testthat::expect_error(
    suppressMessages(
      fetch_frozen_data(
        "toy",
        output_path = out,
        build = "never",
        mirrors = paste0("file://", bad)
      )
    ),
    class = "isoformic_download_reference_error"
  )
  testthat::expect_false(fs::file_exists(fs::path(out, "toy.parquet")))
})
