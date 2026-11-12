# All network access is mocked: `list_remote_dir()` returns canned listings
# and `utils::download.file()` writes a small file.

gencode_human_listing <- c(
  "GRCh38.p14.genome.fa.gz",
  "GRCh38.primary_assembly.genome.fa.gz",
  "gencode.v50.annotation.gff3.gz",
  "gencode.v50.annotation.gtf.gz",
  "gencode.v50.basic.annotation.gff3.gz",
  "gencode.v50.chr_patch_hapl_scaff.annotation.gff3.gz",
  "gencode.v50.chr_patch_hapl_scaff.annotation.gtf.gz",
  "gencode.v50.primary_assembly.annotation.gff3.gz",
  "gencode.v50.pc_transcripts.fa.gz",
  "gencode.v50.transcripts.fa.gz"
)

mock_listing <- function(listings) {
  function(url, ...) {
    for (nm in names(listings)) {
      if (identical(url, nm)) {
        return(listings[[nm]])
      }
    }
    NULL
  }
}

mock_download_ok <- function(url, destfile, ...) {
  writeLines("content", destfile)
  0L
}

local_download_mocks <- function(
  listings,
  download = mock_download_ok,
  envir = parent.frame()
) {
  testthat::local_mocked_bindings(
    has_internet = function(...) TRUE,
    .package = "curl",
    .env = envir
  )
  testthat::local_mocked_bindings(
    list_remote_dir = mock_listing(listings),
    .env = envir
  )
  testthat::local_mocked_bindings(
    download.file = download,
    .package = "utils",
    .env = envir
  )
}


testthat::test_that("download_reference aborts without internet", {
  testthat::local_mocked_bindings(
    has_internet = function(...) FALSE,
    .package = "curl"
  )
  testthat::expect_error(
    download_reference(output_path = withr::local_tempdir()),
    class = "isoformic_download_reference_no_internet_error"
  )
})

testthat::test_that("download_reference discovers GENCODE files from the release listing", {
  local_download_mocks(rlang::list2(
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/" = gencode_human_listing
  ))
  out <- withr::local_tempdir()
  testthat::expect_message(
    res <- download_reference(version = 50, output_path = out),
    "successfully downloaded"
  )
  testthat::expect_equal(
    fs::path_file(res),
    "gencode.v50.chr_patch_hapl_scaff.annotation.gff3.gz"
  )
  testthat::expect_true(fs::file_exists(res))
  testthat::expect_false(fs::file_exists(paste0(res, ".part")))

  testthat::expect_message(
    download_reference(version = "50", output_path = out),
    "already exists"
  )

  res_chr <- download_reference(
    version = "50",
    file_type = "gtf",
    region = "chr",
    output_path = out
  )
  testthat::expect_equal(
    fs::path_file(res_chr),
    "gencode.v50.annotation.gtf.gz"
  )
  res_pri <- download_reference(
    version = "50",
    file_type = "gff",
    region = "primary",
    output_path = out
  )
  testthat::expect_equal(
    fs::path_file(res_pri),
    "gencode.v50.primary_assembly.annotation.gff3.gz"
  )
  res_fa <- download_reference(
    version = "50",
    file_type = "fasta",
    output_path = out
  )
  testthat::expect_equal(fs::path_file(res_fa), "gencode.v50.transcripts.fa.gz")
  res_genome <- download_reference(
    version = "50",
    file_type = "genome_fasta",
    output_path = out
  )
  testthat::expect_equal(fs::path_file(res_genome), "GRCh38.p14.genome.fa.gz")
  res_genome_pri <- download_reference(
    version = "50",
    file_type = "genome_fasta",
    region = "primary",
    output_path = out
  )
  testthat::expect_equal(
    fs::path_file(res_genome_pri),
    "GRCh38.primary_assembly.genome.fa.gz"
  )
})

testthat::test_that("future releases and genome patches work without hardcoded names", {
  listing <- c(
    "GRCh38.p16.genome.fa.gz",
    "gencode.v57.chr_patch_hapl_scaff.annotation.gff3.gz",
    "gencode.v57.transcripts.fa.gz"
  )
  local_download_mocks(rlang::list2(
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_57/" = listing
  ))
  out <- withr::local_tempdir()
  res <- download_reference(
    version = "57",
    file_type = "genome_fasta",
    output_path = out
  )
  testthat::expect_equal(fs::path_file(res), "GRCh38.p16.genome.fa.gz")
  res_gff <- download_reference(version = "v57", output_path = out)
  testthat::expect_match(
    res_gff,
    "gencode.v57.chr_patch_hapl_scaff.annotation.gff3.gz$"
  )
})

testthat::test_that("mouse versions are normalised to M releases", {
  listing <- c(
    "GRCm39.genome.fa.gz",
    "gencode.vM39.chr_patch_hapl_scaff.annotation.gff3.gz"
  )
  local_download_mocks(rlang::list2(
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M39/" = listing
  ))
  out <- withr::local_tempdir()
  res <- download_reference(
    version = "39",
    organism = "mouse",
    output_path = out
  )
  testthat::expect_match(
    res,
    "gencode.vM39.chr_patch_hapl_scaff.annotation.gff3.gz$"
  )
  res2 <- download_reference(
    version = "M39",
    organism = "mouse",
    file_type = "genome_fasta",
    output_path = out
  )
  testthat::expect_match(res2, "GRCm39.genome.fa.gz$")
  testthat::expect_error(
    download_reference(version = "M39", organism = "human", output_path = out),
    class = "isoformic_download_reference_unsupported_version_error"
  )
  testthat::expect_error(
    download_reference(version = "latest", output_path = out),
    class = "isoformic_download_reference_unsupported_version_error"
  )
})

testthat::test_that("missing releases and files give specific errors", {
  local_download_mocks(rlang::list2(
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/" = c(
      "release_47",
      "release_48",
      "release_50",
      "release_50",
      "latest_release"
    ),
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/" = gencode_human_listing
  ))
  out <- withr::local_tempdir()
  err <- testthat::expect_error(
    download_reference(version = "999", output_path = out),
    class = "isoformic_download_reference_release_not_found_error"
  )
  testthat::expect_match(conditionMessage(err), "release_50")
  err2 <- testthat::expect_error(
    download_reference(
      version = "50",
      file_type = "gtf",
      region = "primary",
      output_path = out
    ),
    class = "isoformic_download_reference_file_not_found_error"
  )
  testthat::expect_match(
    conditionMessage(err2),
    "gencode.v50.annotation.gtf.gz"
  )
})

testthat::test_that("failed or empty downloads are reported and leave no partial file", {
  local_download_mocks(
    rlang::list2(
      "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/" = gencode_human_listing
    ),
    download = function(url, destfile, ...) {
      stop("connection reset", call. = FALSE)
    }
  )
  out <- withr::local_tempdir()
  err <- testthat::expect_error(
    download_reference(version = "50", output_path = out),
    class = "isoformic_download_reference_error"
  )
  testthat::expect_match(conditionMessage(err), "connection reset")
  testthat::expect_length(fs::dir_ls(out), 0L)

  local_download_mocks(
    rlang::list2(
      "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/" = gencode_human_listing
    ),
    download = function(url, destfile, ...) {
      file.create(destfile)
      0L
    }
  )
  testthat::expect_error(
    download_reference(version = "50", output_path = out),
    class = "isoformic_download_reference_error"
  )
  testthat::expect_length(fs::dir_ls(out), 0L)
})

testthat::test_that("download_reference resolves T2T-CHM13 (NCBI RefSeq) files", {
  listing <- paste0(
    "GCF_009914755.1_T2T-CHM13v2.0",
    c(
      "_genomic.gff.gz",
      "_genomic.gtf.gz",
      "_rna.fna.gz",
      "_genomic.fna.gz",
      "_protein.faa.gz"
    )
  )
  captured <- new.env(parent = emptyenv())
  captured$value <- NULL
  local_download_mocks(
    rlang::list2(
      "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/" = "GCF_009914755.1_T2T-CHM13v2.0",
      "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0/" = listing
    ),
    download = function(url, destfile, ...) {
      captured$value <- url
      writeLines("x", destfile)
      0L
    }
  )
  out <- withr::local_tempdir()
  res <- download_reference(
    reference = "t2t",
    file_type = "gff",
    output_path = out
  )
  testthat::expect_match(
    captured$value,
    "GCF_009914755.1_T2T-CHM13v2.0_genomic.gff.gz$"
  )
  testthat::expect_match(res, "_genomic.gff.gz$")
  download_reference(reference = "t2t", file_type = "fasta", output_path = out)
  testthat::expect_match(captured$value, "_rna.fna.gz$")
  download_reference(
    reference = "t2t",
    file_type = "genome_fasta",
    output_path = out
  )
  testthat::expect_match(captured$value, "_genomic.fna.gz$")
  testthat::expect_error(
    download_reference(
      reference = "t2t",
      organism = "mouse",
      output_path = out
    ),
    class = "isoformic_download_reference_unsupported_version_error"
  )
  testthat::expect_error(
    download_reference(reference = "t2t", version = "1.1", output_path = out),
    class = "isoformic_download_reference_unsupported_version_error"
  )
})

testthat::test_that("list_remote_dir parses HTML indexes and returns NULL on failure", {
  html <- c(
    "<html><body><a href=\"../\">..</a>",
    "<a href=\"release_50/\">release_50/</a>",
    "<a href=\"gencode.v50.transcripts.fa.gz\">gencode.v50.transcripts.fa.gz</a>",
    "<a href=\"?C=N;O=D\">Name</a></body></html>"
  )
  tmp <- withr::local_tempfile(lines = html, fileext = ".html")
  res <- list_remote_dir(paste0("file://", tmp))
  testthat::expect_equal(res, c("release_50", "gencode.v50.transcripts.fa.gz"))
  testthat::expect_null(list_remote_dir("file:///no/such/dir/index.html"))
})

testthat::test_that("download_reference resolves Ensembl, NCBI accession and MANE files", {
  ens_dir <- "https://ftp.ensembl.org/pub/current_gff3/homo_sapiens/"
  ens_rel <- "https://ftp.ensembl.org/pub/release-116/gtf/homo_sapiens/"
  plants_cdna <- "https://ftp.ebi.ac.uk/ensemblgenomes/pub/plants/current/fasta/arabidopsis_thaliana/cdna/"
  ncbi_parent <- "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/"
  ncbi_dir <- paste0(ncbi_parent, "GCF_000001405.40_GRCh38.p14/")
  mane_cur <- "https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/"
  captured <- new.env(parent = emptyenv())
  captured$value <- NULL
  local_download_mocks(
    rlang::list2(
      !!ens_dir := c(
        "Homo_sapiens.GRCh38.116.abinitio.gff3.gz",
        "Homo_sapiens.GRCh38.116.chr.gff3.gz",
        "Homo_sapiens.GRCh38.116.chr_patch_hapl_scaff.gff3.gz",
        "Homo_sapiens.GRCh38.116.chromosome.1.gff3.gz",
        "Homo_sapiens.GRCh38.116.gff3.gz"
      ),
      !!ens_rel := c(
        "Homo_sapiens.GRCh38.116.chr.gtf.gz",
        "Homo_sapiens.GRCh38.116.gtf.gz"
      ),
      !!plants_cdna := c(
        "Arabidopsis_thaliana.TAIR10.cdna.abinitio.fa.gz",
        "Arabidopsis_thaliana.TAIR10.cdna.all.fa.gz"
      ),
      !!ncbi_parent := c(
        "GCF_000001405.39_GRCh38.p13",
        "GCF_000001405.40_GRCh38.p14"
      ),
      !!ncbi_dir := paste0(
        "GCF_000001405.40_GRCh38.p14",
        c("_genomic.gff.gz", "_rna.fna.gz", "_genomic.fna.gz")
      ),
      !!mane_cur := c(
        "MANE.GRCh38.v1.5.ensembl_genomic.gff.gz",
        "MANE.GRCh38.v1.5.refseq_genomic.gff.gz",
        "MANE.GRCh38.v1.5.ensembl_rna.fna.gz"
      )
    ),
    download = function(url, destfile, ...) {
      captured$value <- url
      writeLines("x", destfile)
      0L
    }
  )
  out <- withr::local_tempdir()
  download_reference(reference = "ensembl", output_path = out)
  testthat::expect_match(
    captured$value,
    "Homo_sapiens.GRCh38.116.chr_patch_hapl_scaff.gff3.gz$"
  )
  download_reference(
    reference = "ensembl",
    region = "primary",
    output_path = out
  )
  testthat::expect_match(captured$value, "/Homo_sapiens.GRCh38.116.gff3.gz$")
  download_reference(reference = "ensembl", region = "chr", output_path = out)
  testthat::expect_match(captured$value, "116.chr.gff3.gz$")
  download_reference(
    reference = "ensembl",
    version = 116,
    file_type = "gtf",
    region = "primary",
    output_path = out
  )
  testthat::expect_match(
    captured$value,
    "release-116/gtf/homo_sapiens/Homo_sapiens.GRCh38.116.gtf.gz$"
  )
  download_reference(
    reference = "ensembl",
    division = "plants",
    species = "arabidopsis_thaliana",
    file_type = "fasta",
    output_path = out
  )
  testthat::expect_match(
    captured$value,
    "Arabidopsis_thaliana.TAIR10.cdna.all.fa.gz$"
  )
  download_reference(
    reference = "ncbi",
    accession = "GCF_000001405.40",
    output_path = out
  )
  testthat::expect_match(
    captured$value,
    "GCF_000001405.40_GRCh38.p14_genomic.gff.gz$"
  )
  download_reference(
    reference = "ncbi",
    accession = "gcf_000001405.40",
    file_type = "fasta",
    output_path = out
  )
  testthat::expect_match(captured$value, "_rna.fna.gz$")
  download_reference(reference = "mane", output_path = out)
  testthat::expect_match(
    captured$value,
    "MANE.GRCh38.v1.5.ensembl_genomic.gff.gz$"
  )
  download_reference(
    reference = "mane",
    mane_source = "refseq",
    output_path = out
  )
  testthat::expect_match(captured$value, "refseq_genomic.gff.gz$")
  testthat::expect_error(
    download_reference(
      reference = "mane",
      file_type = "genome_fasta",
      output_path = out
    ),
    class = "isoformic_download_reference_file_not_found_error"
  )
  testthat::expect_error(
    download_reference(
      reference = "ncbi",
      accession = "GCF_123",
      output_path = out
    ),
    class = "isoformic_download_reference_unsupported_version_error"
  )
  testthat::expect_error(
    download_reference(
      reference = "ncbi",
      accession = "GCF_000001405.99",
      output_path = out
    ),
    class = "isoformic_download_reference_release_not_found_error"
  )
  testthat::expect_error(
    download_reference(
      reference = "ensembl",
      species = "nope_nope",
      output_path = out
    ),
    class = "isoformic_download_reference_release_not_found_error"
  )
  testthat::expect_error(
    download_reference(
      version = "50",
      file_type = "ncrna_fasta",
      output_path = out
    ),
    class = "isoformic_download_reference_file_not_found_error"
  )
})

testthat::test_that("download_reference_url and download_hprc_annotation cache arbitrary files", {
  captured <- new.env(parent = emptyenv())
  captured$value <- character()
  local_download_mocks(
    rlang::list2(),
    download = function(url, destfile, ...) {
      captured$value <- c(captured$value, url)
      writeLines("x", destfile)
      0L
    }
  )
  out <- withr::local_tempdir()
  res <- download_reference_url(
    "https://example.org/annot/chm13v2.0_RefSeq_Liftoff_v5.2.gff3.gz?x=1",
    output_path = out
  )
  testthat::expect_equal(
    fs::path_file(res),
    "chm13v2.0_RefSeq_Liftoff_v5.2.gff3.gz"
  )
  testthat::expect_error(
    download_reference_url("not a url", output_path = out),
    class = "isoformic_download_reference_error"
  )

  files <- download_hprc_annotation("HG00438", output_path = out)
  testthat::expect_named(files, c("hap1", "hap2"))
  testthat::expect_match(captured$value[2], "HG00438.1.chm13.gff3.gz$")
  testthat::expect_match(captured$value[3], "HG00438.2.chm13.gff3.gz$")
  one <- download_hprc_annotation(
    "HG002",
    haplotype = 2,
    assembly = "hg38",
    output_path = out
  )
  testthat::expect_match(captured$value[4], "HG002.2.hg38.gff3.gz$")
  testthat::expect_named(one, "hap2")
  testthat::expect_error(
    download_hprc_annotation("HG002", haplotype = 3, output_path = out),
    class = "isoformic_invalid_haplotypes"
  )
})

testthat::test_that("download_reference resolves FANTOM CAT files and companion tables", {
  lv3 <- "https://fantom.gsc.riken.jp/5/suppl/Hon_et_al_2016/data/assembly/lv3_robust/"
  captured <- new.env(parent = emptyenv())
  captured$value <- character()
  local_download_mocks(
    rlang::list2(
      !!lv3 := c("FANTOM_CAT.lv3_robust.gtf.gz", "FANTOM_CAT.lv3_robust.bed.gz")
    ),
    download = function(url, destfile, ...) {
      captured$value <- c(captured$value, url)
      writeLines("x", destfile)
      0L
    }
  )
  out <- withr::local_tempdir()
  res <- download_reference(
    reference = "fantom",
    file_type = "gtf",
    output_path = out
  )
  testthat::expect_match(res, "FANTOM_CAT.lv3_robust.gtf.gz$")
  testthat::expect_error(
    download_reference(
      reference = "fantom",
      file_type = "gff",
      output_path = out
    ),
    class = "isoformic_download_reference_file_not_found_error"
  )
  testthat::expect_error(
    download_reference(
      reference = "fantom",
      version = "lv9",
      file_type = "gtf",
      output_path = out
    ),
    class = "isoformic_download_reference_unsupported_version_error"
  )
  tables <- download_fantom_cat_tables("lv3_robust", output_path = out)
  testthat::expect_named(tables, c("id_mapping", "gene_info", "cage_clusters"))
  testthat::expect_match(captured$value[2], "info_table.ID_mapping.tsv.gz$")
  testthat::expect_match(captured$value[4], "CAGE_cluster.bed.gz$")
})
