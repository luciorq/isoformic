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
    suppressMessages(download_reference(output_path = withr::local_tempdir())),
    class = "isoformic_download_reference_no_internet_error"
  )
})

testthat::test_that("download_reference discovers GENCODE files from the release listing", {
  local_download_mocks(rlang::list2(
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/" = gencode_human_listing
  ))
  out <- withr::local_tempdir()
  # The "Downloading" message comes first and is not under test
  suppressMessages(testthat::expect_message(
    res <- download_reference(version = 50, output_path = out),
    "successfully downloaded"
  ))
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

  res_chr <- suppressMessages(download_reference(
    version = "50",
    file_type = "gtf",
    region = "chr",
    output_path = out
  ))
  testthat::expect_equal(
    fs::path_file(res_chr),
    "gencode.v50.annotation.gtf.gz"
  )
  res_pri <- suppressMessages(download_reference(
    version = "50",
    file_type = "gff",
    region = "primary",
    output_path = out
  ))
  testthat::expect_equal(
    fs::path_file(res_pri),
    "gencode.v50.primary_assembly.annotation.gff3.gz"
  )
  res_fa <- suppressMessages(download_reference(
    version = "50",
    file_type = "fasta",
    output_path = out
  ))
  testthat::expect_equal(fs::path_file(res_fa), "gencode.v50.transcripts.fa.gz")
  res_genome <- suppressMessages(download_reference(
    version = "50",
    file_type = "genome_fasta",
    output_path = out
  ))
  testthat::expect_equal(fs::path_file(res_genome), "GRCh38.p14.genome.fa.gz")
  res_genome_pri <- suppressMessages(download_reference(
    version = "50",
    file_type = "genome_fasta",
    region = "primary",
    output_path = out
  ))
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
  res <- suppressMessages(download_reference(
    version = "57",
    file_type = "genome_fasta",
    output_path = out
  ))
  testthat::expect_equal(fs::path_file(res), "GRCh38.p16.genome.fa.gz")
  res_gff <- suppressMessages(download_reference(
    version = "v57",
    output_path = out
  ))
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
  res <- suppressMessages(download_reference(
    version = "39",
    organism = "mouse",
    output_path = out
  ))
  testthat::expect_match(
    res,
    "gencode.vM39.chr_patch_hapl_scaff.annotation.gff3.gz$"
  )
  res2 <- suppressMessages(download_reference(
    version = "M39",
    organism = "mouse",
    file_type = "genome_fasta",
    output_path = out
  ))
  testthat::expect_match(res2, "GRCm39.genome.fa.gz$")
  testthat::expect_error(
    suppressMessages(download_reference(
      version = "M39",
      organism = "human",
      output_path = out
    )),
    class = "isoformic_download_reference_unsupported_version_error"
  )
  testthat::expect_error(
    suppressMessages(download_reference(version = "latest", output_path = out)),
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
    suppressMessages(download_reference(version = "999", output_path = out)),
    class = "isoformic_download_reference_release_not_found_error"
  )
  testthat::expect_match(conditionMessage(err), "release_50")
  err2 <- testthat::expect_error(
    suppressMessages(download_reference(
      version = "50",
      file_type = "gtf",
      region = "primary",
      output_path = out
    )),
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
    suppressMessages(download_reference(version = "50", output_path = out)),
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
    suppressMessages(download_reference(version = "50", output_path = out)),
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
  res <- suppressMessages(download_reference(
    reference = "t2t",
    file_type = "gff",
    output_path = out
  ))
  testthat::expect_match(
    captured$value,
    "GCF_009914755.1_T2T-CHM13v2.0_genomic.gff.gz$"
  )
  testthat::expect_match(res, "_genomic.gff.gz$")
  suppressMessages(download_reference(
    reference = "t2t",
    file_type = "fasta",
    output_path = out
  ))
  testthat::expect_match(captured$value, "_rna.fna.gz$")
  suppressMessages(download_reference(
    reference = "t2t",
    file_type = "genome_fasta",
    output_path = out
  ))
  testthat::expect_match(captured$value, "_genomic.fna.gz$")
  testthat::expect_error(
    suppressMessages(download_reference(
      reference = "t2t",
      organism = "mouse",
      output_path = out
    )),
    class = "isoformic_download_reference_unsupported_version_error"
  )
  testthat::expect_error(
    suppressMessages(download_reference(
      reference = "t2t",
      version = "1.1",
      output_path = out
    )),
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
  suppressMessages(download_reference(reference = "ensembl", output_path = out))
  testthat::expect_match(
    captured$value,
    "Homo_sapiens.GRCh38.116.chr_patch_hapl_scaff.gff3.gz$"
  )
  suppressMessages(download_reference(
    reference = "ensembl",
    region = "primary",
    output_path = out
  ))
  testthat::expect_match(captured$value, "/Homo_sapiens.GRCh38.116.gff3.gz$")
  suppressMessages(download_reference(
    reference = "ensembl",
    region = "chr",
    output_path = out
  ))
  testthat::expect_match(captured$value, "116.chr.gff3.gz$")
  suppressMessages(download_reference(
    reference = "ensembl",
    version = 116,
    file_type = "gtf",
    region = "primary",
    output_path = out
  ))
  testthat::expect_match(
    captured$value,
    "release-116/gtf/homo_sapiens/Homo_sapiens.GRCh38.116.gtf.gz$"
  )
  suppressMessages(download_reference(
    reference = "ensembl",
    division = "plants",
    species = "arabidopsis_thaliana",
    file_type = "fasta",
    output_path = out
  ))
  testthat::expect_match(
    captured$value,
    "Arabidopsis_thaliana.TAIR10.cdna.all.fa.gz$"
  )
  suppressMessages(download_reference(
    reference = "ncbi",
    accession = "GCF_000001405.40",
    output_path = out
  ))
  testthat::expect_match(
    captured$value,
    "GCF_000001405.40_GRCh38.p14_genomic.gff.gz$"
  )
  suppressMessages(download_reference(
    reference = "ncbi",
    accession = "gcf_000001405.40",
    file_type = "fasta",
    output_path = out
  ))
  testthat::expect_match(captured$value, "_rna.fna.gz$")
  suppressMessages(download_reference(reference = "mane", output_path = out))
  testthat::expect_match(
    captured$value,
    "MANE.GRCh38.v1.5.ensembl_genomic.gff.gz$"
  )
  suppressMessages(download_reference(
    reference = "mane",
    mane_source = "refseq",
    output_path = out
  ))
  testthat::expect_match(captured$value, "refseq_genomic.gff.gz$")
  testthat::expect_error(
    suppressMessages(download_reference(
      reference = "mane",
      file_type = "genome_fasta",
      output_path = out
    )),
    class = "isoformic_download_reference_file_not_found_error"
  )
  testthat::expect_error(
    suppressMessages(download_reference(
      reference = "ncbi",
      accession = "GCF_123",
      output_path = out
    )),
    class = "isoformic_download_reference_unsupported_version_error"
  )
  testthat::expect_error(
    suppressMessages(download_reference(
      reference = "ncbi",
      accession = "GCF_000001405.99",
      output_path = out
    )),
    class = "isoformic_download_reference_release_not_found_error"
  )
  testthat::expect_error(
    suppressMessages(download_reference(
      reference = "ensembl",
      species = "nope_nope",
      output_path = out
    )),
    class = "isoformic_download_reference_release_not_found_error"
  )
  testthat::expect_error(
    suppressMessages(download_reference(
      version = "50",
      file_type = "ncrna_fasta",
      output_path = out
    )),
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
  res <- suppressMessages(download_reference_url(
    "https://example.org/annot/chm13v2.0_RefSeq_Liftoff_v5.2.gff3.gz?x=1",
    output_path = out
  ))
  testthat::expect_equal(
    fs::path_file(res),
    "chm13v2.0_RefSeq_Liftoff_v5.2.gff3.gz"
  )
  testthat::expect_error(
    suppressMessages(download_reference_url("not a url", output_path = out)),
    class = "isoformic_download_reference_error"
  )

  files <- suppressMessages(download_hprc_annotation(
    "HG00438",
    output_path = out
  ))
  testthat::expect_named(files, c("hap1", "hap2"))
  testthat::expect_match(captured$value[2], "HG00438.1.chm13.gff3.gz$")
  testthat::expect_match(captured$value[3], "HG00438.2.chm13.gff3.gz$")
  one <- suppressMessages(download_hprc_annotation(
    "HG002",
    haplotype = 2,
    assembly = "hg38",
    output_path = out
  ))
  testthat::expect_match(captured$value[4], "HG002.2.hg38.gff3.gz$")
  testthat::expect_named(one, "hap2")
  testthat::expect_error(
    suppressMessages(download_hprc_annotation(
      "HG002",
      haplotype = 3,
      output_path = out
    )),
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
  res <- suppressMessages(download_reference(
    reference = "fantom",
    file_type = "gtf",
    output_path = out
  ))
  testthat::expect_match(res, "FANTOM_CAT.lv3_robust.gtf.gz$")
  testthat::expect_error(
    suppressMessages(download_reference(
      reference = "fantom",
      file_type = "gff",
      output_path = out
    )),
    class = "isoformic_download_reference_file_not_found_error"
  )
  testthat::expect_error(
    suppressMessages(download_reference(
      reference = "fantom",
      version = "lv9",
      file_type = "gtf",
      output_path = out
    )),
    class = "isoformic_download_reference_unsupported_version_error"
  )
  tables <- suppressMessages(download_fantom_cat_tables(
    "lv3_robust",
    output_path = out
  ))
  testthat::expect_named(tables, c("id_mapping", "gene_info", "cage_clusters"))
  testthat::expect_match(captured$value[2], "info_table.ID_mapping.tsv.gz$")
  testthat::expect_match(captured$value[4], "CAGE_cluster.bed.gz$")
})

# Cache first and cached fallback ---------------------------------------------

# Mocks a listing that always fails and counts its calls; any download fails.
local_failed_listing <- function(internet = TRUE, envir = parent.frame()) {
  calls <- new.env(parent = emptyenv())
  calls$urls <- character()
  testthat::local_mocked_bindings(
    has_internet = function(...) internet,
    .package = "curl",
    .env = envir
  )
  testthat::local_mocked_bindings(
    list_remote_dir = function(url, ...) {
      calls$urls <- c(calls$urls, url)
      NULL
    },
    .env = envir
  )
  testthat::local_mocked_bindings(
    download.file = function(url, ...) {
      stop("no download expected", call. = FALSE)
    },
    .package = "utils",
    .env = envir
  )
  calls
}

local_cache_dir <- function(files, envir = parent.frame()) {
  dir <- withr::local_tempdir(.local_envir = envir)
  for (f in files) {
    writeLines("cached", fs::path(dir, f))
  }
  dir
}

testthat::test_that("GENCODE files named by the release are used from the cache without listing", {
  calls <- local_failed_listing()
  out <- local_cache_dir(c(
    "gencode.v50.chr_patch_hapl_scaff.annotation.gff3.gz",
    "gencode.v50.annotation.gtf.gz",
    "gencode.v50.transcripts.fa.gz",
    "gencode.v50.polyAs.gff3.gz",
    "gencode.vM39.chr_patch_hapl_scaff.annotation.gff3.gz"
  ))
  testthat::expect_message(
    res <- download_reference(version = "50", output_path = out),
    "already exists"
  )
  testthat::expect_equal(
    as.character(res),
    as.character(fs::path(
      out,
      "gencode.v50.chr_patch_hapl_scaff.annotation.gff3.gz"
    ))
  )
  res_gtf <- suppressMessages(download_reference(
    version = 50,
    file_type = "gtf",
    region = "chr",
    output_path = out
  ))
  testthat::expect_equal(
    fs::path_file(res_gtf),
    "gencode.v50.annotation.gtf.gz"
  )
  res_fa <- suppressMessages(download_reference(
    file_type = "fasta",
    output_path = out
  ))
  testthat::expect_equal(fs::path_file(res_fa), "gencode.v50.transcripts.fa.gz")
  res_polya <- suppressMessages(download_reference(
    version = "v50",
    file_type = "polya",
    output_path = out
  ))
  testthat::expect_equal(fs::path_file(res_polya), "gencode.v50.polyAs.gff3.gz")
  res_mouse <- suppressMessages(download_reference(
    version = "39",
    organism = "mouse",
    output_path = out
  ))
  testthat::expect_match(res_mouse, "gencode.vM39.chr_patch_hapl_scaff")
  testthat::expect_length(calls$urls, 0L)

  # A file that is not cached still needs the listing
  testthat::expect_error(
    suppressMessages(download_reference(
      version = "50",
      region = "chr",
      output_path = out
    )),
    class = "isoformic_download_reference_release_not_found_error"
  )
  testthat::expect_true(length(calls$urls) > 0L)
})

testthat::test_that("cached files are used offline and missing ones report no connection", {
  calls <- local_failed_listing(internet = FALSE)
  out <- local_cache_dir("gencode.v50.chr_patch_hapl_scaff.annotation.gff3.gz")
  res <- suppressMessages(download_reference(version = "50", output_path = out))
  testthat::expect_true(fs::file_exists(res))
  testthat::expect_error(
    suppressMessages(download_reference(
      version = "50",
      file_type = "gtf",
      output_path = out
    )),
    class = "isoformic_download_reference_no_internet_error"
  )
  testthat::expect_length(calls$urls, 0L)
})

testthat::test_that("a failed GENCODE listing falls back to the one cached genome of the organism", {
  local_failed_listing()
  out <- local_cache_dir(c(
    "GRCh38.p14.genome.fa.gz",
    "GRCh38.primary_assembly.genome.fa.gz",
    "GRCm39.genome.fa.gz"
  ))
  suppressMessages(testthat::expect_message(
    res <- download_reference(
      version = "50",
      file_type = "genome_fasta",
      output_path = out
    ),
    "Could not list[\\s\\S]*Using the cached",
    perl = TRUE
  ))
  testthat::expect_equal(fs::path_file(res), "GRCh38.p14.genome.fa.gz")
  res_pri <- suppressMessages(download_reference(
    version = "50",
    file_type = "genome_fasta",
    region = "primary",
    output_path = out
  ))
  testthat::expect_equal(
    fs::path_file(res_pri),
    "GRCh38.primary_assembly.genome.fa.gz"
  )
  res_mouse <- suppressMessages(download_reference(
    version = "M39",
    organism = "mouse",
    file_type = "genome_fasta",
    output_path = out
  ))
  testthat::expect_equal(fs::path_file(res_mouse), "GRCm39.genome.fa.gz")

  # Two cached patches are ambiguous: the listing error is kept
  writeLines("cached", fs::path(out, "GRCh38.p13.genome.fa.gz"))
  testthat::expect_error(
    suppressMessages(download_reference(
      version = "50",
      file_type = "genome_fasta",
      output_path = out
    )),
    class = "isoformic_download_reference_release_not_found_error"
  )
})

testthat::test_that("a failed Ensembl listing falls back to a cached annotation of an explicit release", {
  local_failed_listing()
  out <- local_cache_dir(c(
    "Homo_sapiens.GRCh38.116.gff3.gz",
    "Homo_sapiens.GRCh38.116.chr_patch_hapl_scaff.gff3.gz",
    "Homo_sapiens.GRCh38.116.chr.gtf.gz",
    "Mus_musculus.GRCm39.116.gff3.gz",
    "Homo_sapiens.GRCh38.cdna.all.fa.gz"
  ))
  suppressMessages(testthat::expect_message(
    res <- download_reference(
      reference = "ensembl",
      version = 116,
      output_path = out
    ),
    "Could not list[\\s\\S]*Using the cached",
    perl = TRUE
  ))
  testthat::expect_equal(
    fs::path_file(res),
    "Homo_sapiens.GRCh38.116.chr_patch_hapl_scaff.gff3.gz"
  )
  res_pri <- suppressMessages(download_reference(
    reference = "ensembl",
    version = "116",
    region = "primary",
    output_path = out
  ))
  testthat::expect_equal(
    fs::path_file(res_pri),
    "Homo_sapiens.GRCh38.116.gff3.gz"
  )
  res_gtf <- suppressMessages(download_reference(
    reference = "ensembl",
    version = "116",
    file_type = "gtf",
    region = "chr",
    output_path = out
  ))
  testthat::expect_equal(
    fs::path_file(res_gtf),
    "Homo_sapiens.GRCh38.116.chr.gtf.gz"
  )
  res_mouse <- suppressMessages(download_reference(
    reference = "ensembl",
    organism = "mouse",
    version = "116",
    region = "primary",
    output_path = out
  ))
  testthat::expect_equal(
    fs::path_file(res_mouse),
    "Mus_musculus.GRCm39.116.gff3.gz"
  )
  # Offline, a plain annotation may be the region = "primary" file of a
  # species with patches, so it never stands in for region = "all"
  testthat::expect_error(
    suppressMessages(download_reference(
      reference = "ensembl",
      organism = "mouse",
      version = "116",
      output_path = out
    )),
    class = "isoformic_download_reference_release_not_found_error"
  )
  # The current release always needs the listing
  testthat::expect_error(
    suppressMessages(download_reference(
      reference = "ensembl",
      output_path = out
    )),
    class = "isoformic_download_reference_release_not_found_error"
  )
  # Transcript FASTA names carry no release, so they have no fallback
  testthat::expect_error(
    suppressMessages(download_reference(
      reference = "ensembl",
      version = "116",
      file_type = "fasta",
      output_path = out
    )),
    class = "isoformic_download_reference_release_not_found_error"
  )
  testthat::expect_error(
    suppressMessages(download_reference(
      reference = "ensembl",
      version = "115",
      output_path = out
    )),
    class = "isoformic_download_reference_release_not_found_error"
  )
})

testthat::test_that("NCBI, MANE and FANTOM CAT use cached files when the listing fails", {
  calls <- local_failed_listing()
  out <- local_cache_dir(c(
    "GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz",
    "GCF_009914755.1_T2T-CHM13v2.0_cds_from_genomic.fna.gz",
    "GCF_009914755.1_T2T-CHM13v2.0_rna_from_genomic.fna.gz",
    "MANE.GRCh38.v1.5.ensembl_genomic.gff.gz",
    "FANTOM_CAT.lv3_robust.gtf.gz"
  ))
  suppressMessages(testthat::expect_message(
    res <- download_reference(
      reference = "t2t",
      file_type = "genome_fasta",
      output_path = out
    ),
    "Could not list[\\s\\S]*Using the cached",
    perl = TRUE
  ))
  testthat::expect_equal(
    fs::path_file(res),
    "GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz"
  )
  testthat::expect_error(
    suppressMessages(download_reference(
      reference = "t2t",
      file_type = "gff",
      output_path = out
    )),
    class = "isoformic_download_reference_release_not_found_error"
  )

  calls$urls <- character()
  res_mane <- suppressMessages(download_reference(
    reference = "mane",
    version = "v1.5",
    output_path = out
  ))
  testthat::expect_equal(
    fs::path_file(res_mane),
    "MANE.GRCh38.v1.5.ensembl_genomic.gff.gz"
  )
  res_fantom <- suppressMessages(download_reference(
    reference = "fantom",
    file_type = "gtf",
    output_path = out
  ))
  testthat::expect_equal(
    fs::path_file(res_fantom),
    "FANTOM_CAT.lv3_robust.gtf.gz"
  )
  testthat::expect_length(calls$urls, 0L)
  # The current MANE release always needs the listing
  testthat::expect_error(
    suppressMessages(download_reference(reference = "mane", output_path = out)),
    class = "isoformic_download_reference_release_not_found_error"
  )
})

testthat::test_that("the Ensembl region 'all' fallback needs a cached chr_patch_hapl_scaff file", {
  local_failed_listing()
  out <- local_cache_dir(c(
    "Homo_sapiens.GRCh38.115.gff3.gz",
    "Homo_sapiens.GRCh38.115.gtf.gz"
  ))
  for (type in c("gff", "gtf")) {
    testthat::expect_error(
      resolve_ensembl_file(
        "115",
        "homo_sapiens",
        "vertebrates",
        type,
        "all",
        out
      ),
      class = "isoformic_download_reference_release_not_found_error"
    )
  }
  suppressMessages(testthat::expect_message(
    res <- resolve_ensembl_file(
      "115",
      "homo_sapiens",
      "vertebrates",
      "gff",
      "primary",
      out
    ),
    "Could not list[\\s\\S]*Using the cached",
    perl = TRUE
  ))
  testthat::expect_equal(res$file, "Homo_sapiens.GRCh38.115.gff3.gz")
})

testthat::test_that("the GENCODE genome fallback only uses the assembly of the release", {
  local_failed_listing()
  out <- local_cache_dir(c("GRCh38.p14.genome.fa.gz", "GRCm39.genome.fa.gz"))
  testthat::expect_identical(gencode_assembly("19", "human"), "GRCh37")
  testthat::expect_identical(gencode_assembly("20", "human"), "GRCh38")
  testthat::expect_identical(gencode_assembly("M25", "mouse"), "GRCm38")
  testthat::expect_identical(gencode_assembly("M26", "mouse"), "GRCm39")
  testthat::expect_error(
    suppressMessages(download_reference(
      version = "19",
      file_type = "genome_fasta",
      output_path = out
    )),
    class = "isoformic_download_reference_release_not_found_error"
  )
  testthat::expect_error(
    suppressMessages(download_reference(
      version = "M20",
      organism = "mouse",
      file_type = "genome_fasta",
      output_path = out
    )),
    class = "isoformic_download_reference_release_not_found_error"
  )
  # A GRCh37 genome cached for release 19 is not used for release 50
  fs::file_delete(fs::path(out, "GRCh38.p14.genome.fa.gz"))
  writeLines("cached", fs::path(out, "GRCh37.p13.genome.fa.gz"))
  testthat::expect_error(
    suppressMessages(download_reference(
      version = "50",
      file_type = "genome_fasta",
      output_path = out
    )),
    class = "isoformic_download_reference_release_not_found_error"
  )
  suppressMessages(testthat::expect_message(
    res <- download_reference(
      version = "19",
      file_type = "genome_fasta",
      output_path = out
    ),
    "Could not list[\\s\\S]*Using the cached",
    perl = TRUE
  ))
  testthat::expect_equal(fs::path_file(res), "GRCh37.p13.genome.fa.gz")
})

testthat::test_that("the GENCODE genome fallback is not used for a release missing from the listing of releases", {
  testthat::local_mocked_bindings(
    has_internet = function(...) TRUE,
    .package = "curl"
  )
  testthat::local_mocked_bindings(
    list_remote_dir = function(url, ...) {
      if (grepl("/Gencode_human/$", url)) {
        c("release_49", "release_50", "latest_release")
      }
    }
  )
  out <- local_cache_dir("GRCh38.p14.genome.fa.gz")
  testthat::expect_error(
    suppressMessages(download_reference(
      version = "99",
      file_type = "genome_fasta",
      output_path = out
    )),
    class = "isoformic_download_reference_release_not_found_error"
  )
  suppressMessages(testthat::expect_message(
    res <- resolve_gencode_file("50", "human", "genome_fasta", "all", out),
    "Could not list[\\s\\S]*Using the cached",
    perl = TRUE
  ))
  testthat::expect_equal(res$file, "GRCh38.p14.genome.fa.gz")
})

testthat::test_that("cached symbolic links to files are used when the listing fails", {
  testthat::skip_on_os("windows")
  local_failed_listing()
  store <- local_cache_dir(c(
    "GRCh38.p14.genome.fa.gz",
    "Homo_sapiens.GRCh38.115.chr_patch_hapl_scaff.gff3.gz"
  ))
  out <- local_cache_dir(character())
  for (f in fs::path_file(fs::dir_ls(store))) {
    fs::link_create(fs::path(store, f), fs::path(out, f))
  }
  # A dangling link is not a cached file
  fs::link_create(
    fs::path(store, "missing.fa.gz"),
    fs::path(out, "GRCh38.p13.genome.fa.gz")
  )
  suppressMessages(testthat::expect_message(
    res <- resolve_gencode_file("50", "human", "genome_fasta", "all", out),
    "Could not list[\\s\\S]*Using the cached",
    perl = TRUE
  ))
  testthat::expect_equal(res$file, "GRCh38.p14.genome.fa.gz")
  suppressMessages(testthat::expect_message(
    res_ens <- resolve_ensembl_file(
      "115",
      "homo_sapiens",
      "vertebrates",
      "gff",
      "all",
      out
    ),
    "Could not list[\\s\\S]*Using the cached",
    perl = TRUE
  ))
  testthat::expect_equal(
    res_ens$file,
    "Homo_sapiens.GRCh38.115.chr_patch_hapl_scaff.gff3.gz"
  )
})

testthat::test_that("the NCBI fallback URL includes the assembly directory", {
  local_failed_listing()
  out <- local_cache_dir("GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz")
  suppressMessages(testthat::expect_message(
    res <- resolve_ncbi_file("GCF_009914755.1", "genome_fasta", out),
    "Could not list[\\s\\S]*Using the cached",
    perl = TRUE
  ))
  testthat::expect_equal(
    res$url,
    paste0(
      "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/",
      "GCF_009914755.1_T2T-CHM13v2.0/",
      "GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz"
    )
  )
})
