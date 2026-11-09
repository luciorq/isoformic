# Frozen reference datasets --------------------------------------------------
#
# External evidence tables used by the productivity and polyadenylation
# analyses. Each dataset is pinned to a source version, rebuilt from the
# source files by a `build_frozen_*()` function (so the processing is code,
# not a one-off script) and distributed as a Parquet file attached to a
# GitHub release of this repository. `fetch_frozen_data()` downloads the
# release file, or rebuilds it from the sources when it cannot, and verifies
# the MD5 recorded in the registry. See data-raw/README.md for the update
# procedure.

#' Frozen Reference Datasets
#'
#' @description
#' `isoformic` ships evidence tables derived from public resources, frozen at
#' a given version so that results are reproducible:
#'
#' * `appris`: APPRIS principal/alternative isoform labels and method
#'   scores (functional residues, structure, conservation, domains) per
#'   transcript, from the APPRIS release files.
#' * `appris_domains`: Pfam domain genomic coordinates per transcript
#'   (APPRIS SPADE).
#' * `riboseq_orfs`: GENCODE Ribo-seq ORF catalog, phase 1 and phase 2,
#'   with exon blocks and compatible transcripts.
#' * `nmd_exons`: NMD-inducing cassette exons validated by NMD inhibition
#'   and nanopore sequencing (Karousis et al. 2021, Genome Biology,
#'   Supplementary Table 1).
#' * `polyasite`: PolyASite 2.0 poly(A) site atlas (GRCh38) with site class
#'   and poly(A) signals.
#' * `gtex_transcript_tpm`: GTEx v8 median transcript TPM per tissue (RSEM).
#'
#' `frozen_datasets()` lists them with version, source URLs and the MD5 of
#' the frozen file. `fetch_frozen_data()` returns the local path of a
#' dataset, downloading the frozen file from the first mirror that has it
#' (`frozen_mirrors()`: any URL in `options(isoformic.frozen_mirrors = )`,
#' the package's GitHub release, then the Zenodo snapshot of the processed
#' files) or, when that is not possible (`build = "auto"`) or asked for
#' (`build = "always"`), rebuilding it from the source files with the
#' matching `build_frozen_*()` function. `read_frozen_data()` fetches and
#' reads it. `bundle_frozen_data()` assembles the snapshot that is deposited
#' on Zenodo: every frozen file, the coding potential model, an `MD5SUMS`
#' file and a README stating versions and sources.
#'
#' @param name Dataset name, one of `names(frozen_datasets())`.
#' @param output_path Directory holding the frozen files; `":cache:"`
#'   (default) uses `get_isoformic_cache("frozen")`.
#' @param build `"auto"` (download, else build from sources), `"never"` or
#'   `"always"`.
#' @param verify Logical, check the MD5 of the file against the registry
#'   (default `TRUE`; a mismatch is a warning, since a rebuilt file from a
#'   newer source is still usable).
#' @param release_tag Tag of the GitHub release hosting the frozen files.
#' @param mirrors Character vector of base URLs tried in order for the
#'   frozen file; default [frozen_mirrors()].
#' @param ... Reserved for future use; must be empty.
#'
#' @returns `frozen_datasets()` returns a tibble; `fetch_frozen_data()` the
#'   file path, invisibly; `read_frozen_data()` a tibble;
#'   `frozen_mirrors()` a character vector; `bundle_frozen_data()` the
#'   bundle directory, invisibly.
#'
#' @examples
#' frozen_datasets()[, c("name", "version", "file")]
#' \dontrun{
#' appris <- read_frozen_data("appris")
#' head(appris)
#' }
#'
#' @name frozen_data
#' @rdname frozen_data
#' @export
frozen_datasets <- function() {
  reg <- frozen_registry()
  tibble::tibble(
    name = names(reg),
    title = vapply(reg, `[[`, character(1L), "title"),
    version = vapply(reg, `[[`, character(1L), "version"),
    file = vapply(reg, `[[`, character(1L), "file"),
    md5 = vapply(reg, `[[`, character(1L), "md5"),
    sources = vapply(
      reg,
      function(r) paste(r$sources, collapse = " "),
      character(1L)
    ),
    builder = vapply(reg, `[[`, character(1L), "builder")
  )
}

frozen_release_base <- function(release_tag) {
  paste0(
    "https://github.com/luciorq/isoformic/releases/download/",
    release_tag,
    "/"
  )
}

#' Zenodo record holding the snapshot of the frozen files
#'
#' Set to the numeric record id once the deposit described in
#' data-raw/README.md is published (a versioned DOI resolves to the record
#' URL `https://zenodo.org/records/<id>`); `NA` while there is none.
#'
#' @keywords internal
#' @noRd
frozen_zenodo_record <- function() {
  NA_character_
}

#' @rdname frozen_data
#' @param zenodo_record Zenodo record id (or `https://zenodo.org/records/<id>`
#'   URL) of a deposit holding the frozen files; defaults to the record
#'   published with the package, `NA` when there is none.
#' @export
frozen_mirrors <- function(
  release_tag = "frozen-data",
  zenodo_record = frozen_zenodo_record()
) {
  mirrors <- getOption("isoformic.frozen_mirrors", character())
  mirrors <- c(mirrors, frozen_release_base(release_tag))
  if (!is.na(zenodo_record) && nzchar(zenodo_record)) {
    id <- sub("^.*/records?/([0-9]+).*$", "\\1", as.character(zenodo_record))
    mirrors <- c(mirrors, paste0("https://zenodo.org/records/", id, "/files/"))
  }
  unique(mirrors)
}

frozen_mirror_url <- function(base, file) {
  url <- paste0(sub("/*$", "/", base), file)
  if (grepl("^https://zenodo\\.org/", url)) {
    url <- paste0(url, "?download=1")
  }
  url
}

#' Registry of the frozen datasets
#'
#' `md5` is the checksum of the frozen Parquet file; update it (and
#' `version`) in the same commit that publishes a rebuilt file, see
#' data-raw/README.md. `NA` means not yet published.
#'
#' @keywords internal
#' @noRd
frozen_registry <- function() {
  appris_base <- "https://apprisws.bioinfo.cnio.es/pub/releases/2024_10.v49/datafiles/homo_sapiens/e113v49/"
  riboseq_base <- "https://ftp.ebi.ac.uk/pub/databases/gencode/riboseq_orfs/"
  list(
    appris = list(
      title = "APPRIS principal isoforms and method scores (human, GENCODE 49)",
      version = "e113v49",
      file = "appris_e113v49.parquet",
      md5 = "fc97d69393054b95c60f2a112d8861fe",
      sources = paste0(appris_base, "appris_data.appris.txt"),
      builder = "build_frozen_appris"
    ),
    appris_domains = list(
      title = "APPRIS SPADE Pfam domain coordinates (human, GENCODE 49)",
      version = "e113v49",
      file = "appris_domains_e113v49.parquet",
      md5 = "ee0da6fa024f40340865930f6580c058",
      sources = paste0(appris_base, "appris_method.spade.gtf.gz"),
      builder = "build_frozen_appris_domains"
    ),
    riboseq_orfs = list(
      title = "GENCODE Ribo-seq ORFs, phase 1 (v35) and phase 2 (v45)",
      version = "phase1_phase2",
      file = "gencode_riboseq_orfs_phase1_phase2.parquet",
      md5 = "2e749f6411597892f45282ac524d39a2",
      sources = c(
        paste0(riboseq_base, "data/Ribo-seq_ORFs.bed"),
        paste0(riboseq_base, "phase2/Ribo-seq_ORFs.comprehensive.bed")
      ),
      builder = "build_frozen_riboseq_orfs"
    ),
    nmd_exons = list(
      title = "NMD-inducing exons validated by NMD inhibition (Karousis et al. 2021)",
      version = "karousis2021",
      file = "nmd_exons_karousis2021.parquet",
      md5 = "729c344fd40f8d83948ad1ebdf4d80f1",
      sources = paste0(
        "https://static-content.springer.com/esm/art%3A10.1186%2Fs13059-021-02439-3/",
        "MediaObjects/13059_2021_2439_MOESM2_ESM.xlsx"
      ),
      builder = "build_frozen_nmd_exons"
    ),
    polyasite = list(
      title = "PolyASite 2.0 poly(A) site atlas (GRCh38.96)",
      version = "2.0_GRCh38.96",
      file = "polyasite_2.0_GRCh38.96.parquet",
      md5 = "c9ece480e5a38fb801adfce207df05bb",
      sources = "https://polyasite.unibas.ch/download/atlas/2.0/GRCh38.96/atlas.clusters.2.0.GRCh38.96.bed.gz",
      builder = "build_frozen_polyasite"
    ),
    gtex_transcript_tpm = list(
      title = "GTEx v8 median transcript TPM per tissue",
      version = "v8",
      file = "gtex_v8_transcript_tpm_median_by_tissue.parquet",
      md5 = "6b755081fabf50d4d9f7f35ce5e80237",
      sources = c(
        paste0(
          "https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/",
          "GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct.gz"
        ),
        paste0(
          "https://storage.googleapis.com/adult-gtex/annotations/v8/metadata-files/",
          "GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"
        )
      ),
      builder = "build_frozen_gtex"
    )
  )
}

frozen_dataset_entry <- function(name) {
  reg <- frozen_registry()
  if (!rlang::is_string(name) || !name %in% names(reg)) {
    cli::cli_abort(
      message = c(
        x = "Unknown frozen dataset {.val {name}}.",
        i = "Available: {.val {names(reg)}}."
      ),
      class = "isoformic_unknown_frozen_dataset"
    )
  }
  reg[[name]]
}

frozen_dir <- function(output_path = ":cache:") {
  if (identical(output_path, ":cache:")) {
    output_path <- get_isoformic_cache("frozen")
  }
  fs::dir_create(output_path)
  as.character(output_path)
}

#' @rdname frozen_data
#' @export
fetch_frozen_data <- function(
  name,
  output_path = ":cache:",
  build = c("auto", "never", "always"),
  verify = TRUE,
  release_tag = "frozen-data",
  mirrors = frozen_mirrors(release_tag),
  ...
) {
  rlang::check_dots_empty()
  build <- rlang::arg_match(build)
  entry <- frozen_dataset_entry(name)
  dir <- frozen_dir(output_path)
  target <- fs::path(dir, entry$file)
  if (isTRUE(fs::file_exists(target)) && !identical(build, "always")) {
    verify_frozen_md5(target, entry, verify)
    return(invisible(as.character(target)))
  }
  if (!identical(build, "always")) {
    urls <- vapply(mirrors, frozen_mirror_url, character(1L), file = entry$file)
    for (url in urls) {
      downloaded <- rlang::try_fetch(
        {
          download_reference_file(url, as.character(target))
          TRUE
        },
        error = function(cnd) FALSE
      )
      if (isTRUE(downloaded)) {
        # A download that fails its checksum is corrupt or stale: discard it
        # and try the next mirror instead of caching it
        if (isTRUE(verify_frozen_md5(target, entry, verify, warn = FALSE))) {
          return(invisible(as.character(target)))
        }
        fs::file_delete(target)
        cli::cli_inform(c(
          "!" = "Discarded {.url {url}}: its MD5 does not match the registry."
        ))
      }
    }
    if (identical(build, "never")) {
      cli::cli_abort(
        message = c(
          x = "Could not download the frozen dataset {.val {name}} from any mirror.",
          i = "Tried: {.url {urls}}.",
          i = "Use {.code build = \"auto\"} to rebuild it from its sources."
        ),
        class = "isoformic_download_reference_error"
      )
    }
    cli::cli_inform(c(
      i = "Frozen file not available for download; building {.val {name}} from its sources."
    ))
  }
  builder <- get(entry$builder, mode = "function")
  path <- builder(output_path = dir)
  verify_frozen_md5(path, entry, verify)
  invisible(as.character(path))
}

#' @rdname frozen_data
#' @export
read_frozen_data <- function(
  name,
  output_path = ":cache:",
  build = c("auto", "never", "always"),
  ...
) {
  rlang::check_dots_empty()
  path <- fetch_frozen_data(name, output_path = output_path, build = build)
  read_parquet(path)
}

#' @rdname frozen_data
#' @param output_dir Directory receiving the bundle (created if needed).
#' @param include_coding_model Logical, also copy the coding potential model
#'   files (default `TRUE`).
#' @export
bundle_frozen_data <- function(
  output_dir,
  output_path = ":cache:",
  build = c("auto", "never", "always"),
  include_coding_model = TRUE,
  verify = TRUE
) {
  check_string(output_dir)
  build <- rlang::arg_match(build)
  fs::dir_create(output_dir)
  reg <- frozen_datasets()
  files <- vapply(
    reg$name,
    function(name) {
      fetch_frozen_data(
        name,
        output_path = output_path,
        build = build,
        verify = verify
      )
    },
    character(1L)
  )
  copied <- vapply(
    files,
    function(f) {
      dest <- fs::path(output_dir, fs::path_file(f))
      fs::file_copy(f, dest, overwrite = TRUE)
      as.character(dest)
    },
    character(1L)
  )
  model_files <- character()
  if (isTRUE(include_coding_model)) {
    model_dir <- system.file("extdata", "coding_model", package = "isoformic")
    model_files <- as.character(fs::dir_ls(model_dir, glob = "*.parquet"))
    for (f in model_files) {
      fs::file_copy(f, fs::path(output_dir, fs::path_file(f)), overwrite = TRUE)
    }
  }
  all_files <- fs::path_file(c(copied, model_files))
  md5 <- unname(tools::md5sum(as.character(fs::path(output_dir, all_files))))
  writeLines(paste(md5, all_files), fs::path(output_dir, "MD5SUMS"))
  version <- as.character(utils::packageVersion("isoformic"))
  readme <- c(
    "# isoformic frozen evidence datasets",
    "",
    paste0(
      "Snapshot of the processed reference tables used by the isoformic R package (version ",
      version,
      ", https://github.com/luciorq/isoformic), built on ",
      format(Sys.Date()),
      " with the `build_frozen_*()` functions of the package. ",
      "Each file is a Parquet table; `MD5SUMS` lists their checksums. ",
      "The package verifies these checksums (`frozen_datasets()$md5`) and can rebuild any file ",
      "from the original sources (`data-raw/README.md` in the repository)."
    ),
    "",
    "| file | dataset | version | sources |",
    "|---|---|---|---|",
    paste0(
      "| `",
      reg$file,
      "` | ",
      reg$name,
      " | ",
      reg$version,
      " | ",
      gsub(" ", " ; ", reg$sources, fixed = TRUE),
      " |"
    ),
    if (length(model_files)) {
      c(
        paste0(
          "| `",
          fs::path_file(model_files),
          "` | coding potential model | GENCODE 49 | GENCODE pc_transcripts.fa.gz and transcripts.fa.gz |"
        )
      )
    },
    "",
    "Use with the package:",
    "",
    "```r",
    "options(isoformic.frozen_mirrors = \"https://zenodo.org/records/<record>/files/\")",
    "fetch_frozen_data(\"appris\")",
    "```"
  )
  writeLines(readme, fs::path(output_dir, "README.md"))
  cli::cli_inform(c(
    v = "Bundled {length(all_files)} file{?s} into {.path {output_dir}} (see MD5SUMS and README.md)."
  ))
  invisible(as.character(output_dir))
}

verify_frozen_md5 <- function(path, entry, verify, warn = TRUE) {
  if (!isTRUE(verify) || is.na(entry$md5)) {
    return(invisible(TRUE))
  }
  md5 <- unname(tools::md5sum(as.character(path)))
  if (!identical(md5, entry$md5)) {
    if (!isTRUE(warn)) {
      return(invisible(FALSE))
    }
    cli::cli_warn(c(
      "!" = paste0(
        "{.path {path}} has MD5 {.val {md5}}, the registry expects ",
        "{.val {entry$md5}} (version {.val {entry$version}})."
      ),
      i = paste0(
        "The file was built from newer sources or is corrupt; ",
        "rebuild with {.code fetch_frozen_data(build = \"always\")}."
      )
    ))
    return(invisible(FALSE))
  }
  invisible(TRUE)
}

#' Download the source files of a frozen dataset into the sources cache
#' @keywords internal
#' @noRd
frozen_sources <- function(name, sources_dir = NULL) {
  entry <- frozen_dataset_entry(name)
  sources_dir <- sources_dir %||% get_isoformic_cache("sources", name)
  fs::dir_create(sources_dir)
  vapply(
    entry$sources,
    function(url) {
      download_reference_file(
        url,
        as.character(fs::path(
          sources_dir,
          fs::path_file(sub("[?#].*$", "", url))
        ))
      )
    },
    character(1L),
    USE.NAMES = FALSE
  )
}

# Builders -------------------------------------------------------------------

#' Build the Frozen Datasets from Their Sources
#'
#' @description
#' Each function downloads the source files of one frozen dataset (see
#' [frozen_datasets()]) into the sources cache, converts them to the Parquet
#' table used by the package and returns the file path. They are what
#' [fetch_frozen_data()] runs when a frozen file cannot be downloaded, and
#' what `data-raw/build_frozen_data.R` runs to publish a new version. The
#' resulting tables:
#'
#' * `build_frozen_appris()`: one row per transcript with
#'   `appris_annotation` (`PRINCIPAL:1` to `ALTERNATIVE:2`, `MINOR`),
#'   `appris_score`, the method scores (`firestar_residues`,
#'   `matador_structure`, `corsair_conservation`, `spade_domain_score`,
#'   `thump_tm_helices`, `crash_signal`, `trifid_score`), `protein_length`,
#'   `ccds_id` and `mane` label.
#' * `build_frozen_appris_domains()`: one row per Pfam domain per transcript
#'   with its genomic span (`seqid`, `start_pos`, `end_pos`, `strand`),
#'   `pfam_acc`, `pfam_name`, `evalue`, `pep_start`, `pep_end`.
#' * `build_frozen_riboseq_orfs()`: one row per Ribo-seq ORF with `phase`,
#'   genomic span, `block_starts`/`block_ends` (1-based, `;`-separated),
#'   `orf_type`, `gene_id`, `gene_name`, `transcript_ids` (compatible
#'   transcripts, phase 1 only), `replicated` and `studies`.
#' * `build_frozen_nmd_exons()`: one row per validated NMD-inducing exon with
#'   `seqid`, `start_pos`, `end_pos`, `strand`, the flanking exons,
#'   `nmd_transcript_ids`, `non_nmd_transcript_ids`, `gene_id`, `gene_name`.
#'   Reading the supplementary spreadsheet needs the `readxl` package.
#' * `build_frozen_polyasite()`: one row per poly(A) site cluster with
#'   `seqid`, `start_pos`, `end_pos`, `cleavage_pos`, `strand`, `tpm`,
#'   `frac_samples`, `n_protocols`, `site_class` (TE terminal exon, EX exon,
#'   IN intron, DS downstream, AE/AI/AU antisense, IG intergenic) and
#'   `signals`.
#' * `build_frozen_gtex()`: one row per transcript and tissue with
#'   `median_tpm` and `n_samples` (GTEx v8, RSEM transcript TPM; identifiers
#'   are GENCODE v26, stored without version).
#'
#' @param output_path Directory receiving the Parquet file (default: the
#'   frozen-data cache).
#' @param sources_dir Directory holding (or receiving) the source files;
#'   default: `get_isoformic_cache("sources", <name>)`.
#'
#' @returns The path of the written Parquet file, invisibly.
#'
#' @examples
#' \dontrun{
#' build_frozen_riboseq_orfs()
#' }
#'
#' @name build_frozen_data
#' @rdname build_frozen_data
NULL

frozen_target <- function(name, output_path) {
  entry <- frozen_dataset_entry(name)
  as.character(fs::path(frozen_dir(output_path), entry$file))
}

report_frozen <- function(path, n) {
  cli::cli_inform(c(
    v = "Wrote {n} row{?s} to {.path {path}} (MD5 {.val {unname(tools::md5sum(path))}})."
  ))
  invisible(path)
}

#' @rdname build_frozen_data
#' @export
build_frozen_appris <- function(output_path = ":cache:", sources_dir = NULL) {
  src <- frozen_sources("appris", sources_dir)
  raw <- read_delimited(src[[1L]], null_strings = c("", "-", "NA"))
  col <- function(pattern) {
    i <- grep(pattern, colnames(raw), ignore.case = TRUE)[1L]
    if (is.na(i)) rep(NA, nrow(raw)) else raw[[i]]
  }
  num <- function(pattern) suppressWarnings(as.numeric(col(pattern)))
  res <- tibble::tibble(
    transcript_id = strip_id_version(as.character(col("^Transcript ID"))),
    gene_id = strip_id_version(as.character(col("^Ensembl Gene ID|^Gene ID"))),
    gene_name = as.character(col("^Gene name")),
    translation_id = as.character(col("^Translation ID")),
    transcript_type = as.character(col("^Transcript type")),
    ccds_id = as.character(col("^CCDS")),
    protein_length = num("^Protein length"),
    firestar_residues = num("firestar"),
    matador_structure = num("Matador"),
    corsair_conservation = num("CORSAIR"),
    spade_domain_score = num("SPADE"),
    thump_tm_helices = num("THUMP"),
    crash_signal = as.character(col("CRASH")),
    trifid_score = num("Trifid"),
    appris_score = num("^APPRIS score"),
    appris_annotation = as.character(col("^APPRIS Annotation")),
    mane = as.character(col("^MANE"))
  )
  res <- res[!is.na(res$transcript_id), ]
  res$appris_principal <- grepl("^PRINCIPAL", res$appris_annotation)
  path <- frozen_target("appris", output_path)
  write_parquet(res, path)
  report_frozen(path, nrow(res))
}

#' @rdname build_frozen_data
#' @export
build_frozen_appris_domains <- function(
  output_path = ":cache:",
  sources_dir = NULL
) {
  src <- frozen_sources("appris_domains", sources_dir)
  gtf <- read_tabular_lines(src[[1L]], n_cols = 9L)
  attrs <- gtf[[9L]]
  get_attr <- function(key) {
    m <- regmatches(attrs, regexpr(paste0(key, ' "[^"]*"'), attrs))
    out <- rep(NA_character_, length(attrs))
    out[grepl(paste0(key, ' "'), attrs)] <- sub(
      paste0("^", key, " \"([^\"]*)\"$"),
      "\\1",
      m
    )
    out
  }
  note <- get_attr("note")
  note_field <- function(key) {
    m <- regmatches(note, regexpr(paste0(key, ":[^,]*"), note))
    out <- rep(NA_character_, length(note))
    out[grepl(paste0(key, ":"), note)] <- sub(paste0("^", key, ":"), "", m)
    out
  }
  res <- tibble::tibble(
    transcript_id = strip_id_version(get_attr("transcript_id")),
    gene_id = strip_id_version(get_attr("gene_id")),
    seqid = ifelse(
      grepl("^chr", gtf[[1L]]),
      gtf[[1L]],
      paste0("chr", gtf[[1L]])
    ),
    start_pos = as.integer(gtf[[4L]]),
    end_pos = as.integer(gtf[[5L]]),
    strand = gtf[[7L]],
    score = suppressWarnings(as.numeric(gtf[[6L]])),
    pfam_acc = note_field("hmm_acc"),
    pfam_name = note_field("hmm_name"),
    evalue = suppressWarnings(as.numeric(note_field("evalue"))),
    pep_start = as.integer(note_field("pep_start")),
    pep_end = as.integer(note_field("pep_end"))
  )
  res <- dplyr::distinct(res)
  path <- frozen_target("appris_domains", output_path)
  write_parquet(res, path)
  report_frozen(path, nrow(res))
}

#' Read a tab-separated file without header as character columns (DuckDB)
#' @keywords internal
#' @noRd
read_tabular_lines <- function(path, n_cols) {
  path <- as.character(fs::path_abs(path))
  skip <- count_comment_lines(path)
  con <- duckdb_connect()
  cols <- paste0("'c", seq_len(n_cols), "': 'VARCHAR'", collapse = ", ")
  sql <- glue::glue_sql(
    "SELECT * FROM read_csv({path}, delim = '\t', header = false, quote = '', skip = {skip},
       columns = {{{DBI::SQL(cols)}}}, null_padding = true, strict_mode = false,
       max_line_size = 10000000)",
    .con = con
  )
  tibble::as_tibble(DBI::dbGetQuery(con, sql))
}

#' Number of leading `#` comment lines of a (possibly gzipped) text file
#' @keywords internal
#' @noRd
count_comment_lines <- function(path, max_lines = 10000L) {
  con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    gzfile(path, open = "rt")
  } else {
    file(path, open = "rt")
  }
  withr::defer(close(con))
  lines <- readLines(con, n = max_lines, warn = FALSE)
  is_comment <- startsWith(lines, "#")
  if (all(is_comment)) {
    return(length(lines))
  }
  which(!is_comment)[[1L]] - 1L
}

#' @rdname build_frozen_data
#' @export
build_frozen_riboseq_orfs <- function(
  output_path = ":cache:",
  sources_dir = NULL
) {
  src <- frozen_sources("riboseq_orfs", sources_dir)
  parse_bed <- function(path, phase_no) {
    bed <- read_tabular_lines(path, n_cols = if (phase_no == 1L) 26L else 12L)
    chrom_start <- as.integer(bed[[2L]])
    sizes <- strsplit(sub(",$", "", bed[[11L]]), ",", fixed = TRUE)
    offsets <- strsplit(sub(",$", "", bed[[12L]]), ",", fixed = TRUE)
    block_start <- lapply(seq_along(sizes), function(i) {
      chrom_start[[i]] + as.integer(offsets[[i]]) + 1L
    })
    block_end <- lapply(seq_along(sizes), function(i) {
      block_start[[i]] + as.integer(sizes[[i]]) - 1L
    })
    name <- bed[[4L]]
    tibble::tibble(
      orf_id = if (phase_no == 1L) name else sub(":.*$", "", name),
      phase = phase_no,
      seqid = bed[[1L]],
      start_pos = chrom_start + 1L,
      end_pos = as.integer(bed[[3L]]),
      strand = bed[[6L]],
      n_blocks = as.integer(bed[[10L]]),
      block_starts = vapply(block_start, paste, character(1L), collapse = ";"),
      block_ends = vapply(block_end, paste, character(1L), collapse = ";"),
      orf_type = if (phase_no == 1L) bed[[17L]] else NA_character_,
      gene_id = if (phase_no == 1L) {
        strip_id_version(bed[[18L]])
      } else {
        NA_character_
      },
      gene_name = if (phase_no == 1L) {
        bed[[19L]]
      } else {
        sub("^[^:]*:([^:]*):.*$", "\\1", name)
      },
      transcript_biotype = if (phase_no == 1L) {
        bed[[21L]]
      } else {
        sub("^.*:", "", name)
      },
      transcript_ids = if (phase_no == 1L) {
        vapply(
          strsplit(bed[[23L]], ",", fixed = TRUE),
          function(x) paste(strip_id_version(x), collapse = ";"),
          character(1L)
        )
      } else {
        NA_character_
      },
      replicated = if (phase_no == 1L) bed[[25L]] == "yes" else NA,
      studies = if (phase_no == 1L) bed[[26L]] else NA_character_,
      protein = if (phase_no == 1L) bed[[22L]] else NA_character_
    )
  }
  res <- dplyr::bind_rows(parse_bed(src[[1L]], 1L), parse_bed(src[[2L]], 2L))
  path <- frozen_target("riboseq_orfs", output_path)
  write_parquet(res, path)
  report_frozen(path, nrow(res))
}

#' @rdname build_frozen_data
#' @export
build_frozen_nmd_exons <- function(
  output_path = ":cache:",
  sources_dir = NULL
) {
  rlang::check_installed(
    "readxl",
    reason = "to read the supplementary spreadsheet."
  )
  src <- frozen_sources("nmd_exons", sources_dir)
  raw <- readxl::read_excel(src[[1L]], sheet = 1L)
  # `GeneID` (StringTie) would otherwise collide with `gene_id` (Ensembl).
  names(raw)[names(raw) == "GeneID"] <- "stringtie_gene_id"
  names(raw) <- to_snake_case(names(raw))
  ids <- function(x) {
    vapply(
      strsplit(as.character(x), "[,; ]+"),
      function(v) paste(strip_id_version(v[nzchar(v)]), collapse = ";"),
      character(1L)
    )
  }
  res <- tibble::tibble(
    seqid = as.character(raw$chr),
    start_pos = as.integer(raw$exon_start_0base) + 1L,
    end_pos = as.integer(raw$exon_end),
    strand = as.character(raw$strand),
    upstream_exon_start = as.integer(raw$upstream_es) + 1L,
    upstream_exon_end = as.integer(raw$upstream_ee),
    downstream_exon_start = as.integer(raw$downstream_es) + 1L,
    downstream_exon_end = as.integer(raw$downstream_ee),
    nmd_transcript_ids = ids(raw$nmd_isoforms),
    non_nmd_transcript_ids = ids(raw$non_nmd_isoforms),
    gene_id = strip_id_version(as.character(raw$gene_id)),
    gene_name = as.character(raw$gene_name),
    source = "Karousis et al. 2021 Genome Biology, Supplementary Table 1"
  )
  res <- dplyr::distinct(res)
  path <- frozen_target("nmd_exons", output_path)
  write_parquet(res, path)
  report_frozen(path, nrow(res))
}

#' @rdname build_frozen_data
#' @export
build_frozen_polyasite <- function(
  output_path = ":cache:",
  sources_dir = NULL
) {
  src <- frozen_sources("polyasite", sources_dir)
  bed <- read_tabular_lines(src[[1L]], n_cols = 11L)
  res <- tibble::tibble(
    site_id = bed[[4L]],
    seqid = ifelse(
      grepl("^chr", bed[[1L]]),
      bed[[1L]],
      paste0("chr", bed[[1L]])
    ),
    start_pos = as.integer(bed[[2L]]) + 1L,
    end_pos = as.integer(bed[[3L]]),
    cleavage_pos = as.integer(sub("^[^:]*:([0-9]+):.*$", "\\1", bed[[4L]])),
    strand = bed[[6L]],
    tpm = suppressWarnings(as.numeric(bed[[5L]])),
    frac_samples = suppressWarnings(as.numeric(bed[[7L]])),
    n_protocols = as.integer(bed[[8L]]),
    site_class = bed[[10L]],
    signals = ifelse(bed[[11L]] %in% c("NA", ""), NA_character_, bed[[11L]])
  )
  path <- frozen_target("polyasite", output_path)
  write_parquet(res, path)
  report_frozen(path, nrow(res))
}

#' @rdname build_frozen_data
#' @export
build_frozen_gtex <- function(output_path = ":cache:", sources_dir = NULL) {
  src <- frozen_sources("gtex_transcript_tpm", sources_dir)
  gct <- as.character(fs::path_abs(src[[1L]]))
  attributes <- read_delimited(src[[2L]], null_strings = "NA")
  attributes <- attributes[
    attributes$SMAFRZE %in% "RNASEQ",
    c("SAMPID", "SMTSD")
  ]
  con <- duckdb_connect()
  header <- DBI::dbGetQuery(
    con,
    glue::glue_sql(
      "SELECT * FROM read_csv({gct}, delim = '\t', header = true, skip = 2, sample_size = 1) LIMIT 0",
      .con = con
    )
  )
  all_samples <- setdiff(colnames(header), c("transcript_id", "gene_id"))
  tissue <- attributes$SMTSD[match(all_samples, attributes$SAMPID)]
  keep <- !is.na(tissue)
  samples <- all_samples[keep]
  tissue <- tissue[keep]
  by_tissue <- split(samples, tissue)
  exprs <- vapply(
    names(by_tissue),
    function(t) {
      cols <- paste(
        DBI::dbQuoteIdentifier(con, by_tissue[[t]]),
        collapse = ", "
      )
      paste0(
        "list_aggregate([",
        cols,
        "], 'median') AS ",
        DBI::dbQuoteIdentifier(con, t)
      )
    },
    character(1L)
  )
  types <- paste0(
    "'transcript_id': 'VARCHAR', 'gene_id': 'VARCHAR', ",
    paste0("'", all_samples, "': 'DOUBLE'", collapse = ", ")
  )
  cli::cli_inform(c(
    i = paste0(
      "Computing median TPM of {length(samples)} samples in {length(by_tissue)} tissues ",
      "(this reads the 3.6 GB GCT file once)."
    )
  ))
  wide <- DBI::dbGetQuery(
    con,
    glue::glue_sql(
      "SELECT transcript_id, gene_id, {DBI::SQL(paste(exprs, collapse = ', '))}
       FROM read_csv({gct}, delim = '\t', header = true, skip = 2, quote = '',
                     columns = {{{DBI::SQL(types)}}}, parallel = false)",
      .con = con
    )
  )
  n_samples <- lengths(by_tissue)
  res <- wide |>
    tidyr_pivot_longer(
      id_cols = c("transcript_id", "gene_id"),
      names_to = "tissue",
      values_to = "median_tpm"
    )
  res$n_samples <- unname(n_samples[res$tissue])
  res$transcript_id <- strip_id_version(res$transcript_id)
  res$gene_id <- strip_id_version(res$gene_id)
  path <- frozen_target("gtex_transcript_tpm", output_path)
  write_parquet(res, path)
  report_frozen(path, nrow(res))
}

#' Minimal wide-to-long reshape without tidyr
#' @keywords internal
#' @noRd
tidyr_pivot_longer <- function(data, id_cols, names_to, values_to) {
  value_cols <- setdiff(colnames(data), id_cols)
  out <- lapply(value_cols, function(col) {
    d <- data[, id_cols, drop = FALSE]
    d[[names_to]] <- col
    d[[values_to]] <- data[[col]]
    d
  })
  tibble::as_tibble(dplyr::bind_rows(out))
}
