#' Download Reference Annotations and Sequences
#'
#' @description
#' Downloads annotation (GFF3/GTF) and sequence (transcript or genome FASTA)
#' files from the main public providers:
#'
#' * `"gencode"` (default): any GENCODE release for human (GRCh38) or mouse
#'   (GRCm39), including future releases.
#' * `"ensembl"`: any Ensembl release, for any species of any division
#'   (`vertebrates`, `plants`, `metazoa`, `fungi`, `protists`).
#' * `"ncbi"`: any NCBI assembly by accession (for example
#'   `GCF_009914755.1` for T2T-CHM13v2.0, `GCF_000001405.40` for GRCh38.p14).
#' * `"t2t"`: shorthand for `reference = "ncbi", accession = "GCF_009914755.1"`
#'   (telomere-to-telomere T2T-CHM13v2.0).
#' * `"mane"`: MANE (Matched Annotation from NCBI and EMBL-EBI) releases, the
#'   canonical transcript per protein-coding gene, in Ensembl or RefSeq
#'   coordinates and identifiers.
#' * `"fantom"`: the FANTOM CAT transcriptome assembly (Hon et al. 2017;
#'   GRCh37/hg19) at one of its levels (`version = "lv3_robust"` by
#'   default; `"lv1_raw"`, `"lv2_permissive"`, `"lv4_stringent"`), GTF
#'   only. FANTOM6 analyses use the lv3 robust catalog. The companion ID
#'   mapping, gene information and CAGE cluster tables come from
#'   [download_fantom_cat_tables()].
#'
#' Instead of assuming file names, the function lists the release directory
#' on the server and picks the matching file, so new releases and new genome
#' patch versions (for example `GRCh38.p15`) work without a package update.
#' Files already present in `output_path` are not downloaded again. Any other
#' file reachable by URL (for example the T2T consortium Liftoff GFF3 files)
#' can be cached the same way with [download_reference_url()], and HPRC
#' per-haplotype annotations with [download_hprc_annotation()].
#'
#' @param version Release version. GENCODE: the release number (`"50"`,
#'   `50`, or `"M39"` / `"39"` for mouse; default `"50"`). Ensembl: the
#'   release number (`"116"`); `NULL` (default) uses the current release.
#'   MANE: `"1.5"`; `NULL` uses the current release. NCBI/T2T: ignored (the
#'   accession identifies the assembly version).
#' @param reference One of `"gencode"`, `"ensembl"`, `"ncbi"`, `"t2t"`,
#'   `"mane"`, `"fantom"`.
#' @param organism `"human"` (default) or `"mouse"`. For Ensembl it sets the
#'   default `species` (`homo_sapiens` / `mus_musculus`).
#' @param file_type `"gff"` (default), `"gtf"`, `"fasta"` (transcript
#'   sequences; for Ensembl the `cdna.all` file), `"ncrna_fasta"` (Ensembl
#'   non-coding transcript sequences), `"genome_fasta"` (genome sequence;
#'   not provided by MANE) or `"polya"` (GENCODE manually annotated poly(A)
#'   sites and signals, `polyAs.gff3.gz`; see [polya_sites()]).
#' @param region Annotation flavour: `"all"` (default; reference chromosomes,
#'   patches, haplotypes and scaffolds), `"chr"` (reference chromosomes only)
#'   or `"primary"` (primary assembly). For `file_type = "genome_fasta"`,
#'   `"primary"` selects the primary assembly genome and the other values the
#'   full patched / toplevel genome. Ignored for NCBI and MANE.
#' @param species Ensembl species name in lower snake case
#'   (`"homo_sapiens"`, `"arabidopsis_thaliana"`). Defaults from `organism`.
#' @param division Ensembl division: `"vertebrates"` (default, the main
#'   Ensembl site), `"plants"`, `"metazoa"`, `"fungi"` or `"protists"`.
#' @param accession NCBI assembly accession (`GCF_...` or `GCA_...`),
#'   required for `reference = "ncbi"`.
#' @param mane_source Identifier space of MANE files: `"ensembl"` (default)
#'   or `"refseq"`.
#' @param output_path Directory where the file is saved. Defaults to
#'   `":cache:"`, the package cache (see [get_isoformic_cache()]).
#' @param timeout_limit Maximum time in seconds for the download. Takes
#'   precedence over `options("timeout")`. Defaults to `3600`.
#' @param method Method used by [utils::download.file()]. Defaults to
#'   `"auto"`.
#'
#' @returns The full path to the downloaded (or already present) file.
#'
#' @section Errors:
#' The function fails early with a specific condition class:
#' * `isoformic_download_reference_no_internet_error`: no connection.
#' * `isoformic_download_reference_release_not_found_error`: the release
#'   (or species, assembly) directory does not exist; for GENCODE the
#'   message lists the latest available releases.
#' * `isoformic_download_reference_file_not_found_error`: the release exists
#'   but has no file of the requested type; the message lists the candidate
#'   files found.
#' * `isoformic_download_reference_error`: the transfer failed or produced
#'   an empty file (partial files are removed).
#'
#' @section Matching annotation and sequences:
#' Salmon indices must be built from the transcript FASTA of the same
#' provider and release as the annotation directory, so that transcript IDs
#' match: GENCODE FASTA headers are `|`-separated (handled by
#' `ignore_after_bar`), Ensembl `cdna.all` headers carry the versioned
#' transcript ID first, NCBI `rna.fna` headers carry the RefSeq accession,
#' and MANE provides both identifier spaces. For lifted annotations (T2T
#' consortium Liftoff, HPRC CAT) extract the FASTA from the assembly with
#' `gffread -w transcripts.fa -g assembly.fa annotation.gff3`.
#'
#' @examples
#' \dontrun{
#' # GENCODE
#' gff_file <- download_reference(version = "50", file_type = "gff")
#' gtf_file <- download_reference(version = "50", file_type = "gtf", region = "chr")
#' gff_file_mouse <- download_reference(version = "M39", organism = "mouse")
#' genome_file <- download_reference(version = "50", file_type = "genome_fasta")
#'
#' # Ensembl: current human release, Arabidopsis from Ensembl Plants
#' ens_gff <- download_reference(reference = "ensembl")
#' ens_cdna <- download_reference(reference = "ensembl", file_type = "fasta")
#' ath_gff <- download_reference(
#'   reference = "ensembl", division = "plants", species = "arabidopsis_thaliana"
#' )
#'
#' # NCBI assemblies by accession, and the T2T-CHM13v2.0 shorthand
#' t2t_gff <- download_reference(reference = "t2t", file_type = "gff")
#' t2t_fasta <- download_reference(reference = "t2t", file_type = "fasta")
#' grch38_gff <- download_reference(reference = "ncbi", accession = "GCF_000001405.40")
#'
#' # MANE canonical transcripts (Ensembl identifiers)
#' mane_gff <- download_reference(reference = "mane", file_type = "gff")
#' c(gff_file, gtf_file, gff_file_mouse, genome_file, ens_gff, ens_cdna,
#'   ath_gff, t2t_gff, t2t_fasta, grch38_gff, mane_gff)
#' }
#'
#' @export
download_reference <- function(
  version = NULL,
  reference = c("gencode", "ensembl", "ncbi", "t2t", "mane", "fantom"),
  organism = c("human", "mouse"),
  file_type = c("gff", "gtf", "fasta", "ncrna_fasta", "genome_fasta", "polya"),
  region = c("all", "chr", "primary"),
  species = NULL,
  division = c("vertebrates", "plants", "metazoa", "fungi", "protists"),
  accession = NULL,
  mane_source = c("ensembl", "refseq"),
  output_path = ":cache:",
  timeout_limit = 3600,
  method = "auto"
) {
  reference <- tolower(reference)
  reference <- rlang::arg_match(
    reference,
    values = c("gencode", "ensembl", "ncbi", "t2t", "mane", "fantom")
  )
  organism <- tolower(organism)
  organism <- rlang::arg_match(organism)
  file_type <- tolower(file_type)
  file_type <- rlang::arg_match(file_type)
  region <- tolower(region)
  region <- rlang::arg_match(region)
  division <- tolower(division)
  division <- rlang::arg_match(division)
  mane_source <- tolower(mane_source)
  mane_source <- rlang::arg_match(mane_source)
  if (identical(file_type, "polya") && !identical(reference, "gencode")) {
    cli::cli_abort(
      message = c(
        x = "Poly(A) site annotations ({.val polya}) are only distributed by GENCODE.",
        i = "Use {.code reference = \"gencode\"}, or the PolyASite atlas ({.code fetch_frozen_data(\"polyasite\")})."
      ),
      class = "isoformic_download_reference_file_not_found_error"
    )
  }
  check_internet()

  if (identical(output_path, ":cache:")) {
    output_path <- get_isoformic_cache()
  }
  check_string(output_path)

  target <- switch(
    reference,
    gencode = resolve_gencode_file(version, organism, file_type, region),
    ensembl = resolve_ensembl_file(
      version = version,
      species = species %||%
        switch(organism, human = "homo_sapiens", mouse = "mus_musculus"),
      division = division,
      file_type = file_type,
      region = region
    ),
    ncbi = resolve_ncbi_file(accession, file_type),
    t2t = {
      check_t2t_request(version, organism)
      resolve_ncbi_file("GCF_009914755.1", file_type)
    },
    mane = resolve_mane_file(version, mane_source, file_type),
    fantom = resolve_fantom_file(version, file_type)
  )
  download_reference_file(
    download_url = target$url,
    full_output_path = as.character(fs::path(output_path, target$file)),
    timeout_limit = timeout_limit,
    method = method
  )
}

#' Download the FANTOM CAT Companion Tables
#'
#' Fetches, for one FANTOM CAT level, the ID mapping table (transcript to
#' CAGE cluster), the gene information table and the CAGE cluster BED file,
#' which [add_fantom_cat_tables()] joins into an annotation directory.
#'
#' @param level FANTOM CAT level: `"lv3_robust"` (default), `"lv1_raw"`,
#'   `"lv2_permissive"` or `"lv4_stringent"`.
#' @param tables Which tables to download.
#' @inheritParams download_reference
#'
#' @returns A named character vector of file paths (`id_mapping`,
#'   `gene_info`, `cage_clusters`).
#'
#' @examples
#' \dontrun{
#' tables <- download_fantom_cat_tables("lv3_robust")
#' tables
#' }
#'
#' @export
download_fantom_cat_tables <- function(
  level = c("lv3_robust", "lv4_stringent", "lv2_permissive", "lv1_raw"),
  tables = c("id_mapping", "gene_info", "cage_clusters"),
  output_path = ":cache:",
  timeout_limit = 3600,
  method = "auto"
) {
  level <- rlang::arg_match(level)
  tables <- rlang::arg_match(tables, multiple = TRUE)
  base <- fantom_cat_base_url(level)
  suffix <- c(
    id_mapping = ".info_table.ID_mapping.tsv.gz",
    gene_info = ".info_table.gene.tsv.gz",
    cage_clusters = ".CAGE_cluster.bed.gz"
  )
  files <- vapply(
    tables,
    function(t) {
      download_reference_url(
        url = paste0(base, "FANTOM_CAT.", level, suffix[[t]]),
        output_path = output_path,
        timeout_limit = timeout_limit,
        method = method
      )
    },
    character(1L)
  )
  names(files) <- tables
  files
}

fantom_cat_base_url <- function(level) {
  paste0(
    "https://fantom.gsc.riken.jp/5/suppl/Hon_et_al_2016/data/assembly/",
    level,
    "/"
  )
}

#' Resolve a FANTOM CAT GTF file
#' @keywords internal
#' @noRd
resolve_fantom_file <- function(version, file_type) {
  level <- version %||% "lv3_robust"
  levels <- c("lv1_raw", "lv2_permissive", "lv3_robust", "lv4_stringent")
  if (!level %in% levels) {
    cli::cli_abort(
      message = c(
        x = "FANTOM CAT level {.val {level}} is not valid.",
        i = "Use one of {.val {levels}} as {.arg version}."
      ),
      class = "isoformic_download_reference_unsupported_version_error"
    )
  }
  if (!identical(file_type, "gtf")) {
    cli::cli_abort(
      message = c(
        x = "FANTOM CAT is distributed as GTF only; use {.code file_type = \"gtf\"}.",
        i = "Transcript sequences must be extracted from the hg19 genome with gffread."
      ),
      class = "isoformic_download_reference_file_not_found_error"
    )
  }
  dir_url <- fantom_cat_base_url(level)
  files <- list_remote_dir(dir_url)
  if (is.null(files)) {
    abort_release_not_found(paste0("FANTOM CAT ", level), dir_url)
  }
  pick_remote_file(
    files = files,
    pattern = paste0("^FANTOM_CAT\\.", level, "\\.gtf\\.gz$"),
    kind_pattern = kind_pattern("gtf"),
    what = paste0("FANTOM CAT ", level, " GTF"),
    dir_url = dir_url
  )
}

#' Download Any Reference File by URL
#'
#' Caches a file from an arbitrary URL in the same way as
#' [download_reference()] (no re-download when present, partial files
#' removed on failure). Useful for annotations without a provider layout,
#' such as the T2T consortium Liftoff files or CHESS releases.
#'
#' @param url URL of the file.
#' @param output_path Directory where the file is saved; `":cache:"` (the
#'   default) is the package cache.
#' @param file_name File name to use; defaults to the last part of the URL.
#' @inheritParams download_reference
#'
#' @returns The full path to the file.
#'
#' @examples
#' \dontrun{
#' liftoff <- download_reference_url(paste0(
#'   "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/",
#'   "assemblies/annotation/chm13v2.0_RefSeq_Liftoff_v5.2.gff3.gz"
#' ))
#' annot <- prepare_isoformic_annotation(liftoff)
#' annot
#' }
#'
#' @export
download_reference_url <- function(
  url,
  output_path = ":cache:",
  file_name = NULL,
  timeout_limit = 3600,
  method = "auto"
) {
  check_string(url)
  if (!grepl("^(https?|ftp)://", url)) {
    cli::cli_abort(
      message = c(
        x = "{.arg url} must start with http://, https:// or ftp://."
      ),
      class = "isoformic_download_reference_error"
    )
  }
  check_internet()
  if (identical(output_path, ":cache:")) {
    output_path <- get_isoformic_cache()
  }
  file_name <- file_name %||% fs::path_file(sub("[?#].*$", "", url))
  check_string(file_name)
  download_reference_file(
    download_url = url,
    full_output_path = as.character(fs::path(output_path, file_name)),
    timeout_limit = timeout_limit,
    method = method
  )
}

#' Download HPRC Per-Haplotype CAT Annotations
#'
#' Downloads the Human Pangenome Reference Consortium year-1 CAT/Liftoff
#' gene annotations of an assembled sample, one GFF3 per haplotype, lifted
#' from GENCODE v38 onto the sample's own assembly in CHM13 or GRCh38
#' coordinates. The files keep GENCODE identifiers in `source_gene` /
#' `source_transcript`, so [merge_annotations()] and [haplotype_summary()]
#' can compare haplotypes.
#'
#' @param sample HPRC sample ID, for example `"HG00438"` or `"HG002"`.
#' @param haplotype Haplotype numbers to download, `1`, `2` or `c(1, 2)`
#'   (default).
#' @param assembly Coordinate space of the CAT run: `"chm13"` (default) or
#'   `"hg38"`.
#' @inheritParams download_reference
#'
#' @returns A named character vector of file paths (names `hap1`, `hap2`),
#'   ready for [merge_annotations()].
#'
#' @examples
#' \dontrun{
#' files <- download_hprc_annotation("HG00438")
#' merged <- merge_annotations(files, annot_name = "HG00438_chm13")
#' merged
#' }
#'
#' @export
download_hprc_annotation <- function(
  sample,
  haplotype = c(1L, 2L),
  assembly = c("chm13", "hg38"),
  output_path = ":cache:",
  timeout_limit = 3600,
  method = "auto"
) {
  check_string(sample)
  assembly <- rlang::arg_match(assembly)
  haplotype <- as.integer(haplotype)
  if (!length(haplotype) || !all(haplotype %in% c(1L, 2L))) {
    cli::cli_abort(
      message = c(x = "{.arg haplotype} must be 1, 2 or c(1, 2)."),
      class = "isoformic_invalid_haplotypes"
    )
  }
  files <- vapply(
    haplotype,
    function(h) {
      download_reference_url(
        url = hprc_cat_url(sample, h, assembly),
        output_path = output_path,
        timeout_limit = timeout_limit,
        method = method
      )
    },
    character(1L)
  )
  names(files) <- paste0("hap", haplotype)
  files
}

hprc_cat_url <- function(sample, haplotype, assembly) {
  paste0(
    "https://s3-us-west-2.amazonaws.com/human-pangenomics/submissions/",
    "FC7E9302-5A35-4A4F-9D4A-B3E957118E36--Y1_CAT_ANNOTATION_GENCODE38/",
    sample,
    ".",
    haplotype,
    ".",
    assembly,
    ".gff3.gz"
  )
}

#' Validate the T2T-CHM13 shorthand (human only, version 2.0)
#' @keywords internal
#' @noRd
check_t2t_request <- function(version, organism) {
  if (!identical(organism, "human")) {
    cli::cli_abort(
      message = c(
        x = "T2T-CHM13 references are only available for {.val human}."
      ),
      class = "isoformic_download_reference_unsupported_version_error"
    )
  }
  if (!is.null(version) && !as.character(version) %in% c("2.0", "2", "v2.0")) {
    cli::cli_abort(
      message = c(
        x = "T2T-CHM13 version {.val {version}} is not supported.",
        i = paste0(
          "Only T2T-CHM13v2.0 ({.val GCF_009914755.1}) is available; ",
          "use {.code reference = \"ncbi\"} with another accession."
        )
      ),
      class = "isoformic_download_reference_unsupported_version_error"
    )
  }
  invisible(TRUE)
}

# Resolution of remote files --------------------------------------------------

gencode_base_url <- function(organism) {
  switch(
    organism,
    human = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/",
    mouse = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/"
  )
}

#' Normalise a GENCODE version string (`50`, `"50"`, `"M39"`, `"39"`)
#' @keywords internal
#' @noRd
normalize_gencode_version <- function(version, organism) {
  if (is.null(version)) {
    version <- if (identical(organism, "mouse")) "M39" else "50"
  }
  version <- toupper(trimws(as.character(version)))
  version <- sub("^(RELEASE_|V)", "", version)
  if (!grepl("^M?[0-9]+$", version)) {
    cli::cli_abort(
      message = c(
        x = "{.arg version} must be a GENCODE release number such as {.val 50} or {.val M39}.",
        i = "Received: {.val {version}}."
      ),
      class = "isoformic_download_reference_unsupported_version_error"
    )
  }
  if (identical(organism, "mouse") && !startsWith(version, "M")) {
    version <- paste0("M", version)
  }
  if (identical(organism, "human") && startsWith(version, "M")) {
    cli::cli_abort(
      message = c(
        x = "Version {.val {version}} looks like a mouse release but {.arg organism} is {.val human}.",
        i = "Use {.code organism = \"mouse\"} for GENCODE M releases."
      ),
      class = "isoformic_download_reference_unsupported_version_error"
    )
  }
  version
}

#' Abort because a remote directory could not be listed
#' @keywords internal
#' @noRd
abort_release_not_found <- function(what, url, releases = character()) {
  cli::cli_abort(
    message = c(
      x = "{what} was not found at {.url {url}}.",
      i = if (length(releases)) {
        "Latest available releases: {.val {utils::tail(releases, 5)}}."
      } else {
        "Check the name/version and the connection."
      }
    ),
    class = "isoformic_download_reference_release_not_found_error"
  )
}

#' Pick the file matching `pattern` in a listing or abort with the candidates
#' @keywords internal
#' @noRd
pick_remote_file <- function(files, pattern, kind_pattern, what, dir_url) {
  matches <- grep(pattern, files, value = TRUE)
  if (!length(matches)) {
    candidates <- grep(kind_pattern, files, value = TRUE)
    cli::cli_abort(
      message = c(
        x = "No file for {what} found in {.url {dir_url}}.",
        i = "Expected a file matching {.code {pattern}}.",
        i = if (length(candidates)) {
          "Files of that kind in the directory: {.val {candidates}}."
        } else {
          "The directory has no file of that kind."
        }
      ),
      class = "isoformic_download_reference_file_not_found_error"
    )
  }
  file <- matches[[1L]]
  list(url = paste0(dir_url, file), file = file)
}

kind_pattern <- function(file_type) {
  switch(
    file_type,
    gff = "\\.gff3?\\.gz$",
    gtf = "\\.gtf\\.gz$",
    fasta = "(transcripts|cdna|rna)[^/]*\\.(fa|fna)\\.gz$",
    ncrna_fasta = "ncrna[^/]*\\.fa\\.gz$",
    genome_fasta = "(genome|dna|genomic)[^/]*\\.(fa|fna)\\.gz$",
    polya = "polyAs\\.gff3\\.gz$"
  )
}

#' Find the GENCODE file for a release by listing the release directory
#' @keywords internal
#' @noRd
resolve_gencode_file <- function(version, organism, file_type, region) {
  version <- normalize_gencode_version(version, organism)
  if (identical(file_type, "ncrna_fasta")) {
    cli::cli_abort(
      message = c(
        x = "GENCODE has no separate ncRNA FASTA; {.val fasta} contains all transcripts.",
        i = "Use {.code file_type = \"fasta\"}."
      ),
      class = "isoformic_download_reference_file_not_found_error"
    )
  }
  base_url <- gencode_base_url(organism)
  release_url <- paste0(base_url, "release_", version, "/")
  files <- list_remote_dir(release_url)
  if (is.null(files)) {
    available <- list_remote_dir(base_url)
    releases <- sort_releases(grep(
      "^release_",
      available %||% character(),
      value = TRUE
    ))
    abort_release_not_found(
      paste0("GENCODE release ", version, " for ", organism),
      release_url,
      releases
    )
  }
  pick_remote_file(
    files = files,
    pattern = gencode_file_pattern(version, file_type, region),
    kind_pattern = kind_pattern(file_type),
    what = paste0("GENCODE ", version, " ", file_type, " (", region, ")"),
    dir_url = release_url
  )
}

#' Regular expression of the expected GENCODE file name
#' @keywords internal
#' @noRd
gencode_file_pattern <- function(version, file_type, region) {
  v <- gsub(".", "\\.", version, fixed = TRUE)
  annot_prefix <- switch(
    region,
    all = "chr_patch_hapl_scaff\\.",
    chr = "",
    primary = "primary_assembly\\."
  )
  switch(
    file_type,
    gff = paste0(
      "^gencode\\.v",
      v,
      "\\.",
      annot_prefix,
      "annotation\\.gff3\\.gz$"
    ),
    gtf = paste0(
      "^gencode\\.v",
      v,
      "\\.",
      annot_prefix,
      "annotation\\.gtf\\.gz$"
    ),
    fasta = paste0("^gencode\\.v", v, "\\.transcripts\\.fa\\.gz$"),
    polya = paste0("^gencode\\.v", v, "\\.polyAs\\.gff3\\.gz$"),
    genome_fasta = if (identical(region, "primary")) {
      "^[A-Za-z0-9]+\\.primary_assembly\\.genome\\.fa\\.gz$"
    } else {
      "^[A-Za-z0-9]+(\\.p[0-9]+)?\\.genome\\.fa\\.gz$"
    }
  )
}

#' Sort GENCODE release directory names by release number
#'
#' Legacy names such as `release_3d` carry a letter suffix; those sort by
#' number first and suffix second.
#'
#' @keywords internal
#' @noRd
sort_releases <- function(x) {
  x <- grep("^release_M?[0-9]+[a-z]?$", x, value = TRUE)
  num <- as.numeric(sub("^release_M?([0-9]+)[a-z]?$", "\\1", x))
  suffix <- sub("^release_M?[0-9]+", "", x)
  x[order(num, suffix)]
}

#' Find an Ensembl file (any division, species and release)
#' @keywords internal
#' @noRd
resolve_ensembl_file <- function(
  version,
  species,
  division,
  file_type,
  region
) {
  check_string(species)
  species <- tolower(species)
  kind_dir <- switch(
    file_type,
    gff = "gff3",
    gtf = "gtf",
    fasta = "fasta",
    ncrna_fasta = "fasta",
    genome_fasta = "fasta"
  )
  sub_dir <- switch(
    file_type,
    fasta = "cdna/",
    ncrna_fasta = "ncrna/",
    genome_fasta = "dna/",
    ""
  )
  if (identical(division, "vertebrates")) {
    base <- "https://ftp.ensembl.org/pub/"
    release_dir <- if (is.null(version)) {
      paste0("current_", kind_dir, "/")
    } else {
      paste0("release-", as.character(version), "/", kind_dir, "/")
    }
  } else {
    base <- paste0("https://ftp.ebi.ac.uk/ensemblgenomes/pub/", division, "/")
    release_dir <- if (is.null(version)) {
      paste0("current/", kind_dir, "/")
    } else {
      paste0("release-", as.character(version), "/", kind_dir, "/")
    }
  }
  dir_url <- paste0(base, release_dir, species, "/", sub_dir)
  files <- list_remote_dir(dir_url)
  if (is.null(files)) {
    abort_release_not_found(
      paste0(
        "Ensembl ",
        division,
        " ",
        if (is.null(version)) "current release" else paste("release", version),
        " for ",
        species
      ),
      dir_url
    )
  }
  pattern <- ensembl_file_pattern(file_type, region)
  files <- files[
    !grepl(
      "\\.chromosome\\.|\\.abinitio\\.|\\.primary_assembly\\.gff|\\.dna_(rm|sm)\\.",
      files
    )
  ]
  res <- pick_remote_file(
    files = files,
    pattern = pattern,
    kind_pattern = kind_pattern(file_type),
    what = paste0("Ensembl ", species, " ", file_type, " (", region, ")"),
    dir_url = dir_url
  )
  # Prefer the full chr_patch_hapl_scaff file for `region = "all"` when present
  if (identical(region, "all") && file_type %in% c("gff", "gtf")) {
    full <- grep(
      "chr_patch_hapl_scaff",
      grep(pattern, files, value = TRUE),
      value = TRUE,
      fixed = TRUE
    )
    if (length(full)) {
      res <- list(url = paste0(dir_url, full[[1L]]), file = full[[1L]])
    }
  }
  res
}

#' File name pattern of an Ensembl file type and region
#' @keywords internal
#' @noRd
ensembl_file_pattern <- function(file_type, region) {
  species_re <- "^[A-Z][a-z0-9_]+\\.[^/]+?"
  switch(
    file_type,
    gff = switch(
      region,
      all = paste0(
        species_re,
        "\\.[0-9]+(\\.chr_patch_hapl_scaff)?\\.gff3\\.gz$"
      ),
      chr = paste0(species_re, "\\.[0-9]+\\.chr\\.gff3\\.gz$"),
      primary = paste0(species_re, "\\.[0-9]+\\.gff3\\.gz$")
    ),
    gtf = switch(
      region,
      all = paste0(
        species_re,
        "\\.[0-9]+(\\.chr_patch_hapl_scaff)?\\.gtf\\.gz$"
      ),
      chr = paste0(species_re, "\\.[0-9]+\\.chr\\.gtf\\.gz$"),
      primary = paste0(species_re, "\\.[0-9]+\\.gtf\\.gz$")
    ),
    fasta = paste0(species_re, "\\.cdna\\.all\\.fa\\.gz$"),
    ncrna_fasta = paste0(species_re, "\\.ncrna\\.fa\\.gz$"),
    genome_fasta = if (identical(region, "primary")) {
      paste0(species_re, "\\.dna\\.primary_assembly\\.fa\\.gz$")
    } else {
      paste0(species_re, "\\.dna\\.toplevel\\.fa\\.gz$")
    }
  )
}

#' Find a file of an NCBI assembly by accession
#' @keywords internal
#' @noRd
resolve_ncbi_file <- function(accession, file_type) {
  check_string(accession)
  accession <- toupper(trimws(accession))
  if (!grepl("^GC[AF]_[0-9]{9}\\.[0-9]+$", accession)) {
    cli::cli_abort(
      message = c(
        x = "{.arg accession} must be an NCBI assembly accession such as {.val GCF_009914755.1}."
      ),
      class = "isoformic_download_reference_unsupported_version_error"
    )
  }
  if (identical(file_type, "ncrna_fasta")) {
    cli::cli_abort(
      message = c(
        x = "NCBI assemblies have no separate ncRNA FASTA; {.val fasta} ({.file _rna.fna.gz}) contains all transcripts."
      ),
      class = "isoformic_download_reference_file_not_found_error"
    )
  }
  digits <- sub("^GC[AF]_([0-9]{9})\\..*$", "\\1", accession)
  prefix <- sub("_.*$", "", accession)
  parent_url <- paste0(
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/",
    prefix,
    "/",
    substr(digits, 1, 3),
    "/",
    substr(digits, 4, 6),
    "/",
    substr(digits, 7, 9),
    "/"
  )
  dirs <- list_remote_dir(parent_url)
  assembly_dir <- grep(
    paste0("^", gsub(".", "\\.", accession, fixed = TRUE), "_"),
    dirs %||% character(),
    value = TRUE
  )
  if (!length(assembly_dir)) {
    abort_release_not_found(
      paste0("NCBI assembly ", accession),
      parent_url,
      grep("^GC", dirs %||% character(), value = TRUE)
    )
  }
  assembly_dir <- assembly_dir[[1L]]
  dir_url <- paste0(parent_url, assembly_dir, "/")
  files <- list_remote_dir(dir_url)
  if (is.null(files)) {
    abort_release_not_found(paste0("NCBI assembly ", accession), dir_url)
  }
  suffix <- switch(
    file_type,
    gff = "_genomic\\.gff\\.gz$",
    gtf = "_genomic\\.gtf\\.gz$",
    fasta = "_rna\\.fna\\.gz$",
    genome_fasta = "_genomic\\.fna\\.gz$"
  )
  pick_remote_file(
    files = files,
    pattern = paste0("^", gsub(".", "\\.", assembly_dir, fixed = TRUE), suffix),
    kind_pattern = kind_pattern(file_type),
    what = paste0("NCBI ", accession, " ", file_type),
    dir_url = dir_url
  )
}

#' Find a MANE release file
#' @keywords internal
#' @noRd
resolve_mane_file <- function(version, mane_source, file_type) {
  if (identical(file_type, "genome_fasta")) {
    cli::cli_abort(
      message = c(
        x = "MANE does not distribute genome sequences.",
        i = "Use {.code reference = \"gencode\"} or {.code reference = \"ncbi\"} for the GRCh38 genome."
      ),
      class = "isoformic_download_reference_file_not_found_error"
    )
  }
  if (identical(file_type, "ncrna_fasta")) {
    cli::cli_abort(
      message = c(
        x = "MANE only covers protein-coding transcripts; use {.code file_type = \"fasta\"}."
      ),
      class = "isoformic_download_reference_file_not_found_error"
    )
  }
  base <- "https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/"
  release_dir <- if (is.null(version)) {
    "current/"
  } else {
    paste0("release_", sub("^v", "", as.character(version)), "/")
  }
  dir_url <- paste0(base, release_dir)
  files <- list_remote_dir(dir_url)
  if (is.null(files)) {
    available <- list_remote_dir(base)
    releases <- grep("^release_", available %||% character(), value = TRUE)
    abort_release_not_found(
      paste0(
        "MANE ",
        if (is.null(version)) "current release" else paste("release", version)
      ),
      dir_url,
      releases
    )
  }
  suffix <- switch(
    file_type,
    gff = "_genomic\\.gff\\.gz$",
    gtf = "_genomic\\.gtf\\.gz$",
    fasta = "_rna\\.fna\\.gz$"
  )
  pick_remote_file(
    files = files,
    pattern = paste0("^MANE\\.GRCh38\\.v[0-9.]+\\.", mane_source, suffix),
    kind_pattern = kind_pattern(file_type),
    what = paste0("MANE ", mane_source, " ", file_type),
    dir_url = dir_url
  )
}

# Network helpers -------------------------------------------------------------

#' Abort when there is no internet connection (only when `curl` is installed)
#' @keywords internal
#' @noRd
check_internet <- function() {
  if (
    requireNamespace("curl", quietly = TRUE) && !isTRUE(curl::has_internet())
  ) {
    cli::cli_abort(
      message = c(x = "No internet connection available."),
      class = "isoformic_download_reference_no_internet_error"
    )
  }
  invisible(TRUE)
}

#' List the entries of a remote directory served as an HTML index
#'
#' Returns the file and directory names (without trailing slashes) found in
#' the `href` attributes of the index page, or `NULL` when the directory
#' cannot be read (for example a 404 for a release that does not exist).
#'
#' @keywords internal
#' @noRd
list_remote_dir <- function(url, timeout = 60) {
  html <- rlang::try_fetch(
    withr::with_options(
      list(timeout = max(timeout, getOption("timeout", 60))),
      {
        con <- url(url, open = "rb")
        withr::defer(close(con))
        readLines(con, warn = FALSE, encoding = "UTF-8")
      }
    ),
    error = function(cnd) NULL,
    warning = function(cnd) NULL
  )
  if (is.null(html)) {
    return(NULL)
  }
  hrefs <- regmatches(html, gregexpr('href="[^"]+"', html))
  hrefs <- unlist(hrefs, use.names = FALSE)
  hrefs <- sub('^href="', "", sub('"$', "", hrefs))
  hrefs <- hrefs[!grepl("^(\\?|/|\\.\\./|#|http)", hrefs)]
  hrefs <- sub("/$", "", hrefs)
  unique(hrefs[nzchar(hrefs)])
}

#' Download a reference file unless it is already cached
#'
#' Downloads to a temporary `.part` file next to the destination and renames
#' it on success, so an interrupted transfer never leaves a truncated file
#' behind that would be mistaken for a complete download.
#'
#' @keywords internal
#' @noRd
download_reference_file <- function(
  download_url,
  full_output_path,
  timeout_limit = 3600,
  method = "auto"
) {
  fs::dir_create(fs::path_dir(full_output_path))
  if (isTRUE(fs::file_exists(full_output_path))) {
    cli::cli_inform(c(i = "{.path {full_output_path}} already exists."))
    return(full_output_path)
  }
  part_path <- paste0(full_output_path, ".part")
  withr::defer(if (fs::file_exists(part_path)) fs::file_delete(part_path))
  cli::cli_inform(c(i = "Downloading {.url {download_url}}."))
  status <- withr::with_options(
    list(timeout = max(unlist(options("timeout")), timeout_limit)),
    rlang::try_fetch(
      utils::download.file(
        url = download_url,
        destfile = part_path,
        method = method,
        mode = "wb",
        quiet = !rlang::is_interactive()
      ),
      error = function(cnd) {
        cli::cli_abort(
          message = c(
            x = "Failed to download {.url {download_url}}.",
            i = "Original error: {conditionMessage(cnd)}"
          ),
          class = "isoformic_download_reference_error",
          parent = cnd
        )
      }
    )
  )
  size <- if (fs::file_exists(part_path)) {
    as.numeric(fs::file_size(part_path))
  } else {
    0
  }
  if (!identical(as.integer(status), 0L) || !isTRUE(size > 0)) {
    cli::cli_abort(
      message = c(
        x = "Download of {.url {download_url}} did not complete (status {status}, {size} bytes).",
        i = "Check the URL in a browser; the file may have moved."
      ),
      class = "isoformic_download_reference_error"
    )
  }
  fs::file_move(part_path, full_output_path)
  cli::cli_inform(c(v = "{.path {full_output_path}} successfully downloaded."))
  full_output_path
}
