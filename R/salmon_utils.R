#' Create the conda environment that holds Salmon
#'
#' Installs Salmon from Bioconda into a conda environment managed by
#' `condathis`. Salmon 2 is a Rust rewrite with a new index format: indices
#' built with Salmon 1.x cannot be read, and vice versa, so keep one
#' environment per Salmon major version (the final C++ release, 1.12.0, is
#' available as the `salmon-cpp` package).
#'
#' @param env_name Name of the conda environment (default `"salmon-env"`).
#' @param version Salmon version to install (default `"2.8.0"`).
#' @param force Logical, recreate the environment if it exists.
#'
#' @keywords internal
#' @noRd
create_salmon_env <- function(
  env_name = "salmon-env",
  version = "2.8.0",
  force = FALSE
) {
  rlang::check_installed("condathis")
  condathis::create_env(
    paste0("bioconda::salmon==", version),
    env_name = env_name,
    additional_channels = "bioconda",
    verbose = "output",
    overwrite = force
  )
}

check_env_installed <- function(env_name = "salmon-env") {
  rlang::check_installed("condathis")
  rlang::check_required(env_name)
  if (isFALSE(condathis::env_exists(env_name))) {
    rlang::try_fetch(
      expr = {
        create_salmon_env(force = FALSE, env_name = env_name)
      },
      error = function(e) {
        cli::cli_abort(
          c(
            "The conda environment '{env_name}' could not be created or installed.",
            "i" = "Please try running {.code create_salmon_env()} manually.",
            "i" = "Original error message: {e$message}"
          ),
          class = "isoformic_salmon_env_create_error"
        )
      }
    )
  }
  return(invisible(TRUE))
}

#' Build Salmon Index
#'
#' Create a Salmon index from a reference transcriptome FASTA file.
#'
#' @details
#' The index must be built with the same Salmon major version used for
#' `salmon_quant()`: Salmon 2 (the Rust rewrite) cannot read indices built by
#' Salmon 1.x and rejects them with an error. Rebuild the index after
#' upgrading.
#'
#' Transcripts with identical sequences are all kept (`--keepDuplicates`),
#' so that every annotated transcript is quantified.
#'
#' With `decoy_fasta`, the index is decoy-aware, as recommended by the Salmon
#' authors: the function writes, in a temporary directory, a `gentrome` FASTA
#' (the transcripts followed by the decoy sequences) and a `decoys.txt` file
#' with the names of the decoy sequences (the first word of each FASTA
#' header), passes them to `salmon index` as `--transcripts` and `--decoys`,
#' and removes them afterwards. Reads that align better to a decoy than to a
#' transcript are then discarded instead of being assigned to a transcript.
#' The `gentrome` is gzipped when `decoy_fasta` is (the decoy file is
#' appended byte by byte, and concatenated gzip streams are valid gzip), so a
#' gzipped genome is not decompressed to disk; the transcript FASTA, plain or
#' gzipped, is re-encoded to match. A decoy compressed with xz or bzip2 is
#' decompressed into a plain `gentrome`. Decoy names that are also
#' transcript names are an error.
#'
#' @param fasta_path Path to the reference transcriptome FASTA file.
#' @param index_path Directory path to save the Salmon index (default is
#' "salmon_index").
#' @param kmer_len K-mer length for the index (default is 31).
#' @param num_threads Number of threads to use (default is 2).
#' @param env_name Name of the conda environment with Salmon installed
#' (default is "salmon-env").
#' @param is_gencode Logical indicating if the FASTA is from GENCODE
#' (default is `FALSE`); the transcript names are then truncated at the
#' first `|` so that they match the annotation.
#' @param decoy_fasta Optional path to a FASTA file, plain or compressed
#' with gzip, xz or bzip2, with
#' the decoy sequences (default is `NULL`, an index without decoys). Usually
#' the genome assembly the transcripts come from (for GENCODE, the
#' `GRCh38.primary_assembly.genome.fa.gz` of the same release). A full human
#' genome makes the index several times larger and slower to build (tens of
#' GB of memory and disk instead of a few), and needs free space in
#' [tempdir()] for a copy of the genome while the index is built.
#' @param clip_poly_a Logical indicating whether to clip poly-A tails
#' (default is `TRUE`).
#'
#' @returns processx style output list.
#'
#' @export
salmon_index <- function(
  fasta_path,
  index_path = "salmon_index",
  kmer_len = 31,
  num_threads = 2,
  env_name = "salmon-env",
  is_gencode = FALSE,
  decoy_fasta = NULL,
  clip_poly_a = TRUE
) {
  rlang::check_installed("condathis")
  rlang::check_required(fasta_path)
  if (!is.null(decoy_fasta)) {
    check_fasta_file(fasta_path)
    check_fasta_file(decoy_fasta)
  }
  check_env_installed(env_name)

  if (isTRUE(is_gencode)) {
    gencode_arg <- "--gencode"
  } else {
    gencode_arg <- NULL
  }
  decoy_args <- NULL
  if (!is.null(decoy_fasta)) {
    decoy_dir <- tempfile("salmon_decoys_")
    dir.create(decoy_dir)
    on.exit(unlink(decoy_dir, recursive = TRUE), add = TRUE)
    decoy_files <- write_salmon_gentrome(fasta_path, decoy_fasta, decoy_dir)
    fasta_path <- decoy_files$gentrome
    decoy_args <- c("--decoys", decoy_files$decoys)
  }

  if (isFALSE(clip_poly_a)) {
    no_clip_arg <- "--no-clip"
  } else {
    no_clip_arg <- NULL
  }
  condathis::run(
    "salmon",
    "index",
    "--transcripts",
    fasta_path,
    "--kmerLen",
    kmer_len,
    "--index",
    index_path,
    gencode_arg,
    "--keepDuplicates",
    "--threads",
    num_threads,
    decoy_args,
    no_clip_arg,
    env_name = env_name,
    verbose = "output"
  )
}

#' Write the gentrome FASTA and the decoy names for a decoy-aware index
#'
#' @param fasta_path Path to the transcript FASTA (plain or compressed).
#' @param decoy_fasta Path to the decoy FASTA (plain, or compressed with gzip,
#'   xz or bzip2).
#' @param out_dir Existing directory for the output files.
#'
#' @returns A list with the paths `gentrome` and `decoys`.
#'
#' @keywords internal
#' @noRd
write_salmon_gentrome <- function(fasta_path, decoy_fasta, out_dir) {
  check_fasta_file(fasta_path)
  check_fasta_file(decoy_fasta)
  decoy_is_gz <- is_gzip_file(decoy_fasta)
  # A plain or gzipped FASTA decoy with LF line endings is appended byte by
  # byte (fast, and concatenated gzip streams are valid gzip); any other file
  # (xz, bzip2, or plain or gzipped text with CRLF line endings or not
  # starting with a header) is decompressed and written line by line, because
  # appending its raw bytes would corrupt the gentrome.
  decoy_append <- is_plain_fasta_file(decoy_fasta)
  gentrome <- file.path(
    out_dir,
    if (decoy_is_gz) "gentrome.fa.gz" else "gentrome.fa"
  )
  decoys <- file.path(out_dir, "decoys.txt")

  # Transcripts first, re-encoded to the compression of the decoy file;
  # writeLines() also guarantees a final newline before the decoys.
  out_con <- if (decoy_is_gz) {
    gzfile(gentrome, open = "wb", compression = 1L)
  } else {
    file(gentrome, open = "wb")
  }
  write_ok <- FALSE
  on.exit(if (!write_ok) close(out_con), add = TRUE)
  found <- new.env(parent = emptyenv())
  found$tx <- character()
  found$decoy <- character()
  for_each_line_chunk(fasta_path, function(lines) {
    writeLines(lines, out_con)
    found$tx <- c(found$tx, fasta_header_names(lines))
  })
  if (!decoy_append) {
    for_each_line_chunk(decoy_fasta, function(lines) {
      writeLines(sub("\r$", "", lines), out_con)
    })
  }
  write_ok <- TRUE
  close(out_con)

  # A plain decoy file without a final newline still ends the gentrome
  if (decoy_append && !file.append(gentrome, decoy_fasta)) {
    cli::cli_abort(
      "Could not append {.path {decoy_fasta}} to {.path {gentrome}}.",
      class = "isoformic_gentrome_error"
    )
  }

  # Decoy names: first word of each header, without '>'
  name_con <- file(decoys, open = "wb")
  names_ok <- FALSE
  on.exit(if (!names_ok) close(name_con), add = TRUE)
  for_each_line_chunk(decoy_fasta, function(lines) {
    ids <- fasta_header_names(lines)
    if (length(ids) > 0L) {
      writeLines(ids, name_con)
      found$decoy <- c(found$decoy, ids)
    }
  })
  names_ok <- TRUE
  close(name_con)
  if (length(found$decoy) == 0L) {
    cli::cli_abort(
      "{.path {decoy_fasta}} has no FASTA header ({.code >}).",
      class = "isoformic_gentrome_error"
    )
  }
  # Salmon cannot tell a decoy from a transcript of the same name (GENCODE
  # names are compared up to the first '|', as with `is_gencode = TRUE`)
  clash <- intersect(
    sub("[|].*$", "", found$decoy),
    sub("[|].*$", "", found$tx)
  )
  if (length(clash) > 0L) {
    cli::cli_abort(
      c(
        x = "Decoy sequence names are also transcript names: {.val {utils::head(clash, 5L)}}.",
        i = "Rename the decoy sequences in {.path {decoy_fasta}}."
      ),
      class = "isoformic_gentrome_error"
    )
  }
  list(gentrome = gentrome, decoys = decoys)
}

#' First word of each FASTA header line, without `>` (and without a CR)
#' @keywords internal
#' @noRd
fasta_header_names <- function(lines) {
  headers <- lines[startsWith(lines, ">")]
  sub("^>[[:space:]]*([^[:space:]]*).*$", "\\1", headers)
}

#' Is a file a plain or gzipped FASTA file with LF line endings?
#'
#' Checks that the (decompressed, for a gzip file) content starts with a
#' header line and has no carriage return anywhere: the whole file is read,
#' in chunks of `chunk_size` bytes, because a CRLF part after the first
#' bytes would otherwise be appended raw to the gentrome. Files compressed in
#' other formats do not start with `>` and return `FALSE`.
#' @keywords internal
#' @noRd
is_plain_fasta_file <- function(path, chunk_size = 16777216L) {
  con <- if (is_gzip_file(path)) {
    gzfile(path, open = "rb")
  } else {
    file(path, open = "rb")
  }
  on.exit(close(con), add = TRUE)
  cr <- as.raw(0x0d)
  bytes <- readBin(con, what = "raw", n = chunk_size)
  if (length(bytes) == 0L || !identical(bytes[[1L]], charToRaw(">"))) {
    return(FALSE)
  }
  repeat {
    if (length(grepRaw(cr, bytes, fixed = TRUE))) {
      return(FALSE)
    }
    bytes <- readBin(con, what = "raw", n = chunk_size)
    if (length(bytes) == 0L) {
      return(TRUE)
    }
  }
}

check_fasta_file <- function(path, arg = rlang::caller_arg(path)) {
  if (!rlang::is_string(path) || !file.exists(path) || dir.exists(path)) {
    cli::cli_abort(
      "{.arg {arg}} must be the path to an existing FASTA file.",
      class = "isoformic_file_not_found"
    )
  }
  invisible(TRUE)
}

for_each_line_chunk <- function(path, fun, chunk_size = 1000000L) {
  # gzfile() reads plain and gzipped files alike
  con <- gzfile(path, open = "rb")
  on.exit(close(con), add = TRUE)
  repeat {
    lines <- readLines(con, n = chunk_size, warn = FALSE)
    if (length(lines) == 0L) {
      break
    }
    fun(lines)
  }
  invisible(NULL)
}

is_gzip_file <- function(path) {
  magic <- readBin(path, what = "raw", n = 2L)
  length(magic) == 2L && identical(magic, as.raw(c(0x1f, 0x8b)))
}

#' Run Salmon Quantification
#'
#' Perform transcript quantification using Salmon's
#' selective-alignment-based mode from raw RNA-seq reads.
#'
#' @details
#' The command detects the library type (`--libType A`), enables the
#' sequence, GC and positional bias corrections,
#' allows dovetailing mates and soft-clipped overhangs, writes the
#' equivalence classes with their weights (`--dumpEq --dumpEqWeights`) and
#' draws `num_gibbs` Gibbs samples (`aux_info/bootstrap/`) so that
#' inferential replicates are available for `swish`. Selective alignment and
#' the variational Bayes optimiser are Salmon defaults and are not passed
#' explicitly; Salmon 2 quantifies deterministically (identical output for
#' any thread count). Single-end reads are passed with `--unmatedReads`.
#'
#' @param input_r1 Path to the FASTQ file for read 1 (or single-end reads).
#' @param input_r2 Optional path to the FASTQ file for read 2
#' (paired-end reads).
#' @param index_path Path to the Salmon index directory (default is
#' "salmon_index").
#' @param output_dir Directory to save the quantification output
#' (default is "quant_output").
#' @param num_threads Number of threads to use (default is 8).
#' @param num_gibbs Number of Gibbs samples (inferential replicates) drawn per
#' sample (default is 30). They quantify how well each transcript abundance is
#' determined by the reads and are used by [add_uncertainty()] and by `swish`.
#' The genome-wide evaluation of the uncertainty evidence article used 30 and
#' did not test more draws for its calls. For isoforms near a cutoff (few
#' unique fragments or reads, retained introns, lncRNAs), more draws (for
#' example 100) are expected to make the `uncertain` flag more stable, since
#' InfRV from more draws has a smaller sampling error, but this was not
#' tested; they cost time and disk space. TreeTerminus trees of lncRNAs often
#' differed between a 30-draw run and a separate 100-draw run, a comparison
#' that does not separate the number of draws from run-to-run randomness (see
#' the quantification stability article).
#' Set it to `0` to skip them.
#' @param min_score_fraction Minimum score fraction for alignments
#' (default is "0.65").
#' @param env_name Name of the conda environment with Salmon installed
#' (default is "salmon-env").
#'
#' @returns processx style output list.
#'
#' @export
salmon_quant <- function(
  input_r1,
  input_r2 = NULL,
  index_path = "salmon_index",
  output_dir = "quant_output",
  num_threads = 8,
  num_gibbs = 30,
  min_score_fraction = "0.65",
  env_name = "salmon-env"
) {
  rlang::check_installed("condathis")
  check_env_installed(env_name)
  rlang::check_required(input_r1)
  if (isFALSE(fs::dir_exists(index_path))) {
    cli::cli_abort(
      c(
        "The Salmon index directory '{index_path}' does not exist.",
        "i" = "Please run salmon_index() first to create the index."
      ),
      class = "isoformic_salmon_missing_index_error"
    )
  }

  gibbs_args <- if (isTRUE(num_gibbs > 0)) {
    c("--numGibbsSamples", num_gibbs)
  } else {
    NULL
  }
  if (is.null(input_r2)) {
    read_args <- c("--unmatedReads", input_r1)
  } else {
    read_args <- c("--mates1", input_r1, "--mates2", input_r2)
  }
  condathis::run(
    "salmon",
    "quant",
    "--libType",
    "A",
    "--index",
    index_path,
    read_args,
    "--output",
    output_dir,
    "--seqBias",
    "--gcBias",
    "--posBias",
    "--threads",
    num_threads,
    "--minScoreFraction",
    min_score_fraction,
    "--allowDovetail",
    "--softclipOverhangs",
    "--dumpEq",
    "--dumpEqWeights",
    gibbs_args,
    env_name = env_name,
    verbose = "output"
  )
}
