# Transcript sequences -------------------------------------------------------
#
# Sequence-dependent features (ORF projection, NMD rules beyond the 50-nt rule,
# uORFs, coding potential) need the spliced sequence of every transcript. They
# are stored once per annotation directory in `annot_data_sequences.parquet`
# and read back with the other annotation tables.

#' Extract Transcript Sequences into an Annotation Directory
#'
#' @description
#' Adds the spliced sequence of every transcript of an annotation directory
#' as `annot_data_sequences.parquet`, read by [transcript_orfs()],
#' [coding_potential()] and the sequence-aware rules of
#' [transcript_productivity()]. Sequences come either from a **genome FASTA**
#' (`genome_fasta`; exons are cut out of each chromosome, reverse complemented
#' on the minus strand and concatenated in transcript order, so novel
#' transcripts from long-read or assembly annotations are covered) or from a
#' **transcript FASTA** (`transcript_fasta`, for example the GENCODE
#' `transcripts.fa.gz`; sequences are matched by the first token of the
#' header, with and without the identifier version).
#'
#' The genome is streamed in chunks (plain or gzip-compressed FASTA, no
#' index needed) and only one chromosome is held in memory at a time; a
#' human genome takes a few minutes. Chromosome names are matched to the
#' annotation `seqid` directly, or after adding or removing a `chr` prefix.
#' Soft-masked (lower-case) bases are upper-cased.
#'
#' @param x An annotation directory path or an [IsoformicExperiment].
#' @param genome_fasta Path to a genome FASTA file (`.fa`, `.fa.gz`), see
#'   `download_reference(file_type = "genome_fasta")`.
#' @param transcript_fasta Path to a transcript FASTA file, see
#'   `download_reference(file_type = "fasta")`.
#' @param overwrite Logical, replace existing sequences (default `FALSE`).
#' @param chunk_bytes Number of bytes read per chunk while streaming (default 64 MB).
#'
#' @returns The path of the written Parquet file, invisibly. The table has
#'   columns `transcript_id`, `sequence` and `seq_length`.
#'
#' @details
#' The sequences table is optional: annotation directories without it keep
#' working, and [has_sequences()] tells whether it is present. Transcripts
#' whose chromosome is absent from the FASTA (or whose identifier is absent
#' from a transcript FASTA) are reported and left out.
#'
#' @examples
#' gff <- system.file("extdata", "example.gff3", package = "isoformic")
#' fasta <- system.file("extdata", "example_transcripts.fa.gz", package = "isoformic")
#' annot <- prepare_isoformic_annotation(gff, output_path = tempdir(), annot_name = "seq_demo")
#' extract_transcript_sequences(annot, transcript_fasta = fasta, overwrite = TRUE)
#' has_sequences(annot)
#' head(transcript_sequences(annot)[, c("transcript_id", "seq_length")])
#'
#' \dontrun{
#' # From a genome FASTA (covers novel transcripts of any annotation)
#' genome <- download_reference(version = "49", file_type = "genome_fasta", region = "primary")
#' extract_transcript_sequences(annot, genome_fasta = genome)
#' }
#'
#' @seealso [transcript_orfs()], [coding_potential()]
#' @name transcript_sequences
#' @rdname transcript_sequences
#' @export
extract_transcript_sequences <- function(
  x,
  genome_fasta = NULL,
  transcript_fasta = NULL,
  overwrite = FALSE,
  chunk_bytes = 64e6
) {
  annot_dir <- annotation_dir_of(x)
  out_file <- as.character(annot_file_path(annot_dir, "sequences"))
  if (isTRUE(fs::file_exists(out_file)) && !isTRUE(overwrite)) {
    cli::cli_inform(c(
      i = "{.path {out_file}} already exists; use {.code overwrite = TRUE} to replace it."
    ))
    return(invisible(out_file))
  }
  if (is.null(genome_fasta) == is.null(transcript_fasta)) {
    cli::cli_abort(
      message = c(
        x = "Supply exactly one of {.arg genome_fasta} or {.arg transcript_fasta}."
      ),
      class = "isoformic_sequence_source"
    )
  }
  transcripts <- read_annotation_table(
    annot_dir,
    "transcripts",
    columns = c("transcript_id", "seqid", "strand")
  )
  if (!is.null(genome_fasta)) {
    assert_file_exists(genome_fasta)
    exons <- read_annotation_table(
      annot_dir,
      "exons",
      columns = c("transcript_id", "seqid", "start_pos", "end_pos", "strand")
    )
    res <- sequences_from_genome(exons, genome_fasta, chunk_bytes = chunk_bytes)
    source_desc <- paste0("genome:", fs::path_file(genome_fasta))
  } else {
    assert_file_exists(transcript_fasta)
    res <- sequences_from_tx_fasta(
      transcript_fasta,
      transcripts$transcript_id,
      chunk_bytes = chunk_bytes
    )
    source_desc <- paste0("transcripts:", fs::path_file(transcript_fasta))
  }
  missing <- setdiff(transcripts$transcript_id, res$transcript_id)
  if (length(missing)) {
    missing_seqids <- unique(transcripts$seqid[
      transcripts$transcript_id %in% missing
    ])
    cli::cli_warn(c(
      "!" = "No sequence found for {length(missing)} of {nrow(transcripts)} transcript{?s}.",
      i = "Sequence regions concerned: {.val {utils::head(missing_seqids, 10)}}."
    ))
  }
  res <- res[order(match(res$transcript_id, transcripts$transcript_id)), ]
  res$seq_length <- nchar(res$sequence)
  write_parquet(res, out_file)
  meta <- kv_to_list(read_annotation_table(annot_dir, "metadata"))
  meta$sequences_source <- source_desc
  meta$sequences_n <- as.character(nrow(res))
  write_parquet(list_to_kv(meta), annot_file_path(annot_dir, "metadata"))
  cli::cli_inform(c(
    v = "Wrote {nrow(res)} transcript sequence{?s} to {.path {out_file}}."
  ))
  invisible(out_file)
}

#' @rdname transcript_sequences
#' @export
has_sequences <- function(x) {
  annot_dir <- annotation_dir_of(x)
  isTRUE(fs::file_exists(annot_file_path(annot_dir, "sequences")))
}

#' @rdname transcript_sequences
#' @param transcript_id Optional transcript identifiers to read; `NULL` reads
#'   all sequences.
#' @export
transcript_sequences <- function(x, transcript_id = NULL) {
  annot_dir <- annotation_dir_of(x)
  seqs <- read_annotation_table(annot_dir, "sequences")
  if (!is.null(transcript_id)) {
    seqs <- seqs[seqs$transcript_id %in% transcript_id, ]
  }
  seqs
}

#' Annotation directory of an object or path
#' @keywords internal
#' @noRd
annotation_dir_of <- function(x) {
  annot_dir <- if (S7::S7_inherits(x, IsoformicExperiment)) x@annot_path else x
  if (!is_annotation_dir(annot_dir)) {
    cli::cli_abort(
      message = c(
        x = "{.arg x} must be an annotation directory or an {.cls IsoformicExperiment}."
      ),
      class = "isoformic_annot_path_dont_exist"
    )
  }
  as.character(annot_dir)
}

#' Open a FASTA file (plain or gzip) for binary streaming
#' @keywords internal
#' @noRd
open_fasta <- function(path) {
  path <- as.character(fs::path_abs(path))
  if (grepl("\\.(gz|bgz)$", path, ignore.case = TRUE)) {
    gzfile(path, open = "rb")
  } else {
    file(path, open = "rb")
  }
}

#' Stream a FASTA file and call `on_record(name, sequence)` for every record
#'
#' The file is read in large byte chunks (not line by line, which is slow in
#' R for a genome). `keep(name)` decides whether the sequence of a record is
#' accumulated; records not kept are skipped without storing their sequence.
#'
#' @keywords internal
#' @noRd
stream_fasta <- function(path, keep, on_record, chunk_bytes = 64e6) {
  con <- open_fasta(path)
  withr::defer(close(con))
  current <- NULL
  current_kept <- FALSE
  # The sequence pieces of the current record, in an environment so that the
  # helper below can append to them
  acc <- new.env(parent = emptyenv())
  acc$pieces <- list()
  carry <- ""
  flush <- function() {
    if (!is.null(current) && current_kept) {
      on_record(
        current,
        paste(unlist(acc$pieces, use.names = FALSE), collapse = "")
      )
    }
  }
  add_sequence <- function(text) {
    if (current_kept && nzchar(text)) {
      acc$pieces[[length(acc$pieces) + 1L]] <- gsub(
        "\n",
        "",
        text,
        fixed = TRUE
      )
    }
  }
  repeat {
    chunk <- readChar(con, nchars = chunk_bytes, useBytes = TRUE)
    if (!length(chunk) || !nzchar(chunk)) {
      break
    }
    text <- paste0(carry, chunk)
    # Keep the trailing incomplete line for the next chunk
    newlines <- gregexpr("\n", text, fixed = TRUE)[[1L]]
    if (identical(as.integer(newlines[[1L]]), -1L)) {
      carry <- text
      next
    }
    last_nl <- newlines[[length(newlines)]]
    carry <- substr(text, last_nl + 1L, nchar(text))
    text <- substr(text, 1L, last_nl)
    headers <- gregexpr("(?m)^>[^\n]*", text, perl = TRUE)[[1L]]
    if (identical(as.integer(headers[[1L]]), -1L)) {
      add_sequence(text)
      next
    }
    h_start <- as.integer(headers)
    h_len <- attr(headers, "match.length")
    pos <- 1L
    for (i in seq_along(h_start)) {
      if (h_start[[i]] > pos) {
        add_sequence(substr(text, pos, h_start[[i]] - 1L))
      }
      flush()
      header <- substring(
        text,
        h_start[[i]] + 1L,
        h_start[[i]] + h_len[[i]] - 1L
      )
      current <- sub("^(\\S+).*$", "\\1", header)
      current_kept <- isTRUE(keep(current))
      acc$pieces <- list()
      pos <- h_start[[i]] + h_len[[i]] + 1L
    }
    if (pos <= nchar(text)) {
      add_sequence(substr(text, pos, nchar(text)))
    }
  }
  if (nzchar(carry)) {
    if (startsWith(carry, ">")) {
      flush()
      current <- sub("^>(\\S+).*$", "\\1", carry)
      current_kept <- isTRUE(keep(current))
      acc$pieces <- list()
    } else {
      add_sequence(carry)
    }
  }
  flush()
  invisible(NULL)
}

#' Resolve a FASTA record name to an annotation seqid (with/without `chr`)
#' @keywords internal
#' @noRd
match_seqid <- function(name, seqids) {
  candidates <- c(name, paste0("chr", name), sub("^chr", "", name))
  hit <- candidates[candidates %in% seqids]
  if (length(hit)) hit[[1L]] else NA_character_
}

#' Build transcript sequences from a genome FASTA and an exon table
#' @keywords internal
#' @noRd
sequences_from_genome <- function(exons, genome_fasta, chunk_bytes = 64e6) {
  seqids <- unique(exons$seqid)
  results <- new.env(parent = emptyenv())
  on_record <- function(name, chr_seq) {
    seqid <- match_seqid(name, seqids)
    ex <- exons[exons$seqid == seqid, ]
    ex$seq <- substring(chr_seq, ex$start_pos, ex$end_pos)
    assign(seqid, assemble_transcripts(ex), envir = results)
  }
  stream_fasta(
    genome_fasta,
    keep = function(name) !is.na(match_seqid(name, seqids)),
    on_record = on_record,
    chunk_bytes = chunk_bytes
  )
  found <- ls(results)
  if (!length(found)) {
    cli::cli_abort(
      message = c(
        x = "None of the annotation sequence regions ({.val {utils::head(seqids, 5)}}) were found in the FASTA.",
        i = "File: {.path {genome_fasta}}."
      ),
      class = "isoformic_sequence_source"
    )
  }
  dplyr::bind_rows(mget(found, envir = results))
}

#' Concatenate exon sequences into transcript sequences (DuckDB)
#'
#' Exons are ordered along the transcript (descending genomic position on the
#' minus strand), reverse complemented when needed and joined.
#'
#' @keywords internal
#' @noRd
assemble_transcripts <- function(ex) {
  if (!nrow(ex)) {
    return(tibble::tibble(transcript_id = character(), sequence = character()))
  }
  ex <- as.data.frame(ex[, c("transcript_id", "strand", "start_pos", "seq")])
  con <- duckdb_connect()
  duckdb::duckdb_register(con, "ex", ex)
  res <- DBI::dbGetQuery(
    con,
    "SELECT transcript_id,
            string_agg(
              CASE WHEN strand = '-'
                   THEN reverse(translate(upper(seq), 'ACGTRYKMBDHVN', 'TGCAYRMKVHDBN'))
                   ELSE upper(seq) END,
              ''
              ORDER BY CASE WHEN strand = '-' THEN -start_pos ELSE start_pos END
            ) AS sequence
     FROM ex
     GROUP BY transcript_id"
  )
  tibble::as_tibble(res)
}

#' Read transcript sequences from a transcript FASTA
#' @keywords internal
#' @noRd
sequences_from_tx_fasta <- function(
  transcript_fasta,
  transcript_ids,
  chunk_bytes = 64e6
) {
  wanted <- unique(transcript_ids)
  # Matching is lenient: any trailing `.N` is ignored on both sides. Hashed
  # lookups keep the per-record cost constant for large catalogs.
  unversion <- function(x) sub("\\.[0-9]+$", "", x)
  exact <- list2env(
    as.list(stats::setNames(wanted, wanted)),
    parent = emptyenv()
  )
  lenient <- list2env(
    as.list(stats::setNames(wanted, unversion(wanted))),
    parent = emptyenv()
  )
  resolve_id <- function(name) {
    id <- sub(
      "|",
      "\\1",
      regmatches(name, regexpr("^[^|]+", name)),
      fixed = TRUE
    )
    hit <- get0(id, envir = exact, inherits = FALSE)
    if (!is.null(hit)) {
      return(hit)
    }
    hit <- get0(unversion(id), envir = lenient, inherits = FALSE)
    if (is.null(hit)) NA_character_ else hit
  }
  found <- new.env(parent = emptyenv())
  on_record <- function(name, sequence) {
    id <- resolve_id(name)
    if (!exists(id, envir = found, inherits = FALSE)) {
      assign(id, toupper(sequence), envir = found)
    }
  }
  stream_fasta(
    transcript_fasta,
    keep = function(name) !is.na(resolve_id(name)),
    on_record = on_record,
    chunk_bytes = chunk_bytes
  )
  ids <- ls(found)
  tibble::tibble(
    transcript_id = ids,
    sequence = unlist(mget(ids, envir = found), use.names = FALSE) %||%
      character()
  )
}

#' Strip a trailing `.N` version from Ensembl/GENCODE and RefSeq identifiers
#'
#' Other identifiers (StringTie `MSTRG.1.2`, FANTOM, custom) are returned
#' unchanged, since their last dotted field is not a version.
#'
#' @keywords internal
#' @noRd
strip_id_version <- function(x) {
  versioned <- grepl(
    "^(ENS[A-Z]*[0-9]{6,}|[NX][MRP]_[0-9]+)\\.[0-9]+(_PAR_Y)?$",
    x
  )
  x[versioned] <- sub("\\.[0-9]+(_PAR_Y)?$", "", x[versioned])
  x
}

# Sequence utilities ---------------------------------------------------------

#' Reverse complement of DNA strings (IUPAC aware)
#' @keywords internal
#' @noRd
reverse_complement <- function(x) {
  comp <- chartr("ACGTRYKMBDHVNacgtrykmbdhvn", "TGCAYRMKVHDBNtgcayrmkvhdbn", x)
  vapply(
    strsplit(comp, "", fixed = TRUE),
    function(ch) paste(rev(ch), collapse = ""),
    character(1L)
  )
}

#' Translate DNA (standard code) into protein; `*` marks stops
#' @keywords internal
#' @noRd
translate_dna <- function(x) {
  n <- nchar(x)
  if (n < 3L) {
    return("")
  }
  starts <- seq.int(1L, n - 2L, by = 3L)
  codons <- substring(x, starts, starts + 2L)
  aa <- codon_table()[codons]
  aa[is.na(aa)] <- "X"
  paste(aa, collapse = "")
}

codon_table <- function() {
  bases <- c("T", "C", "A", "G")
  # Order matches the standard genetic code table laid out TCAG x TCAG x TCAG
  codons <- paste0(
    rep(bases, each = 16L),
    rep(rep(bases, each = 4L), times = 4L),
    rep(bases, times = 16L)
  )
  aa <- strsplit(
    "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    ""
  )[[1L]]
  stats::setNames(aa, codons)
}
