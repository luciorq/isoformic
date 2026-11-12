# Splice-site strength --------------------------------------------------------
#
# The sequences around the donor and acceptor site of every intron are cut out
# of the genome once per annotation directory
# (`annot_data_splice_sites.parquet`). Position weight matrices are built from
# the introns of the canonical transcripts of the same annotation, so no
# external model is needed.

#' Splice-Site Strength from the Annotation and the Genome
#'
#' @description
#' Scores the donor (5') and acceptor (3') splice site of every intron with
#' position weight matrices built from the annotation itself, and reports the
#' weakest site of every transcript. A weak site is one way for an isoform to
#' be rare or to be a splicing error; the score is evidence to read next to
#' [transcript_productivity()], not a classifier.
#'
#' * `extract_splice_sites(x, genome_fasta)` cuts the sequences around the
#'   splice sites out of a genome FASTA and stores them in the annotation
#'   directory. [extract_transcript_sequences()] does the same when it is
#'   given a genome FASTA, so this is only needed for directories whose
#'   sequences came from a transcript FASTA.
#' * `splice_site_strength(x)` returns one row per transcript, or one row per
#'   intron with `by = "intron"`.
#'
#' @details
#' The donor window holds the last 3 exonic and the first 6 intronic bases,
#' the acceptor window the last 20 intronic and the first 3 exonic bases, in
#' the direction of transcription (the windows of MaxEntScan, Yeo and Burge
#' 2004). The matrices are trained on the unique introns of the canonical
#' transcripts (`is_canonical`, see [prepare_isoformic_annotation()]), or on
#' all introns when the annotation marks fewer than 100 canonical introns.
#' A score is the sum over the window of `log2(frequency / background)`, with
#' the background taken from the base composition of the training windows
#' and half a pseudocount per base. Positive means closer to the consensus.
#'
#' Limits. The matrices assume that positions are independent, so they are
#' weaker than MaxEntScan or a sequence model, and they describe the major
#' (GT-AG) spliceosome: U12 introns and other non-canonical sites get low
#' scores whatever their real efficiency (`n_noncanonical` counts the introns
#' that are not GT-AG). Canonical transcripts train the matrices, so their
#' own scores are slightly optimistic. Windows that run off the sequence or
#' hold a base other than A, C, G or T are not scored.
#'
#' @param x An annotation directory path or an [IsoformicExperiment].
#' @param self An annotation directory path or an [IsoformicExperiment]. An
#'   experiment reports its own transcripts.
#' @param genome_fasta Path to a genome FASTA file (`.fa`, `.fa.gz`).
#' @param overwrite Logical, replace an existing table (default `FALSE`).
#' @param chunk_bytes Number of bytes read per chunk while streaming.
#' @param gene_name Optional gene names or IDs to keep.
#' @param by `"transcript"` (default) or `"intron"`.
#' @param ... Unused.
#'
#' @returns
#' * `extract_splice_sites()`: the path of the written Parquet file,
#'   invisibly. The table has `seqid`, `strand`, `intron_start`, `intron_end`
#'   (1-based, inclusive), `donor_seq` and `acceptor_seq`.
#' * `splice_site_strength()`: a tibble. By transcript: `transcript_id`,
#'   `gene_id`, `gene_name`, `n_introns`, `n_scored`, `weakest_donor`,
#'   `weakest_acceptor`, `mean_donor`, `mean_acceptor` and `n_noncanonical`;
#'   transcripts without introns have `NA` scores. By intron: one row per
#'   transcript and intron with the coordinates, the two window sequences,
#'   `donor_score`, `acceptor_score` and `canonical_dinucleotides`.
#'
#' @references
#' Yeo, G., Burge, C.B. (2004). Maximum entropy modeling of short sequence
#' motifs with applications to RNA splicing signals. *Journal of
#' Computational Biology* 11(2-3), 377-394.
#'
#' @examples
#' gff <- system.file("extdata", "example.gff3", package = "isoformic")
#' annot <- prepare_isoformic_annotation(
#'   gff,
#'   output_path = tempdir(),
#'   annot_name = "splice_demo"
#' )
#' has_splice_sites(annot)
#' \dontrun{
#' genome <- download_reference(file_type = "genome_fasta", region = "primary")
#' extract_splice_sites(annot, genome_fasta = genome)
#' splice_site_strength(annot)
#' }
#'
#' @seealso [extract_transcript_sequences()], [transcript_productivity()]
#' @name splice_sites
#' @rdname splice_sites
#' @export
extract_splice_sites <- function(
  x,
  genome_fasta,
  overwrite = FALSE,
  chunk_bytes = 64e6
) {
  annot_dir <- annotation_dir_of(x)
  out_file <- as.character(annot_file_path(annot_dir, "splice_sites"))
  if (isTRUE(fs::file_exists(out_file)) && !isTRUE(overwrite)) {
    cli::cli_inform(c(
      i = "{.path {out_file}} already exists; use {.code overwrite = TRUE} to replace it."
    ))
    return(invisible(out_file))
  }
  assert_file_exists(genome_fasta)
  exons <- read_annotation_table(
    annot_dir,
    "exons",
    columns = c("transcript_id", "seqid", "start_pos", "end_pos", "strand")
  )
  res <- sequences_from_genome(
    exons,
    genome_fasta,
    chunk_bytes = chunk_bytes,
    transcripts = FALSE
  )
  write_splice_sites(annot_dir, attr(res, "splice_sites"), genome_fasta)
}

#' @rdname splice_sites
#' @export
has_splice_sites <- function(x) {
  annot_dir <- annotation_dir_of(x)
  isTRUE(fs::file_exists(annot_file_path(annot_dir, "splice_sites")))
}

#' @rdname splice_sites
#' @export
splice_site_strength <- S7::new_generic(
  "splice_site_strength",
  "self",
  function(self, gene_name = NULL, by = c("transcript", "intron"), ...) {
    S7::S7_dispatch()
  }
)

S7::method(splice_site_strength, IsoformicExperiment) <- function(
  self,
  gene_name = NULL,
  by = c("transcript", "intron"),
  ...
) {
  by <- rlang::arg_match(by)
  res <- compute_splice_site_strength(self@annot_path, by)
  res <- res[res$transcript_id %in% row_names(self), ]
  filter_genes(res, gene_name)
}

S7::method(splice_site_strength, S7::class_character) <- function(
  self,
  gene_name = NULL,
  by = c("transcript", "intron"),
  ...
) {
  by <- rlang::arg_match(by)
  filter_genes(
    compute_splice_site_strength(annotation_dir_of(self), by),
    gene_name
  )
}

# Internals --------------------------------------------------------------------

splice_window <- function() {
  list(
    donor_exon = 3L,
    donor_intron = 6L,
    acceptor_intron = 20L,
    acceptor_exon = 3L
  )
}

#' Transcript x intron table from an exon table
#' @keywords internal
#' @noRd
transcript_introns <- function(exons) {
  ex <- exons[
    order(exons$transcript_id, exons$start_pos),
    c("transcript_id", "seqid", "strand", "start_pos", "end_pos")
  ]
  n <- nrow(ex)
  if (n < 2L) {
    return(tibble::tibble(
      transcript_id = character(),
      seqid = character(),
      strand = character(),
      intron_start = integer(),
      intron_end = integer()
    ))
  }
  same <- ex$transcript_id[-n] == ex$transcript_id[-1L]
  res <- tibble::tibble(
    transcript_id = ex$transcript_id[-n][same],
    seqid = ex$seqid[-n][same],
    strand = ex$strand[-n][same],
    intron_start = ex$end_pos[-n][same] + 1L,
    intron_end = ex$start_pos[-1L][same] - 1L
  )
  res[res$intron_end >= res$intron_start, ]
}

#' Unique introns of an exon table
#' @keywords internal
#' @noRd
unique_introns <- function(exons) {
  introns <- transcript_introns(exons)
  introns$transcript_id <- NULL
  dplyr::distinct(introns)
}

#' Donor and acceptor windows of the introns of one sequence region
#'
#' Windows are in the direction of transcription. A window that runs off the
#' sequence is `NA`.
#' @keywords internal
#' @noRd
splice_site_windows <- function(introns, chr_seq) {
  w <- splice_window()
  n_chr <- nchar(chr_seq)
  plus <- !introns$strand %in% "-"
  cut_out <- function(first, last) {
    ok <- first >= 1L & last <= n_chr
    out <- rep(NA_character_, length(first))
    if (!any(ok)) {
      # `substring()` refuses empty positions
      return(out)
    }
    out[ok] <- toupper(substring(chr_seq, first[ok], last[ok]))
    minus <- ok & !plus
    out[minus] <- reverse_complement(out[minus])
    out
  }
  introns$donor_seq <- cut_out(
    ifelse(
      plus,
      introns$intron_start - w$donor_exon,
      introns$intron_end - w$donor_intron + 1L
    ),
    ifelse(
      plus,
      introns$intron_start + w$donor_intron - 1L,
      introns$intron_end + w$donor_exon
    )
  )
  introns$acceptor_seq <- cut_out(
    ifelse(
      plus,
      introns$intron_end - w$acceptor_intron + 1L,
      introns$intron_start - w$acceptor_exon
    ),
    ifelse(
      plus,
      introns$intron_end + w$acceptor_exon,
      introns$intron_start + w$acceptor_intron - 1L
    )
  )
  tibble::as_tibble(introns)
}

#' Write the splice-site table of an annotation directory
#' @keywords internal
#' @noRd
write_splice_sites <- function(annot_dir, sites, genome_fasta) {
  out_file <- as.character(annot_file_path(annot_dir, "splice_sites"))
  sites <- sites[order(sites$seqid, sites$intron_start, sites$intron_end), ]
  write_parquet(sites, out_file)
  meta <- kv_to_list(read_annotation_table(annot_dir, "metadata"))
  meta$splice_sites_source <- paste0("genome:", fs::path_file(genome_fasta))
  meta$splice_sites_n <- as.character(nrow(sites))
  write_parquet(list_to_kv(meta), annot_file_path(annot_dir, "metadata"))
  cli::cli_inform(c(
    v = "Wrote the splice sites of {nrow(sites)} intron{?s} to {.path {out_file}}."
  ))
  invisible(out_file)
}

#' Log-odds position weight matrix of equally long sequences
#'
#' Rows are A, C, G, T; columns the positions. Half a pseudocount per base.
#' @keywords internal
#' @noRd
position_weight_matrix <- function(seqs, background) {
  bases <- c("A", "C", "G", "T")
  width <- nchar(seqs[[1L]])
  counts <- vapply(
    seq_len(width),
    function(j) tabulate(match(substring(seqs, j, j), bases), 4L),
    numeric(4L)
  )
  freq <- (counts + 0.5) / (length(seqs) + 2)
  pwm <- log2(freq / background)
  rownames(pwm) <- bases
  pwm
}

#' Score sequences with a position weight matrix; `NA` when not scorable
#' @keywords internal
#' @noRd
score_pwm <- function(seqs, pwm) {
  width <- ncol(pwm)
  score <- numeric(length(seqs))
  for (j in seq_len(width)) {
    score <- score + pwm[match(substring(seqs, j, j), rownames(pwm)), j]
  }
  score[is.na(seqs) | nchar(seqs) != width] <- NA_real_
  unname(score)
}

#' Sequences that can train or be scored by a matrix
#' @keywords internal
#' @noRd
scorable <- function(seqs, width) {
  !is.na(seqs) & nchar(seqs) == width & !grepl("[^ACGT]", seqs)
}

#' Matrices trained on the introns of canonical transcripts
#' @keywords internal
#' @noRd
train_splice_matrices <- function(sites, training, min_training = 100L) {
  w <- splice_window()
  d_width <- w$donor_exon + w$donor_intron
  a_width <- w$acceptor_intron + w$acceptor_exon
  usable <- scorable(sites$donor_seq, d_width) &
    scorable(sites$acceptor_seq, a_width)
  set <- "canonical"
  if (sum(training & usable) < min_training) {
    training <- rep(TRUE, nrow(sites))
    set <- "all"
  }
  train <- sites[training & usable, ]
  if (!nrow(train)) {
    cli::cli_abort(
      message = c(
        x = "No intron has splice-site sequences that can be scored."
      ),
      class = "isoformic_no_splice_sites"
    )
  }
  letters <- unlist(
    strsplit(c(train$donor_seq, train$acceptor_seq), "", fixed = TRUE),
    use.names = FALSE
  )
  background <- tabulate(match(letters, c("A", "C", "G", "T")), 4L)
  background <- background / sum(background)
  list(
    donor = position_weight_matrix(train$donor_seq, background),
    acceptor = position_weight_matrix(train$acceptor_seq, background),
    background = background,
    n_training = nrow(train),
    training_set = set
  )
}

#' Scores per intron or per transcript for an annotation directory
#' @keywords internal
#' @noRd
compute_splice_site_strength <- function(annot_dir, by) {
  w <- splice_window()
  sites <- read_annotation_table(annot_dir, "splice_sites")
  exons <- read_annotation_table(
    annot_dir,
    "exons",
    columns = c("transcript_id", "seqid", "start_pos", "end_pos", "strand")
  )
  tx <- read_annotation_table(annot_dir, "transcripts")
  genes <- read_annotation_table(
    annot_dir,
    "genes",
    columns = c("gene_id", "gene_name")
  )
  introns <- transcript_introns(exons)
  key <- function(t) {
    paste(t$seqid, t$strand, t$intron_start, t$intron_end, sep = ":")
  }
  site_key <- key(sites)
  canonical_tx <- if ("is_canonical" %in% colnames(tx)) {
    tx$transcript_id[tx$is_canonical %in% TRUE]
  } else {
    character()
  }
  training <- site_key %in%
    key(introns[introns$transcript_id %in% canonical_tx, ])
  model <- train_splice_matrices(sites, training)
  sites$donor_score <- score_pwm(sites$donor_seq, model$donor)
  sites$acceptor_score <- score_pwm(sites$acceptor_seq, model$acceptor)
  sites$canonical_dinucleotides <- substr(
    sites$donor_seq,
    w$donor_exon + 1L,
    w$donor_exon + 2L
  ) ==
    "GT" &
    substr(sites$acceptor_seq, w$acceptor_intron - 1L, w$acceptor_intron) ==
      "AG"
  idx <- match(key(introns), site_key)
  introns <- dplyr::bind_cols(
    introns,
    sites[
      idx,
      c(
        "donor_seq",
        "acceptor_seq",
        "donor_score",
        "acceptor_score",
        "canonical_dinucleotides"
      )
    ]
  )
  gene_of <- tx$gene_id[match(introns$transcript_id, tx$transcript_id)]
  if (identical(by, "intron")) {
    res <- tibble::add_column(
      introns,
      gene_id = gene_of,
      gene_name = genes$gene_name[match(gene_of, genes$gene_id)],
      .after = "transcript_id"
    )
  } else {
    scored <- !is.na(introns$donor_score) & !is.na(introns$acceptor_score)
    per_tx <- introns |>
      dplyr::mutate(scored = scored) |>
      dplyr::summarise(
        n_introns = dplyr::n(),
        n_scored = sum(.data$scored),
        weakest_donor = min_or_na(.data$donor_score),
        weakest_acceptor = min_or_na(.data$acceptor_score),
        mean_donor = mean(.data$donor_score, na.rm = TRUE),
        mean_acceptor = mean(.data$acceptor_score, na.rm = TRUE),
        n_noncanonical = sum(!.data$canonical_dinucleotides, na.rm = TRUE),
        .by = "transcript_id"
      )
    res <- tibble::tibble(
      transcript_id = tx$transcript_id,
      gene_id = tx$gene_id,
      gene_name = genes$gene_name[match(tx$gene_id, genes$gene_id)]
    ) |>
      dplyr::left_join(per_tx, by = "transcript_id")
    res$n_introns[is.na(res$n_introns)] <- 0L
    res$n_scored[is.na(res$n_scored)] <- 0L
    res$n_noncanonical[is.na(res$n_noncanonical)] <- 0L
    res$mean_donor[is.nan(res$mean_donor)] <- NA_real_
    res$mean_acceptor[is.nan(res$mean_acceptor)] <- NA_real_
  }
  attr(res, "n_training") <- model$n_training
  attr(res, "training_set") <- model$training_set
  res
}

min_or_na <- function(x) {
  x <- x[!is.na(x)]
  if (length(x)) min(x) else NA_real_
}
