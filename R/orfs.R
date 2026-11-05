# ORF projection and sequence-based NMD features -----------------------------
#
# Annotation-only productivity calls depend on the CDS the annotator wrote.
# Novel transcripts (StringTie, SQANTI, FANTOM CAT, HPRC CAT) mostly have no
# CDS. With transcript sequences (see extract_transcript_sequences()) the
# canonical start codon of each gene is projected onto every isoform and
# translated through, so ORFs, premature stops, frameshifts and uORFs are
# computed rather than read.

#' Open Reading Frames of Transcripts from Their Sequences
#'
#' @description
#' Computes, for every assayed transcript, an open reading frame (ORF) and the
#' sequence-level features that decide whether it is productive. Requires
#' transcript sequences in the annotation directory (see
#' [extract_transcript_sequences()]). The ORF of each transcript comes from
#' the first applicable source (`orf_source`):
#'
#' 1. `"annotation"`: a complete annotated CDS. Its stop codon is added when
#'    the annotation excludes it (as GENCODE does) and the ORF is re-read
#'    from the annotated start codon to the first in-frame stop, so a CDS
#'    that does not match the sequence is corrected
#'    (`annotation_cds_agrees = FALSE`). Transcripts tagged `seleno`,
#'    `ribosomal_slippage` or `non_ATG_start` keep the annotated CDS.
#' 2. `"projected"`: the start codon of the gene's canonical isoform lies in
#'    an exon of the transcript and reads `ATG` there; translation proceeds
#'    from that position to the first in-frame stop (the ORFanage approach).
#'    This is what rescues novel isoforms without CDS.
#' 3. `"annotation_partial"`: an annotated CDS flagged incomplete
#'    (`cds_start_NF` / `cds_end_NF`) whose start could not be projected; the
#'    annotated span is kept, and when its 3' end is known the NMD rule still
#'    applies (GENCODE annotates many NMD transcripts this way).
#' 4. `"longest"`: none of the above; the longest `ATG`-to-stop ORF of at
#'    least `min_orf_length` nt is used and `coding_probability` tells
#'    whether it is believable.
#' 5. `"none"`: no ORF of at least `min_orf_length` nt.
#'
#' From the ORF the function derives the stop codon (`NA` when translation
#' runs into the transcript end, a non-stop decay substrate), the UTR
#' lengths, uORFs in the 5' UTR (`n_uorf`, complete upstream ORFs;
#' `n_overlapping_uorf`, out-of-frame upstream ORFs overlapping the main
#' start), the genomic coordinates of start and stop, whether start and stop
#' are shared with the canonical isoform and whether the reading frame at the
#' stop codon is in frame with the canonical ORF (`frame_vs_canonical`), and
#' the NMD class:
#'
#' * `"ptc"`: the stop lies more than `nmd_rule` nt upstream of the last
#'   exon-exon junction (the 50-nt rule), NMD predicted.
#' * `"ptc_start_proximal"`: same, but the ORF is shorter than
#'   `start_proximal` nt; still called unproductive, with the note that such
#'   transcripts can escape NMD by translation re-initiation.
#' * `"non_stop"`: no in-frame stop codon.
#' * `"long_3utr"`: no premature stop, but a 3' UTR longer than `long_utr3`
#'   nt (default `Inf`, so the flag is off unless asked for).
#' * `"none"`: no NMD feature.
#' * `NA`: no stop codon in a transcript whose 3' end is annotated as
#'   incomplete (`cds_end_NF` / `mRNA_end_NF`, column `end_incomplete`), so
#'   nothing can be said.
#'
#' `coding_potential()` scores every transcript's longest ORF with the
#' Fickett TESTCODE statistic and an in-frame hexamer usage score, combined
#' by a logistic model trained on GENCODE coding versus lncRNA transcripts
#' (see `data-raw/build_coding_model.R`); `coding_probability` above
#' `coding_cutoff` (0.364, the CPAT human cutoff) is called coding.
#'
#' @param self An [IsoformicExperiment] or an annotation directory path.
#' @param gene_name Optional genes (names or IDs) to restrict the result to.
#' @param min_orf_length Minimum ORF length in nt, stop codon included
#'   (default 90, i.e. 29 amino acids).
#' @param nmd_rule Distance (nt) from the stop codon to the last junction
#'   above which NMD is predicted; default 50.
#' @param start_proximal ORF length (nt) below which a premature stop is
#'   classed `ptc_start_proximal` instead of `ptc`; default 150.
#' @param long_utr3 3' UTR length (nt) above which `long_3utr` is flagged;
#'   default `Inf`.
#' @param ... Passed on to methods.
#'
#' @returns A tibble with one row per transcript.
#'
#' @examples
#' iso <- example_isoformic()
#' orfs <- transcript_orfs(iso, gene_name = "RBPJ")
#' orfs[, c("transcript_name", "orf_source", "orf_length", "nmd_class", "coding_probability")]
#'
#' @seealso [transcript_productivity()], [transcript_features()]
#' @name transcript_orfs
#' @rdname transcript_orfs
#' @export
transcript_orfs <- S7::new_generic(
  "transcript_orfs",
  "self",
  function(
    self,
    gene_name = NULL,
    min_orf_length = 90L,
    nmd_rule = 50,
    start_proximal = 150,
    long_utr3 = Inf,
    ...
  ) {
    S7::S7_dispatch()
  }
)

S7::method(transcript_orfs, IsoformicExperiment) <- function(
  self,
  gene_name = NULL,
  min_orf_length = 90L,
  nmd_rule = 50,
  start_proximal = 150,
  long_utr3 = Inf,
  ...
) {
  signature <- paste(
    annot_full_signature(self),
    annot_signature(self@annot_path, "sequences"),
    rlang::hash(row_names(self)),
    min_orf_length,
    nmd_rule,
    start_proximal,
    long_utr3
  )
  res <- object_cache_fetch(self, "transcript_orfs", signature, function() {
    compute_transcript_orfs(
      transcripts = row_data(self),
      exons = self@annot_data_exons,
      sequences = transcript_sequences(self, transcript_id = row_names(self)),
      min_orf_length = min_orf_length,
      nmd_rule = nmd_rule,
      start_proximal = start_proximal,
      long_utr3 = long_utr3
    )
  })
  filter_genes(res, gene_name)
}

S7::method(transcript_orfs, S7::class_character) <- function(
  self,
  gene_name = NULL,
  min_orf_length = 90L,
  nmd_rule = 50,
  start_proximal = 150,
  long_utr3 = Inf,
  ...
) {
  annot <- prepare_annotation(self)
  transcripts <- annot$transcripts |>
    dplyr::left_join(annot$genes[, c("gene_id", "gene_name")], by = "gene_id")
  res <- compute_transcript_orfs(
    transcripts = transcripts,
    exons = annot$exons,
    sequences = transcript_sequences(self),
    min_orf_length = min_orf_length,
    nmd_rule = nmd_rule,
    start_proximal = start_proximal,
    long_utr3 = long_utr3
  )
  filter_genes(res, gene_name)
}

stop_codons <- function() {
  c("TAA", "TAG", "TGA")
}

#' Compute ORFs and sequence features for a transcript table
#' @keywords internal
#' @noRd
compute_transcript_orfs <- function(
  transcripts,
  exons,
  sequences,
  min_orf_length = 90L,
  nmd_rule = 50,
  start_proximal = 150,
  long_utr3 = Inf
) {
  tx <- transcripts
  n <- nrow(tx)
  for (col in c("cds_start", "cds_end", "cds_length")) {
    if (!col %in% colnames(tx)) {
      tx[[col]] <- NA_integer_
    }
  }
  if (!"strand" %in% colnames(tx)) {
    strand_by_tx <- exons[
      !duplicated(exons$transcript_id),
      c("transcript_id", "strand")
    ]
    tx <- dplyr::left_join(tx, strand_by_tx, by = "transcript_id")
  }
  seq <- sequences$sequence[match(tx$transcript_id, sequences$transcript_id)]
  seq[is.na(seq)] <- ""
  seq_len <- nchar(seq)
  has_seq <- seq_len > 0L

  ex <- exon_cum_table(exons, tx)
  tx_len <- stats::setNames(
    as.vector(rowsum(ex$ex_len, ex$transcript_id)),
    sort(unique(ex$transcript_id))
  )
  tx_length <- unname(tx_len[tx$transcript_id])
  is_last <- !duplicated(ex$transcript_id, fromLast = TRUE)
  last_len <- stats::setNames(ex$ex_len[is_last], ex$transcript_id[is_last])
  last_junction <- tx_length - unname(last_len[tx$transcript_id])
  n_exons <- unname(stats::setNames(
    ex$idx_in_tx[is_last],
    ex$transcript_id[is_last]
  )[
    tx$transcript_id
  ])
  length_mismatch <- has_seq & !is.na(tx_length) & seq_len != tx_length
  if (any(length_mismatch)) {
    cli::cli_warn(c(
      "!" = paste0(
        "{sum(length_mismatch)} transcript{?s} {?has/have} a sequence length ",
        "differing from the exon length; {?its/their} ORF{?s} {?is/are} skipped."
      )
    ))
    seq[length_mismatch] <- ""
    has_seq <- nzchar(seq)
  }

  # 1. Annotated complete CDS -----------------------------------------------
  tags <- if ("tags" %in% colnames(tx)) tx$tags else rep(NA_character_, n)
  tags[is.na(tags)] <- ""
  has_cds <- !is.na(tx$cds_length) & tx$cds_length > 0
  cds_complete <- has_cds & !grepl("cds_(start|end)_NF", tags)
  minus <- tx$strand == "-"
  cds_first_g <- ifelse(minus, tx$cds_end, tx$cds_start)
  cds_last_g <- ifelse(minus, tx$cds_start, tx$cds_end)
  cds_tx <- genomic_to_transcript_coord(
    ex,
    tx$transcript_id,
    cds_first_g,
    cds_last_g
  )
  orf_start <- rep(NA_real_, n)
  orf_end <- rep(NA_real_, n)
  orf_source <- rep("none", n)
  ann_ok <- has_seq & cds_complete & !is.na(cds_tx$first) & !is.na(cds_tx$last)
  orf_start[ann_ok] <- cds_tx$first[ann_ok]
  orf_end[ann_ok] <- cds_tx$last[ann_ok]
  orf_source[ann_ok] <- "annotation"
  # GENCODE excludes the stop codon from the CDS: extend when the next codon is a stop
  codon_at <- function(s, p) substring(s, p, p + 2L)
  ends_with_stop <- ann_ok & codon_at(seq, orf_end - 2) %in% stop_codons()
  next_is_stop <- ann_ok &
    !ends_with_stop &
    codon_at(seq, orf_end + 1) %in% stop_codons()
  orf_end[next_is_stop] <- orf_end[next_is_stop] + 3
  # The sequence is the ground truth: translate from the annotated start and
  # take the first in-frame stop, unless the annotation documents a reason
  # for a non-standard ORF (selenocysteine, ribosomal slippage, non-ATG start)
  stops <- find_stop_positions(seq)
  trusted <- grepl("seleno|ribosomal_slippage|non_ATG_start", tags)
  retranslate <- ann_ok & !trusted & codon_at(seq, orf_start) == "ATG"
  ann_end <- first_inframe_stop(stops, ifelse(retranslate, orf_start, NA_real_))
  annotation_cds_agrees <- rep(NA, n)
  annotation_cds_agrees[ann_ok] <- TRUE
  disagree <- retranslate & !is.na(ann_end) & ann_end != orf_end
  annotation_cds_agrees[disagree] <- FALSE
  orf_end[disagree] <- ann_end[disagree]

  # 2. Projection of the canonical start codon -----------------------------
  canonical <- canonical_transcripts(tx)
  canonical_id <- unname(canonical[tx$gene_id])
  can_idx <- match(canonical_id, tx$transcript_id)
  can_start_g <- cds_first_g[can_idx]
  can_start_g[!(cds_complete[can_idx] %in% TRUE)] <- NA_real_
  need <- has_seq & orf_source == "none" & !is.na(can_start_g)
  proj_pos <- rep(NA_real_, n)
  proj_pos[need] <- genomic_to_transcript_pos(
    ex,
    tx$transcript_id[need],
    can_start_g[need]
  )
  proj_ok <- !is.na(proj_pos) & codon_at(seq, proj_pos) == "ATG"
  proj_end <- first_inframe_stop(stops, proj_pos)
  proj_end_eff <- ifelse(
    is.na(proj_end),
    proj_pos + floor((seq_len - proj_pos + 1) / 3) * 3 - 1,
    proj_end
  )
  proj_use <- proj_ok & (proj_end_eff - proj_pos + 1) >= min_orf_length
  orf_start[proj_use] <- proj_pos[proj_use]
  orf_end[proj_use] <- proj_end_eff[proj_use]
  orf_source[proj_use] <- "projected"

  # 3. Incomplete annotated CDS (cds_start_NF / cds_end_NF): the annotated
  #    span is kept as a partial ORF; a known stop codon still allows the NMD rule
  partial <- has_seq &
    orf_source == "none" &
    has_cds &
    !cds_complete &
    !is.na(cds_tx$first) &
    !is.na(cds_tx$last)
  orf_start[partial] <- cds_tx$first[partial]
  orf_end[partial] <- cds_tx$last[partial]
  orf_source[partial] <- "annotation_partial"
  partial_next_stop <- partial &
    !(codon_at(seq, orf_end - 2) %in% stop_codons()) &
    codon_at(seq, orf_end + 1) %in% stop_codons()
  orf_end[partial_next_stop] <- orf_end[partial_next_stop] + 3

  # 4. Longest ORF (also the basis of the coding potential) ----------------
  longest <- longest_orfs(seq, stops, min_orf_length = min_orf_length)
  use_longest <- has_seq & orf_source == "none" & !is.na(longest$start)
  orf_start[use_longest] <- longest$start[use_longest]
  orf_end[use_longest] <- longest$end[use_longest]
  orf_source[use_longest] <- "longest"

  # Derived features --------------------------------------------------------
  has_orf <- orf_source != "none"
  stop_codon <- ifelse(has_orf, codon_at(seq, orf_end - 2), NA_character_)
  stop_codon[!stop_codon %in% stop_codons()] <- NA_character_
  start_codon <- ifelse(has_orf, codon_at(seq, orf_start), NA_character_)
  orf_length <- ifelse(has_orf, orf_end - orf_start + 1, NA_real_)
  protein_length <- ifelse(
    has_orf,
    floor((orf_length - ifelse(is.na(stop_codon), 0, 3)) / 3),
    NA_real_
  )
  utr5_length <- ifelse(has_orf, orf_start - 1, NA_real_)
  utr3_length <- ifelse(has_orf, tx_length - orf_end, NA_real_)
  stop_to_last_junction <- ifelse(
    has_orf & n_exons > 1L & !is.na(stop_codon),
    last_junction - orf_end,
    NA_real_
  )
  nmd_class <- rep(NA_character_, n)
  nmd_class[has_orf] <- "none"
  nmd_class[has_orf & is.na(stop_codon)] <- "non_stop"
  # A missing stop codon in a transcript whose 3' end is annotated as
  # incomplete (cds_end_NF / mRNA_end_NF) is not evidence of non-stop decay
  end_incomplete <- grepl("cds_end_NF|mRNA_end_NF", tags)
  nmd_class[has_orf & is.na(stop_codon) & end_incomplete] <- NA_character_
  ptc <- has_orf &
    !is.na(stop_to_last_junction) &
    stop_to_last_junction > nmd_rule
  nmd_class[ptc & orf_length < start_proximal] <- "ptc_start_proximal"
  nmd_class[ptc & orf_length >= start_proximal] <- "ptc"
  nmd_class[
    has_orf &
      nmd_class == "none" &
      !is.na(utr3_length) &
      utr3_length > long_utr3
  ] <- "long_3utr"

  uorfs <- count_uorfs(seq, stops, orf_start)

  orf_start_g <- transcript_to_genomic_coord(ex, tx$transcript_id, orf_start)
  stop_first <- ifelse(is.na(stop_codon), NA_real_, orf_end - 2)
  orf_stop_g <- transcript_to_genomic_coord(ex, tx$transcript_id, stop_first)
  orf_end_g <- transcript_to_genomic_coord(ex, tx$transcript_id, orf_end)
  can_orf_start_g <- orf_start_g[can_idx]
  can_stop_g <- orf_stop_g[can_idx]
  start_same <- ifelse(
    has_orf & !is.na(can_orf_start_g),
    orf_start_g == can_orf_start_g,
    NA
  )
  stop_same <- ifelse(
    has_orf & !is.na(can_stop_g) & !is.na(orf_stop_g),
    orf_stop_g == can_stop_g,
    NA
  )
  # Frame of this ORF's stop codon within the canonical transcript
  can_stop_tx <- genomic_to_transcript_pos(ex, canonical_id, orf_stop_g)
  can_orf_start_tx <- orf_start[can_idx]
  frame_vs_canonical <- rep(NA_character_, n)
  ok_frame <- has_orf & !is.na(can_stop_tx) & !is.na(can_orf_start_tx)
  frame_vs_canonical[ok_frame] <- ifelse(
    (can_stop_tx[ok_frame] - can_orf_start_tx[ok_frame]) %% 3 == 0,
    "in_frame",
    "frameshift"
  )
  novel <- has_orf & !is.na(can_orf_start_tx) & is.na(can_stop_tx)
  frame_vs_canonical[novel] <- "novel_sequence"
  frame_vs_canonical[
    has_orf & (tx$transcript_id == canonical_id) %in% TRUE
  ] <- "in_frame"

  potential <- coding_potential_scores(seq, longest)

  tibble::tibble(
    transcript_id = tx$transcript_id,
    transcript_name = if ("transcript_name" %in% colnames(tx)) {
      tx$transcript_name
    } else {
      tx$transcript_id
    },
    gene_id = tx$gene_id,
    gene_name = if ("gene_name" %in% colnames(tx)) tx$gene_name else tx$gene_id,
    has_sequence = has_seq,
    orf_source = orf_source,
    annotation_cds_agrees = annotation_cds_agrees,
    orf_tx_start = orf_start,
    orf_tx_end = orf_end,
    orf_length = orf_length,
    protein_length = protein_length,
    start_codon = start_codon,
    stop_codon = stop_codon,
    orf_genomic_start = orf_start_g,
    orf_genomic_end = orf_end_g,
    utr5_length = utr5_length,
    utr3_length = utr3_length,
    stop_to_last_junction = stop_to_last_junction,
    nmd_class = nmd_class,
    end_incomplete = end_incomplete,
    n_uorf = uorfs$n_uorf,
    n_overlapping_uorf = uorfs$n_overlapping_uorf,
    canonical_transcript_id = canonical_id,
    start_same_as_canonical = start_same,
    stop_same_as_canonical = stop_same,
    frame_vs_canonical = frame_vs_canonical,
    longest_orf_length = longest$length,
    orf_coverage = potential$orf_coverage,
    fickett_score = potential$fickett_score,
    hexamer_score = potential$hexamer_score,
    coding_probability = potential$coding_probability
  )
}

#' Positions of all stop codons (any frame) in each sequence
#'
#' Stop codons never overlap each other, so a single non-overlapping regex
#' scan finds them all.
#'
#' @keywords internal
#' @noRd
find_stop_positions <- function(seq) {
  res <- gregexpr("TAA|TAG|TGA", seq, perl = TRUE)
  lapply(res, function(m) {
    m <- as.integer(m)
    m[m > 0L]
  })
}

#' First stop codon in frame with and downstream of `start` (position of its
#' first base plus 2, i.e. the last base of the ORF), `NA` when none
#' @keywords internal
#' @noRd
first_inframe_stop <- function(stops, start) {
  out <- rep(NA_real_, length(start))
  idx <- which(!is.na(start))
  for (i in idx) {
    s <- stops[[i]]
    s <- s[s > start[[i]] & (s - start[[i]]) %% 3 == 0]
    if (length(s)) {
      out[[i]] <- s[[1L]] + 2
    }
  }
  out
}

#' Longest ATG-to-stop ORF per sequence
#'
#' Returns `start`, `end` (last base of the stop codon, or of the last full
#' codon when the ORF is open-ended) and `length`; `NA` when no ORF reaches
#' `min_orf_length`. Open-ended ORFs (no in-frame stop) are considered.
#'
#' @keywords internal
#' @noRd
longest_orfs <- function(seq, stops, min_orf_length = 90L) {
  n <- length(seq)
  out <- list(
    start = rep(NA_real_, n),
    end = rep(NA_real_, n),
    length = rep(NA_real_, n)
  )
  atgs <- gregexpr("ATG", seq, fixed = TRUE)
  seq_len <- nchar(seq)
  for (i in seq_len(n)) {
    a <- as.integer(atgs[[i]])
    a <- a[a > 0L]
    if (!length(a)) {
      next
    }
    s <- stops[[i]]
    best_len <- 0
    best_start <- NA_real_
    best_end <- NA_real_
    for (f in 0:2) {
      af <- a[a %% 3 == f]
      if (!length(af)) {
        next
      }
      sf <- s[s %% 3 == f]
      # Next in-frame stop after each ATG; open-ended when none
      pos <- findInterval(af, sf) + 1L
      nxt <- ifelse(pos <= length(sf), sf[pmin(pos, length(sf))], NA_real_)
      end <- ifelse(
        is.na(nxt),
        af + floor((seq_len[[i]] - af + 1) / 3) * 3 - 1,
        nxt + 2
      )
      # Keep only the first ATG of each ORF (same stop => same ORF)
      keep <- !duplicated(end)
      len <- end[keep] - af[keep] + 1
      j <- which.max(len)
      if (length(j) && len[[j]] > best_len) {
        best_len <- len[[j]]
        best_start <- af[keep][[j]]
        best_end <- end[keep][[j]]
      }
    }
    if (best_len >= min_orf_length) {
      out$start[[i]] <- best_start
      out$end[[i]] <- best_end
      out$length[[i]] <- best_len
    }
  }
  out
}

#' Count uORFs in the 5' UTR upstream of `orf_start`
#'
#' A uORF is an ATG upstream of the main start with an in-frame stop that
#' ends before the main start; an overlapping uORF is an out-of-frame
#' upstream ATG whose first in-frame stop lies at or after the main start.
#' In-frame upstream ATGs without an intervening stop are N-terminal
#' extensions and are not counted.
#'
#' @keywords internal
#' @noRd
count_uorfs <- function(seq, stops, orf_start, min_uorf_length = 9L) {
  n <- length(seq)
  n_uorf <- rep(NA_integer_, n)
  n_over <- rep(NA_integer_, n)
  idx <- which(!is.na(orf_start) & orf_start > 1)
  utr <- substring(seq[idx], 1L, orf_start[idx] - 1)
  atgs <- gregexpr("ATG", utr, fixed = TRUE)
  for (k in seq_along(idx)) {
    i <- idx[[k]]
    a <- as.integer(atgs[[k]])
    a <- a[a > 0L]
    if (!length(a)) {
      n_uorf[[i]] <- 0L
      n_over[[i]] <- 0L
      next
    }
    s <- stops[[i]]
    p <- orf_start[[i]]
    stop_end <- vapply(
      a,
      function(x) {
        sx <- s[s > x & (s - x) %% 3 == 0]
        if (length(sx)) sx[[1L]] + 2 else NA_real_
      },
      numeric(1L)
    )
    complete <- !is.na(stop_end) &
      stop_end < p &
      (stop_end - a + 1) >= min_uorf_length
    in_frame <- (p - a) %% 3 == 0
    overlapping <- !in_frame & (is.na(stop_end) | stop_end >= p)
    n_uorf[[i]] <- sum(complete)
    n_over[[i]] <- sum(overlapping)
  }
  list(n_uorf = n_uorf, n_overlapping_uorf = n_over)
}

# Coding potential -----------------------------------------------------------

#' @rdname transcript_orfs
#' @param coding_cutoff Coding probability above which a transcript is called
#'   coding (default 0.364).
#' @export
coding_potential <- S7::new_generic(
  "coding_potential",
  "self",
  function(
    self,
    gene_name = NULL,
    min_orf_length = 90L,
    coding_cutoff = 0.364,
    ...
  ) {
    S7::S7_dispatch()
  }
)

S7::method(coding_potential, IsoformicExperiment) <- function(
  self,
  gene_name = NULL,
  min_orf_length = 90L,
  coding_cutoff = 0.364,
  ...
) {
  orfs <- transcript_orfs(
    self,
    gene_name = gene_name,
    min_orf_length = min_orf_length
  )
  select_coding_cols(orfs, coding_cutoff)
}

S7::method(coding_potential, S7::class_character) <- function(
  self,
  gene_name = NULL,
  min_orf_length = 90L,
  coding_cutoff = 0.364,
  ...
) {
  orfs <- transcript_orfs(
    self,
    gene_name = gene_name,
    min_orf_length = min_orf_length
  )
  select_coding_cols(orfs, coding_cutoff)
}

select_coding_cols <- function(orfs, coding_cutoff) {
  res <- orfs[, c(
    "transcript_id",
    "transcript_name",
    "gene_id",
    "gene_name",
    "longest_orf_length",
    "orf_coverage",
    "fickett_score",
    "hexamer_score",
    "coding_probability"
  )]
  res$is_coding <- !is.na(res$coding_probability) &
    res$coding_probability > coding_cutoff
  res
}

#' Coding potential scores of sequences given their longest ORFs
#' @keywords internal
#' @noRd
coding_potential_scores <- function(seq, longest) {
  seq_len <- nchar(seq)
  raw <- coding_potential_scores_raw(seq, longest)
  fickett <- raw$fickett_score
  hexamer <- raw$hexamer_score
  coverage <- raw$orf_coverage
  orf_len <- ifelse(is.na(longest$length), 0, longest$length)
  model <- coding_model()
  eta <- model$intercept +
    model$log_orf_length * log1p(orf_len) +
    model$orf_coverage * coverage +
    model$fickett_score * fickett +
    model$hexamer_score * hexamer
  prob <- 1 / (1 + exp(-eta))
  prob[seq_len == 0L] <- NA_real_
  fickett[seq_len == 0L] <- NA_real_
  hexamer[seq_len == 0L] <- NA_real_
  coverage[seq_len == 0L] <- NA_real_
  list(
    orf_coverage = coverage,
    fickett_score = fickett,
    hexamer_score = hexamer,
    coding_probability = prob
  )
}

#' Fickett TESTCODE score (Fickett 1982), as implemented in CPAT
#'
#' Combines the position asymmetry of each base over the three codon
#' positions with the base composition, through the published lookup tables.
#'
#' @keywords internal
#' @noRd
fickett_score <- function(seq) {
  position_prob <- list(
    A = c(0.94, 0.68, 0.84, 0.93, 0.58, 0.68, 0.45, 0.34, 0.20, 0.22),
    C = c(0.80, 0.70, 0.70, 0.81, 0.66, 0.48, 0.51, 0.33, 0.30, 0.23),
    G = c(0.90, 0.88, 0.74, 0.64, 0.53, 0.48, 0.27, 0.16, 0.08, 0.08),
    T = c(0.97, 0.97, 0.91, 0.68, 0.69, 0.44, 0.54, 0.20, 0.09, 0.09)
  )
  position_weight <- c(A = 0.26, C = 0.18, G = 0.31, T = 0.33)
  position_para <- c(1.9, 1.8, 1.7, 1.6, 1.5, 1.4, 1.3, 1.2, 1.1, 0.0)
  content_prob <- list(
    A = c(0.28, 0.49, 0.44, 0.55, 0.62, 0.49, 0.54, 0.52, 0.36, 0.51),
    C = c(0.82, 0.64, 0.51, 0.64, 0.59, 0.59, 0.56, 0.49, 0.26, 0.36),
    G = c(0.40, 0.54, 0.47, 0.64, 0.64, 0.73, 0.41, 0.41, 0.33, 0.29),
    T = c(0.28, 0.24, 0.39, 0.40, 0.55, 0.75, 0.56, 0.69, 0.51, 0.58)
  )
  content_weight <- c(A = 0.11, C = 0.12, G = 0.15, T = 0.14)
  content_para <- c(0.33, 0.31, 0.29, 0.27, 0.25, 0.23, 0.21, 0.17, 0.17, 0)
  bases <- c("A", "C", "G", "T")
  lookup <- function(value, para, prob) {
    i <- which(value >= para)[1L]
    if (is.na(i)) prob[[length(prob)]] else prob[[i]]
  }
  vapply(
    seq,
    function(s) {
      n <- nchar(s)
      if (n < 3L) {
        return(0)
      }
      ch <- strsplit(s, "", fixed = TRUE)[[1L]]
      frame <- (seq_along(ch) - 1L) %% 3L
      score <- 0
      for (b in bases) {
        counts <- vapply(0:2, function(f) sum(ch[frame == f] == b), numeric(1L))
        asym <- max(counts) / (min(counts) + 1)
        comp <- sum(ch == b) / n
        score <- score +
          lookup(asym, position_para, position_prob[[b]]) *
            position_weight[[b]] +
          lookup(comp, content_para, content_prob[[b]]) * content_weight[[b]]
      }
      score
    },
    numeric(1L),
    USE.NAMES = FALSE
  )
}

#' Mean in-frame hexamer log-odds (coding vs non-coding) of ORF sequences
#' @keywords internal
#' @noRd
hexamer_score <- function(orf_seq) {
  table <- hexamer_table()
  vapply(
    orf_seq,
    function(s) {
      n <- nchar(s)
      if (n < 6L) {
        return(0)
      }
      starts <- seq.int(1L, n - 5L, by = 3L)
      hex <- substring(s, starts, starts + 5L)
      v <- table[hex]
      v <- v[!is.na(v)]
      if (length(v)) mean(v) else 0
    },
    numeric(1L),
    USE.NAMES = FALSE
  )
}

#' Frozen hexamer log-odds table (shipped in inst/extdata/coding_model)
#' @keywords internal
#' @noRd
hexamer_table <- function() {
  cache_fetch(annotation_cache, "hexamer_table", "v1", function() {
    path <- system.file(
      "extdata",
      "coding_model",
      "hexamer_logodds.parquet",
      package = "isoformic"
    )
    tbl <- read_parquet(path)
    stats::setNames(tbl$log_odds, tbl$hexamer)
  })
}

#' Frozen logistic model coefficients (shipped in inst/extdata/coding_model)
#' @keywords internal
#' @noRd
coding_model <- function() {
  cache_fetch(annotation_cache, "coding_model", "v1", function() {
    path <- system.file(
      "extdata",
      "coding_model",
      "coding_model_coefficients.parquet",
      package = "isoformic"
    )
    tbl <- read_parquet(path)
    as.list(stats::setNames(tbl$estimate, tbl$term))
  })
}

#' Train the Coding Potential Model
#'
#' Builds the hexamer log-odds table and the logistic regression used by
#' [coding_potential()] from a set of coding ORF sequences and non-coding
#' transcript sequences, and writes both as Parquet files. The shipped model
#' was trained on GENCODE human protein-coding CDS versus lncRNA transcripts
#' by `data-raw/build_coding_model.R`; rerun that script to refresh it for a
#' new GENCODE release or another organism.
#'
#' @param coding_cds Character vector of CDS sequences (start to stop codon),
#'   used to build the hexamer table.
#' @param coding Character vector of full-length coding transcript
#'   sequences (mRNAs), the positive class of the logistic model.
#' @param noncoding Character vector of non-coding transcript sequences, the
#'   negative class (and the hexamer background).
#' @param output_dir Directory receiving `hexamer_logodds.parquet` and
#'   `coding_model_coefficients.parquet`.
#' @param min_orf_length Minimum ORF length used when scoring the training
#'   sequences.
#' @param pseudocount Pseudocount added to hexamer counts.
#'
#' @returns A list with the hexamer table, the fitted `glm` and the training
#'   accuracy, invisibly.
#'
#' @examples
#' \dontrun{
#' build_coding_model(cds_sequences, mrna_sequences, lncrna_sequences, output_dir = "coding_model")
#' }
#'
#' @export
build_coding_model <- function(
  coding_cds,
  coding,
  noncoding,
  output_dir,
  min_orf_length = 90L,
  pseudocount = 1
) {
  check_string(output_dir)
  fs::dir_create(output_dir)
  coding_cds <- toupper(coding_cds)
  coding <- toupper(coding)
  noncoding <- toupper(noncoding)
  hex_counts <- function(seqs, step) {
    counts <- integer()
    for (s in seqs) {
      n <- nchar(s)
      if (n < 6L) {
        next
      }
      starts <- seq.int(1L, n - 5L, by = step)
      hex <- substring(s, starts, starts + 5L)
      tab <- table(hex)
      counts[names(tab)] <- ifelse(
        is.na(counts[names(tab)]),
        0L,
        counts[names(tab)]
      ) +
        as.integer(tab)
    }
    counts
  }
  all_hex <- do.call(
    paste0,
    expand.grid(rep(list(c("A", "C", "G", "T")), 6L), stringsAsFactors = FALSE)
  )
  cod <- hex_counts(coding_cds, step = 3L)
  non <- hex_counts(noncoding, step = 1L)
  cod_v <- ifelse(is.na(cod[all_hex]), 0, cod[all_hex]) + pseudocount
  non_v <- ifelse(is.na(non[all_hex]), 0, non[all_hex]) + pseudocount
  log_odds <- log(cod_v / sum(cod_v)) - log(non_v / sum(non_v))
  hex_tbl <- tibble::tibble(hexamer = all_hex, log_odds = unname(log_odds))
  write_parquet(hex_tbl, fs::path(output_dir, "hexamer_logodds.parquet"))
  assign(
    "hexamer_table",
    list(
      signature = "v1",
      value = stats::setNames(hex_tbl$log_odds, hex_tbl$hexamer)
    ),
    envir = annotation_cache
  )
  features <- function(seqs) {
    stops <- find_stop_positions(seqs)
    longest <- longest_orfs(seqs, stops, min_orf_length = min_orf_length)
    sc <- coding_potential_scores_raw(seqs, longest)
    tibble::tibble(
      log_orf_length = log1p(ifelse(is.na(longest$length), 0, longest$length)),
      orf_coverage = sc$orf_coverage,
      fickett_score = sc$fickett_score,
      hexamer_score = sc$hexamer_score
    )
  }
  train <- dplyr::bind_rows(features(coding), features(noncoding))
  train$coding <- rep(c(1L, 0L), c(length(coding), length(noncoding)))
  fit <- stats::glm(
    coding ~ log_orf_length + orf_coverage + fickett_score + hexamer_score,
    data = train,
    family = stats::binomial()
  )
  coefs <- stats::coef(fit)
  names(coefs)[names(coefs) == "(Intercept)"] <- "intercept"
  coef_tbl <- tibble::tibble(term = names(coefs), estimate = unname(coefs))
  write_parquet(
    coef_tbl,
    fs::path(output_dir, "coding_model_coefficients.parquet")
  )
  assign(
    "coding_model",
    list(signature = "v1", value = as.list(coefs)),
    envir = annotation_cache
  )
  pred <- stats::predict(fit, type = "response") > 0.5
  accuracy <- mean(pred == (train$coding == 1L))
  cli::cli_inform(c(
    v = "Coding model written to {.path {output_dir}} (training accuracy {round(accuracy, 3)})."
  ))
  invisible(list(hexamers = hex_tbl, model = fit, accuracy = accuracy))
}

#' Raw coding potential features without the logistic model
#' @keywords internal
#' @noRd
coding_potential_scores_raw <- function(seq, longest) {
  seq_len <- nchar(seq)
  has_orf <- !is.na(longest$start)
  orf_seq <- rep("", length(seq))
  orf_seq[has_orf] <- substring(
    seq[has_orf],
    longest$start[has_orf],
    longest$end[has_orf]
  )
  list(
    orf_coverage = ifelse(has_orf, longest$length / seq_len, 0),
    fickett_score = fickett_score(ifelse(has_orf, orf_seq, seq)),
    hexamer_score = hexamer_score(orf_seq)
  )
}
