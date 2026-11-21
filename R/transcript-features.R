# Per-transcript structural and coding features ------------------------------

#' Structural and Coding Features of Transcripts
#'
#' @description
#' Computes, for every assayed transcript, features that go beyond the
#' biotype label, from the exon and CDS coordinates of the annotation:
#'
#' * **Ends**: `tss_pos`, `tes_pos` (genomic positions of the transcription
#'   start and end), `tss_group` and `tes_group` (integer labels of
#'   alternative promoters / 3' ends within the gene; starts or ends closer
#'   than `end_window` bp share a group).
#' * **ORF**: `has_cds`, `cds_complete` (`FALSE` when the annotation flags
#'   an incomplete CDS with `cds_start_NF`/`cds_end_NF`), `orf_length`
#'   (nt), `utr5_length`, `utr3_length`.
#' * **NMD rule**: `stop_to_last_junction`, the distance in transcript
#'   coordinates from the stop codon to the last exon-exon junction, and
#'   `nmd_predicted` (`TRUE` when that distance exceeds `nmd_rule` nt, the
#'   classical 50-55 nt rule).
#' * **Comparison with the canonical isoform**: `canonical_transcript_id`
#'   (the MANE Select / Ensembl canonical transcript of the gene, else the
#'   transcript with the longest CDS, else the longest transcript),
#'   `orf_fraction_of_canonical` and `cds_differs_from_canonical`.
#' * **Splicing events relative to the canonical isoform** (see
#'   [splicing_events()]): `n_skipped_exons`, `n_novel_exons`,
#'   `n_retained_introns`, `n_alt_5ss`, `n_alt_3ss`, `alt_first_exon`,
#'   `alt_last_exon`.
#'
#' `splicing_events()` returns the event counts alone; both are cached per
#' object. With an annotation directory path every transcript of the
#' annotation is described.
#'
#' @param self An [IsoformicExperiment] or an annotation directory path.
#' @param gene_name Optional genes (names or IDs) to restrict the result to.
#' @param nmd_rule Distance in nucleotides from the stop codon to the last
#'   junction above which NMD is predicted (default 50).
#' @param end_window Distance in bp within which transcription starts (or
#'   ends) are grouped as the same promoter (or 3' end); default 50.
#' @param ... Not used; reserved for methods.
#'
#' @returns A tibble with one row per transcript. `transcript_features()`:
#'   `transcript_id`, `transcript_name`, `transcript_type`, `gene_id`,
#'   `gene_name`, `strand`, `n_exons`, `tx_length`, `is_canonical` and the
#'   feature columns listed above. `splicing_events()`: `transcript_id`,
#'   `transcript_name`, `gene_id`, `gene_name`, `canonical_transcript_id`
#'   and the event columns.
#'
#' @examples
#' iso <- example_isoformic()
#' feats <- transcript_features(iso, gene_name = "RBPJ")
#' feats[, c("transcript_name", "orf_length", "nmd_predicted", "tss_group", "n_skipped_exons")]
#' splicing_events(iso, gene_name = "EGFR")
#'
#' @name transcript_features
#' @rdname transcript_features
NULL

#' @rdname transcript_features
#' @export
transcript_features <- S7::new_generic(
  "transcript_features",
  "self",
  function(self, gene_name = NULL, nmd_rule = 50, end_window = 50, ...) {
    S7::S7_dispatch()
  }
)

S7::method(transcript_features, IsoformicExperiment) <- function(
  self,
  gene_name = NULL,
  nmd_rule = 50,
  end_window = 50,
  ...
) {
  signature <- paste(
    annot_full_signature(self),
    rlang::hash(row_names(self)),
    nmd_rule,
    end_window
  )
  res <- object_cache_fetch(self, "transcript_features", signature, function() {
    compute_transcript_features(
      transcripts = row_data(self),
      exons = self@annot_data_exons,
      nmd_rule = nmd_rule,
      end_window = end_window
    )
  })
  filter_genes(res, gene_name)
}

S7::method(transcript_features, S7::class_character) <- function(
  self,
  gene_name = NULL,
  nmd_rule = 50,
  end_window = 50,
  ...
) {
  # A GFF3/GTF file is prepared first, as in prepare_annotation()
  annot_dir <- if (is_annotation_dir(self)) {
    self
  } else {
    prepare_isoformic_annotation(input_path = self)
  }
  annot <- annotation_of_genes(annot_dir, gene_name)
  res <- compute_transcript_features(
    annot$transcripts,
    annot$exons(),
    nmd_rule,
    end_window
  )
  filter_genes(res, gene_name)
}

#' Transcripts and exons of selected genes of an annotation directory
#'
#' The features and ORFs of a transcript depend only on the transcripts of
#' its gene (the canonical isoform is chosen within the gene), so a call
#' restricted to some genes reads only their transcripts and exons, selected
#' in DuckDB, instead of the whole annotation. The genes are matched as in
#' `filter_genes()`: by gene name or gene ID. When no gene matches, the
#' transcripts of the first gene of the annotation are returned, so that the
#' result filtered afterwards is an empty table with the usual columns.
#'
#' @returns A list with `transcripts` (with `gene_name` from the genes
#'   table), `annotation` (the same transcripts without `gene_name`) and
#'   `exons`, a function reading the exons of those transcripts.
#' @keywords internal
#' @noRd
annotation_of_genes <- function(annot_dir, gene_name = NULL) {
  genes <- read_annotation_table(
    annot_dir,
    "genes",
    columns = c("gene_id", "gene_name")
  )
  if (is.null(gene_name) || anyNA(gene_name)) {
    transcripts <- read_annotation_table(annot_dir, "transcripts")
    exons <- function() read_annotation_table(annot_dir, "exons")
  } else {
    gene_name <- as.character(gene_name)
    wanted <- unique(c(
      gene_name,
      genes$gene_id[genes$gene_name %in% gene_name]
    ))
    transcripts <- read_annotation_table(
      annot_dir,
      "transcripts",
      filter = list(gene_id = wanted)
    )
    if (!nrow(transcripts)) {
      first <- read_parquet(
        annot_file_path(annot_dir, "transcripts"),
        columns = "gene_id",
        n_max = 1L
      )
      transcripts <- read_annotation_table(
        annot_dir,
        "transcripts",
        filter = list(gene_id = utils::head(first$gene_id, 1L))
      )
    }
    ids <- transcripts$transcript_id
    exons <- function() {
      read_annotation_table(
        annot_dir,
        "exons",
        filter = list(transcript_id = ids)
      )
    }
  }
  list(
    transcripts = dplyr::left_join(transcripts, genes, by = "gene_id"),
    annotation = transcripts,
    exons = exons
  )
}

#' @rdname transcript_features
#' @export
splicing_events <- S7::new_generic(
  "splicing_events",
  "self",
  function(self, gene_name = NULL, ...) {
    S7::S7_dispatch()
  }
)

S7::method(splicing_events, IsoformicExperiment) <- function(
  self,
  gene_name = NULL,
  ...
) {
  feats <- transcript_features(self, gene_name = gene_name)
  select_event_cols(feats)
}

S7::method(splicing_events, S7::class_character) <- function(
  self,
  gene_name = NULL,
  ...
) {
  feats <- transcript_features(self, gene_name = gene_name)
  select_event_cols(feats)
}

select_event_cols <- function(feats) {
  dplyr::select(
    feats,
    dplyr::any_of(c(
      "transcript_id",
      "transcript_name",
      "gene_id",
      "gene_name",
      "canonical_transcript_id",
      "n_skipped_exons",
      "n_novel_exons",
      "n_retained_introns",
      "n_alt_5ss",
      "n_alt_3ss",
      "alt_first_exon",
      "alt_last_exon"
    ))
  )
}

#' Compute transcript features from transcript and exon tables
#' @keywords internal
#' @noRd
compute_transcript_features <- function(
  transcripts,
  exons,
  nmd_rule = 50,
  end_window = 50
) {
  tx <- transcripts
  # Annotation directories created before CDS aggregation lack these columns
  for (col in c("cds_start", "cds_end", "cds_length")) {
    if (!col %in% colnames(tx)) {
      tx[[col]] <- NA_integer_
    }
  }
  needed <- c("transcript_id", "gene_id", "strand")
  if (!all(needed %in% colnames(tx))) {
    strand_by_tx <- exons |>
      dplyr::distinct(.data$transcript_id, .keep_all = TRUE) |>
      dplyr::select(dplyr::any_of(c("transcript_id", "strand")))
    tx <- dplyr::left_join(tx, strand_by_tx, by = "transcript_id")
  }
  exons <- exon_cum_table(exons, tx)
  ex_strand <- exons$strand
  ex_len <- exons$ex_len
  tx_id <- exons$transcript_id
  cum_before <- exons$cum_before
  idx_in_tx <- exons$idx_in_tx
  n_ex_vec <- as.vector(rowsum(rep(1L, length(tx_id)), tx_id))
  names(n_ex_vec) <- sort(unique(tx_id))
  n_ex <- unname(n_ex_vec[tx_id])
  tx_len_vec <- as.vector(rowsum(ex_len, tx_id))
  names(tx_len_vec) <- sort(unique(tx_id))
  is_last <- idx_in_tx == n_ex
  last_len <- stats::setNames(ex_len[is_last], tx_id[is_last])

  # Genomic -> transcript coordinate of the CDS ends (DuckDB containment join)
  minus <- tx$strand == "-"
  orf_first <- ifelse(minus, tx$cds_end, tx$cds_start) # 5' end of the ORF
  orf_last <- ifelse(minus, tx$cds_start, tx$cds_end) # 3' end of the ORF
  to_tx_coord <- genomic_to_transcript_coord(
    exon_tbl = data.frame(
      transcript_id = tx_id,
      start_pos = exons$start_pos,
      end_pos = exons$end_pos,
      strand = ex_strand,
      cum_before = cum_before,
      stringsAsFactors = FALSE
    ),
    transcript_id = tx$transcript_id,
    first_pos = orf_first,
    last_pos = orf_last
  )
  orf_first_tx <- to_tx_coord$first
  orf_last_tx <- to_tx_coord$last

  tx_length <- unname(tx_len_vec[tx$transcript_id])
  last_exon_len <- unname(last_len[tx$transcript_id])
  n_exons <- unname(n_ex_vec[tx$transcript_id])
  has_cds <- !is.na(tx$cds_length) & tx$cds_length > 0
  tags <- if ("tags" %in% colnames(tx)) {
    tx$tags
  } else {
    rep(NA_character_, nrow(tx))
  }
  tags[is.na(tags)] <- ""
  cds_complete <- has_cds & !grepl("cds_(start|end)_NF", tags)
  last_junction_tx <- tx_length - last_exon_len
  stop_to_last_junction <- ifelse(
    has_cds & n_exons > 1L,
    last_junction_tx - orf_last_tx,
    NA_real_
  )
  nmd_predicted <- has_cds &
    !is.na(stop_to_last_junction) &
    stop_to_last_junction > nmd_rule

  is_first <- idx_in_tx == 1L
  first_exon_pos <- ifelse(
    ex_strand[is_first] == "-",
    exons$end_pos[is_first],
    exons$start_pos[is_first]
  )
  names(first_exon_pos) <- tx_id[is_first]
  last_exon_pos <- ifelse(
    ex_strand[is_last] == "-",
    exons$start_pos[is_last],
    exons$end_pos[is_last]
  )
  names(last_exon_pos) <- tx_id[is_last]
  tss_pos <- unname(first_exon_pos[tx$transcript_id])
  tes_pos <- unname(last_exon_pos[tx$transcript_id])

  canonical <- canonical_transcripts(tx)
  canonical_id <- canonical[tx$gene_id]
  canonical_cds <- tx$cds_length[match(canonical_id, tx$transcript_id)]
  canonical_start <- tx$cds_start[match(canonical_id, tx$transcript_id)]
  canonical_end <- tx$cds_end[match(canonical_id, tx$transcript_id)]
  orf_fraction <- ifelse(
    has_cds & !is.na(canonical_cds) & canonical_cds > 0,
    tx$cds_length / canonical_cds,
    NA_real_
  )
  cds_differs <- ifelse(
    has_cds & !is.na(canonical_start),
    tx$cds_start != canonical_start |
      tx$cds_end != canonical_end |
      tx$cds_length != canonical_cds,
    NA
  )

  events <- compute_splicing_events(tx, exons, canonical)

  res <- tibble::tibble(
    transcript_id = tx$transcript_id,
    transcript_name = if ("transcript_name" %in% colnames(tx)) {
      tx$transcript_name
    } else {
      tx$transcript_id
    },
    transcript_type = if ("transcript_type" %in% colnames(tx)) {
      tx$transcript_type
    } else {
      NA_character_
    },
    gene_id = tx$gene_id,
    gene_name = if ("gene_name" %in% colnames(tx)) tx$gene_name else tx$gene_id,
    strand = tx$strand,
    n_exons = n_exons,
    tx_length = tx_length,
    tss_pos = tss_pos,
    tes_pos = tes_pos,
    tss_group = end_groups(tx$gene_id, tss_pos, end_window),
    tes_group = end_groups(tx$gene_id, tes_pos, end_window),
    has_cds = has_cds,
    cds_complete = cds_complete,
    orf_length = ifelse(has_cds, tx$cds_length, NA_integer_),
    utr5_length = ifelse(has_cds, orf_first_tx - 1, NA_real_),
    utr3_length = ifelse(has_cds, tx_length - orf_last_tx, NA_real_),
    stop_to_last_junction = stop_to_last_junction,
    nmd_predicted = nmd_predicted,
    canonical_transcript_id = unname(canonical_id),
    is_canonical = tx$transcript_id == unname(canonical_id),
    orf_fraction_of_canonical = orf_fraction,
    cds_differs_from_canonical = cds_differs
  )
  dplyr::left_join(res, events, by = "transcript_id")
}

#' Exons of `tx` ordered along each transcript with cumulative offsets
#'
#' Returns the exon rows (transcript_id, start_pos, end_pos, strand) sorted by
#' transcript and by position along the strand, with `ex_len`, `cum_before`
#' (transcript length upstream of the exon) and `idx_in_tx` (1-based exon
#' rank from the TSS). This is the table all coordinate conversions use.
#'
#' @keywords internal
#' @noRd
exon_cum_table <- function(exons, tx) {
  exons <- exons[exons$transcript_id %in% tx$transcript_id, ]
  ex_strand <- tx$strand[match(exons$transcript_id, tx$transcript_id)]
  ord <- order(
    exons$transcript_id,
    ifelse(ex_strand == "-", -exons$start_pos, exons$start_pos)
  )
  exons <- exons[ord, ]
  exons$strand <- ex_strand[ord]
  ex_len <- exons$end_pos - exons$start_pos + 1L
  tx_id <- exons$transcript_id
  first_idx <- match(tx_id, tx_id)
  cs <- cumsum(as.numeric(ex_len))
  exons$ex_len <- ex_len
  exons$cum_before <- cs - ex_len - (cs - ex_len)[first_idx]
  exons$idx_in_tx <- seq_along(tx_id) - first_idx + 1L
  exons
}

#' Convert genomic positions to transcript coordinates, one position per row
#'
#' Unlike `genomic_to_transcript_coord()`, `transcript_id` may repeat (for
#' example the canonical transcript of every isoform of a gene); results are
#' matched by row index.
#'
#' @keywords internal
#' @noRd
genomic_to_transcript_pos <- function(
  exon_tbl,
  transcript_id,
  pos,
  con = NULL
) {
  out <- rep(NA_real_, length(transcript_id))
  pos_tbl <- data.frame(
    idx = seq_along(transcript_id),
    transcript_id = transcript_id,
    pos = pos,
    stringsAsFactors = FALSE
  )
  pos_tbl <- pos_tbl[!is.na(pos_tbl$pos) & !is.na(pos_tbl$transcript_id), ]
  if (!nrow(pos_tbl)) {
    return(out)
  }
  res <- query_registered(
    con,
    list(pos = pos_tbl, ex = as.data.frame(exon_tbl)),
    "SELECT p.idx,
            min(e.cum_before
                + CASE WHEN e.strand = '-' THEN e.end_pos - p.pos ELSE p.pos - e.start_pos END
                + 1) AS tx_pos
     FROM pos p JOIN ex e
       ON e.transcript_id = p.transcript_id AND p.pos BETWEEN e.start_pos AND e.end_pos
     GROUP BY p.idx"
  )
  out[res$idx] <- res$tx_pos
  out
}

#' Convert transcript coordinates (1-based from the TSS) to genomic positions
#' @keywords internal
#' @noRd
transcript_to_genomic_coord <- function(
  exon_tbl,
  transcript_id,
  pos,
  con = NULL
) {
  out <- rep(NA_real_, length(transcript_id))
  pos_tbl <- data.frame(
    idx = seq_along(transcript_id),
    transcript_id = transcript_id,
    pos = pos,
    stringsAsFactors = FALSE
  )
  pos_tbl <- pos_tbl[!is.na(pos_tbl$pos) & !is.na(pos_tbl$transcript_id), ]
  if (!nrow(pos_tbl)) {
    return(out)
  }
  res <- query_registered(
    con,
    list(pos = pos_tbl, ex = as.data.frame(exon_tbl)),
    "SELECT p.idx,
            min(CASE WHEN e.strand = '-'
                     THEN e.end_pos - (p.pos - e.cum_before - 1)
                     ELSE e.start_pos + (p.pos - e.cum_before - 1) END) AS gpos
     FROM pos p JOIN ex e
       ON e.transcript_id = p.transcript_id
      AND p.pos > e.cum_before AND p.pos <= e.cum_before + e.ex_len
     GROUP BY p.idx"
  )
  out[res$idx] <- res$gpos
  out
}

#' Canonical transcript per gene: tagged canonical, else longest CDS, else longest
#' @keywords internal
#' @noRd
canonical_transcripts <- function(tx) {
  tagged <- if ("is_canonical" %in% colnames(tx)) {
    tx$is_canonical %in% TRUE
  } else {
    rep(FALSE, nrow(tx))
  }
  cds <- if ("cds_length" %in% colnames(tx)) {
    ifelse(is.na(tx$cds_length), 0, tx$cds_length)
  } else {
    rep(0, nrow(tx))
  }
  len <- if ("tx_length" %in% colnames(tx)) tx$tx_length else rep(0, nrow(tx))
  ord <- order(tx$gene_id, -as.integer(tagged), -cds, -len, tx$transcript_id)
  first <- ord[!duplicated(tx$gene_id[ord])]
  stats::setNames(tx$transcript_id[first], tx$gene_id[first])
}

#' Group positions within a gene that lie within `window` bp of each other
#' @keywords internal
#' @noRd
end_groups <- function(gene_id, pos, window) {
  out <- rep(NA_integer_, length(pos))
  ok <- !is.na(pos)
  ord <- order(gene_id[ok], pos[ok])
  g <- gene_id[ok][ord]
  p <- pos[ok][ord]
  new_group <- c(
    TRUE,
    g[-1] != g[-length(g)] | (p[-1] - p[-length(p)]) > window
  )
  grp <- cumsum(new_group)
  # renumber within gene: subtract the group number of the gene's first position
  grp <- grp - grp[match(g, g)] + 1L
  out[ok][ord] <- as.integer(grp)
  out
}

#' Splicing events of each transcript relative to its gene's canonical isoform
#'
#' Overlap join in DuckDB between the exons of each transcript and the exons
#' of the canonical isoform of its gene. The keys are integers (the row of
#' the transcript in `tx`, and the row of its canonical transcript), and the
#' rank of each exon, the exon count and the span of each transcript are
#' computed in R, so DuckDB hashes no strings and runs no window function.
#' The overlap join runs once: novel exons (transcript exons without an
#' overlapping canonical exon) and skipped exons (canonical exons inside the
#' transcript span without an overlapping transcript exon) are anti joins
#' with its result.
#'
#' @keywords internal
#' @noRd
compute_splicing_events <- function(tx, exons, canonical_id) {
  counts <- c(
    "n_skipped_exons",
    "n_novel_exons",
    "n_retained_introns",
    "n_alt_5ss",
    "n_alt_3ss"
  )
  flags <- c("alt_first_exon", "alt_last_exon")
  out <- tibble::tibble(transcript_id = tx$transcript_id)
  for (col in counts) {
    out[[col]] <- rep(0L, nrow(tx))
  }
  for (col in flags) {
    out[[col]] <- rep(FALSE, nrow(tx))
  }
  if (!nrow(exons)) {
    return(out)
  }
  # Integer keys: the first row of each transcript ID in `tx`, as a string
  # join on transcript_id would match it
  tid <- match(exons$transcript_id, tx$transcript_id, incomparables = NA)
  canon_tid <- match(
    unname(canonical_id),
    tx$transcript_id,
    incomparables = NA
  )
  gene_canon <- canon_tid[match(
    tx$gene_id[tid],
    names(canonical_id),
    incomparables = NA
  )]
  is_cex <- !is.na(tid) & tid %in% canon_tid[!is.na(canon_tid)]
  is_tex <- !is.na(tid) & !is.na(gene_canon) & tid != gene_canon
  strand <- match(tx$strand[tid], c("+", "-"))
  strand <- c(1L, -1L)[strand]
  strand[is.na(strand)] <- 0L
  # Exons of the selected rows sorted by transcript and start, with the rank
  # of each exon from the 5'-most one, the exon count and the span of the
  # transcript
  exon_rank <- function(keep) {
    rows <- which(keep)
    rows <- rows[order(tid[rows], exons$start_pos[rows])]
    id <- tid[rows]
    start <- exons$start_pos[rows]
    end <- exons$end_pos[rows]
    first <- match(id, id)
    by_end <- order(id, -end)
    last_end <- end[by_end][!duplicated(id[by_end])]
    data.frame(
      row = rows,
      tid = id,
      rk = seq_along(id) - first + 1L,
      n = tabulate(id, nbins = nrow(tx))[id],
      start_pos = start,
      end_pos = end,
      tx_start = start[first],
      tx_end = last_end[match(id, unique(id[by_end]))]
    )
  }
  cex <- exon_rank(is_cex)[, c("tid", "rk", "n", "start_pos", "end_pos")]
  tex <- exon_rank(is_tex)
  tex$ctid <- gene_canon[tex$row]
  tex$strand <- strand[tex$row]
  tex$row <- NULL
  res <- query_registered(
    con = NULL,
    tables = list(tex = tex, cex = cex),
    sql = "WITH ov AS MATERIALIZED (
       -- transcript exons overlapping canonical exons
       SELECT t.tid, t.rk AS t_rk, t.n AS t_n, t.start_pos AS t_start,
              t.end_pos AS t_end, c.rk AS c_rk, c.n AS c_n,
              c.start_pos AS c_start, c.end_pos AS c_end, t.strand
       FROM tex t JOIN cex c
         ON c.tid = t.ctid
        AND t.start_pos <= c.end_pos AND t.end_pos >= c.start_pos
     ),
     novel AS (
       SELECT t.tid, count(*) AS n_novel_exons
       FROM tex t ANTI JOIN ov o ON o.tid = t.tid AND o.t_rk = t.rk
       GROUP BY t.tid
     ),
     span AS (
       SELECT DISTINCT tid, ctid, tx_start, tx_end FROM tex
     ),
     skipped AS (
       -- canonical exons inside the transcript span that it does not overlap
       SELECT s.tid, count(DISTINCT c.rk) AS n_skipped_exons
       FROM span s JOIN cex c
         ON c.tid = s.ctid
        AND c.start_pos >= s.tx_start AND c.end_pos <= s.tx_end
       ANTI JOIN ov o ON o.tid = s.tid AND o.c_rk = c.rk
       GROUP BY s.tid
     ),
     retained AS (
       -- a transcript exon spanning two consecutive canonical exons
       SELECT o1.tid, count(DISTINCT o1.c_rk) AS n_retained_introns
       FROM ov o1 JOIN ov o2
         ON o1.tid = o2.tid AND o1.t_rk = o2.t_rk AND o2.c_rk = o1.c_rk + 1
       GROUP BY o1.tid
     ),
     alt AS (
       SELECT tid,
         sum(CASE WHEN (strand = 1 AND t_rk > 1 AND t_start <> c_start)
                    OR (strand = -1 AND t_rk < t_n AND t_end <> c_end)
                  THEN 1 ELSE 0 END) AS n_alt_3ss,
         sum(CASE WHEN (strand = 1 AND t_rk < t_n AND t_end <> c_end)
                    OR (strand = -1 AND t_rk > 1 AND t_start <> c_start)
                  THEN 1 ELSE 0 END) AS n_alt_5ss,
         max(CASE WHEN (strand = 1 AND t_rk = 1 AND c_rk = 1)
                    OR (strand = -1 AND t_rk = t_n AND c_rk = c_n)
                  THEN 1 ELSE 0 END) AS first_shared,
         max(CASE WHEN (strand = 1 AND t_rk = t_n AND c_rk = c_n)
                    OR (strand = -1 AND t_rk = 1 AND c_rk = 1)
                  THEN 1 ELSE 0 END) AS last_shared
       FROM ov GROUP BY tid
     )
     SELECT t.tid,
            coalesce(s.n_skipped_exons, 0)::INTEGER AS n_skipped_exons,
            coalesce(nv.n_novel_exons, 0)::INTEGER AS n_novel_exons,
            coalesce(r.n_retained_introns, 0)::INTEGER AS n_retained_introns,
            coalesce(a.n_alt_5ss, 0)::INTEGER AS n_alt_5ss,
            coalesce(a.n_alt_3ss, 0)::INTEGER AS n_alt_3ss,
            coalesce(a.first_shared, 0) = 0 AS alt_first_exon,
            coalesce(a.last_shared, 0) = 0 AS alt_last_exon
     FROM (SELECT DISTINCT tid FROM tex) t
     LEFT JOIN skipped s USING (tid)
     LEFT JOIN novel nv USING (tid)
     LEFT JOIN retained r USING (tid)
     LEFT JOIN alt a USING (tid)"
  )
  # Rows with a duplicated transcript ID get the events of its first row;
  # canonical transcripts keep no events
  hit <- match(match(tx$transcript_id, tx$transcript_id), res$tid)
  keep <- !is.na(hit) & !(tx$transcript_id %in% unname(canonical_id))
  for (col in c(counts, flags)) {
    out[[col]][keep] <- res[[col]][hit[keep]]
  }
  out
}

# Annotation columns from external tables -----------------------------------

#' Attach Extra Columns to the Annotation Tables
#'
#' @description
#' Joins columns of an external per-transcript or per-gene table into the
#' annotation directory (persisted in the Parquet tables, so every object
#' using the annotation sees them and caches refresh automatically). This is
#' the hook for evidence that the package does not compute itself: FANTOM
#' CAGE clusters, GTEx dominant isoforms per tissue, Ribo-seq translation
#' calls, NMD-inhibition responses, protein domain retention, long-read
#' full-length counts. A column named `productivity_evidence` holding one of
#' the classes of [productivity_palette()] (`"productive"`,
#' `"altered_productive"`, `"unproductive"`, `"noncoding"`, `"uncertain"`)
#' overrides the call of [transcript_productivity()]; other values are
#' ignored.
#'
#' @param x An annotation directory path or an [IsoformicExperiment].
#' @param table A data frame, or a path to a tabular file readable by
#'   [read_tabular()].
#' @param level `"transcript"` (default) or `"gene"`.
#' @param by Name of the key column in `table`; defaults to
#'   `transcript_id`/`gene_id`, or the first column.
#' @param prefix Optional prefix added to the new column names.
#'
#' @details
#' Only the first row of each key of `table` is used. An annotation column
#' with the name of a new column is replaced, and rows of the annotation
#' without a match get `NA`.
#'
#' @returns The annotation directory path, invisibly.
#'
#' @examples
#' gff <- system.file("extdata", "example.gff3.gz", package = "isoformic")
#' annot <- prepare_isoformic_annotation(gff, output_path = tempdir(), annot_name = "evidence_demo")
#' tx <- read_parquet(file.path(annot, "annot_data_transcripts.parquet"))
#' evidence <- data.frame(
#'   transcript_id = tx$transcript_id[1:2],
#'   riboseq_translated = c(TRUE, FALSE),
#'   productivity_evidence = c("productive", "unproductive")
#' )
#' add_annotation_columns(annot, evidence)
#' tx <- read_parquet(file.path(annot, "annot_data_transcripts.parquet"))
#' head(tx[, c("transcript_id", "riboseq_translated")], 3)
#'
#' @export
add_annotation_columns <- function(
  x,
  table,
  level = c("transcript", "gene"),
  by = NULL,
  prefix = NULL
) {
  level <- rlang::arg_match(level)
  annot_dir <- if (S7::S7_inherits(x, IsoformicExperiment)) x@annot_path else x
  if (!is_annotation_dir(annot_dir)) {
    cli::cli_abort(
      message = c(
        x = "{.arg x} must be an annotation directory or an {.cls IsoformicExperiment}."
      ),
      class = "isoformic_annot_path_dont_exist"
    )
  }
  if (!is.data.frame(table)) {
    table <- read_tabular(table)
  }
  table <- tibble::as_tibble(table)
  key <- paste0(level, "_id")
  by <- by %||% if (key %in% colnames(table)) key else colnames(table)[1L]
  if (!by %in% colnames(table)) {
    cli::cli_abort(
      message = c(x = "Column {.field {by}} not found in {.arg table}."),
      class = "isoformic_invalid_evidence"
    )
  }
  table[[by]] <- as.character(table[[by]])
  table <- dplyr::distinct(table, .data[[by]], .keep_all = TRUE)
  new_cols <- setdiff(colnames(table), by)
  if (!is.null(prefix)) {
    names(table)[names(table) %in% new_cols] <- paste0(prefix, new_cols)
    new_cols <- paste0(prefix, new_cols)
  }
  type <- if (identical(level, "transcript")) "transcripts" else "genes"
  annot <- read_annotation_table(annot_dir, type)
  n_hit <- sum(annot[[key]] %in% table[[by]])
  if (!n_hit) {
    cli::cli_abort(
      message = c(
        x = "No {level} of the annotation matches the {.field {by}} column of {.arg table}."
      ),
      class = "isoformic_invalid_evidence"
    )
  }
  annot <- annot[, setdiff(colnames(annot), new_cols)]
  names(table)[names(table) == by] <- key
  annot <- dplyr::left_join(annot, table, by = key)
  write_parquet(annot, annot_file_path(annot_dir, type))
  cli::cli_inform(c(
    v = paste0(
      "Added {length(new_cols)} column{?s} ({.field {new_cols}}) ",
      "for {n_hit} of {nrow(annot)} {type} in {.path {annot_dir}}."
    )
  ))
  invisible(annot_dir)
}

#' Attach FANTOM CAT Tables to an Annotation
#'
#' FANTOM CAT distributes, next to the GTF, an ID mapping table
#' (`geneID`, `transcriptID`, `CAGEClusterID`) and a gene information table
#' (`geneClass`, `geneCategory`, `TIEScore_tier`, `splicingType`,
#' `directionality`, transposable element associations, ...). This function
#' joins them into the annotation directory: the CAGE cluster (the
#' CAGE-defined promoter of each transcript) becomes `cage_cluster_id`, with
#' its coordinates parsed into `cage_cluster_seqid`, `cage_cluster_start`,
#' `cage_cluster_end`, and the gene information columns are added to the
#' gene table in snake case (except the location, strand and the attributes
#' the GTF already provides, such as `gene_class` and `gene_category`). See
#' [download_fantom_cat_tables()].
#'
#' @param x An annotation directory path or an [IsoformicExperiment].
#' @param id_mapping Path (or data frame) of the ID mapping table.
#' @param gene_info Optional path (or data frame) of the gene info table.
#'
#' @returns The annotation directory path, invisibly.
#'
#' @examples
#' \dontrun{
#' gtf <- download_reference(
#'   reference = "fantom",
#'   version = "lv3_robust",
#'   file_type = "gtf"
#' )
#' tables <- download_fantom_cat_tables("lv3_robust")
#' annot <- prepare_isoformic_annotation(gtf)
#' add_fantom_cat_tables(annot, tables[["id_mapping"]], tables[["gene_info"]])
#' }
#'
#' @export
add_fantom_cat_tables <- function(x, id_mapping, gene_info = NULL) {
  annot_dir <- if (S7::S7_inherits(x, IsoformicExperiment)) x@annot_path else x
  mapping <- if (is.data.frame(id_mapping)) {
    tibble::as_tibble(id_mapping)
  } else {
    read_tabular(id_mapping)
  }
  required <- c("transcriptID", "CAGEClusterID")
  if (!all(required %in% colnames(mapping))) {
    cli::cli_abort(
      message = c(
        x = "The FANTOM ID mapping table needs columns {.field {required}}."
      ),
      class = "isoformic_invalid_evidence"
    )
  }
  cluster <- as.character(mapping$CAGEClusterID)
  parts <- regmatches(
    cluster,
    regexec("^([^:]+):([0-9]+)\\.\\.([0-9]+),([+-])$", cluster)
  )
  get <- function(i) {
    vapply(
      parts,
      function(p) if (length(p) == 5L) p[[i]] else NA_character_,
      character(1L)
    )
  }
  tx_table <- tibble::tibble(
    transcript_id = as.character(mapping$transcriptID),
    cage_cluster_id = cluster,
    cage_cluster_seqid = get(2L),
    cage_cluster_start = as.integer(get(3L)),
    cage_cluster_end = as.integer(get(4L))
  )
  add_annotation_columns(annot_dir, tx_table, level = "transcript")
  if (!is.null(gene_info)) {
    info <- if (is.data.frame(gene_info)) {
      tibble::as_tibble(gene_info)
    } else {
      read_tabular(gene_info)
    }
    if (!"geneID" %in% colnames(info)) {
      cli::cli_abort(
        message = c(
          x = "The FANTOM gene info table needs a {.field geneID} column."
        ),
        class = "isoformic_invalid_evidence"
      )
    }
    names(info) <- to_snake_case(names(info))
    names(info)[names(info) == "gene_id"] <- "gene_id"
    drop <- c(
      "loc",
      "strand",
      "gene_name",
      "gene_class",
      "gene_category",
      "dhs_type",
      "gene_type",
      "coding_status",
      "cumulative_support"
    )
    info <- info[, setdiff(colnames(info), drop), drop = FALSE]
    add_annotation_columns(annot_dir, info, level = "gene", by = "gene_id")
  }
  invisible(annot_dir)
}

to_snake_case <- function(x) {
  x <- gsub("([A-Z]+)([A-Z][a-z])", "\\1_\\2", x)
  x <- gsub("([a-z0-9])([A-Z])", "\\1_\\2", x)
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- tolower(x)
  gsub("_+", "_", x)
}

#' Convert genomic positions to transcript coordinates (1-based from the TSS)
#'
#' Containment join in DuckDB between positions and the exons of the same
#' transcript; `cum_before` is the transcript length upstream of each exon.
#'
#' @keywords internal
#' @noRd
genomic_to_transcript_coord <- function(
  exon_tbl,
  transcript_id,
  first_pos,
  last_pos,
  con = NULL
) {
  pos_tbl <- data.frame(
    transcript_id = c(transcript_id, transcript_id),
    which = rep(c("first", "last"), each = length(transcript_id)),
    pos = c(first_pos, last_pos),
    stringsAsFactors = FALSE
  )
  pos_tbl <- pos_tbl[!is.na(pos_tbl$pos), ]
  out <- list(
    first = rep(NA_real_, length(transcript_id)),
    last = rep(NA_real_, length(transcript_id))
  )
  if (!nrow(pos_tbl)) {
    return(out)
  }
  res <- query_registered(
    con,
    list(pos = pos_tbl, ex = exon_tbl),
    "SELECT p.transcript_id, p.which,
            min(e.cum_before
                + CASE WHEN e.strand = '-' THEN e.end_pos - p.pos ELSE p.pos - e.start_pos END
                + 1) AS tx_pos
     FROM pos p JOIN ex e
       ON e.transcript_id = p.transcript_id AND p.pos BETWEEN e.start_pos AND e.end_pos
     GROUP BY p.transcript_id, p.which"
  )
  for (w in c("first", "last")) {
    sub <- res[res$which == w, ]
    out[[w]][match(sub$transcript_id, transcript_id)] <- sub$tx_pos
  }
  out
}
