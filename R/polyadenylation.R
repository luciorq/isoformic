# Alternative polyadenylation ------------------------------------------------
#
# Alternative cleavage and polyadenylation (APA) changes the 3' end of a
# transcript without touching its ORF: tandem 3' UTR shortening/lengthening
# removes or adds miRNA and RBP sites, intronic and alternative-last-exon
# APA changes the protein. isoformic already groups transcript ends
# (tes_group); this module classifies the events, validates the ends against
# poly(A) signals in the sequence and poly(A) site atlases, and quantifies
# proximal/distal usage per sample and group.

#' Poly(A) Sites and Alternative Polyadenylation of Transcripts
#'
#' @description
#' `polya_sites()` describes the 3' end of every assayed transcript:
#'
#' * `tes_pos`, `tes_group` (alternative 3' ends of the gene, from
#'   [transcript_features()]), `tes_rank` (1 = most proximal 3' end along
#'   the direction of transcription), `n_tes_groups`, `is_proximal`,
#'   `is_distal`.
#' * `apa_type` relative to the canonical isoform: `"same_end"`,
#'   `"tandem_utr"` (different 3' end, same stop codon: 3' UTR shortening or
#'   lengthening), `"alternative_last_exon"` (different last exon),
#'   `"intronic"` (3' end inside an intron of the canonical isoform, typically
#'   a composite terminal exon) or `"other"`.
#' * `utr3_length` (from [transcript_orfs()] when sequences are available,
#'   else from the annotated CDS).
#' * `polya_signal` and `polya_signal_distance`: the poly(A) signal hexamer
#'   (`AATAAA`, `ATTAAA` or one of the weaker variants) found in the last
#'   `signal_window` nt of the transcript sequence, and its distance to the
#'   3' end (canonical signals sit 10-30 nt upstream); `NA` without
#'   sequences.
#' * With `polya_atlas`, the nearest atlas site within `window` bp of the
#'   3' end: `pas_site_id`, `pas_class`, `pas_tpm`, `pas_distance`,
#'   `pas_supported`. The atlas is the frozen PolyASite 2.0 table
#'   (`polya_atlas = "polyasite"`, see [frozen_data]), a data frame with
#'   `seqid`, `cleavage_pos`, `strand` (and optionally `site_id`,
#'   `site_class`, `tpm`), or a GENCODE `polyAs.gff3.gz` file
#'   (`download_reference(file_type = "polya")`).
#'
#' `apa_usage()` quantifies per gene, for genes with at least two 3' end
#' groups: the share of expression ending at the distal-most 3' end
#' (`distal_fraction`) and the abundance-weighted 3' UTR length
#' (`weighted_utr3_length`), per sample or per group with the change
#' between the last and the first group level. `test_apa()` tests the change
#' in distal usage between two groups with the quasi-binomial test used by
#' [test_productive_fraction()], on counts.
#'
#' @param self An [IsoformicExperiment].
#' @param gene_name Optional genes (names or IDs) to restrict the result to.
#' @param polya_atlas `NULL` (no atlas), `"polyasite"`, a data frame or a
#'   file path (Parquet table or GENCODE poly(A) GFF3).
#' @param window Distance in bp within which an atlas site supports a 3' end
#'   (default 50).
#' @param signal_window Length of the 3' terminal sequence searched for a
#'   poly(A) signal (default 50 nt).
#' @param group Name of a `col_data` column defining groups; `NULL` for
#'   per-sample values.
#' @param levels Levels of `group` to use, in this order; the change is
#'   `last - first`. Samples of other groups are ignored, as if the object
#'   had been subset to these samples, so every returned value comes from
#'   the samples kept. `NULL` takes all the levels of `group`. The tests
#'   (`test_productive_fraction()`, `test_apa()`) need exactly two levels.
#' @param assay Transcript assay to use (`"abundance"` for usage,
#'   `"counts"` for the test).
#' @param min_expr Minimum assay value for an isoform to count (default 0).
#' @param as_matrix Logical, return the genes x samples matrix of distal
#'   fractions instead of the summary table.
#' @param ... Passed on to methods.
#'
#' @returns Tibbles with one row per transcript (`polya_sites()`) or per gene
#'   (`apa_usage()`, `test_apa()`); a matrix for `as_matrix = TRUE`.
#'
#' @examples
#' iso <- example_isoformic()
#' pa <- polya_sites(iso, gene_name = "RBPJ")
#' pa[, c("transcript_name", "tes_rank", "apa_type", "utr3_length", "polya_signal")]
#' head(apa_usage(iso, group = "condition"))
#' head(test_apa(iso, group = "condition"))
#'
#' @name polyadenylation
#' @rdname polyadenylation
NULL

polya_signal_hexamers <- function() {
  c(
    "AATAAA",
    "ATTAAA",
    "AGTAAA",
    "TATAAA",
    "CATAAA",
    "GATAAA",
    "AATATA",
    "AATACA",
    "AATAGA",
    "AAAAAG",
    "ACTAAA",
    "AAGAAA"
  )
}

#' @rdname polyadenylation
#' @export
polya_sites <- S7::new_generic(
  "polya_sites",
  "self",
  function(
    self,
    gene_name = NULL,
    polya_atlas = NULL,
    window = 50L,
    signal_window = 50L,
    ...
  ) {
    S7::S7_dispatch()
  }
)

S7::method(polya_sites, IsoformicExperiment) <- function(
  self,
  gene_name = NULL,
  polya_atlas = NULL,
  window = 50L,
  signal_window = 50L,
  ...
) {
  feats <- transcript_features(self)
  orfs <- if (has_sequences(self)) transcript_orfs(self) else NULL
  seqs <- if (has_sequences(self)) {
    transcript_sequences(self, transcript_id = feats$transcript_id)
  } else {
    NULL
  }
  res <- compute_polya_sites(
    feats = feats,
    exons = self@annot_data_exons,
    orfs = orfs,
    seqs = seqs,
    signal_window = signal_window
  )
  atlas <- resolve_polya_atlas(polya_atlas)
  if (!is.null(atlas)) {
    res <- join_polya_atlas(res, atlas, window = window)
  }
  filter_genes(res, gene_name)
}

#' @keywords internal
#' @noRd
compute_polya_sites <- function(feats, exons, orfs, seqs, signal_window = 50L) {
  n <- nrow(feats)
  minus <- feats$strand == "-"
  seqid <- exons$seqid[match(feats$transcript_id, exons$transcript_id)]
  # Rank of each 3' end group along the direction of transcription
  direction <- ifelse(minus, -feats$tes_pos, feats$tes_pos)
  grp <- tibble::tibble(
    gene_id = feats$gene_id,
    tes_group = feats$tes_group,
    direction = direction
  ) |>
    dplyr::summarise(
      direction = min(.data$direction),
      .by = c("gene_id", "tes_group")
    ) |>
    dplyr::arrange(.data$gene_id, .data$direction) |>
    dplyr::mutate(rank = dplyr::row_number(), .by = "gene_id")
  key <- paste(feats$gene_id, feats$tes_group)
  rank <- grp$rank[match(key, paste(grp$gene_id, grp$tes_group))]
  n_groups <- as.integer(stats::ave(rank, feats$gene_id, FUN = max))
  can_idx <- match(feats$canonical_transcript_id, feats$transcript_id)
  can_tes_group <- feats$tes_group[can_idx]
  if (!is.null(orfs)) {
    o <- orfs[match(feats$transcript_id, orfs$transcript_id), ]
    # 3' UTR lengths are only meaningful for ORFs anchored on the annotation
    utr3 <- ifelse(
      o$orf_source %in% c("annotation", "projected"),
      o$utr3_length,
      NA_real_
    )
    same_stop <- o$stop_same_as_canonical
  } else {
    utr3 <- feats$utr3_length
    same_stop <- !(feats$cds_differs_from_canonical %in% TRUE) & feats$has_cds
  }
  inside_canonical_intron <- tes_in_canonical_intron(feats, exons)
  apa_type <- rep("other", n)
  apa_type[feats$tes_group == can_tes_group] <- "same_end"
  diff_end <- feats$tes_group != can_tes_group
  apa_type[
    diff_end & same_stop %in% TRUE & !(feats$alt_last_exon %in% TRUE)
  ] <- "tandem_utr"
  apa_type[diff_end & feats$alt_last_exon %in% TRUE] <- "alternative_last_exon"
  apa_type[diff_end & inside_canonical_intron] <- "intronic"
  apa_type[is.na(can_tes_group)] <- NA_character_
  signal <- rep(NA_character_, n)
  signal_distance <- rep(NA_integer_, n)
  if (!is.null(seqs)) {
    s <- seqs$sequence[match(feats$transcript_id, seqs$transcript_id)]
    tail_seq <- substring(s, pmax(1L, nchar(s) - signal_window + 1L), nchar(s))
    found <- find_polya_signal(tail_seq)
    signal <- found$signal
    signal_distance <- found$distance
  }
  tibble::tibble(
    transcript_id = feats$transcript_id,
    transcript_name = feats$transcript_name,
    gene_id = feats$gene_id,
    gene_name = feats$gene_name,
    seqid = seqid,
    strand = feats$strand,
    tes_pos = feats$tes_pos,
    tes_group = feats$tes_group,
    tes_rank = rank,
    n_tes_groups = n_groups,
    is_proximal = rank == 1L & n_groups > 1L,
    is_distal = rank == n_groups & n_groups > 1L,
    canonical_transcript_id = feats$canonical_transcript_id,
    apa_type = apa_type,
    utr3_length = utr3,
    polya_signal = signal,
    polya_signal_distance = signal_distance
  )
}

#' Is the 3' end of each transcript inside an intron of its canonical isoform?
#' @keywords internal
#' @noRd
tes_in_canonical_intron <- function(feats, exons) {
  can <- unique(feats$canonical_transcript_id)
  ex <- exons[
    exons$transcript_id %in% can,
    c("transcript_id", "start_pos", "end_pos")
  ]
  if (!nrow(ex)) {
    return(rep(FALSE, nrow(feats)))
  }
  span <- ex |>
    dplyr::summarise(
      g_start = min(.data$start_pos),
      g_end = max(.data$end_pos),
      .by = "transcript_id"
    )
  pos <- data.frame(
    idx = seq_len(nrow(feats)),
    canonical_transcript_id = feats$canonical_transcript_id,
    tes_pos = feats$tes_pos,
    stringsAsFactors = FALSE
  )
  con <- duckdb_connect()
  duckdb::duckdb_register(con, "pos", pos)
  duckdb::duckdb_register(con, "ex", as.data.frame(ex))
  duckdb::duckdb_register(con, "span", as.data.frame(span))
  res <- DBI::dbGetQuery(
    con,
    "SELECT p.idx
     FROM pos p
     JOIN span s ON s.transcript_id = p.canonical_transcript_id
     WHERE p.tes_pos BETWEEN s.g_start AND s.g_end
       AND NOT EXISTS (
         SELECT 1 FROM ex e
         WHERE e.transcript_id = p.canonical_transcript_id
           AND p.tes_pos BETWEEN e.start_pos AND e.end_pos
       )"
  )
  out <- rep(FALSE, nrow(feats))
  out[res$idx] <- TRUE
  out
}

#' Strongest poly(A) signal hexamer in 3' terminal sequences
#'
#' Hexamers are searched in order of strength; the match closest to the 3'
#' end of the strongest hexamer present is reported with its distance (nt
#' between the signal start and the transcript end).
#'
#' @keywords internal
#' @noRd
find_polya_signal <- function(tail_seq) {
  n <- length(tail_seq)
  signal <- rep(NA_character_, n)
  distance <- rep(NA_integer_, n)
  todo <- !is.na(tail_seq) & nchar(tail_seq) >= 6L
  for (hex in polya_signal_hexamers()) {
    if (!any(todo)) {
      break
    }
    pos <- gregexpr(hex, tail_seq[todo], fixed = TRUE)
    last <- vapply(pos, function(m) max(as.integer(m)), integer(1L))
    hit <- last > 0L
    idx <- which(todo)[hit]
    signal[idx] <- hex
    distance[idx] <- nchar(tail_seq[idx]) - last[hit]
    todo[idx] <- FALSE
  }
  list(signal = signal, distance = distance)
}

#' Resolve the `polya_atlas` argument into a table (or NULL)
#' @keywords internal
#' @noRd
resolve_polya_atlas <- function(polya_atlas) {
  if (is.null(polya_atlas)) {
    return(NULL)
  }
  if (is.data.frame(polya_atlas)) {
    atlas <- tibble::as_tibble(polya_atlas)
  } else if (identical(polya_atlas, "polyasite")) {
    atlas <- read_frozen_data("polyasite")
  } else if (grepl("\\.gff3?(\\.gz)?$", polya_atlas, ignore.case = TRUE)) {
    atlas <- read_polya_gff(polya_atlas)
  } else {
    atlas <- read_parquet(polya_atlas)
  }
  needed <- c("seqid", "cleavage_pos", "strand")
  if (!all(needed %in% colnames(atlas))) {
    cli::cli_abort(
      message = c(
        x = "The poly(A) atlas needs columns {.field {needed}}."
      ),
      class = "isoformic_invalid_evidence"
    )
  }
  if (!"site_id" %in% colnames(atlas)) {
    atlas$site_id <- paste0(
      atlas$seqid,
      ":",
      atlas$cleavage_pos,
      ":",
      atlas$strand
    )
  }
  if (!"site_class" %in% colnames(atlas)) {
    atlas$site_class <- NA_character_
  }
  if (!"tpm" %in% colnames(atlas)) {
    atlas$tpm <- NA_real_
  }
  atlas
}

#' @rdname polyadenylation
#' @param path Path to a GENCODE `polyAs.gff3.gz` file.
#' @export
read_polya_gff <- function(path) {
  assert_file_exists(path)
  gff <- read_tabular_lines(path, n_cols = 9L)
  gff <- gff[gff[[3L]] %in% c("polyA_site", "polyA_signal", "pseudo_polyA"), ]
  strand <- gff[[7L]]
  start <- as.integer(gff[[4L]])
  end <- as.integer(gff[[5L]])
  tibble::tibble(
    site_id = sub("^ID=([^;]+).*$", "\\1", gff[[9L]]),
    seqid = gff[[1L]],
    start_pos = start,
    end_pos = end,
    cleavage_pos = ifelse(strand == "-", start, end),
    strand = strand,
    tpm = NA_real_,
    site_class = gff[[3L]],
    signals = NA_character_
  )
}

#' Join the nearest atlas site within `window` of each 3' end
#' @keywords internal
#' @noRd
join_polya_atlas <- function(res, atlas, window = 50L) {
  tes <- data.frame(
    idx = seq_len(nrow(res)),
    seqid = as.character(res$seqid),
    strand = res$strand,
    tes_pos = res$tes_pos,
    stringsAsFactors = FALSE
  )
  atl <- as.data.frame(atlas[, c(
    "site_id",
    "seqid",
    "cleavage_pos",
    "strand",
    "site_class",
    "tpm"
  )])
  atl$seqid <- as.character(atl$seqid)
  if (!any(atl$seqid %in% tes$seqid)) {
    atl$seqid <- if (any(grepl("^chr", tes$seqid))) {
      paste0("chr", sub("^chr", "", atl$seqid))
    } else {
      sub("^chr", "", atl$seqid)
    }
  }
  con <- duckdb_connect()
  duckdb::duckdb_register(con, "tes", tes)
  duckdb::duckdb_register(con, "atl", atl)
  hits <- DBI::dbGetQuery(
    con,
    glue::glue_sql(
      "SELECT idx, site_id, site_class, tpm, distance FROM (
         SELECT t.idx, a.site_id, a.site_class, a.tpm,
                abs(a.cleavage_pos - t.tes_pos) AS distance,
                row_number() OVER (PARTITION BY t.idx ORDER BY abs(a.cleavage_pos - t.tes_pos), a.tpm DESC) AS rn
         FROM tes t JOIN atl a
           ON a.seqid = t.seqid AND a.strand = t.strand
          AND abs(a.cleavage_pos - t.tes_pos) <= {as.integer(window)}
       ) WHERE rn = 1",
      .con = con
    )
  )
  res$pas_site_id <- NA_character_
  res$pas_class <- NA_character_
  res$pas_tpm <- NA_real_
  res$pas_distance <- NA_integer_
  res$pas_site_id[hits$idx] <- hits$site_id
  res$pas_class[hits$idx] <- hits$site_class
  res$pas_tpm[hits$idx] <- hits$tpm
  res$pas_distance[hits$idx] <- as.integer(hits$distance)
  res$pas_supported <- !is.na(res$pas_site_id)
  res
}

#' @rdname polyadenylation
#' @export
apa_usage <- S7::new_generic(
  "apa_usage",
  "self",
  function(
    self,
    group = NULL,
    levels = NULL,
    assay = "abundance",
    min_expr = 0,
    as_matrix = FALSE,
    gene_name = NULL,
    ...
  ) {
    S7::S7_dispatch()
  }
)

S7::method(apa_usage, IsoformicExperiment) <- function(
  self,
  group = NULL,
  levels = NULL,
  assay = "abundance",
  min_expr = 0,
  as_matrix = FALSE,
  gene_name = NULL,
  ...
) {
  assay <- resolve_assay_name(self, assay)
  pa <- polya_sites(self)
  tx_mat <- assay(self, assay)
  grp <- group_samples(self, group, colnames(tx_mat), levels = levels)
  groups <- grp$groups
  tx_mat <- tx_mat[, grp$samples, drop = FALSE]
  tx_mat[tx_mat < min_expr] <- 0
  idx <- match(rownames(tx_mat), pa$transcript_id)
  gene_ids <- pa$gene_id[idx]
  multi <- pa$n_tes_groups[idx] > 1L
  is_distal <- pa$is_distal[idx] %in% TRUE
  utr3 <- pa$utr3_length[idx]
  total <- rowsum(tx_mat, gene_ids)
  distal <- rowsum(tx_mat * is_distal, gene_ids)
  frac <- distal / total
  frac[!is.finite(frac)] <- NA_real_
  has_utr <- !is.na(utr3)
  w_utr <- rowsum(tx_mat * ifelse(has_utr, utr3, 0), gene_ids) /
    rowsum(tx_mat * has_utr, gene_ids)
  w_utr[!is.finite(w_utr)] <- NA_real_
  keep <- rownames(frac) %in% unique(gene_ids[multi])
  frac <- frac[keep, , drop = FALSE]
  w_utr <- w_utr[keep, , drop = FALSE]
  if (isTRUE(as_matrix)) {
    return(frac)
  }
  genes <- self@annot_data_genes[, c("gene_id", "gene_name")]
  types <- pa |>
    dplyr::filter(
      .data$gene_id %in% rownames(frac),
      !.data$apa_type %in% "same_end"
    ) |>
    dplyr::summarise(
      apa_types = paste(sort(unique(.data$apa_type)), collapse = ";"),
      .by = "gene_id"
    )
  res <- tibble::tibble(
    gene_id = rownames(frac) %||% character(),
    gene_name = genes$gene_name[match(rownames(frac), genes$gene_id)],
    n_tes_groups = pa$n_tes_groups[match(rownames(frac), pa$gene_id)],
    apa_types = types$apa_types[match(rownames(frac), types$gene_id)],
    mean_gene_expr = as.vector(rowMeans(total[rownames(frac), , drop = FALSE])),
    distal_fraction = as.vector(rowMeans(frac, na.rm = TRUE)),
    weighted_utr3_length = as.vector(rowMeans(w_utr, na.rm = TRUE))
  )
  res$distal_fraction[is.nan(res$distal_fraction)] <- NA_real_
  res$weighted_utr3_length[is.nan(res$weighted_utr3_length)] <- NA_real_
  if (!is.null(group)) {
    cols <- c(
      group_fraction_columns(frac, groups, "distal_fraction"),
      group_fraction_columns(w_utr, groups, "weighted_utr3_length")
    )
    for (nm in names(cols)) {
      res[[nm]] <- cols[[nm]]
    }
    if ("delta_distal_fraction" %in% names(cols)) {
      res <- dplyr::arrange(res, dplyr::desc(abs(.data$delta_distal_fraction)))
    }
  }
  filter_genes(res, gene_name)
}

#' @rdname polyadenylation
#' @export
test_apa <- S7::new_generic(
  "test_apa",
  "self",
  function(
    self,
    group,
    levels = NULL,
    assay = "counts",
    min_expr = 0,
    gene_name = NULL,
    ...
  ) {
    S7::S7_dispatch()
  }
)

S7::method(test_apa, IsoformicExperiment) <- function(
  self,
  group,
  levels = NULL,
  assay = "counts",
  min_expr = 0,
  gene_name = NULL,
  ...
) {
  rlang::check_required(group)
  assay <- resolve_assay_name(self, assay)
  pa <- polya_sites(self)
  tx_mat <- assay(self, assay)
  tx_mat[tx_mat < min_expr] <- 0
  idx <- match(rownames(tx_mat), pa$transcript_id)
  gene_ids <- pa$gene_id[idx]
  multi <- unique(gene_ids[pa$n_tes_groups[idx] > 1L])
  is_distal <- pa$is_distal[idx] %in% TRUE
  total <- rowsum(tx_mat, gene_ids)
  distal <- rowsum(tx_mat * is_distal, gene_ids)
  keep <- rownames(total) %in% multi
  total <- total[keep, , drop = FALSE]
  distal <- distal[keep, , drop = FALSE]
  setup <- two_group_setup(self, group, colnames(total), levels)
  total <- total[, setup$keep, drop = FALSE]
  distal <- distal[, setup$keep, drop = FALSE]
  test <- two_group_proportion_test(distal, total, setup$groups)
  genes <- self@annot_data_genes[, c("gene_id", "gene_name")]
  res <- dplyr::bind_cols(
    tibble::tibble(
      gene_id = rownames(total) %||% character(),
      gene_name = genes$gene_name[match(rownames(total), genes$gene_id)]
    ),
    test
  )
  names(res) <- sub("^proportion_", "distal_fraction_", names(res))
  names(res)[names(res) == "delta"] <- "delta_distal_fraction"
  res <- dplyr::arrange(res, .data$pvalue)
  filter_genes(res, gene_name)
}
