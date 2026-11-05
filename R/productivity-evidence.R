# External evidence for productivity -----------------------------------------
#
# Joins of the frozen datasets (see frozen-data.R) into the annotation
# directory, through add_annotation_columns(), plus the comparisons that need
# the object (domain retention, GTEx reference). transcript_productivity()
# reads the columns these functions add.

#' Attach External Productivity Evidence to an Annotation
#'
#' @description
#' These functions add evidence columns to the transcript table of the
#' annotation directory of `x` from the frozen datasets ([frozen_data]), so
#' that [transcript_productivity()] can use observed translation and decay
#' instead of rules alone. Identifiers are matched without version.
#'
#' * `add_appris()` adds `appris_annotation` (`PRINCIPAL:1` ...
#'   `ALTERNATIVE:2`, `MINOR`), `appris_principal`, `appris_score` and the
#'   method scores (`firestar_residues`, `matador_structure`,
#'   `corsair_conservation`, `spade_domain_score`, `thump_tm_helices`,
#'   `trifid_score`). Principal isoforms are productive evidence.
#' * `add_riboseq_orfs()` adds `riboseq_n_orfs` (Ribo-seq ORFs compatible
#'   with the transcript's exon structure), `riboseq_orf_types`,
#'   `riboseq_main_orf` (`TRUE` when a translated ORF coincides with the
#'   transcript's own ORF, from [transcript_orfs()] when sequences are
#'   available or else from the annotated CDS) and `riboseq_uorf`. A
#'   translated main ORF is productive evidence. Phase 1 ORFs list their
#'   compatible transcripts; all ORFs are also matched by exon-block
#'   compatibility (every block inside an exon, consecutive blocks in
#'   consecutive exons, internal block boundaries on exon boundaries).
#' * `add_nmd_exons()` adds `nmd_validated_exon`, `TRUE` when the transcript
#'   contains an exon with the exact coordinates of an NMD-inducing cassette
#'   exon validated by NMD inhibition (Karousis et al. 2021), or is one of
#'   the listed NMD isoforms. This is unproductive evidence.
#'
#' Each function returns the annotation directory path invisibly and
#' invalidates the object caches of `x` (the next call to a method rereads
#' the annotation).
#'
#' @param x An annotation directory path or an [IsoformicExperiment].
#' @param table A data frame with the dataset, or a path to its Parquet
#'   file; `NULL` (default) fetches the frozen dataset with
#'   [fetch_frozen_data()].
#'
#' @returns The annotation directory path, invisibly.
#'
#' @examples
#' \dontrun{
#' iso <- example_isoformic()
#' add_appris(iso)
#' add_riboseq_orfs(iso)
#' add_nmd_exons(iso)
#' prod <- transcript_productivity(iso, gene_name = "RBPJ")
#' prod[, c("transcript_name", "productivity", "reason", "source")]
#' }
#'
#' @name productivity_evidence
#' @rdname productivity_evidence
NULL

evidence_table <- function(table, name) {
  if (is.null(table)) {
    return(read_frozen_data(name))
  }
  if (is.data.frame(table)) {
    return(tibble::as_tibble(table))
  }
  read_parquet(table)
}

#' @rdname productivity_evidence
#' @export
add_appris <- function(x, table = NULL) {
  annot_dir <- annotation_dir_of(x)
  appris <- evidence_table(table, "appris")
  tx <- read_annotation_table(
    annot_dir,
    "transcripts",
    columns = "transcript_id"
  )
  idx <- match(strip_id_version(tx$transcript_id), appris$transcript_id)
  cols <- c(
    "appris_annotation",
    "appris_principal",
    "appris_score",
    "firestar_residues",
    "matador_structure",
    "corsair_conservation",
    "spade_domain_score",
    "thump_tm_helices",
    "trifid_score"
  )
  cols <- intersect(cols, colnames(appris))
  res <- tibble::tibble(transcript_id = tx$transcript_id[!is.na(idx)])
  for (col in cols) {
    res[[col]] <- appris[[col]][idx[!is.na(idx)]]
  }
  add_evidence(x, annot_dir, res, "APPRIS")
}

#' @rdname productivity_evidence
#' @export
add_riboseq_orfs <- function(x, table = NULL) {
  annot_dir <- annotation_dir_of(x)
  orfs <- evidence_table(table, "riboseq_orfs")
  tx <- read_annotation_table(
    annot_dir,
    "transcripts",
    columns = c("transcript_id", "seqid", "strand", "cds_start", "cds_end")
  )
  exons <- read_annotation_table(
    annot_dir,
    "exons",
    columns = c("transcript_id", "seqid", "start_pos", "end_pos", "strand")
  )
  # Phase 1: compatible transcripts are listed
  listed <- orfs[!is.na(orfs$transcript_ids) & nzchar(orfs$transcript_ids), ]
  listed_ids <- strsplit(listed$transcript_ids, ";", fixed = TRUE)
  pairs_listed <- tibble::tibble(
    orf_id = rep(listed$orf_id, lengths(listed_ids)),
    transcript_key = as.character(unlist(listed_ids)) %||% character()
  )
  tx_key <- strip_id_version(tx$transcript_id)
  pairs_listed$transcript_id <- tx$transcript_id[match(
    pairs_listed$transcript_key,
    tx_key
  )]
  pairs_listed <- pairs_listed[
    !is.na(pairs_listed$transcript_id),
    c("orf_id", "transcript_id")
  ]
  # All ORFs: exon-block compatibility
  pairs_blocks <- compatible_orf_transcripts(orfs, exons)
  pairs <- dplyr::distinct(dplyr::bind_rows(pairs_listed, pairs_blocks))
  if (!nrow(pairs)) {
    cli::cli_abort(
      message = c(
        x = "No Ribo-seq ORF is compatible with a transcript of the annotation."
      ),
      class = "isoformic_invalid_evidence"
    )
  }
  pairs <- dplyr::left_join(
    pairs,
    orfs[, c("orf_id", "orf_type", "start_pos", "end_pos", "strand")],
    by = "orf_id"
  )
  # Main ORF of each transcript: sequence-derived when available, else the CDS
  main <- main_orf_span(annot_dir, tx)
  pairs <- dplyr::left_join(pairs, main, by = "transcript_id")
  pairs$is_main <- !is.na(pairs$main_start) &
    pairs$start_pos == pairs$main_start &
    pairs$end_pos == pairs$main_end
  res <- pairs |>
    dplyr::summarise(
      riboseq_n_orfs = dplyr::n(),
      riboseq_orf_types = paste(
        sort(unique(ifelse(
          is.na(.data$orf_type),
          "unclassified",
          .data$orf_type
        ))),
        collapse = ";"
      ),
      riboseq_main_orf = any(.data$is_main),
      riboseq_uorf = any(.data$orf_type %in% c("uORF", "uoORF")),
      .by = "transcript_id"
    )
  add_evidence(x, annot_dir, res, "Ribo-seq ORF")
}

#' Genomic span of each transcript's main ORF (sequence-derived or annotated)
#' @keywords internal
#' @noRd
main_orf_span <- function(annot_dir, tx) {
  if (has_sequences(annot_dir)) {
    orfs <- transcript_orfs(annot_dir)
    orfs <- orfs[orfs$orf_source %in% c("annotation", "projected"), ]
    return(tibble::tibble(
      transcript_id = orfs$transcript_id,
      main_start = pmin(orfs$orf_genomic_start, orfs$orf_genomic_end),
      main_end = pmax(orfs$orf_genomic_start, orfs$orf_genomic_end)
    ))
  }
  tibble::tibble(
    transcript_id = tx$transcript_id,
    main_start = as.numeric(tx$cds_start),
    main_end = as.numeric(tx$cds_end)
  )
}

#' Pairs (orf_id, transcript_id) whose exon blocks fit the transcript exons
#'
#' Every block must lie inside an exon of the transcript on the same strand
#' and chromosome; block boundaries that are internal to the ORF must
#' coincide with exon boundaries (so the ORF uses the transcript's splice
#' junctions). Single-block ORFs need only containment.
#'
#' @keywords internal
#' @noRd
compatible_orf_transcripts <- function(orfs, exons) {
  if (!nrow(orfs) || !nrow(exons)) {
    return(tibble::tibble(orf_id = character(), transcript_id = character()))
  }
  starts <- strsplit(orfs$block_starts, ";", fixed = TRUE)
  ends <- strsplit(orfs$block_ends, ";", fixed = TRUE)
  n_blocks <- lengths(starts)
  blocks <- data.frame(
    orf_id = rep(orfs$orf_id, n_blocks),
    seqid = rep(orfs$seqid, n_blocks),
    strand = rep(orfs$strand, n_blocks),
    block_no = sequence(n_blocks),
    n_blocks = rep(n_blocks, n_blocks),
    b_start = as.numeric(unlist(starts)),
    b_end = as.numeric(unlist(ends)),
    stringsAsFactors = FALSE
  )
  ex <- as.data.frame(exons)
  ex$seqid <- as.character(ex$seqid)
  # Match chromosome naming (chr prefix) to the annotation
  if (!any(blocks$seqid %in% ex$seqid)) {
    alt <- if (any(grepl("^chr", ex$seqid))) {
      paste0("chr", sub("^chr", "", blocks$seqid))
    } else {
      sub("^chr", "", blocks$seqid)
    }
    blocks$seqid <- alt
  }
  con <- duckdb_connect()
  duckdb::duckdb_register(con, "blk", blocks)
  duckdb::duckdb_register(con, "ex", ex)
  res <- DBI::dbGetQuery(
    con,
    "WITH ranked AS (
       SELECT *, row_number() OVER (PARTITION BY transcript_id ORDER BY start_pos) AS exon_rank
       FROM ex
     ),
     hits AS (
       SELECT b.orf_id, e.transcript_id, b.n_blocks, e.exon_rank,
              CASE WHEN (b.block_no > 1 AND b.b_start <> e.start_pos)
                     OR (b.block_no < b.n_blocks AND b.b_end <> e.end_pos)
                   THEN 0 ELSE 1 END AS ok
       FROM blk b JOIN ranked e
         ON e.seqid = b.seqid AND e.strand = b.strand
        AND b.b_start >= e.start_pos AND b.b_end <= e.end_pos
     )
     SELECT orf_id, transcript_id
     FROM hits
     GROUP BY orf_id, transcript_id, n_blocks
     HAVING count(*) = n_blocks AND min(ok) = 1
        AND max(exon_rank) - min(exon_rank) = n_blocks - 1"
  )
  tibble::as_tibble(res)
}

#' @rdname productivity_evidence
#' @export
add_nmd_exons <- function(x, table = NULL) {
  annot_dir <- annotation_dir_of(x)
  nmd <- evidence_table(table, "nmd_exons")
  exons <- read_annotation_table(
    annot_dir,
    "exons",
    columns = c("transcript_id", "seqid", "start_pos", "end_pos", "strand")
  )
  tx <- read_annotation_table(
    annot_dir,
    "transcripts",
    columns = "transcript_id"
  )
  nmd_seqid <- nmd$seqid
  if (!any(nmd_seqid %in% exons$seqid)) {
    nmd_seqid <- if (any(grepl("^chr", exons$seqid))) {
      paste0("chr", sub("^chr", "", nmd_seqid))
    } else {
      sub("^chr", "", nmd_seqid)
    }
  }
  key_nmd <- paste(nmd_seqid, nmd$start_pos, nmd$end_pos, nmd$strand)
  key_ex <- paste(exons$seqid, exons$start_pos, exons$end_pos, exons$strand)
  by_exon <- unique(exons$transcript_id[key_ex %in% key_nmd])
  listed <- unique(unlist(strsplit(nmd$nmd_transcript_ids, ";", fixed = TRUE)))
  by_id <- tx$transcript_id[strip_id_version(tx$transcript_id) %in% listed]
  hits <- union(by_exon, by_id)
  if (!length(hits)) {
    cli::cli_abort(
      message = c(
        x = "No transcript of the annotation contains a validated NMD exon."
      ),
      class = "isoformic_invalid_evidence"
    )
  }
  res <- tibble::tibble(transcript_id = hits, nmd_validated_exon = TRUE)
  add_evidence(x, annot_dir, res, "validated NMD exon")
}

#' Write evidence columns and drop the derived caches of the object
#' @keywords internal
#' @noRd
add_evidence <- function(x, annot_dir, res, what) {
  if (!nrow(res)) {
    cli::cli_abort(
      message = c(
        x = "No transcript of the annotation matches the {what} table."
      ),
      class = "isoformic_invalid_evidence"
    )
  }
  add_annotation_columns(annot_dir, res, level = "transcript")
  if (S7::S7_inherits(x, IsoformicExperiment)) {
    clear_cache(x)
  }
  invisible(annot_dir)
}

# Domain retention -------------------------------------------------------------

#' Protein Domain Retention of Isoforms
#'
#' @description
#' Compares the Pfam domains of every isoform with those of the canonical
#' isoform of its gene, using the APPRIS SPADE domain coordinates (frozen
#' dataset `appris_domains`, [frozen_data]). For isoforms present in APPRIS
#' the domain sets are compared directly. For other isoforms (novel
#' transcripts) the genomic span of each canonical domain is projected onto
#' the isoform: a domain is `retained` when the isoform covers the same
#' exonic bases over the domain span and its ORF (see [transcript_orfs()],
#' or the annotated CDS) extends over it without a frameshift, `partial`
#' when only part of the span is covered, and `lost` otherwise.
#'
#' An `altered_productive` isoform that keeps all domains of the canonical
#' protein is a milder event than one that loses a domain; the
#' `domain_retention_fraction` column grades this and is included in
#' [transcript_productivity()] output when the dataset is available.
#'
#' @param self An [IsoformicExperiment].
#' @param gene_name Optional genes (names or IDs) to restrict the result to.
#' @param domains Data frame or Parquet path of the domain table; `NULL`
#'   fetches the frozen `appris_domains` dataset.
#' @param ... Passed on to methods.
#'
#' @returns A tibble with one row per transcript: `n_domains_canonical`,
#'   `n_retained`, `n_partial`, `n_lost`, `domains_lost` (Pfam names) and
#'   `domain_retention_fraction` (`NA` when the canonical isoform has no
#'   domain).
#'
#' @examples
#' \dontrun{
#' iso <- example_isoformic()
#' domain_retention(iso, gene_name = "RBPJ")
#' }
#'
#' @export
domain_retention <- S7::new_generic(
  "domain_retention",
  "self",
  function(self, gene_name = NULL, domains = NULL, ...) {
    S7::S7_dispatch()
  }
)

S7::method(domain_retention, IsoformicExperiment) <- function(
  self,
  gene_name = NULL,
  domains = NULL,
  ...
) {
  dom <- evidence_table(domains, "appris_domains")
  feats <- transcript_features(self)
  res <- compute_domain_retention(
    feats = feats,
    exons = self@annot_data_exons,
    orfs = if (has_sequences(self)) transcript_orfs(self) else NULL,
    transcripts = row_data(self),
    domains = dom
  )
  filter_genes(res, gene_name)
}

#' @keywords internal
#' @noRd
compute_domain_retention <- function(feats, exons, orfs, transcripts, domains) {
  key <- strip_id_version(feats$transcript_id)
  canonical_key <- strip_id_version(feats$canonical_transcript_id)
  dom <- domains[domains$transcript_id %in% c(key, canonical_key), ]
  dom$label <- paste0(dom$pfam_name, "@", dom$start_pos, "-", dom$end_pos)
  dom_by_tx <- split(dom$label, dom$transcript_id)
  in_appris <- key %in% domains$transcript_id
  # ORF genomic span of each transcript (sequence-derived or annotated CDS)
  if (!is.null(orfs)) {
    o <- orfs[match(feats$transcript_id, orfs$transcript_id), ]
    orf_lo <- pmin(o$orf_genomic_start, o$orf_genomic_end)
    orf_hi <- pmax(o$orf_genomic_start, o$orf_genomic_end)
    frameshift <- o$frame_vs_canonical %in% "frameshift"
  } else {
    t <- transcripts[match(feats$transcript_id, transcripts$transcript_id), ]
    orf_lo <- as.numeric(t$cds_start)
    orf_hi <- as.numeric(t$cds_end)
    frameshift <- rep(FALSE, nrow(feats))
  }
  # Exonic coverage of each canonical domain span by each isoform (DuckDB)
  can_dom <- dom[dom$transcript_id %in% canonical_key, ]
  cov <- domain_coverage(feats, canonical_key, can_dom, exons)
  n <- nrow(feats)
  n_can <- integer(n)
  n_ret <- integer(n)
  n_part <- integer(n)
  n_lost <- integer(n)
  lost_names <- rep(NA_character_, n)
  for (i in seq_len(n)) {
    can_labels <- dom_by_tx[[canonical_key[[i]]]]
    if (is.null(can_labels)) {
      n_can[[i]] <- 0L
      next
    }
    can_labels <- unique(can_labels)
    n_can[[i]] <- length(can_labels)
    if (in_appris[[i]]) {
      own <- unique(dom_by_tx[[key[[i]]]])
      retained <- can_labels %in% own
      n_ret[[i]] <- sum(retained)
      n_lost[[i]] <- sum(!retained)
      lost <- can_labels[!retained]
    } else {
      c_i <- cov[cov$idx == i, ]
      status <- vapply(
        can_labels,
        function(lab) {
          r <- c_i[c_i$label == lab, ]
          if (!nrow(r) || r$covered == 0) {
            return("lost")
          }
          full <- r$covered >= r$canonical_covered
          within_orf <- !is.na(orf_lo[[i]]) &&
            orf_lo[[i]] <= r$start_pos &&
            orf_hi[[i]] >= r$end_pos
          if (full && within_orf && !frameshift[[i]]) "retained" else "partial"
        },
        character(1L)
      )
      n_ret[[i]] <- sum(status == "retained")
      n_part[[i]] <- sum(status == "partial")
      n_lost[[i]] <- sum(status == "lost")
      lost <- can_labels[status == "lost"]
    }
    if (length(lost)) {
      lost_names[[i]] <- paste(sub("@.*$", "", lost), collapse = ";")
    }
  }
  tibble::tibble(
    transcript_id = feats$transcript_id,
    transcript_name = feats$transcript_name,
    gene_id = feats$gene_id,
    gene_name = feats$gene_name,
    canonical_transcript_id = feats$canonical_transcript_id,
    in_domain_table = in_appris,
    n_domains_canonical = n_can,
    n_retained = n_ret,
    n_partial = n_part,
    n_lost = n_lost,
    domains_lost = lost_names,
    domain_retention_fraction = ifelse(n_can > 0L, n_ret / n_can, NA_real_)
  )
}

#' Exonic bases of each isoform and of its canonical within each domain span
#' @keywords internal
#' @noRd
domain_coverage <- function(feats, canonical_key, can_dom, exons) {
  empty <- tibble::tibble(
    idx = integer(),
    label = character(),
    start_pos = numeric(),
    end_pos = numeric(),
    covered = numeric(),
    canonical_covered = numeric()
  )
  if (!nrow(can_dom)) {
    return(empty)
  }
  pairs <- data.frame(
    idx = seq_len(nrow(feats)),
    transcript_id = feats$transcript_id,
    canonical_transcript_id = feats$canonical_transcript_id,
    canonical_key = canonical_key,
    stringsAsFactors = FALSE
  )
  dom <- as.data.frame(can_dom[, c(
    "transcript_id",
    "label",
    "seqid",
    "start_pos",
    "end_pos"
  )])
  names(dom)[1L] <- "canonical_key"
  ex <- as.data.frame(exons[, c(
    "transcript_id",
    "seqid",
    "start_pos",
    "end_pos"
  )])
  ex$seqid <- as.character(ex$seqid)
  if (!any(dom$seqid %in% ex$seqid)) {
    dom$seqid <- if (any(grepl("^chr", ex$seqid))) {
      paste0("chr", sub("^chr", "", dom$seqid))
    } else {
      sub("^chr", "", dom$seqid)
    }
  }
  con <- duckdb_connect()
  duckdb::duckdb_register(con, "pairs", pairs)
  duckdb::duckdb_register(con, "dom", dom)
  duckdb::duckdb_register(con, "ex", ex)
  res <- DBI::dbGetQuery(
    con,
    "WITH spans AS (
       SELECT p.idx, p.transcript_id, p.canonical_transcript_id, d.label, d.seqid,
              d.start_pos AS d_start, d.end_pos AS d_end
       FROM pairs p JOIN dom d ON d.canonical_key = p.canonical_key
     ),
     cov AS (
       SELECT s.idx, s.label, s.d_start, s.d_end,
              coalesce(sum(CASE WHEN e.transcript_id IS NULL THEN 0
                           ELSE least(e.end_pos, s.d_end) - greatest(e.start_pos, s.d_start) + 1 END), 0) AS covered
       FROM spans s LEFT JOIN ex e
         ON e.transcript_id = s.transcript_id AND e.seqid = s.seqid
        AND e.start_pos <= s.d_end AND e.end_pos >= s.d_start
       GROUP BY s.idx, s.label, s.d_start, s.d_end
     ),
     can_cov AS (
       SELECT s.idx, s.label,
              coalesce(sum(CASE WHEN e.transcript_id IS NULL THEN 0
                           ELSE least(e.end_pos, s.d_end) - greatest(e.start_pos, s.d_start) + 1 END),
                       0) AS canonical_covered
       FROM spans s LEFT JOIN ex e
         ON e.transcript_id = s.canonical_transcript_id AND e.seqid = s.seqid
        AND e.start_pos <= s.d_end AND e.end_pos >= s.d_start
       GROUP BY s.idx, s.label
     )
     SELECT c.idx, c.label, c.d_start AS start_pos, c.d_end AS end_pos,
            c.covered, k.canonical_covered
     FROM cov c JOIN can_cov k ON k.idx = c.idx AND k.label = c.label"
  )
  tibble::as_tibble(res)
}

# GTEx reference ---------------------------------------------------------------

#' Tissue Reference Values from GTEx
#'
#' @description
#' Computes, from the GTEx v8 median transcript TPM of one tissue (frozen
#' dataset `gtex_transcript_tpm`, [frozen_data]) and the productivity calls
#' of `x`, the reference state of every gene in healthy tissue: its
#' productive fraction, its dominant isoform and the dominant isoform's
#' productivity. [productive_fraction()] accepts the result (or a tissue
#' name) as `reference` and reports the difference between the samples and
#' the tissue baseline, which separates a disease-specific loss of
#' productive isoforms from a tissue's normal state.
#'
#' Transcript identifiers are matched without version; GTEx v8 uses GENCODE
#' v26, so transcripts annotated later have no reference value.
#'
#' @param self An [IsoformicExperiment].
#' @param tissue GTEx tissue name (`SMTSD`, for example `"Whole Blood"`,
#'   `"Brain - Cortex"`); see `gtex_tissues()`.
#' @param gtex Data frame or Parquet path of the GTEx table; `NULL` fetches
#'   the frozen dataset.
#' @param include_altered Logical, count `altered_productive` isoforms as
#'   productive (default `TRUE`).
#' @param ... Passed on to methods.
#'
#' @returns A tibble with one row per gene of `x`: `reference_tissue`,
#'   `reference_gene_tpm`, `reference_productive_fraction`,
#'   `reference_dominant_transcript_id`, `reference_dominant_productivity`,
#'   `n_reference_transcripts`.
#'
#' @examples
#' \dontrun{
#' iso <- example_isoformic()
#' ref <- gtex_reference(iso, tissue = "Whole Blood")
#' productive_fraction(iso, group = "condition", reference = ref)
#' }
#'
#' @name gtex_reference
#' @rdname gtex_reference
#' @export
gtex_reference <- S7::new_generic(
  "gtex_reference",
  "self",
  function(self, tissue, gtex = NULL, include_altered = TRUE, ...) {
    S7::S7_dispatch()
  }
)

S7::method(gtex_reference, IsoformicExperiment) <- function(
  self,
  tissue,
  gtex = NULL,
  include_altered = TRUE,
  ...
) {
  rlang::check_required(tissue)
  check_string(tissue)
  gt <- evidence_table(gtex, "gtex_transcript_tpm")
  if (!tissue %in% gt$tissue) {
    cli::cli_abort(
      message = c(
        x = "Tissue {.val {tissue}} not found in the GTEx table.",
        i = "Available tissues: {.val {sort(unique(gt$tissue))}}."
      ),
      class = "isoformic_unknown_tissue"
    )
  }
  gt <- gt[gt$tissue == tissue, ]
  prod <- transcript_productivity(self)
  key <- strip_id_version(prod$transcript_id)
  tpm <- gt$median_tpm[match(key, gt$transcript_id)]
  productive_classes <- c(
    "productive",
    if (isTRUE(include_altered)) "altered_productive"
  )
  d <- tibble::tibble(
    gene_id = prod$gene_id,
    gene_name = prod$gene_name,
    transcript_id = prod$transcript_id,
    productivity = as.character(prod$productivity),
    tpm = tpm
  )
  d <- d[!is.na(d$tpm), ]
  res <- d |>
    dplyr::arrange(dplyr::desc(.data$tpm)) |>
    dplyr::summarise(
      reference_gene_tpm = sum(.data$tpm),
      reference_productive_fraction = ifelse(
        sum(.data$tpm) > 0,
        sum(.data$tpm[.data$productivity %in% productive_classes]) /
          sum(.data$tpm),
        NA_real_
      ),
      reference_dominant_transcript_id = dplyr::first(.data$transcript_id),
      reference_dominant_productivity = dplyr::first(.data$productivity),
      n_reference_transcripts = dplyr::n(),
      .by = c("gene_id", "gene_name")
    )
  res$reference_productive_fraction[
    res$gene_id %in% prod$gene_id[prod$productivity %in% "noncoding"]
  ] <- NA_real_
  res$reference_tissue <- tissue
  dplyr::relocate(res, "reference_tissue", .after = "gene_name")
}

#' @rdname gtex_reference
#' @export
gtex_tissues <- function(gtex = NULL) {
  gt <- evidence_table(gtex, "gtex_transcript_tpm")
  tbl <- dplyr::distinct(gt[, c("tissue", "n_samples")])
  dplyr::arrange(tbl, .data$tissue)
}
