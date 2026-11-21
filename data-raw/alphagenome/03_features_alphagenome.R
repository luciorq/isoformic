# AlphaGenome pilot benchmark, step 3: AlphaGenome features
#
# Sourced by 03_evaluate.R when <work_dir>/scores/ of step 2 exists. Uses
# `work_dir`, `target_dir`, `tx`, `introns` (per transcript and intron, with
# `novel`) and `min_or()` of that script, and defines `ag` (one row per
# transcript), `scored` (gene IDs) and `diagnostics`.
#
# Only tracks with imputation_mode == "observed" are used: the AlphaGenome 2
# documentation says not to score splicing on imputed tracks.

score_dir <- fs::path(work_dir, "scores")
read_scores <- function(sub) {
  dplyr::bind_rows(lapply(
    fs::dir_ls(fs::path(score_dir, sub), glob = "*.parquet"),
    read_parquet
  ))
}
tracks <- read_parquet(fs::path(score_dir, "tracks.parquet")) |>
  dplyr::rename(track_strand = "strand")
genes <- read_parquet(fs::path(target_dir, "genes.parquet"))
site_targets <- read_parquet(fs::path(target_dir, "sites.parquet"))
ag_junctions <- read_scores("junctions")
scored <- unique(ag_junctions$gene_id)
track_key <- c("output", "condition", "track")

# Splice sites. `pos` is the first or last base of the intron; the signal is
# on the neighbouring exonic base: offset -1 at the left end of the intron and
# +1 at the right end (checked below on the data).
ag_sites <- read_scores("sites") |>
  dplyr::inner_join(
    dplyr::distinct(
      site_targets,
      .data$gene_id,
      .data$strand,
      .data$site,
      .data$pos
    ),
    by = c("gene_id", "site", "pos")
  ) |>
  dplyr::inner_join(
    tracks[, c(
      track_key,
      "track_strand",
      "splice_site_type",
      "imputation_mode"
    )],
    by = track_key
  ) |>
  dplyr::filter(
    .data$track_strand == .data$strand,
    .data$imputation_mode == "observed",
    .data$output == "splice_site_usage" | .data$splice_site_type == .data$site
  ) |>
  dplyr::mutate(
    right_end = (.data$site == "acceptor") == (.data$strand == "+"),
    expected = dplyr::if_else(.data$right_end, 1L, -1L)
  )
offset_check <- ag_sites |>
  dplyr::summarise(
    value = mean(.data$value),
    .by = c("output", "right_end", "offset", "expected")
  ) |>
  dplyr::slice_max(.data$value, n = 1L, by = c("output", "right_end"))
stopifnot(offset_check$offset == offset_check$expected)
site_scores <- ag_sites |>
  dplyr::filter(.data$offset == .data$expected) |>
  dplyr::summarise(
    value = mean(.data$value),
    .by = c("gene_id", "site", "pos", "output")
  ) |>
  tidyr::pivot_wider(names_from = "output", values_from = "value")

# Junctions: 0-based start, so intron_start = junction_start + 1
observed_junction <- tracks |>
  dplyr::filter(
    .data$output == "splice_junctions",
    .data$imputation_mode == "observed"
  )
n_junction_tracks <- nrow(observed_junction)
junction_scores <- ag_junctions |>
  dplyr::semi_join(observed_junction, by = c("condition", "track")) |>
  dplyr::mutate(
    relative = .data$value / max(.data$value),
    .by = c("gene_id", "condition", "track")
  ) |>
  dplyr::summarise(
    ag_junction = sum(.data$relative) / n_junction_tracks,
    .by = c("gene_id", "junction_start", "junction_end")
  ) |>
  dplyr::transmute(
    .data$gene_id,
    intron_start = .data$junction_start + 1L,
    intron_end = .data$junction_end,
    .data$ag_junction
  )

site_columns <- function(site, prefix) {
  out <- site_scores[site_scores$site == site, ]
  out <- out[, c("gene_id", "pos", "splice_sites", "splice_site_usage")]
  colnames(out) <- c("gene_id", paste0(prefix, c("_pos", "_prob", "_usage")))
  out
}
plus <- introns$strand == "+"
ag_introns <- introns |>
  dplyr::mutate(
    donor_pos = ifelse(plus, .data$intron_start, .data$intron_end),
    acceptor_pos = ifelse(plus, .data$intron_end, .data$intron_start)
  ) |>
  dplyr::filter(.data$gene_id %in% scored) |>
  dplyr::left_join(
    junction_scores,
    by = c("gene_id", "intron_start", "intron_end")
  ) |>
  dplyr::left_join(
    site_columns("donor", "donor"),
    by = c("gene_id", "donor_pos")
  ) |>
  dplyr::left_join(
    site_columns("acceptor", "acceptor"),
    by = c("gene_id", "acceptor_pos")
  ) |>
  dplyr::mutate(
    predicted = !is.na(.data$ag_junction),
    ag_junction = dplyr::coalesce(.data$ag_junction, 0),
    ag_site = pmin(.data$donor_prob, .data$acceptor_prob),
    ag_usage = pmin(.data$donor_usage, .data$acceptor_usage)
  )
unique_introns <- dplyr::distinct(
  ag_introns,
  .data$gene_id,
  .data$intron_start,
  .data$intron_end,
  .keep_all = TRUE
)
mid <- function(x) stats::median(x, na.rm = TRUE)
diagnostics <- tibble::tibble(
  genes_scored = length(scored),
  n_unique_introns = nrow(unique_introns),
  canonical_predicted = mean(unique_introns$predicted[!unique_introns$novel]),
  novel_predicted = mean(unique_introns$predicted[unique_introns$novel]),
  site_canonical = mid(unique_introns$ag_site[!unique_introns$novel]),
  site_novel = mid(unique_introns$ag_site[unique_introns$novel]),
  usage_canonical = mid(unique_introns$ag_usage[!unique_introns$novel]),
  usage_novel = mid(unique_introns$ag_usage[unique_introns$novel]),
  junction_tracks = n_junction_tracks
)

ag <- ag_introns |>
  dplyr::summarise(
    ag_weakest_junction = min_or(.data$ag_junction, TRUE, 0),
    ag_novel_junction = min_or(
      .data$ag_junction,
      .data$novel,
      mid(ag_introns$ag_junction)
    ),
    ag_weakest_site = min_or(.data$ag_site, TRUE, mid(ag_introns$ag_site)),
    ag_novel_site = min_or(
      .data$ag_site,
      .data$novel,
      mid(ag_introns$ag_site)
    ),
    ag_weakest_usage = min_or(.data$ag_usage, TRUE, mid(ag_introns$ag_usage)),
    ag_novel_usage = min_or(
      .data$ag_usage,
      .data$novel,
      mid(ag_introns$ag_usage)
    ),
    .by = "transcript_id"
  )

# Intron retention: predicted coverage of a canonical intron that lies inside
# an exon of the transcript, relative to the canonical exons of the gene
coverage <- read_scores("regions") |>
  dplyr::inner_join(genes[, c("gene_id", "strand")], by = "gene_id") |>
  dplyr::inner_join(
    tracks[, c(track_key, "track_strand", "imputation_mode")],
    by = track_key
  ) |>
  dplyr::filter(
    .data$track_strand == .data$strand,
    .data$imputation_mode == "observed"
  ) |>
  dplyr::summarise(
    coverage = mean(.data$value),
    .by = c("gene_id", "region", "start_pos", "end_pos")
  )
tx_exons <- read_parquet(fs::path(target_dir, "exons.parquet")) |>
  dplyr::filter(.data$gene_id %in% scored)
canonical_ids <- tx$transcript_id[tx$is_canonical %in% TRUE]
exon_level <- tx_exons |>
  dplyr::filter(.data$transcript_id %in% canonical_ids) |>
  dplyr::inner_join(
    dplyr::filter(coverage, .data$region == "exon"),
    by = c("gene_id", "start_pos", "end_pos")
  ) |>
  dplyr::summarise(exon_coverage = mean(.data$coverage), .by = "gene_id")
floor_coverage <- 1e-4
intron_coverage <- coverage[coverage$region == "intron", ]
intron_coverage <- tibble::tibble(
  gene_id = intron_coverage$gene_id,
  intron_start = intron_coverage$start_pos,
  intron_end = intron_coverage$end_pos,
  intron_coverage = intron_coverage$coverage
)
intron_level <- unique_introns |>
  dplyr::filter(!.data$novel) |>
  dplyr::select("gene_id", "intron_start", "intron_end") |>
  dplyr::inner_join(
    intron_coverage,
    by = c("gene_id", "intron_start", "intron_end")
  ) |>
  dplyr::inner_join(exon_level, by = "gene_id") |>
  dplyr::mutate(
    retention = log2(
      (.data$intron_coverage + floor_coverage) /
        (.data$exon_coverage + floor_coverage)
    )
  )
retained <- tx_exons |>
  dplyr::inner_join(
    intron_level,
    by = dplyr::join_by(
      gene_id,
      start_pos <= intron_start,
      end_pos >= intron_end
    )
  ) |>
  dplyr::summarise(ag_retention = max(.data$retention), .by = "transcript_id")
ag <- ag |>
  dplyr::full_join(retained, by = "transcript_id") |>
  dplyr::mutate(
    has_contained = as.integer(!is.na(.data$ag_retention)),
    ag_retention = dplyr::coalesce(
      .data$ag_retention,
      mid(intron_level$retention)
    )
  )
