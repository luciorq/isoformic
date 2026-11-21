# Splice-site strength by transcript group on a full GENCODE annotation
#
# The numbers of backlog item A4 in
# .github/devdocs/ALPHAGENOME-BENCHMARK.md: median weakest donor and acceptor
# per transcript biotype, per-intron medians of introns of canonical
# transcripts against introns seen only in NMD transcripts, and the AUCs
# between them. The splice sites are extracted again with the current code
# into a copy of the annotation directory (the cached one is not touched) and
# compared with the table of the cached directory.
#
# Usage (from the package root; streams the genome, so on a shared machine
# run it under a memory cap and a lock for heavy jobs):
#   ISOFORMIC_SPLICE_WORK_DIR=/path/to/workdir \
#     pixi run --frozen Rscript data-raw/alphagenome/splice_site_groups.R
# Optional environment variables:
#   ISOFORMIC_SPLICE_GENCODE  GENCODE release (default: the package default)
#   ISOFORMIC_SPLICE_EXTRACT  "false" skips the fresh extraction and computes
#                             the tables from the splice sites of the cached
#                             annotation directory (light: no genome is read)
#
# Writes <work_dir>/annotation_copy/ (a copy of the annotation tables with a
# fresh annot_data_splice_sites.parquet; not with ISOFORMIC_SPLICE_EXTRACT =
# false) and <work_dir>/splice_site_groups.rds.

devtools::load_all(".", quiet = TRUE)

work_dir <- Sys.getenv("ISOFORMIC_SPLICE_WORK_DIR", unset = NA)
if (is.na(work_dir) || !nzchar(work_dir)) {
  stop("Set ISOFORMIC_SPLICE_WORK_DIR to a work directory.", call. = FALSE)
}
release <- Sys.getenv("ISOFORMIC_SPLICE_GENCODE", unset = "")
release <- if (nzchar(release)) release else NULL

extract <- !identical(
  tolower(Sys.getenv("ISOFORMIC_SPLICE_EXTRACT", unset = "true")),
  "false"
)
gff <- download_reference(version = release, file_type = "gff")
annot_dir <- prepare_isoformic_annotation(gff)

# 1. Fresh extraction into a copy of the annotation tables -------------------
if (extract) {
  genome <- download_reference(
    version = release,
    file_type = "genome_fasta",
    region = "primary"
  )
  copy_dir <- fs::path(work_dir, "annotation_copy")
  if (fs::dir_exists(copy_dir)) {
    fs::dir_delete(copy_dir)
  }
  fs::dir_create(copy_dir)
  core <- c(
    "annot_data_exons.parquet",
    "annot_data_genes.parquet",
    "annot_data_transcripts.parquet",
    "annot_metadata.parquet",
    "annot_seq_regions.parquet"
  )
  fs::file_copy(fs::path(annot_dir, core), fs::path(copy_dir, core))
  t0 <- Sys.time()
  extract_splice_sites(copy_dir, genome_fasta = genome)
  extract_seconds <- as.numeric(Sys.time() - t0, units = "secs")

  fresh <- as.data.frame(read_parquet(fs::path(
    copy_dir,
    "annot_data_splice_sites.parquet"
  )))
  cached_file <- fs::path(annot_dir, "annot_data_splice_sites.parquet")
  same_as_cached <- NA
  if (fs::file_exists(cached_file)) {
    cached <- as.data.frame(read_parquet(cached_file))
    ord <- function(x) x[do.call(order, unname(as.list(x))), , drop = FALSE]
    same_as_cached <- identical(colnames(fresh), colnames(cached)) &&
      isTRUE(all.equal(ord(fresh), ord(cached), check.attributes = FALSE))
  }
} else {
  stopifnot(has_splice_sites(annot_dir))
  copy_dir <- annot_dir
  fresh <- as.data.frame(read_parquet(fs::path(
    annot_dir,
    "annot_data_splice_sites.parquet"
  )))
  same_as_cached <- NA
  extract_seconds <- NA_real_
}

# 2. Per transcript group ------------------------------------------------------
tx <- read_parquet(fs::path(copy_dir, "annot_data_transcripts.parquet"))
by_t <- splice_site_strength(copy_dir)
by_t$type <- tx$transcript_type[match(by_t$transcript_id, tx$transcript_id)]
groups <- c(
  "protein_coding",
  "nonsense_mediated_decay",
  "retained_intron",
  "lncRNA",
  "protein_coding_CDS_not_defined"
)
group_tbl <- by_t |>
  dplyr::filter(.data$type %in% groups, .data$n_scored > 0) |>
  dplyr::summarise(
    n = dplyr::n(),
    median_introns = stats::median(.data$n_introns),
    weakest_donor = stats::median(.data$weakest_donor),
    weakest_acceptor = stats::median(.data$weakest_acceptor),
    .by = "type"
  )

# 3. Per intron, main chromosomes ---------------------------------------------
by_i <- splice_site_strength(copy_dir, by = "intron")
n_training <- attr(by_i, "n_training")
by_i$type <- tx$transcript_type[match(by_i$transcript_id, tx$transcript_id)]
by_i$canon <- tx$is_canonical[match(by_i$transcript_id, tx$transcript_id)] %in%
  TRUE
introns <- by_i |>
  dplyr::filter(grepl("^chr[0-9XY]+$", .data$seqid)) |>
  dplyr::summarise(
    in_canon = any(.data$canon),
    in_pc = any(.data$type == "protein_coding"),
    in_nmd = any(.data$type == "nonsense_mediated_decay"),
    in_lnc = any(.data$type == "lncRNA"),
    donor = dplyr::first(.data$donor_score),
    acceptor = dplyr::first(.data$acceptor_score),
    .by = c("seqid", "strand", "intron_start", "intron_end")
  ) |>
  dplyr::mutate(
    class = dplyr::case_when(
      .data$in_canon ~ "canonical",
      .data$in_pc ~ "other protein_coding",
      .data$in_nmd ~ "NMD only",
      .data$in_lnc ~ "lncRNA",
      .default = "other"
    )
  )
intron_tbl <- introns |>
  dplyr::filter(!is.na(.data$donor)) |>
  dplyr::summarise(
    n = dplyr::n(),
    donor = stats::median(.data$donor),
    acceptor = stats::median(.data$acceptor),
    .by = "class"
  )
auc <- function(score, lab) {
  r <- rank(score)
  n1 <- as.numeric(sum(lab))
  n0 <- as.numeric(sum(!lab))
  (sum(r[lab]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
x <- introns |>
  dplyr::filter(
    !is.na(.data$donor),
    .data$class %in% c("canonical", "NMD only")
  )
y <- by_t |>
  dplyr::filter(
    .data$type %in% c("protein_coding", "nonsense_mediated_decay"),
    .data$n_scored > 0
  )
auc_tbl <- tibble::tibble(
  comparison = c(
    "canonical vs NMD-only introns, donor",
    "canonical vs NMD-only introns, acceptor",
    "canonical vs NMD-only introns, smaller of the two",
    "protein_coding vs NMD transcripts, weakest donor",
    "protein_coding vs NMD transcripts, weakest acceptor"
  ),
  auc = c(
    auc(x$donor, x$class == "canonical"),
    auc(x$acceptor, x$class == "canonical"),
    auc(pmin(x$donor, x$acceptor), x$class == "canonical"),
    auc(y$weakest_donor, y$type == "protein_coding"),
    auc(y$weakest_acceptor, y$type == "protein_coding")
  )
)

res <- list(
  n_unique_introns = nrow(fresh),
  n_training = n_training,
  same_as_cached = same_as_cached,
  groups = group_tbl,
  introns = intron_tbl,
  auc = auc_tbl,
  extract_seconds = extract_seconds
)
saveRDS(res, fs::path(work_dir, "splice_site_groups.rds"))

cli::cli_h2("Splice-site strength by group")
cli::cli_inform(c(
  i = "Unique introns: {nrow(fresh)}; training introns: {n_training}",
  i = "Same table as the cached annotation directory: {same_as_cached}",
  i = "Extraction took {round(extract_seconds)} s"
))
print(as.data.frame(group_tbl), digits = 3)
print(as.data.frame(intron_tbl), digits = 3)
print(as.data.frame(auc_tbl), digits = 3)
