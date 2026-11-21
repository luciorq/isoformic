# Benchmark of the sequence rules of transcript_productivity() on a full
# GENCODE annotation
#
# The table "Benchmark of the sequence rules" in
# .github/devdocs/ARCHITECTURE.md (backlog item A7 of
# .github/devdocs/ALPHAGENOME-BENCHMARK.md) comes from this script:
# benchmark_productivity() (use_biotype = FALSE, use_evidence = FALSE) over
# every transcript of the annotation that has a sequence from the primary
# assembly genome, and the break-down of the calls that the text below the
# table quotes. Rerun it when the rules, the ORF code or the coding potential
# model (inst/extdata/coding_model/) change.
#
# History: the first v50 table (2026-11-11, file benchmark_v50.rds of the
# AlphaGenome work directory) came from an inline `Rscript -e` call with the
# same steps (benchmark_productivity() and the protein_coding reasons); this
# script replaces it.
#
# Usage (from the package root). The whole annotation is in memory and the
# ORFs of every transcript are computed again when the coding model or the
# package version changed (one R process); on a shared machine run it under
# a memory cap and the heavy-job lock, for example:
#   ISOFORMIC_BENCHMARK_OUT=/path/to/benchmark.rds \
#     <memory-cap wrapper> \
#     pixi run --frozen Rscript data-raw/benchmark_rules.R
# Optional environment variables:
#   ISOFORMIC_BENCHMARK_GENCODE  GENCODE release (default: the package default)
#
# Downloads the GENCODE GFF3 (chr_patch_hapl_scaff) and the primary assembly
# genome when they are not in the isoformic cache, and extracts the transcript
# sequences into the annotation directory when it has none. The ORFs are
# stored in the annotation directory (annot_data_orfs.parquet, a cache keyed
# by the coding model and the package version).
#
# Writes an RDS file with:
#   bm          list(summary, confusion), as benchmark_productivity() returns
#   pc_reasons  reasons of protein_coding transcripts not called productive
#               or altered productive
#   ri_calls    calls of the retained_intron transcripts
#   lof_calls   calls of the protein_coding_LoF transcripts
#   n_with_sequence, n_transcripts, n_canonical, release, elapsed, git_head,
#   git_changes (uncommitted files under R/ and inst/ when the code was
#   loaded), coding_model (coefficients in use)

devtools::load_all(".", quiet = TRUE)
# The code that load_all() read: commit and uncommitted changes of R/, inst/
git_head <- tryCatch(
  system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
  error = function(cnd) NA_character_
)
git_changes <- tryCatch(
  system2("git", c("status", "--porcelain", "--", "R", "inst"), stdout = TRUE),
  error = function(cnd) NA_character_
)

out_file <- Sys.getenv("ISOFORMIC_BENCHMARK_OUT", unset = NA)
if (is.na(out_file) || !nzchar(out_file)) {
  stop("Set ISOFORMIC_BENCHMARK_OUT to the output .rds file.", call. = FALSE)
}
release <- Sys.getenv("ISOFORMIC_BENCHMARK_GENCODE", unset = "")
release <- if (nzchar(release)) release else NULL

gff <- download_reference(version = release, file_type = "gff")
annot_dir <- prepare_isoformic_annotation(gff)
if (!has_sequences(annot_dir)) {
  genome <- download_reference(
    version = release,
    file_type = "genome_fasta",
    region = "primary"
  )
  extract_transcript_sequences(annot_dir, genome_fasta = genome)
}
ids <- transcript_sequences(annot_dir)$transcript_id
n_transcripts <- nrow(read_parquet(fs::path(
  annot_dir,
  "annot_data_transcripts.parquet"
)))
cli::cli_inform(c(
  i = "Annotation: {.path {annot_dir}}",
  i = "{length(ids)} of {n_transcripts} transcripts have a sequence"
))
iso <- IsoformicExperiment(
  experiment_name = "rules_benchmark",
  annot_path = annot_dir,
  assays = list(
    abundance = matrix(
      1,
      nrow = length(ids),
      ncol = 1L,
      dimnames = list(ids, "s1")
    )
  ),
  col_data = data.frame(sample_id = "s1")
)

# The arguments of benchmark_productivity(), so that one classification
# gives both the table and the break-down
t0 <- Sys.time()
rules <- transcript_productivity(
  iso,
  nmd_rule = 50,
  use_biotype = FALSE,
  use_sequence = TRUE,
  use_evidence = FALSE
)
elapsed <- Sys.time() - t0
bm <- benchmark_table(rules)

productivity <- as.character(rules$productivity)
pc <- rules$transcript_type %in%
  "protein_coding" &
  !productivity %in% c("productive", "altered_productive")
pc_reasons <- table(rules$reason[pc])
ri_calls <- table(productivity[rules$transcript_type %in% "retained_intron"])
lof_calls <- table(productivity[
  rules$transcript_type %in% "protein_coding_LoF"
])
res <- list(
  bm = bm,
  pc_reasons = pc_reasons,
  ri_calls = ri_calls,
  lof_calls = lof_calls,
  n_with_sequence = length(ids),
  n_transcripts = n_transcripts,
  n_canonical = sum(rules$is_canonical %in% TRUE),
  release = release %||% "package default",
  elapsed = elapsed,
  git_head = git_head,
  git_changes = git_changes,
  coding_model = coding_model()
)
saveRDS(res, out_file)

cli::cli_h2("Benchmark of the sequence rules")
cli::cli_inform(c(i = "Classification took {format(elapsed)}"))
print(as.data.frame(bm$summary), digits = 4)
keep <- c(
  "protein_coding",
  "nonsense_mediated_decay",
  "retained_intron",
  "protein_coding_CDS_not_defined",
  "protein_coding_LoF",
  "non_stop_decay"
)
confusion <- as.data.frame(bm$confusion)
print(tidyr::pivot_wider(
  confusion[confusion$transcript_type %in% keep, ],
  names_from = "productivity",
  values_from = "n",
  values_fill = 0L
))
print(sort(pc_reasons, decreasing = TRUE))
print(ri_calls)
print(lof_calls)
