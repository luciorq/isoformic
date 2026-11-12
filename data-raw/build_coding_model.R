# Train the coding potential model shipped in inst/extdata/coding_model/.
#
# Training data: GENCODE human transcripts (the release cached by
# download_reference(); update `release` to retrain on a newer one).
#   * hexamer table: CDS sequences (start to stop codon) of `protein_coding`
#                transcripts with a complete CDS, taken from the GENCODE
#                pc_transcripts FASTA header coordinates (`CDS:start-end`),
#                against full `lncRNA` transcript sequences.
#   * logistic model: full mRNA sequences of those coding transcripts versus
#                the lncRNA transcripts, scored on their longest ORF.
# A balanced random sample (`n_per_class`, seed 1) keeps the run short and the
# logistic regression well conditioned; one third is held out to report the
# test accuracy. Run from the package root:
#   Rscript data-raw/build_coding_model.R
# Outputs: inst/extdata/coding_model/hexamer_logodds.parquet and
#          inst/extdata/coding_model/coding_model_coefficients.parquet
# Record the release and accuracy printed at the end in data-raw/README.md.

devtools::load_all(".", quiet = TRUE)
set.seed(1)
release <- "50"
n_per_class <- 15000L

gff <- download_reference(version = release, file_type = "gff")
fasta <- download_reference(version = release, file_type = "fasta")
# The protein-coding transcript FASTA carries the CDS coordinates in its
# headers (`CDS:start-end`, stop codon included); transcripts.fa.gz does not.
pc_fasta <- download_reference_url(
  paste0(
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_",
    release,
    "/gencode.v",
    release,
    ".pc_transcripts.fa.gz"
  )
)
annot <- prepare_isoformic_annotation(gff)
tx <- read_parquet(fs::path(annot, "annot_data_transcripts.parquet"))
coding_ids <- tx$transcript_id[
  tx$transcript_type %in%
    "protein_coding" &
    !is.na(tx$cds_length) &
    !grepl("cds_(start|end)_NF", ifelse(is.na(tx$tags), "", tx$tags))
]
noncoding_ids <- tx$transcript_id[tx$transcript_type %in% "lncRNA"]
coding_ids <- sample(coding_ids, min(n_per_class, length(coding_ids)))
noncoding_ids <- sample(noncoding_ids, min(n_per_class, length(noncoding_ids)))

read_fasta_subset <- function(path, ids) {
  acc <- new.env(parent = emptyenv())
  acc$headers <- character()
  acc$seqs <- character()
  stream_fasta(
    path,
    keep = function(name) {
      strsplit(name, "|", fixed = TRUE)[[1L]][[1L]] %in% ids
    },
    on_record = function(name, sequence) {
      acc$headers[[length(acc$headers) + 1L]] <- name
      acc$seqs[[length(acc$seqs) + 1L]] <- toupper(sequence)
    }
  )
  list(headers = acc$headers, seqs = acc$seqs)
}
pc <- read_fasta_subset(pc_fasta, coding_ids)
cds_idx <- grepl("CDS:[0-9]+-[0-9]+", pc$headers)
cds_field <- regmatches(pc$headers, regexpr("CDS:[0-9]+-[0-9]+", pc$headers))
cds_start <- as.integer(sub("CDS:([0-9]+)-.*", "\\1", cds_field))
cds_end <- as.integer(sub("CDS:[0-9]+-([0-9]+)", "\\1", cds_field))
cds_seqs <- substring(pc$seqs[cds_idx], cds_start, cds_end)
mrna_seqs <- pc$seqs[cds_idx]
# Keep CDS that start with ATG, end with a stop codon and are in frame
ok <- nchar(cds_seqs) %% 3L == 0L &
  nchar(cds_seqs) >= 90L &
  startsWith(cds_seqs, "ATG") &
  substring(cds_seqs, nchar(cds_seqs) - 2L) %in% c("TAA", "TAG", "TGA")
cds_seqs <- cds_seqs[ok]
mrna_seqs <- mrna_seqs[ok]
noncoding_seqs <- read_fasta_subset(fasta, noncoding_ids)$seqs

# Hold out one third of each class to measure accuracy on unseen transcripts
test_c <- seq_along(mrna_seqs) %% 3L == 0L
test_n <- seq_along(noncoding_seqs) %% 3L == 0L
res <- build_coding_model(
  coding_cds = cds_seqs[!test_c],
  coding = mrna_seqs[!test_c],
  noncoding = noncoding_seqs[!test_n],
  output_dir = "inst/extdata/coding_model"
)
score <- function(seqs) {
  stops <- find_stop_positions(seqs)
  longest <- longest_orfs(seqs, stops)
  coding_potential_scores(seqs, longest)$coding_probability
}
p_c <- score(mrna_seqs[test_c])
p_n <- score(noncoding_seqs[test_n])
cat(
  "GENCODE release:",
  release,
  "\ncoding training sequences:",
  sum(!test_c),
  "\nnoncoding training sequences:",
  sum(!test_n),
  "\ntraining accuracy:",
  round(res$accuracy, 4),
  "\nheld-out sensitivity (coding > 0.364):",
  round(mean(p_c > 0.364), 4),
  "\nheld-out specificity (lncRNA <= 0.364):",
  round(mean(p_n <= 0.364), 4),
  "\n"
)
print(summary(res$model)$coefficients)
