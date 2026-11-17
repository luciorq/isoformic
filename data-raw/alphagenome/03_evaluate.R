# AlphaGenome pilot benchmark, step 3: evaluation
#
# The baseline: what the structure rules and the
# splice-site strength of the annotation (splice_site_strength(), no
# AlphaGenome) already tell about productive against unproductive isoforms on
# the pilot genes. When <work_dir>/scores/ of step 2 exists, the AlphaGenome
# features of 03_features_alphagenome.R are added and every model is fitted
# on the scored genes only, so that the rows compare:
#   junction   predicted junction count of an intron relative to the largest
#              of the gene, mean over tracks; 0 when not predicted
#   site       splice-site probability (the lower of donor and acceptor)
#   usage      splice-site usage (the lower of donor and acceptor)
#   each as the weakest over all introns of the transcript and over the
#   introns that the canonical transcript does not have ("novel")
#   ag_retention  for a canonical intron inside an exon of the transcript:
#              log2 of its predicted RNA-seq coverage over the coverage of
#              the canonical exons; the largest per transcript
#
# Truths:
#   karousis        validated NMD against non-NMD transcripts (Karousis 2021)
#   nmd_biotype     nonsense_mediated_decay against protein_coding
#   ri_biotype      retained_intron against protein_coding
# Models: logistic regression, cross-validation grouped by gene, pooled
# out-of-fold AUC. For the Karousis pairs also the share of pairs in which the
# NMD transcript has the lower value of each score.
#
# Usage (from the package root, after 01_targets.R):
#   ISOFORMIC_ALPHAGENOME_DIR=/path/to/workdir Rscript data-raw/alphagenome/03_evaluate.R
# Writes <work_dir>/evaluation/{features,auc,pairs,diagnostics}.parquet.

devtools::load_all(".", quiet = TRUE)

work_dir <- Sys.getenv("ISOFORMIC_ALPHAGENOME_DIR", unset = NA)
if (is.na(work_dir) || !nzchar(work_dir)) {
  stop("Set ISOFORMIC_ALPHAGENOME_DIR to a work directory.", call. = FALSE)
}
release <- Sys.getenv("ISOFORMIC_ALPHAGENOME_GENCODE", unset = "")
release <- if (nzchar(release)) release else NULL
seed <- as.integer(Sys.getenv("ISOFORMIC_ALPHAGENOME_SEED", unset = "1"))
n_folds <- 5L
target_dir <- fs::path(work_dir, "targets")
out_dir <- fs::dir_create(fs::path(work_dir, "evaluation"))

tx <- read_parquet(fs::path(target_dir, "transcripts.parquet"))
pairs <- read_parquet(fs::path(target_dir, "karousis_pairs.parquet"))

# 1. Splice-site strength per intron ------------------------------------------
gff <- download_reference(version = release, file_type = "gff")
annot_dir <- prepare_isoformic_annotation(gff)
if (!has_splice_sites(annot_dir)) {
  genome <- download_reference(
    version = release,
    file_type = "genome_fasta",
    region = "primary"
  )
  extract_splice_sites(annot_dir, genome_fasta = genome)
}
introns <- splice_site_strength(annot_dir, by = "intron") |>
  dplyr::filter(.data$transcript_id %in% tx$transcript_id)

# An intron is "novel" when the canonical transcript of the gene does not
# have it: the sites that an unproductive isoform adds are among these.
canonical <- introns |>
  dplyr::filter(
    .data$transcript_id %in% tx$transcript_id[tx$is_canonical %in% TRUE]
  ) |>
  dplyr::distinct(.data$gene_id, .data$intron_start, .data$intron_end) |>
  dplyr::mutate(in_canonical = TRUE)
introns <- introns |>
  dplyr::left_join(
    canonical,
    by = c("gene_id", "intron_start", "intron_end")
  ) |>
  dplyr::mutate(novel = is.na(.data$in_canonical))

min_or <- function(x, keep, default) {
  x <- x[keep & !is.na(x)]
  if (length(x)) min(x) else default
}
# A transcript without novel introns gets the median score of all introns:
# no evidence in either direction; `n_novel` carries the difference.
mid_donor <- stats::median(introns$donor_score, na.rm = TRUE)
mid_acceptor <- stats::median(introns$acceptor_score, na.rm = TRUE)
pwm <- introns |>
  dplyr::summarise(
    weakest_donor = min_or(.data$donor_score, TRUE, mid_donor),
    weakest_acceptor = min_or(.data$acceptor_score, TRUE, mid_acceptor),
    novel_donor = min_or(.data$donor_score, .data$novel, mid_donor),
    novel_acceptor = min_or(.data$acceptor_score, .data$novel, mid_acceptor),
    n_novel = sum(.data$novel),
    n_introns = dplyr::n(),
    .by = "transcript_id"
  )

features <- tx |>
  dplyr::select(
    "transcript_id",
    "gene_id",
    "transcript_type",
    "is_canonical",
    "karousis",
    "rule_productivity"
  ) |>
  dplyr::left_join(pwm, by = "transcript_id") |>
  dplyr::mutate(
    rule_unproductive = as.integer(.data$rule_productivity == "unproductive"),
    rule_altered = as.integer(.data$rule_productivity == "altered_productive"),
    rule_uncertain = as.integer(.data$rule_productivity == "uncertain"),
    has_novel = as.integer(.data$n_novel > 0L),
    log_introns = log2(.data$n_introns)
  )

has_ag <- fs::file_exists(fs::path(work_dir, "scores", "tracks.parquet"))
if (has_ag) {
  source("data-raw/alphagenome/03_features_alphagenome.R")
  features <- features |>
    dplyr::filter(.data$gene_id %in% scored) |>
    dplyr::left_join(ag, by = "transcript_id")
}

feature_sets <- list(
  rules = c("rule_unproductive", "rule_altered", "rule_uncertain"),
  pwm = c(
    "weakest_donor",
    "weakest_acceptor",
    "novel_donor",
    "novel_acceptor",
    "has_novel"
  ),
  # Control: the weakest site falls with the number of introns
  intron_number = "log_introns"
)
feature_sets$pwm_and_number <- c(feature_sets$pwm, "log_introns")
feature_sets$rules_and_pwm <- c(feature_sets$rules, feature_sets$pwm)
feature_sets$baseline <- c(feature_sets$rules_and_pwm, "log_introns")
if (has_ag) {
  feature_sets$ag_junction <- c("ag_weakest_junction", "ag_novel_junction")
  feature_sets$ag_site <- c(
    "ag_weakest_site",
    "ag_novel_site",
    "ag_weakest_usage",
    "ag_novel_usage"
  )
  feature_sets$ag_retention <- "ag_retention"
  # Structure only: is there a canonical intron inside an exon at all?
  feature_sets$contained_intron <- "has_contained"
  feature_sets$ag <- c(
    feature_sets$ag_junction,
    feature_sets$ag_site,
    "ag_retention",
    "has_novel"
  )
  feature_sets$rules_and_ag <- c(feature_sets$rules, feature_sets$ag)
  feature_sets$baseline_and_contained <- c(
    feature_sets$baseline,
    "has_contained"
  )
  feature_sets$baseline_and_ag <- c(
    feature_sets$baseline_and_contained,
    setdiff(feature_sets$ag, "has_novel")
  )
}

truths <- list(
  karousis = features |>
    dplyr::filter(!is.na(.data$karousis)) |>
    dplyr::mutate(y = as.integer(.data$karousis == "nmd")),
  nmd_biotype = features |>
    dplyr::filter(
      .data$transcript_type %in% c("nonsense_mediated_decay", "protein_coding")
    ) |>
    dplyr::mutate(
      y = as.integer(.data$transcript_type == "nonsense_mediated_decay")
    ),
  ri_biotype = features |>
    dplyr::filter(
      .data$transcript_type %in% c("retained_intron", "protein_coding")
    ) |>
    dplyr::mutate(y = as.integer(.data$transcript_type == "retained_intron"))
)

# 2. Cross-validated AUC -------------------------------------------------------
auc <- function(score, y) {
  r <- rank(score)
  n1 <- sum(y == 1L)
  n0 <- sum(y == 0L)
  (sum(r[y == 1L]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
cv_auc <- function(data, vars) {
  data <- data[stats::complete.cases(data[, vars, drop = FALSE]), ]
  genes <- unique(data$gene_id)
  fold <- sample(rep_len(seq_len(n_folds), length(genes)))
  fold <- fold[match(data$gene_id, genes)]
  pred <- numeric(nrow(data))
  form <- stats::reformulate(vars, response = "y")
  for (f in seq_len(n_folds)) {
    fit <- suppressWarnings(stats::glm(
      form,
      data = data[fold != f, ],
      family = stats::binomial()
    ))
    pred[fold == f] <- suppressWarnings(stats::predict(
      fit,
      newdata = data[fold == f, ]
    ))
  }
  tibble::tibble(
    n = nrow(data),
    n_positive = sum(data$y),
    auc = auc(pred, data$y)
  )
}
set.seed(seed)
auc_tbl <- dplyr::bind_rows(lapply(names(truths), function(truth) {
  dplyr::bind_rows(lapply(names(feature_sets), function(set) {
    dplyr::bind_cols(
      tibble::tibble(truth = truth, features = set),
      cv_auc(truths[[truth]], feature_sets[[set]])
    )
  }))
}))

# 3. Karousis pairs: is the NMD transcript the one with the weaker site? ------
scores <- c(
  "weakest_donor",
  "weakest_acceptor",
  "novel_donor",
  "novel_acceptor"
)
if (has_ag) {
  scores <- c(scores, feature_sets$ag_junction, feature_sets$ag_site)
}
at <- function(ids, column) {
  features[[column]][match(ids, features$transcript_id)]
}
pair_tbl <- dplyr::bind_rows(lapply(scores, function(column) {
  delta <- at(pairs$nmd_tx, column) - at(pairs$non_nmd_tx, column)
  delta <- delta[!is.na(delta)]
  informative <- delta != 0
  tibble::tibble(
    score = column,
    n_pairs = length(delta),
    n_differ = sum(informative),
    nmd_lower = sum(delta < 0),
    share_nmd_lower = sum(delta < 0) / sum(informative),
    p_sign = stats::binom.test(sum(delta < 0), sum(informative))$p.value
  )
}))

write_parquet(features, fs::path(out_dir, "features.parquet"))
write_parquet(auc_tbl, fs::path(out_dir, "auc.parquet"))
write_parquet(pair_tbl, fs::path(out_dir, "pairs.parquet"))

if (has_ag) {
  write_parquet(diagnostics, fs::path(out_dir, "diagnostics.parquet"))
  cli::cli_h2("AlphaGenome diagnostics")
  print(as.data.frame(diagnostics), digits = 3)
}
cli::cli_h2(if (has_ag) "Baseline and AlphaGenome" else "Baseline")
print(as.data.frame(auc_tbl), digits = 3)
print(as.data.frame(pair_tbl), digits = 3)
