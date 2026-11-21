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
# out-of-fold AUC, repeated with several seeds (the fold assignment changes;
# with the first seed the draws are those of the single run of the first
# version of this script, which fitted each feature set on its own complete
# cases). Uncertainty: bootstrap over genes (genes drawn with replacement,
# all transcripts of a drawn gene kept), applied to the out-of-fold
# predictions, so the models are not refitted in each replicate; each
# replicate averages the AUC over the seeds. The differences between feature
# sets in `contrasts` are paired: the same genes in each replicate. For the
# Karousis pairs also the share of pairs in which the NMD transcript has the
# lower value of each score, with a sign test.
#
# Usage (from the package root, after 01_targets.R):
#   ISOFORMIC_ALPHAGENOME_DIR=/path/to/workdir Rscript data-raw/alphagenome/03_evaluate.R
# Optional environment variables:
#   ISOFORMIC_ALPHAGENOME_GENCODE     GENCODE release (default: the package default)
#   ISOFORMIC_ALPHAGENOME_TARGETS     sub-directory with the targets (default
#                                     "targets")
#   ISOFORMIC_ALPHAGENOME_EVALUATION  sub-directory to write into (default
#                                     "evaluation")
#   ISOFORMIC_ALPHAGENOME_SEEDS       comma-separated seeds of the
#                                     cross-validation (default "1,2,3,4,5")
#   ISOFORMIC_ALPHAGENOME_BOOT        bootstrap replicates (default 2000; 0
#                                     skips the bootstrap)
#   ISOFORMIC_ALPHAGENOME_SCORES      "false" leaves the AlphaGenome features
#                                     out even when <work_dir>/scores/ exists:
#                                     the baseline alone, on all target genes
#                                     (default "true")
# The AlphaGenome scores are always read from <work_dir>/scores/.
# Within a truth every feature set is fitted on the same transcripts: those
# with a value for every feature of every set (single-exon transcripts have
# no splice-site score and are left out), so that the differences between
# feature sets compare the same transcripts.
# Writes <work_dir>/<evaluation>/:
#   features.parquet      one row per transcript
#   auc.parquet           pooled out-of-fold AUC per truth, feature set and seed
#   auc_summary.parquet   per truth and feature set: AUC of the first seed,
#                         mean, minimum and maximum over the seeds, bootstrap
#                         95% interval of the seed mean
#   auc_contrasts.parquet differences between feature sets, with the
#                         bootstrap 95% interval and the share of replicates
#                         at or below zero
#   predictions.parquet   out-of-fold linear predictors per truth, feature set
#                         and seed
#   pairs.parquet         Karousis pairs and sign tests
#   diagnostics.parquet   AlphaGenome coverage of the introns (when scored)

devtools::load_all(".", quiet = TRUE)

work_dir <- Sys.getenv("ISOFORMIC_ALPHAGENOME_DIR", unset = NA)
if (is.na(work_dir) || !nzchar(work_dir)) {
  stop("Set ISOFORMIC_ALPHAGENOME_DIR to a work directory.", call. = FALSE)
}
release <- Sys.getenv("ISOFORMIC_ALPHAGENOME_GENCODE", unset = "")
release <- if (nzchar(release)) release else NULL
seeds <- as.integer(strsplit(
  Sys.getenv("ISOFORMIC_ALPHAGENOME_SEEDS", unset = "1,2,3,4,5"),
  ","
)[[1L]])
n_boot <- as.integer(Sys.getenv("ISOFORMIC_ALPHAGENOME_BOOT", unset = "2000"))
boot_seed <- 20261117L
n_folds <- 5L
target_dir <- fs::path(
  work_dir,
  Sys.getenv("ISOFORMIC_ALPHAGENOME_TARGETS", unset = "targets")
)
out_dir <- fs::dir_create(fs::path(
  work_dir,
  Sys.getenv("ISOFORMIC_ALPHAGENOME_EVALUATION", unset = "evaluation")
))

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

has_ag <- fs::file_exists(fs::path(work_dir, "scores", "tracks.parquet")) &&
  !identical(
    tolower(Sys.getenv("ISOFORMIC_ALPHAGENOME_SCORES", unset = "true")),
    "false"
  )
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

truth_karousis <- features |>
  dplyr::filter(!is.na(.data$karousis)) |>
  dplyr::mutate(y = as.integer(.data$karousis == "nmd"))
truth_nmd_biotype <- features |>
  dplyr::filter(
    .data$transcript_type %in% c("nonsense_mediated_decay", "protein_coding")
  ) |>
  dplyr::mutate(
    y = as.integer(.data$transcript_type == "nonsense_mediated_decay")
  )
truth_ri_biotype <- features |>
  dplyr::filter(
    .data$transcript_type %in% c("retained_intron", "protein_coding")
  ) |>
  dplyr::mutate(y = as.integer(.data$transcript_type == "retained_intron"))
truths <- list(
  karousis = truth_karousis,
  nmd_biotype = truth_nmd_biotype,
  ri_biotype = truth_ri_biotype
)
# The same transcripts for every feature set of a truth
all_vars <- unique(unlist(feature_sets, use.names = FALSE))
n_dropped <- vapply(
  truths,
  function(d) sum(!stats::complete.cases(d[, all_vars, drop = FALSE])),
  integer(1L)
)
truths <- lapply(truths, function(d) {
  d[stats::complete.cases(d[, all_vars, drop = FALSE]), ]
})

# 2. Cross-validated AUC -------------------------------------------------------
auc <- function(score, y) {
  r <- rank(score)
  n1 <- sum(y == 1L)
  n0 <- sum(y == 0L)
  (sum(r[y == 1L]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
# Out-of-fold linear predictors; draws the fold assignment from the current
# random state
cv_predict <- function(data, vars) {
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
    transcript_id = data$transcript_id,
    gene_id = data$gene_id,
    y = data$y,
    pred = pred
  )
}
# Same order of draws as the single-seed version: one set.seed() per seed,
# then truths and feature sets in order
pred_tbl <- dplyr::bind_rows(lapply(seeds, function(s) {
  set.seed(s)
  dplyr::bind_rows(lapply(names(truths), function(truth) {
    dplyr::bind_rows(lapply(names(feature_sets), function(set) {
      dplyr::bind_cols(
        tibble::tibble(seed = s, truth = truth, features = set),
        cv_predict(truths[[truth]], feature_sets[[set]])
      )
    }))
  }))
}))
auc_tbl <- pred_tbl |>
  dplyr::summarise(
    n = dplyr::n(),
    n_positive = sum(.data$y),
    auc = auc(.data$pred, .data$y),
    .by = c("truth", "features", "seed")
  )

# Bootstrap over genes. AUC with integer weights (the number of times the
# gene of each transcript was drawn): sort once, then per replicate sum the
# weights of positives and negatives within groups of tied predictions.
weighted_auc_setup <- function(pred, y, gene_idx) {
  o <- order(pred)
  tie <- cumsum(c(TRUE, diff(pred[o]) != 0))
  list(y = y[o], gene_idx = gene_idx[o], tie = tie, n_tie = max(tie))
}
weighted_auc <- function(setup, mult) {
  w <- mult[setup$gene_idx]
  w1 <- w * setup$y
  w0 <- w - w1
  g1 <- tabulate_weights(setup$tie, w1, setup$n_tie)
  g0 <- tabulate_weights(setup$tie, w0, setup$n_tie)
  below <- cumsum(g0) - g0
  sum(g1 * (below + g0 / 2)) / (sum(g1) * sum(g0))
}
tabulate_weights <- function(group, w, n) {
  out <- numeric(n)
  sums <- rowsum(w, group, reorder = FALSE)
  out[as.integer(rownames(sums))] <- sums[, 1L]
  out
}
contrasts <- list(
  c("rules_and_pwm", "rules"),
  c("baseline", "rules_and_pwm"),
  c("baseline_and_contained", "baseline"),
  c("baseline_and_ag", "baseline"),
  c("baseline_and_ag", "baseline_and_contained"),
  c("rules_and_ag", "rules")
)
contrasts <- Filter(function(x) all(x %in% names(feature_sets)), contrasts)
boot_tbl <- NULL
if (n_boot > 0L) {
  set.seed(boot_seed)
  boot_tbl <- dplyr::bind_rows(lapply(names(truths), function(truth) {
    sub <- pred_tbl[pred_tbl$truth == truth, ]
    genes <- unique(truths[[truth]]$gene_id)
    keys <- unique(sub[, c("features", "seed")])
    setups <- lapply(seq_len(nrow(keys)), function(i) {
      rows <- sub$features == keys$features[i] & sub$seed == keys$seed[i]
      weighted_auc_setup(
        sub$pred[rows],
        sub$y[rows],
        match(sub$gene_id[rows], genes)
      )
    })
    # With every gene drawn once the weighted AUC is the plain one
    unit <- vapply(
      setups,
      weighted_auc,
      numeric(1L),
      mult = rep(1L, length(genes))
    )
    plain <- auc_tbl$auc[match(
      paste(truth, keys$features, keys$seed),
      paste(auc_tbl$truth, auc_tbl$features, auc_tbl$seed)
    )]
    stopifnot(isTRUE(all.equal(unit, plain)))
    values <- matrix(NA_real_, n_boot, nrow(keys))
    for (b in seq_len(n_boot)) {
      mult <- tabulate(
        sample.int(length(genes), replace = TRUE),
        nbins = length(genes)
      )
      values[b, ] <- vapply(setups, weighted_auc, numeric(1L), mult = mult)
    }
    # Mean over seeds within each replicate
    set_mean <- vapply(
      names(feature_sets),
      function(set) rowMeans(values[, keys$features == set, drop = FALSE]),
      numeric(n_boot)
    )
    tibble::tibble(
      truth = truth,
      replicate = rep(seq_len(n_boot), times = ncol(set_mean)),
      features = rep(colnames(set_mean), each = n_boot),
      auc = as.vector(set_mean)
    )
  }))
}
auc_summary <- auc_tbl |>
  dplyr::summarise(
    n = .data$n[1L],
    n_positive = .data$n_positive[1L],
    auc_first_seed = .data$auc[.data$seed == seeds[1L]],
    auc_mean = mean(.data$auc),
    auc_min = min(.data$auc),
    auc_max = max(.data$auc),
    .by = c("truth", "features")
  )
contrast_tbl <- NULL
if (!is.null(boot_tbl)) {
  boot_ci <- boot_tbl |>
    dplyr::summarise(
      ci_low = stats::quantile(.data$auc, 0.025, names = FALSE),
      ci_high = stats::quantile(.data$auc, 0.975, names = FALSE),
      .by = c("truth", "features")
    )
  auc_summary <- auc_summary |>
    dplyr::left_join(boot_ci, by = c("truth", "features"))
  contrast_tbl <- dplyr::bind_rows(lapply(contrasts, function(x) {
    dplyr::bind_rows(lapply(names(truths), function(truth) {
      a <- boot_tbl[boot_tbl$truth == truth & boot_tbl$features == x[1L], ]
      b <- boot_tbl[boot_tbl$truth == truth & boot_tbl$features == x[2L], ]
      stopifnot(identical(a$replicate, b$replicate))
      delta <- a$auc - b$auc
      point <- auc_summary[auc_summary$truth == truth, ]
      seed_delta <- auc_tbl$auc[
        auc_tbl$truth == truth & auc_tbl$features == x[1L]
      ] -
        auc_tbl$auc[auc_tbl$truth == truth & auc_tbl$features == x[2L]]
      tibble::tibble(
        truth = truth,
        features = x[1L],
        minus = x[2L],
        delta_mean = point$auc_mean[point$features == x[1L]] -
          point$auc_mean[point$features == x[2L]],
        delta_min = min(seed_delta),
        delta_max = max(seed_delta),
        ci_low = stats::quantile(delta, 0.025, names = FALSE),
        ci_high = stats::quantile(delta, 0.975, names = FALSE),
        share_at_or_below_zero = mean(delta <= 0)
      )
    }))
  }))
}

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
write_parquet(auc_summary, fs::path(out_dir, "auc_summary.parquet"))
write_parquet(pred_tbl, fs::path(out_dir, "predictions.parquet"))
write_parquet(pair_tbl, fs::path(out_dir, "pairs.parquet"))
if (!is.null(contrast_tbl)) {
  write_parquet(contrast_tbl, fs::path(out_dir, "auc_contrasts.parquet"))
}

if (has_ag) {
  write_parquet(diagnostics, fs::path(out_dir, "diagnostics.parquet"))
  cli::cli_h2("AlphaGenome diagnostics")
  print(as.data.frame(diagnostics), digits = 3)
}
cli::cli_h2(if (has_ag) "Baseline and AlphaGenome" else "Baseline")
cli::cli_inform(c(
  i = "Seeds: {seeds}; bootstrap replicates over genes: {n_boot}",
  i = "Transcripts without a value for some feature, left out: {paste(names(n_dropped), n_dropped, sep = ' ', collapse = ', ')}"
))
print(as.data.frame(auc_summary), digits = 3)
if (!is.null(contrast_tbl)) {
  print(as.data.frame(contrast_tbl), digits = 3)
}
print(as.data.frame(pair_tbl), digits = 3)
