# Molecular features of isoforms and the stability of their inferential replicates
#
# Part A computes, from the GENCODE v50 annotation and transcript sequences,
# features that can make an isoform hard to quantify: length, exons, isoforms of
# the gene, GC content, support tags, alternative and NAGNAG splice sites, the
# exonic bases unique to the isoform, the distance to its most similar sibling
# and whether its exons lie inside another isoform's.
#
# Part B measures, on three GSE129540 samples quantified with 100 Gibbs draws
# (and with 30 draws, quants_terminus/), how stable the replicates make each
# isoform: the inferential relative variance with 100 draws, whether two
# independent sets of 30 draws put the isoform on different sides of the 0.5
# cutoff, how much the variance estimate moves between them, whether its
# TreeTerminus tree changes between 30 and 100 draws, and its unique fragments.
#
# Part C summarises both into small tables for the article
# vignettes/articles/isoformic-quantification-stability.qmd.
#
# Usage:
#   ISOFORMIC_KNOCKDOWN_DIR=/path/to/workdir Rscript data-raw/knockdown/isoform_features.R
#
# Needs g100/<run> (Salmon, 100 draws), quants_terminus/<run> (30 draws) and
# tt_single/<run>_<30|100> (TreeTerminus group per sample) for the three runs.
suppressPackageStartupMessages({
  library(isoformic)
  library(dplyr)
})

work_dir <- Sys.getenv("ISOFORMIC_KNOCKDOWN_DIR", unset = NA)
if (is.na(work_dir) || !dir.exists(work_dir)) {
  stop("Set ISOFORMIC_KNOCKDOWN_DIR to the work directory.", call. = FALSE)
}
runs <- c("SRR8870780", "SRR8870781", "SRR8870783")
cache <- file.path(work_dir, "analysis", "isoform_features_A.rds")
dir.create(dirname(cache), recursive = TRUE, showWarnings = FALSE)

# Part A: annotation and sequence features -------------------------------------
annot <- prepare_isoformic_annotation(download_reference(
  version = "50",
  file_type = "gff"
))
tx <- read_parquet(fs::path(annot, "annot_data_transcripts.parquet"))
if (file.exists(cache)) {
  features <- readRDS(cache)
} else {
  exons <- read_parquet(
    fs::path(annot, "annot_data_exons.parquet"),
    columns = c(
      "transcript_id",
      "exon_number",
      "start_pos",
      "end_pos",
      "strand"
    )
  )
  exons$gene_id <- tx$gene_id[match(exons$transcript_id, tx$transcript_id)]
  by_gene <- split(exons, exons$gene_id)
  gene_features <- function(e) {
    ids <- unique(e$transcript_id)
    out <- data.frame(
      transcript_id = ids,
      unique_bp = NA_real_,
      closest_bp = NA_real_,
      contained = FALSE,
      nagnag = FALSE,
      near_acceptor = FALSE,
      near_donor = FALSE
    )
    if (length(ids) < 2L) {
      out$unique_bp <- tx$tx_length[match(ids, tx$transcript_id)]
      return(out)
    }
    # Disjoint exonic segments of the gene and which isoform covers each
    breaks <- sort(unique(c(e$start_pos, e$end_pos + 1L)))
    seg_start <- breaks[-length(breaks)]
    seg_end <- breaks[-1L] - 1L
    covered <- matrix(FALSE, length(ids), length(seg_start))
    row <- match(e$transcript_id, ids)
    for (k in seq_len(nrow(e))) {
      covered[
        row[k],
        seg_start >= e$start_pos[k] & seg_end <= e$end_pos[k]
      ] <- TRUE
    }
    keep <- colSums(covered) > 0L
    covered <- covered[, keep, drop = FALSE]
    len <- (seg_end - seg_start + 1L)[keep]
    m <- covered * 1
    shared <- m %*% (t(m) * len)
    own <- diag(shared)
    dist <- outer(own, own, "+") - 2 * shared
    diag(dist) <- Inf
    out$closest_bp <- apply(dist, 1L, min)
    only <- colSums(covered) == 1L
    out$unique_bp <- as.vector(m[, only, drop = FALSE] %*% len[only])
    diag(shared) <- -1
    out$contained <- apply(shared == own, 1L, any)
    # Splice sites in transcription order: acceptors start exons 2..n,
    # donors end exons 1..n-1; on the minus strand the genomic ends swap
    plus <- e$strand == "+"
    last <- ave(e$exon_number, e$transcript_id, FUN = max)
    acc <- ifelse(plus, e$start_pos, e$end_pos)[e$exon_number > 1L]
    acc_tx <- e$transcript_id[e$exon_number > 1L]
    don <- ifelse(plus, e$end_pos, e$start_pos)[e$exon_number < last]
    don_tx <- e$transcript_id[e$exon_number < last]
    near <- function(pos, owner, lo, hi) {
      hit <- vapply(
        seq_along(pos),
        function(k) {
          d <- abs(pos[k] - pos[owner != owner[k]])
          any(d >= lo & d <= hi)
        },
        logical(1L)
      )
      tapply(hit, owner, any)
    }
    if (length(acc)) {
      na <- near(acc, acc_tx, 3L, 3L)
      out$nagnag <- out$transcript_id %in% names(na)[na]
      nr <- near(acc, acc_tx, 1L, 30L)
      out$near_acceptor <- out$transcript_id %in% names(nr)[nr]
    }
    if (length(don)) {
      nd <- near(don, don_tx, 1L, 30L)
      out$near_donor <- out$transcript_id %in% names(nd)[nd]
    }
    out
  }
  message("structure of ", length(by_gene), " genes")
  structure <- dplyr::bind_rows(lapply(by_gene, gene_features))

  # GC content from the transcript sequences
  fasta <- download_reference(version = "50", file_type = "fasta")
  gc <- new.env(parent = emptyenv())
  gc$id <- character()
  gc$value <- numeric()
  isoformic:::stream_fasta(
    fasta,
    keep = function(name) TRUE,
    on_record = function(name, sequence) {
      gc$id[[length(gc$id) + 1L]] <- sub("\\|.*$", "", name)
      bases <- nchar(sequence)
      gc$value[[length(gc$value) + 1L]] <- (bases -
        nchar(gsub("[GCgc]", "", sequence))) /
        bases
    }
  )
  tags <- ifelse(is.na(tx$tags), "", tx$tags)
  features <- tibble::tibble(
    transcript_id = tx$transcript_id,
    gene_id = tx$gene_id,
    transcript_type = tx$transcript_type,
    tx_length = tx$tx_length,
    n_exons = tx$n_exons,
    mane = grepl("MANE_Select", tags, fixed = TRUE),
    basic = grepl("basic", tags, fixed = TRUE),
    tsl = suppressWarnings(as.integer(sub(
      "^([0-9]).*",
      "\\1",
      tx$transcript_support_level
    ))),
    gc = gc$value[match(tx$transcript_id, gc$id)]
  ) |>
    dplyr::add_count(.data$gene_id, name = "isoforms_in_gene") |>
    dplyr::left_join(structure, by = "transcript_id")
  saveRDS(features, cache)
}
message(nrow(features), " transcripts with features")

# Part B: stability of the replicates, per isoform and sample ------------------
infrv_cols <- function(m, cols) {
  v <- m[, cols, drop = FALSE]
  mu <- rowMeans(v)
  va <- rowSums((v - mu)^2) / (length(cols) - 1L)
  pmax(va - mu, 0) / (mu + 5) + 0.01
}
tree_key <- function(dir) {
  f <- list.files(
    dir,
    pattern = "group_nwk.txt",
    recursive = TRUE,
    full.names = TRUE
  )
  trees <- read_treeterminus(
    f,
    salmon_dir = file.path(work_dir, "quants_terminus", runs[[1L]])
  )
  leaves <- trees[!is.na(trees$transcript_id), ]
  key <- tapply(leaves$transcript_id, leaves$tree_id, function(v) {
    paste(sort(v), collapse = ";")
  })
  stats::setNames(key[as.character(leaves$tree_id)], leaves$transcript_id)
}
per_sample <- lapply(runs, function(r) {
  dir <- c(s = file.path(work_dir, "g100", r))
  quant <- read_tabular(file.path(dir, "quant.sf"))
  ids <- sub("\\|.*$", "", quant$Name)
  expressed <- ids[quant$NumReads >= 5]
  reps <- read_inferential_replicates(dir, transcripts = expressed)
  m <- do.call(cbind, lapply(reps, function(x) x[, 1L]))
  a30 <- infrv_cols(m, 1:30)
  b30 <- infrv_cols(m, 31:60)
  k30 <- tree_key(file.path(work_dir, "tt_single", paste0(r, "_30")))
  k100 <- tree_key(file.path(work_dir, "tt_single", paste0(r, "_100")))
  in_trees <- rownames(m) %in% names(k30) | rownames(m) %in% names(k100)
  same_tree <- rownames(m) %in%
    names(k30) &
    rownames(m) %in% names(k100) &
    k30[rownames(m)] %in% k100[rownames(m)]
  support <- equivalence_support(dir, transcripts = rownames(m))
  tibble::tibble(
    run = r,
    transcript_id = rownames(m),
    reads = quant$NumReads[match(rownames(m), ids)],
    inf_rv = infrv_cols(m, seq_len(ncol(m))),
    rv_disagreement = abs(log((a30 + 0.01) / (b30 + 0.01))),
    cutoff_flip = (a30 >= 0.5) != (b30 >= 0.5),
    tree_changed = ifelse(in_trees, !same_tree, NA),
    unique_fraction = support$unique_fraction[match(
      rownames(m),
      support$transcript_id
    )]
  )
})
stability <- dplyr::bind_rows(per_sample)
# One row per isoform: means over the samples where it is expressed
iso <- stability |>
  dplyr::group_by(.data$transcript_id) |>
  dplyr::summarise(
    samples = dplyr::n(),
    reads = mean(.data$reads),
    inf_rv = mean(.data$inf_rv),
    rv_disagreement = mean(.data$rv_disagreement),
    cutoff_flip = mean(.data$cutoff_flip),
    tree_changed = mean(.data$tree_changed, na.rm = TRUE),
    unique_fraction = mean(.data$unique_fraction, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::left_join(features, by = "transcript_id")
iso$tree_changed[is.nan(iso$tree_changed)] <- NA
iso$unique_fraction[is.nan(iso$unique_fraction)] <- NA
saveRDS(iso, file.path(work_dir, "analysis", "isoform_stability.rds"))
message(
  nrow(iso),
  " isoforms expressed (>= 5 reads) in at least one of ",
  length(runs),
  " samples"
)

# Part C: tables for the article -----------------------------------------------
main_types <- c(
  "protein_coding",
  "nonsense_mediated_decay",
  "retained_intron",
  "protein_coding_CDS_not_defined",
  "lncRNA"
)
iso$type <- ifelse(
  iso$transcript_type %in% main_types,
  iso$transcript_type,
  "other"
)
iso$type <- factor(iso$type, levels = c(main_types, "other"))

by_type <- iso |>
  dplyr::group_by(.data$type) |>
  dplyr::summarise(
    isoforms = dplyr::n(),
    median_reads = stats::median(.data$reads),
    median_length = stats::median(.data$tx_length),
    median_exons = stats::median(.data$n_exons),
    median_isoforms_in_gene = stats::median(.data$isoforms_in_gene),
    median_gc = stats::median(.data$gc, na.rm = TRUE),
    median_unique_bp = stats::median(.data$unique_bp, na.rm = TRUE),
    median_closest_bp = stats::median(.data$closest_bp[is.finite(
      .data$closest_bp
    )]),
    contained = mean(.data$contained, na.rm = TRUE),
    nagnag = mean(.data$nagnag, na.rm = TRUE),
    near_acceptor = mean(.data$near_acceptor, na.rm = TRUE),
    near_donor = mean(.data$near_donor, na.rm = TRUE),
    mane = mean(.data$mane),
    median_unique_fraction = stats::median(.data$unique_fraction, na.rm = TRUE),
    median_inf_rv = stats::median(.data$inf_rv),
    median_rv_disagreement = stats::median(.data$rv_disagreement),
    cutoff_flip = mean(.data$cutoff_flip),
    tree_changed = mean(.data$tree_changed, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::mutate(dplyr::across(
    where(is.numeric) & !c("isoforms"),
    ~ signif(.x, 3)
  )) |>
  dplyr::mutate(type = as.character(.data$type))

bin <- function(x, breaks, labels) {
  cut(x, breaks, labels = labels, right = FALSE, include.lowest = TRUE)
}
binned <- list(
  reads = bin(
    iso$reads,
    c(5, 20, 100, 1000, Inf),
    c("5-20", "20-100", "100-1000", ">1000")
  ),
  tx_length = bin(
    iso$tx_length,
    c(0, 500, 1000, 2000, 4000, Inf),
    c("<500", "500-1k", "1-2k", "2-4k", ">4k")
  ),
  n_exons = bin(
    iso$n_exons,
    c(1, 2, 3, 6, 11, Inf),
    c("1", "2", "3-5", "6-10", ">10")
  ),
  isoforms_in_gene = bin(
    iso$isoforms_in_gene,
    c(1, 2, 4, 11, 31, Inf),
    c("1", "2-3", "4-10", "11-30", ">30")
  ),
  gc = bin(
    iso$gc,
    c(0, 0.4, 0.45, 0.5, 0.55, 1),
    c("<40%", "40-45%", "45-50%", "50-55%", ">55%")
  ),
  unique_bp = bin(
    iso$unique_bp,
    c(0, 1, 50, 200, 1000, Inf),
    c("0", "1-50", "50-200", "200-1000", ">1000")
  ),
  closest_bp = bin(
    ifelse(is.finite(iso$closest_bp), iso$closest_bp, NA),
    c(0, 10, 100, 500, 2000, Inf),
    c("<10", "10-100", "100-500", "500-2000", ">2000")
  ),
  unique_fraction = bin(
    iso$unique_fraction,
    c(0, 0.001, 0.05, 0.25, 0.75, 1),
    c("0", "<5%", "5-25%", "25-75%", ">75%")
  ),
  tsl = factor(ifelse(is.na(iso$tsl), "NA", iso$tsl), levels = c(1:5, "NA")),
  contained = factor(
    iso$contained,
    levels = c(FALSE, TRUE),
    labels = c("no", "yes")
  ),
  nagnag = factor(iso$nagnag, levels = c(FALSE, TRUE), labels = c("no", "yes")),
  near_acceptor = factor(
    iso$near_acceptor,
    levels = c(FALSE, TRUE),
    labels = c("no", "yes")
  ),
  near_donor = factor(
    iso$near_donor,
    levels = c(FALSE, TRUE),
    labels = c("no", "yes")
  ),
  mane = factor(iso$mane, levels = c(FALSE, TRUE), labels = c("no", "yes"))
)
feature_bins <- dplyr::bind_rows(lapply(names(binned), function(f) {
  d <- iso
  d$level <- binned[[f]]
  d <- d[!is.na(d$level), ]
  both <- dplyr::bind_rows(d, dplyr::mutate(d, type = "all"))
  both |>
    dplyr::group_by(type = as.character(.data$type), .data$level) |>
    dplyr::summarise(
      isoforms = dplyr::n(),
      median_inf_rv = stats::median(.data$inf_rv),
      median_rv_disagreement = stats::median(.data$rv_disagreement),
      cutoff_flip = mean(.data$cutoff_flip),
      tree_changed = mean(.data$tree_changed, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      feature = f,
      level_order = as.integer(.data$level),
      level = as.character(.data$level)
    )
})) |>
  dplyr::relocate("feature", "type", "level", "level_order") |>
  dplyr::mutate(dplyr::across(
    c("median_inf_rv", "median_rv_disagreement", "cutoff_flip", "tree_changed"),
    ~ signif(.x, 3)
  ))

# Multivariable models with standardised covariates, so effects are comparable
z <- function(x) (x - mean(x, na.rm = TRUE)) / stats::sd(x, na.rm = TRUE)
md <- iso |>
  dplyr::filter(
    is.finite(.data$closest_bp),
    !is.na(.data$gc),
    !is.na(.data$unique_fraction)
  ) |>
  dplyr::transmute(
    type = .data$type,
    log_reads = z(log10(.data$reads)),
    log_length = z(log10(.data$tx_length)),
    log_exons = z(log10(.data$n_exons)),
    log_isoforms = z(log10(.data$isoforms_in_gene)),
    gc = z(.data$gc),
    log_unique_bp = z(log10(.data$unique_bp + 1)),
    log_closest_bp = z(log10(.data$closest_bp + 1)),
    unique_fraction = z(.data$unique_fraction),
    contained = .data$contained,
    nagnag = .data$nagnag,
    near_acceptor = .data$near_acceptor,
    near_donor = .data$near_donor,
    mane = .data$mane,
    log_inf_rv = log10(.data$inf_rv),
    log_rv_disagreement = log10(.data$rv_disagreement + 0.01),
    tree_changed = .data$tree_changed,
    samples = .data$samples
  )
rhs <- "type + log_reads + log_length + log_exons + log_isoforms + gc + log_unique_bp + log_closest_bp + unique_fraction + contained + nagnag + near_acceptor + near_donor + mane"
fits <- list(
  log_inf_rv = stats::lm(
    stats::as.formula(paste("log_inf_rv ~", rhs)),
    data = md
  ),
  log_rv_disagreement = stats::lm(
    stats::as.formula(paste("log_rv_disagreement ~", rhs)),
    data = md
  ),
  tree_changed = stats::glm(
    stats::as.formula(paste("tree_changed ~", rhs)),
    data = md[!is.na(md$tree_changed), ],
    family = stats::quasibinomial(),
    weights = md$samples[!is.na(md$tree_changed)]
  )
)
model_terms <- dplyr::bind_rows(lapply(names(fits), function(nm) {
  co <- summary(fits[[nm]])$coefficients
  tibble::tibble(
    outcome = nm,
    term = rownames(co),
    estimate = signif(co[, 1L], 3),
    std_error = signif(co[, 2L], 3),
    p_value = signif(co[, 4L], 3)
  )
}))
importance <- dplyr::bind_rows(lapply(names(fits), function(nm) {
  dr <- stats::drop1(fits[[nm]])
  dev <- if ("Deviance" %in% names(dr)) dr$Deviance else dr$`Sum of Sq`
  base <- if ("Deviance" %in% names(dr)) dr$Deviance[[1L]] else dr$RSS[[1L]]
  added <- if ("Deviance" %in% names(dr)) dev - base else dev
  tibble::tibble(outcome = nm, term = rownames(dr), deviance_explained = added)[
    -1L,
  ] |>
    dplyr::mutate(
      share = signif(
        .data$deviance_explained / sum(.data$deviance_explained),
        3
      )
    ) |>
    dplyr::arrange(dplyr::desc(.data$share))
}))
r2 <- tibble::tibble(
  outcome = c("log_inf_rv", "log_rv_disagreement"),
  r_squared = signif(
    c(
      summary(fits$log_inf_rv)$r.squared,
      summary(fits$log_rv_disagreement)$r.squared
    ),
    3
  ),
  isoforms = nrow(md)
)

out <- fs::path("inst", "extdata", "knockdown", "stability")
fs::dir_create(out)
write_parquet(by_type, fs::path(out, "by_type.parquet"))
write_parquet(feature_bins, fs::path(out, "feature_bins.parquet"))
write_parquet(model_terms, fs::path(out, "model_terms.parquet"))
write_parquet(importance, fs::path(out, "importance.parquet"))
write_parquet(r2, fs::path(out, "model_fit.parquet"))
# Tree changes matter for the isoforms a cut would group: split by uncertainty
tree_uncertainty <- iso |>
  dplyr::mutate(uncertain = .data$inf_rv >= 0.5) |>
  dplyr::group_by(.data$type, .data$uncertain) |>
  dplyr::summarise(
    isoforms = dplyr::n(),
    cutoff_flip = signif(mean(.data$cutoff_flip), 3),
    tree_changed = signif(mean(.data$tree_changed, na.rm = TRUE), 3),
    .groups = "drop"
  ) |>
  dplyr::mutate(type = as.character(.data$type))
write_parquet(tree_uncertainty, fs::path(out, "tree_uncertainty.parquet"))
print(as.data.frame(by_type))
print(as.data.frame(importance))
print(as.data.frame(r2))
