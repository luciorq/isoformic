# Placenta store: import, haplotype folding, QC, universe U and strata.
#
# Purpose: shared prerequisite (package S0) of the placenta analyses of the
# uncertainty evidence study. Imports the 25 placenta Salmon runs (GENCODE
# 50, 30 Gibbs draws each), folds patch and haplotype copies onto their
# primary-chromosome transcript, writes sample QC (depth, mapping, library
# type, fetal sex), checks the A/B technical pairs, defines the universe U
# and the fixed strata from the 8 preeclampsia samples only (the null never
# tests them), validates the matrix code path against edgeR::catchSalmon and
# fishpond::computeInfRV, and stores everything as Parquet.
#
# Inputs:
#   prereg.yaml and the package snapshot (01_preregister.R)
#   $ISOFORMIC_V50_DIR/human/placenta/quants/<run>/   Salmon 2 output
#   $ISOFORMIC_V50_DIR/human/ref/gencode.v50.chr_patch_hapl_scaff.annotation.gtf.gz
#   $ISOFORMIC_V50_DIR/human/ref/salmon_index_v50/duplicate_clusters.tsv
#   ~/.cache/R/isoformic/gencode_v50_chr_patch_hapl_scaff_annotation_gff3
#     (annotation directory, for the productivity classes only)
#
# Env vars:
#   ISOFORMIC_UE_OUT      output root (required)
#   ISOFORMIC_UE_THREADS  threads (DuckDB, data.table), default 4
#   ISOFORMIC_UE_PKG      package snapshot to load (default from prereg.yaml)
#   ISOFORMIC_V50_DIR     work directory of the Salmon runs (required)
#   ISOFORMIC_UE_ANNOT    annotation directory for productivity, default the
#                         GENCODE 50 directory in the isoformic cache
#
# Memory (estimate): ~12 transcript x 25 matrices of 654,828 doubles (~1.6
# GB), one sample's draws at a time (654,828 x 30 x 8 B = 157 MB, plus the
# folded copy), long Parquet tables of 16.4 M rows (~1-1.5 GB transient),
# edgeR::catchSalmon on up to 25 samples (~2-3 GB): 6-10 GB peak. Measured:
# 222 min, peak 9.1 GB, most of it identifiability (checkpointed per sample
# in store/identifiability_pe/, so a rerun resumes there).
#
# Run from the repository root, inside a memory cap of about 16 GB and not
# alongside other heavy jobs (see README.md):
#   pixi run --frozen Rscript data-raw/uncertainty_evidence/10_placenta_prepare.R
#
# Outputs ($ISOFORMIC_UE_OUT):
#   store/  tx2gene, fold_map, point_counts(_unfolded), moments(_unfolded),
#           draws/<run>.parquet (folded, FLOAT, units with >= 1 read in >= 3
#           samples), draws_unfolded_involved/<run>.parquet (unfolded draws of
#           the transcripts that folding merged, for the sensitivity
#           analysis), strata_tx, strata_gene, strata_cutpoints, qc,
#           control_units, productivity, eq_folded/<run>/aux_info/eq_classes.txt.gz
#           (Salmon text format, so that the package functions read it)
#   S0/     fold_summary, ab_correlations, ab_check, validation, run_summary,
#           stratum_sizes, session_info, log.txt

source("data-raw/uncertainty_evidence/ue_common.R")
ue <- ue_setup("10_placenta_prepare")
pre <- ue$prereg
t_start <- Sys.time()

v50 <- ue_v50_dir()
annot_dir <- Sys.getenv(
  "ISOFORMIC_UE_ANNOT",
  path.expand(
    "~/.cache/R/isoformic/gencode_v50_chr_patch_hapl_scaff_annotation_gff3"
  )
)
root <- ue_out_root()
# Test mode (code check on 8 samples, written to ISOFORMIC_UE_TEST_OUT)
test_mode <- nzchar(Sys.getenv("ISOFORMIC_UE_TEST_OUT"))
if (test_mode) {
  root <- Sys.getenv("ISOFORMIC_UE_TEST_OUT")
}
store <- file.path(root, "store")
rep_dir <- file.path(root, "S0")
dir.create(store, recursive = TRUE, showWarnings = FALSE)
dir.create(rep_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(rep_dir, "log.txt")
say <- function(...) {
  msg <- paste0(format(Sys.time(), "%H:%M:%S"), " ", paste0(...))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}
say("start; snapshot ", ue$git_sha, " sha256 ", ue$snapshot_sha256)

pe <- pre$data$preeclampsia
ctrl <- pre$data$control
if (test_mode) {
  pe <- pe[1:3]
  ctrl <- c(
    "SRR11498063",
    "SRR11498064",
    "SRR11498065",
    "SRR11498076",
    "SRR11498077"
  )
}
runs <- c(pe, ctrl)
quant_dir <- pre$data$placenta_quants
dirs <- stats::setNames(file.path(quant_dir, runs), runs)
stopifnot(file.exists(file.path(dirs, "quant.sf")))
condition <- stats::setNames(
  rep(c("preeclampsia", "control"), c(length(pe), length(ctrl))),
  runs
)

# 1. Annotation from the GTF (DuckDB) ----------------------------------------
gtf <- pre$data$gtf
con <- DBI::dbConnect(duckdb::duckdb())
DBI::dbExecute(con, sprintf("SET threads = %d", ue_threads()))
gtf_src <- sprintf(
  "read_csv('%s', delim = '\t', header = false, auto_detect = false, comment = '#', quote = '', escape = '', columns = {'c0': 'VARCHAR', 'c1': 'VARCHAR', 'c2': 'VARCHAR', 'c3': 'BIGINT', 'c4': 'BIGINT', 'c5': 'VARCHAR', 'c6': 'VARCHAR', 'c7': 'VARCHAR', 'c8': 'VARCHAR'})",
  gtf
)
attr_sql <- function(key) {
  sprintf("regexp_extract(c8, '%s \"([^\"]+)\"', 1)", key)
}
tx <- DBI::dbGetQuery(
  con,
  sprintf(
    "SELECT c0 AS seqname, c6 AS strand, %s AS transcript_id, %s AS gene_id,
     %s AS gene_name, %s AS gene_type, %s AS transcript_name,
     %s AS transcript_type,
     c8 LIKE '%%tag \"Ensembl_canonical\"%%' AS ensembl_canonical,
     c8 LIKE '%%tag \"MANE_Select\"%%' AS mane_select,
     c8 LIKE '%%tag \"basic\"%%' AS basic
   FROM %s WHERE c2 = 'transcript'",
    attr_sql("transcript_id"),
    attr_sql("gene_id"),
    attr_sql("gene_name"),
    attr_sql("gene_type"),
    attr_sql("transcript_name"),
    attr_sql("transcript_type"),
    gtf_src
  )
)
exw <- DBI::dbGetQuery(
  con,
  sprintf(
    "WITH e AS (SELECT %s AS transcript_id, c6 AS strand, c3 AS s, c4 AS e
              FROM %s WHERE c2 = 'exon')
   SELECT transcript_id, count(*) AS n_exons,
     string_agg(CAST(e - s + 1 AS VARCHAR), ',' ORDER BY CASE WHEN strand = '-' THEN -s ELSE s END) AS exon_widths
   FROM e GROUP BY transcript_id",
    attr_sql("transcript_id"),
    gtf_src
  )
)
DBI::dbDisconnect(con, shutdown = TRUE)
tx <- merge(tx, exw, by = "transcript_id", all.x = TRUE, sort = FALSE)
primary_chroms <- c(paste0("chr", c(1:22, "X", "Y", "M")))
tx$primary <- tx$seqname %in% primary_chroms
say(
  "GTF transcripts: ",
  nrow(tx),
  "; on primary chromosomes: ",
  sum(tx$primary)
)

# Salmon targets (order of quant.sf, identical for all runs; checked below)
q1 <- read_tabular(file.path(dirs[[1]], "quant.sf"))
targets <- q1$Name
n_targets <- length(targets)
stopifnot(targets %in% tx$transcript_id)
say("Salmon targets: ", n_targets)

# 3. Haplotype folding map ---------------------------------------------------
# Salmon already merged identical sequences (duplicate_clusters.tsv); a
# primary transcript removed there is represented by its retained copy.
dup <- utils::read.delim(file.path(
  v50,
  "human/ref/salmon_index_v50/duplicate_clusters.tsv"
))
dup$retained_primary <- tx$primary[match(dup$RetainedRef, tx$transcript_id)]
dup$duplicate_primary <- tx$primary[match(dup$DuplicateRef, tx$transcript_id)]
n_dup_primary_on_alt <- sum(!dup$retained_primary & dup$duplicate_primary)
say(
  "Salmon duplicate refs: ",
  nrow(dup),
  "; primary removed with a non-primary retained: ",
  n_dup_primary_on_alt
)
# The isoformic annotation has no source_transcript_id for GENCODE (checked:
# all NA) and alt-contig copies have their own transcript names, so the
# mapping is by gene_name plus exon structure (exon widths in transcript
# order): exact, then internal exons (>= 3 exons). Ambiguous matches are left
# unmapped.
tx_idx <- match(targets, tx$transcript_id)
ann <- tx[tx_idx, ]
# Salmon duplicate clusters whose retained member is on an alt contig while a
# removed member is on a primary chromosome: the retained target stands for
# the primary transcript (identical sequence), so it takes the primary
# transcript's identity (first primary member, sorted, when several).
dup_alias <- dup[
  !dup$retained_primary & dup$duplicate_primary & dup$RetainedRef %in% targets,
]
dup_alias <- dup_alias[order(dup_alias$RetainedRef, dup_alias$DuplicateRef), ]
dup_alias <- dup_alias[!duplicated(dup_alias$RetainedRef), ]
alias_of <- stats::setNames(rep(NA_character_, n_targets), targets)
alias_of[dup_alias$RetainedRef] <- dup_alias$DuplicateRef
is_alias <- !is.na(alias_of)
ann_eff <- ann
ann_eff[is_alias, ] <- tx[match(alias_of[is_alias], tx$transcript_id), ]
say("Salmon duplicate aliases (alt retained, primary removed): ", sum(is_alias))
prim <- ann_eff[ann_eff$primary, ]
np <- ann_eff[!ann_eff$primary, ]
prim_target <- targets[ann_eff$primary]
np_target <- targets[!ann_eff$primary]
key_exact <- function(d) paste(d$gene_name, d$exon_widths, sep = "|")
key_internal <- function(d) {
  w <- strsplit(d$exon_widths, ",", fixed = TRUE)
  inner <- vapply(
    w,
    function(v) {
      if (length(v) >= 3L) {
        paste(v[-c(1L, length(v))], collapse = ",")
      } else {
        NA_character_
      }
    },
    character(1)
  )
  ifelse(
    is.na(inner),
    NA_character_,
    paste(d$gene_name, lengths(w), inner, sep = "|")
  )
}
unique_match <- function(k_from, k_to) {
  dupk <- unique(k_to[duplicated(k_to)])
  m <- match(k_from, k_to)
  m[k_from %in% dupk] <- NA_integer_
  list(m = m, ambiguous = !is.na(k_from) & k_from %in% dupk)
}
ke <- unique_match(key_exact(np), key_exact(prim))
ki_to <- key_internal(prim)
ki <- unique_match(key_internal(np), ki_to)
ki$m[is.na(key_internal(np))] <- NA_integer_
tier <- ifelse(
  !is.na(ke$m),
  "exact_structure",
  ifelse(
    !is.na(ki$m),
    "internal_exons",
    ifelse(ke$ambiguous | ki$ambiguous, "ambiguous", "unmapped")
  )
)
target_of <- ifelse(
  !is.na(ke$m),
  prim$transcript_id[ke$m],
  ifelse(!is.na(ki$m), prim$transcript_id[ki$m], np$transcript_id)
)
fold_map <- data.frame(
  transcript_id = c(prim_target, np_target),
  folded_id = c(prim$transcript_id, target_of),
  seqname = c(ann$seqname[ann_eff$primary], ann$seqname[!ann_eff$primary]),
  primary = c(rep(TRUE, nrow(prim)), rep(FALSE, nrow(np))),
  fold_tier = c(
    ifelse(is_alias[ann_eff$primary], "salmon_duplicate_alias", "primary"),
    tier
  ),
  stringsAsFactors = FALSE
)
fold_map <- fold_map[match(targets, fold_map$transcript_id), ]
stopifnot(identical(fold_map$transcript_id, targets))
fold_map$same_type <- ann$transcript_type ==
  tx$transcript_type[match(fold_map$folded_id, tx$transcript_id)]
units <- unique(fold_map$folded_id)
fold_idx <- match(fold_map$folded_id, units)
grp_size <- tabulate(fold_idx, length(units))
fold_map$involved <- grp_size[fold_idx] > 1L
# Target that represents each unit (its own row, or the aliased retained copy)
rep_target <- match(
  units,
  ifelse(
    fold_map$fold_tier %in%
      c("primary", "salmon_duplicate_alias", "ambiguous", "unmapped"),
    fold_map$folded_id,
    NA_character_
  )
)
stopifnot(!anyNA(rep_target))
say(
  "non-primary targets: ",
  nrow(np),
  "; tiers: ",
  paste(names(table(tier)), table(tier), collapse = ", ")
)
say("folded units: ", length(units))

# 4. Point counts and QC -----------------------------------------------------
read_meta <- function(d) {
  jsonlite::read_json(file.path(d, "aux_info", "meta_info.json"))
}
counts_u <- matrix(0, n_targets, length(runs), dimnames = list(targets, runs))
efflen_u <- counts_u
tpm_u <- counts_u
qc <- vector("list", length(runs))
for (i in seq_along(runs)) {
  q <- if (i == 1L) q1 else read_tabular(file.path(dirs[[i]], "quant.sf"))
  stopifnot(identical(q$Name, targets))
  counts_u[, i] <- q$NumReads
  efflen_u[, i] <- q$EffectiveLength
  tpm_u[, i] <- q$TPM
  m <- read_meta(dirs[[i]])
  lib <- jsonlite::read_json(file.path(dirs[[i]], "lib_format_counts.json"))
  qc[[i]] <- data.frame(
    sample_id = runs[[i]],
    condition = condition[[i]],
    salmon_version = m$salmon_version,
    samp_type = m$samp_type,
    num_draws = m$num_bootstraps,
    opt_type = m$opt_type,
    em_converged = isTRUE(m$em_converged),
    num_processed = m$num_processed,
    num_mapped = m$num_mapped,
    percent_mapped = m$percent_mapped,
    num_decoy_fragments = m$num_decoy_fragments %||% NA_real_,
    library_type = lib$expected_format,
    detected_library_type = m$detected_library_type %||% NA_character_,
    compatible_fragment_ratio = lib$compatible_fragment_ratio,
    strand_mapping_bias = lib$strand_mapping_bias,
    frag_length_mean = m$frag_length_mean,
    num_eq_classes = m$num_eq_classes,
    index_seq_hash = m$index_seq_hash,
    sum_num_reads = sum(q$NumReads),
    stringsAsFactors = FALSE
  )
}
qc <- do.call(rbind, qc)
stopifnot(length(unique(qc$index_seq_hash)) == 1L, qc$num_draws == 30)
counts_f <- rowsum(counts_u, fold_idx, reorder = TRUE)
rownames(counts_f) <- units
# Effective length of a unit: that of its primary transcript (own one when unmapped)
efflen_f <- efflen_u[rep_target, , drop = FALSE]
rownames(efflen_f) <- units
lib_size <- colSums(counts_f)
cpm_f <- t(t(counts_f) / lib_size) * 1e6
keep_draws <- rowSums(counts_f >= 1) >= 3L
keep_ids <- units[keep_draws]
say("units with >= 1 read in >= 3 samples (draws stored): ", length(keep_ids))

# 2. Annotation of the units --------------------------------------------------
ua <- tx[
  match(units, tx$transcript_id),
  c(
    "transcript_id",
    "gene_id",
    "gene_name",
    "gene_type",
    "transcript_name",
    "transcript_type",
    "seqname",
    "strand",
    "primary",
    "ensembl_canonical",
    "mane_select",
    "basic",
    "n_exons"
  )
]
ua$n_folded <- grp_size
ua$on_primary_chrom <- ua$seqname %in% primary_chroms
iso_per_gene <- table(ua$gene_id[ua$on_primary_chrom])
ua$n_isoforms_gene <- as.integer(iso_per_gene[ua$gene_id])
ua$n_isoforms_gene[!ua$on_primary_chrom] <- NA_integer_
ua$isoform_bin <- cut(
  ua$n_isoforms_gene,
  c(0, 1, 4, 9, Inf),
  labels = c("1", "2-4", "5-9", "10+")
)
ua$isoform_bin <- as.character(ua$isoform_bin)

# Draw loop: moments (unfolded, folded, folded first 15) and draw tables ------
set.seed(pre$seeds$validation_subsets)
val_pool <- which(rowSums(counts_u) > 0)
val_rows <- sort(sample(val_pool, 20000L))
val_draws <- array(
  0,
  c(length(val_rows), length(runs), 30L),
  dimnames = list(targets[val_rows], runs, NULL)
)
mom_u <- list(mean = counts_u * 0, var = counts_u * 0, n = 30L)
mom_f <- list(mean = counts_f * 0, var = counts_f * 0, n = 30L)
mom_f15 <- list(mean = counts_f * 0, var = counts_f * 0, n = 15L)
keep_row <- match(keep_ids, units)
involved_rows <- which(fold_map$involved & keep_draws[fold_idx])
dir.create(file.path(store, "draws"), showWarnings = FALSE)
dir.create(file.path(store, "draws_unfolded_involved"), showWarnings = FALSE)
for (i in seq_along(runs)) {
  t0 <- Sys.time()
  x <- read_bootstrap_dir(dirs[[i]])
  stopifnot(identical(rownames(x$counts), targets), ncol(x$counts) == 30L)
  mu <- replicate_moments(x$counts, 30L)
  mom_u$mean[, i] <- mu$mean
  mom_u$var[, i] <- mu$var
  val_draws[, i, ] <- x$counts[val_rows, ]
  dr <- x$counts[involved_rows, , drop = FALSE]
  ue_write_float(
    data.frame(
      transcript_id = rep(targets[involved_rows], 30L),
      sample_id = runs[[i]],
      draw = rep(seq_len(30L), each = length(involved_rows)),
      count = as.vector(dr)
    ),
    file.path(store, "draws_unfolded_involved", paste0(runs[[i]], ".parquet")),
    "count",
    ue
  )
  f <- rowsum(x$counts, fold_idx, reorder = TRUE)
  rm(x, dr)
  mf <- replicate_moments(f, 30L)
  mom_f$mean[, i] <- mf$mean
  mom_f$var[, i] <- mf$var
  mf <- replicate_moments(f, 15L)
  mom_f15$mean[, i] <- mf$mean
  mom_f15$var[, i] <- mf$var
  fk <- f[keep_row, , drop = FALSE]
  ue_write_float(
    data.frame(
      transcript_id = rep(keep_ids, 30L),
      sample_id = runs[[i]],
      draw = rep(seq_len(30L), each = length(keep_ids)),
      count = as.vector(fk)
    ),
    file.path(store, "draws", paste0(runs[[i]], ".parquet")),
    "count",
    ue
  )
  rm(f, fk, mf)
  invisible(gc())
  say(
    "draws ",
    runs[[i]],
    " ",
    round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1),
    " s; peak RSS ",
    round(ue_peak_rss_gb(), 2),
    " GB"
  )
}
t_import <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
rss_import <- ue_peak_rss_gb()
say(
  "import done in ",
  round(t_import, 1),
  " min; peak RSS ",
  round(rss_import, 2),
  " GB"
)

# Point counts and moments, long tables
long <- function(m, value) {
  data.frame(
    transcript_id = rep(rownames(m), ncol(m)),
    sample_id = rep(colnames(m), each = nrow(m)),
    value = as.vector(m),
    stringsAsFactors = FALSE
  ) |>
    stats::setNames(c("transcript_id", "sample_id", value))
}
pc <- long(counts_f, "count")
pc$effective_length <- as.vector(efflen_f)
ue_write(pc, file.path(store, "point_counts.parquet"), ue)
pc <- long(counts_u, "count")
pc$effective_length <- as.vector(efflen_u)
pc$tpm <- as.vector(tpm_u)
ue_write(pc, file.path(store, "point_counts_unfolded.parquet"), ue)
mt <- long(mom_f$mean, "mean")
mt$var <- as.vector(mom_f$var)
mt$mean_15 <- as.vector(mom_f15$mean)
mt$var_15 <- as.vector(mom_f15$var)
mt$n_draws <- 30L
ue_write(mt, file.path(store, "moments.parquet"), ue)
mt <- long(mom_u$mean, "mean")
mt$var <- as.vector(mom_u$var)
mt$n_draws <- 30L
ue_write(mt, file.path(store, "moments_unfolded.parquet"), ue)
rm(pc, mt)
invisible(gc())

# Folding summary
fs_tab <- as.data.frame(
  table(fold_tier = fold_map$fold_tier),
  stringsAsFactors = FALSE
)
names(fs_tab)[2] <- "n_targets"
fs_tab$n_with_reads <- vapply(
  fs_tab$fold_tier,
  function(t) {
    sum(rowSums(counts_u[fold_map$fold_tier == t, , drop = FALSE]) > 0)
  },
  numeric(1)
)
fs_tab$reads_total <- vapply(
  fs_tab$fold_tier,
  function(t) sum(counts_u[fold_map$fold_tier == t, , drop = FALSE]),
  numeric(1)
)
fs_tab <- rbind(
  fs_tab,
  data.frame(
    fold_tier = c(
      "salmon_duplicates_total",
      "salmon_duplicates_primary_removed_on_nonprimary"
    ),
    n_targets = c(nrow(dup), n_dup_primary_on_alt),
    n_with_reads = NA,
    reads_total = NA
  )
)
ue_write(fs_tab, file.path(rep_dir, "fold_summary.parquet"), ue)
ue_write(fold_map, file.path(store, "fold_map.parquet"), ue)

# 4b. Fetal sex ---------------------------------------------------------------
gene_counts <- rowsum(
  counts_f[ua$on_primary_chrom, ],
  ua$gene_name[ua$on_primary_chrom]
)
gene_cpm <- t(t(gene_counts) / lib_size) * 1e6
y_genes <- c("RPS4Y1", "DDX3Y", "KDM5D", "UTY")
stopifnot(c("XIST", y_genes) %in% rownames(gene_cpm))
qc$xist_log2cpm <- log2(gene_cpm["XIST", ] + 1)
qc$y_log2cpm <- log2(colSums(gene_cpm[y_genes, ]) + 1)
qc$sex_score <- qc$y_log2cpm - qc$xist_log2cpm
set.seed(pre$seeds$sex_kmeans)
km <- stats::kmeans(qc$sex_score, centers = 2L, nstart = 25L)
male_cluster <- which.max(km$centers[, 1])
qc$fetal_sex <- ifelse(km$cluster == male_cluster, "male", "female")
qc$sex_gap <- min(qc$sex_score[qc$fetal_sex == "male"]) -
  max(qc$sex_score[qc$fetal_sex == "female"])
say(
  "fetal sex: ",
  paste(names(table(qc$fetal_sex)), table(qc$fetal_sex), collapse = ", "),
  "; gap ",
  round(qc$sex_gap[1], 2)
)

# 5. A/B technical pairs -------------------------------------------------------
cl <- ctrl
expr_tx <- ua$on_primary_chrom & rowMeans(cpm_f[, cl]) >= 1
lc_tx <- log2(cpm_f[expr_tx, cl] + 0.5)
expr_g <- rowMeans(gene_cpm[, cl]) >= 1
lc_g <- log2(gene_cpm[expr_g, cl] + 0.5)
pairs <- utils::combn(cl, 2)
ab <- lapply(pre$data$technical_pairs, sort)
ab_key <- vapply(ab, paste, character(1), collapse = "_")
pair_tab <- data.frame(
  a = pairs[1, ],
  b = pairs[2, ],
  rho_tx = vapply(
    seq_len(ncol(pairs)),
    function(k) {
      stats::cor(
        lc_tx[, pairs[1, k]],
        lc_tx[, pairs[2, k]],
        method = "spearman"
      )
    },
    numeric(1)
  ),
  rho_gene = vapply(
    seq_len(ncol(pairs)),
    function(k) {
      stats::cor(lc_g[, pairs[1, k]], lc_g[, pairs[2, k]], method = "spearman")
    },
    numeric(1)
  )
)
pair_tab$key <- paste(
  pmin(pair_tab$a, pair_tab$b),
  pmax(pair_tab$a, pair_tab$b),
  sep = "_"
)
pair_tab$is_ab <- pair_tab$key %in% ab_key
other <- pair_tab[!pair_tab$is_ab, ]
ab_check <- do.call(
  rbind,
  lapply(names(ab), function(nm) {
    r <- pair_tab[pair_tab$key == ab_key[[nm]], ]
    depth <- qc$num_mapped[match(ab[[nm]], qc$sample_id)]
    data.frame(
      pair = nm,
      a = ab[[nm]][1],
      b = ab[[nm]][2],
      rho_tx = r$rho_tx,
      rho_gene = r$rho_gene,
      pct_rank_tx = 100 * mean(other$rho_tx < r$rho_tx),
      pct_rank_gene = 100 * mean(other$rho_gene < r$rho_gene),
      other_tx_median = stats::median(other$rho_tx),
      other_tx_max = max(other$rho_tx),
      other_gene_median = stats::median(other$rho_gene),
      other_gene_max = max(other$rho_gene),
      n_other_pairs = nrow(other),
      n_tx = nrow(lc_tx),
      n_genes = nrow(lc_g),
      mapped_a = depth[1],
      mapped_b = depth[2],
      kept_primary = ab[[nm]][which.max(depth)],
      kept_sensitivity = ab[[nm]][which.min(depth)]
    )
  })
)
ue_write(pair_tab, file.path(rep_dir, "ab_correlations.parquet"), ue)
ue_write(ab_check, file.path(rep_dir, "ab_check.parquet"), ue)
units_primary <- setdiff(cl, ab_check$kept_sensitivity)
units_sens <- setdiff(cl, ab_check$kept_primary)
ue_write(
  data.frame(
    sample_id = cl,
    unit_primary = cl %in% units_primary,
    unit_sensitivity = cl %in% units_sens
  ),
  file.path(store, "control_units.parquet"),
  ue
)
qc$control_unit_primary <- qc$sample_id %in% units_primary
qc$control_unit_sensitivity <- qc$sample_id %in% units_sens
ue_write(qc, file.path(store, "qc.parquet"), ue)
say(
  "A/B: ",
  paste(
    ab_check$pair,
    round(ab_check$rho_tx, 4),
    round(ab_check$pct_rank_tx, 1),
    collapse = "; "
  )
)

# Folded equivalence classes (Salmon text format) ------------------------------
write_folded_eq <- function(run) {
  ec <- read_eq_classes(dirs[[run]])
  fid <- match(
    fold_map$folded_id[match(ec$transcripts, fold_map$transcript_id)],
    units
  )
  dt <- data.table::data.table(
    cls = ec$members$class,
    f = fid[ec$members$tid],
    count = ec$members$count
  )
  dt <- dt[!is.na(f)]
  dt <- unique(dt, by = c("cls", "f"))
  data.table::setorder(dt, cls, f)
  cl_tab <- dt[,
    list(key = paste(f - 1L, collapse = "\t"), size = .N, count = count[1L]),
    by = cls
  ]
  cl_tab <- cl_tab[, list(size = size[1L], count = sum(count)), by = key]
  path <- file.path(store, "eq_folded", run, "aux_info")
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  gz <- gzfile(file.path(path, "eq_classes.txt.gz"), "w")
  writeLines(
    c(
      as.character(length(units)),
      as.character(nrow(cl_tab)),
      units,
      paste(
        cl_tab$size,
        cl_tab$key,
        format(cl_tab$count, scientific = FALSE, trim = TRUE),
        sep = "\t"
      )
    ),
    gz
  )
  close(gz)
  file.path(store, "eq_folded", run)
}
eq_dirs <- vapply(runs, write_folded_eq, character(1))
say("folded equivalence classes written")

# 6. Universe U and strata from the PE samples ---------------------------------
cpm_pe <- cpm_f[, pe]
in_u <- ua$on_primary_chrom &
  rowSums(cpm_pe >= pre$universe$cpm_floor) >= pre$universe$min_samples
u_ids <- units[in_u]
say("|U| = ", length(u_ids))
rv <- ue_infrv(mom_f$mean[, pe], mom_f$var[, pe])
rv15 <- ue_infrv(mom_f15$mean[, pe], mom_f15$var[, pe])
mean_rv <- rowMeans(rv)
mean_rv15 <- rowMeans(rv15)
cuts <- stats::quantile(mean_rv[in_u], c(1, 2) / 3, names = FALSE)
cuts15 <- stats::quantile(mean_rv15[in_u], c(1, 2) / 3, names = FALSE)
tert <- function(x, k) {
  as.character(cut(x, c(-Inf, k, Inf), labels = c("low", "mid", "high")))
}
t_supp <- Sys.time()
supp <- equivalence_support(eq_dirs[pe])
supp_u <- equivalence_support(dirs[pe])
say(
  "equivalence support ",
  round(as.numeric(difftime(Sys.time(), t_supp, units = "mins")), 1),
  " min"
)
# Identifiability one sample at a time, checkpointed (the slow step): the
# per-sample tables are summed exactly as equivalence_identifiability() does
# over several samples.
t_id <- Sys.time()
id_dir <- file.path(store, "identifiability_pe")
dir.create(id_dir, showWarnings = FALSE)
id_tabs <- lapply(pe, function(run) {
  f <- file.path(id_dir, paste0(run, ".parquet"))
  if (!file.exists(f)) {
    t0 <- Sys.time()
    one <- equivalence_identifiability(eq_dirs[run])
    ue_write(
      as.data.frame(one)[, c("transcript_id", "n_samples", "n_identifiable")],
      f,
      ue
    )
    say(
      "identifiability ",
      run,
      " ",
      round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1),
      " min; peak RSS ",
      round(ue_peak_rss_gb(), 2),
      " GB"
    )
  }
  read_parquet(f)[, c("transcript_id", "n_samples", "n_identifiable")]
})
id_all <- do.call(rbind, id_tabs)
ident <- data.frame(
  transcript_id = names(tapply(id_all$n_samples, id_all$transcript_id, sum)),
  n_samples = as.integer(tapply(id_all$n_samples, id_all$transcript_id, sum)),
  n_identifiable = as.integer(tapply(
    id_all$n_identifiable,
    id_all$transcript_id,
    sum
  )),
  stringsAsFactors = FALSE
)
ident$identifiable_fraction <- ident$n_identifiable / ident$n_samples
ident$identifiable <- ident$n_identifiable == ident$n_samples
rm(id_all, id_tabs)
t_id <- as.numeric(difftime(Sys.time(), t_id, units = "mins"))
say("identifiability on 8 PE samples ", round(t_id, 1), " min")
mean_cpm_pe <- rowMeans(cpm_pe)
dec_cuts <- stats::quantile(
  mean_cpm_pe[in_u],
  seq(0.1, 0.9, 0.1),
  names = FALSE
)
strata <- data.frame(
  ua,
  in_u = in_u,
  draws_stored = keep_draws,
  mean_cpm_pe = mean_cpm_pe,
  mean_count_pe = rowMeans(counts_f[, pe]),
  mean_infrv_pe = mean_rv,
  infrv_tertile = ifelse(in_u, tert(mean_rv, cuts), NA_character_),
  uncertain = mean_rv >= pre$strata$infrv_uncertain_cutoff,
  mean_infrv_pe_15 = mean_rv15,
  infrv_tertile_15 = ifelse(in_u, tert(mean_rv15, cuts15), NA_character_),
  uncertain_15 = mean_rv15 >= pre$strata$infrv_uncertain_cutoff,
  stringsAsFactors = FALSE
)
m <- match(units, supp$transcript_id)
strata$unique_fraction <- supp$unique_fraction[m]
strata$n_fragments_pe <- supp$n_fragments[m]
strata$zero_unique <- strata$unique_fraction == 0
strata$unique_fraction_unfolded <- supp_u$unique_fraction[match(
  units,
  supp_u$transcript_id
)]
m <- match(units, ident$transcript_id)
strata$identifiable_n_samples <- ident$n_samples[m]
strata$identifiable_fraction <- ident$identifiable_fraction[m]
strata$non_identifiable <- !ident$identifiable[m]
strata$expr_decile <- ifelse(
  in_u,
  findInterval(mean_cpm_pe, dec_cuts) + 1L,
  NA_integer_
)
strata$productivity <- NA_character_

ue_write(strata, file.path(store, "strata_tx.parquet"), ue)
ue_write(
  strata[, c(
    "transcript_id",
    "gene_id",
    "gene_name",
    "transcript_name",
    "transcript_type",
    "seqname",
    "primary",
    "on_primary_chrom",
    "n_isoforms_gene",
    "n_folded"
  )],
  file.path(store, "tx2gene.parquet"),
  ue
)
ue_write(
  data.frame(
    what = c(
      "infrv_tertile",
      "infrv_tertile",
      "infrv_tertile_15",
      "infrv_tertile_15",
      paste0("expr_decile_", 1:9)
    ),
    cut = c(
      "low|mid",
      "mid|high",
      "low|mid",
      "mid|high",
      paste0(1:9, "|", 2:10)
    ),
    value = c(cuts, cuts15, dec_cuts)
  ),
  file.path(store, "strata_cutpoints.parquet"),
  ue
)

# Gene-level strata on U
su <- strata[in_u, ]
gl <- data.frame(
  gene_id = sort(unique(su$gene_id)),
  stringsAsFactors = FALSE
)
gl$gene_name <- su$gene_name[match(gl$gene_id, su$gene_id)]
gl$n_isoforms_u <- as.integer(table(su$gene_id)[gl$gene_id])
gl$max_infrv_pe <- tapply(su$mean_infrv_pe, su$gene_id, max)[gl$gene_id]
gl$any_non_identifiable <- tapply(su$non_identifiable, su$gene_id, function(v) {
  any(v, na.rm = TRUE)
})[gl$gene_id]
gl$any_uncertain <- tapply(su$uncertain, su$gene_id, any)[gl$gene_id]
gl$gene_cpm_pe <- tapply(su$mean_cpm_pe, su$gene_id, sum)[gl$gene_id]
gl$max_infrv_tertile <- tert(
  gl$max_infrv_pe,
  stats::quantile(gl$max_infrv_pe, c(1, 2) / 3, names = FALSE)
)
ue_write(gl, file.path(store, "strata_gene.parquet"), ue)

# Stratum sizes (on U)
count_by <- function(col) {
  tab <- table(su[[col]], useNA = "ifany")
  data.frame(
    stratum = col,
    level = as.character(names(tab)),
    n = as.integer(tab)
  )
}
sizes <- rbind(
  data.frame(stratum = "U", level = "all", n = nrow(su)),
  count_by("infrv_tertile"),
  count_by("uncertain"),
  count_by("infrv_tertile_15"),
  count_by("uncertain_15"),
  count_by("non_identifiable"),
  count_by("zero_unique"),
  count_by("isoform_bin"),
  count_by("expr_decile"),
  data.frame(stratum = "genes_U", level = "all", n = nrow(gl)),
  data.frame(
    stratum = "genes_U_any_non_identifiable",
    level = c("TRUE", "FALSE", "NA"),
    n = c(
      sum(gl$any_non_identifiable %in% TRUE),
      sum(gl$any_non_identifiable %in% FALSE),
      sum(is.na(gl$any_non_identifiable))
    )
  ),
  data.frame(
    stratum = "genes_U_max_infrv_tertile",
    level = names(table(gl$max_infrv_tertile)),
    n = as.integer(table(gl$max_infrv_tertile))
  )
)
ue_write(sizes, file.path(rep_dir, "stratum_sizes.parquet"), ue)
say("strata written")

# 7. Validation of the matrix code path ----------------------------------------
set.seed(pre$seeds$validation_subsets + 1L)
subsets <- list(
  all25 = runs,
  pe8 = pe,
  control_units15 = units_primary,
  random6 = sort(sample(runs, min(6L, length(runs))))
)
rel <- function(a, b) {
  ok <- is.finite(a) & is.finite(b)
  max(abs(a[ok] - b[ok]) / pmax(abs(b[ok]), .Machine$double.eps))
}
val <- list()
for (nm in names(subsets)) {
  s <- subsets[[nm]]
  t0 <- Sys.time()
  cs <- edgeR::catchSalmon(unname(dirs[s]), verbose = FALSE)
  t_cs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  ours <- ue_overdispersion(
    mom_u$mean[, s, drop = FALSE],
    mom_u$var[, s, drop = FALSE],
    30L
  )
  stopifnot(identical(rownames(cs$annotation), ours$transcript_id))
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = stats::setNames(
      lapply(seq_len(30L), function(k) val_draws[, s, k]),
      paste0("infRep", seq_len(30L))
    )
  )
  se <- fishpond::computeInfRV(se, meanVariance = FALSE)
  rv_ours <- rowMeans(ue_infrv(
    mom_u$mean[val_rows, s, drop = FALSE],
    mom_u$var[val_rows, s, drop = FALSE]
  ))
  val[[nm]] <- data.frame(
    subset = nm,
    n_samples = length(s),
    samples = paste(s, collapse = ","),
    n_transcripts_overdispersion = nrow(ours),
    n_with_df = sum(ours$df > 0),
    max_rel_diff_overdispersion = rel(
      ours$overdispersion,
      cs$annotation$Overdispersion
    ),
    prior_ours = attr(ours, "prior"),
    prior_catchSalmon = cs$overdispersion.prior,
    rel_diff_prior = abs(attr(ours, "prior") - cs$overdispersion.prior) /
      cs$overdispersion.prior,
    max_rel_diff_counts = rel(counts_u[, s], cs$counts),
    n_transcripts_infrv = length(val_rows),
    max_rel_diff_infrv = rel(
      rv_ours,
      SummarizedExperiment::rowData(se)$meanInfRV
    ),
    catchSalmon_seconds = t_cs
  )
  say(
    "validation ",
    nm,
    ": overdispersion ",
    signif(val[[nm]]$max_rel_diff_overdispersion, 3),
    ", InfRV ",
    signif(val[[nm]]$max_rel_diff_infrv, 3)
  )
  rm(cs, se)
  invisible(gc())
}
val <- do.call(rbind, val)
ue_write(val, file.path(rep_dir, "validation.parquet"), ue)
# Folded overdispersion of the PE samples, for the record
od_pe <- ue_overdispersion(mom_f$mean[, pe], mom_f$var[, pe], 30L)

# Productivity classes (last: optional, a failure leaves the rest intact) -----
# Units on primary chromosomes with draws stored (>= 1 read in >= 3 samples).
t_prod <- Sys.time()
prod_ok <- tryCatch(
  {
    prod_ids <- units[ua$on_primary_chrom & keep_draws]
    iso <- IsoformicExperiment(
      annot_path = annot_dir,
      assays = list(counts = counts_f[prod_ids, , drop = FALSE])
    )
    prod <- transcript_productivity(iso)
    prod <- as.data.frame(prod)
    prod <- prod[, intersect(
      c(
        "transcript_id",
        "productivity",
        "reason",
        "source",
        "nmd_predicted",
        "coding_probability"
      ),
      colnames(prod)
    )]
    prod$productivity <- as.character(prod$productivity)
    ue_write(prod, file.path(store, "productivity.parquet"), ue)
    strata$productivity <- prod$productivity[match(units, prod$transcript_id)]
    ue_write(strata, file.path(store, "strata_tx.parquet"), ue)
    pt <- table(strata$productivity[strata$in_u], useNA = "ifany")
    ue_write(
      rbind(
        sizes,
        data.frame(
          stratum = "productivity",
          level = as.character(names(pt)),
          n = as.integer(pt)
        )
      ),
      file.path(rep_dir, "stratum_sizes.parquet"),
      ue
    )
    rm(iso)
    TRUE
  },
  error = function(e) {
    say("productivity FAILED: ", conditionMessage(e))
    FALSE
  }
)
t_prod <- as.numeric(difftime(Sys.time(), t_prod, units = "mins"))
say("productivity ", round(t_prod, 1), " min")

# Run summary -------------------------------------------------------------------
t_total <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
summary_tab <- data.frame(
  metric = c(
    "wall_minutes_total",
    "wall_minutes_import",
    "peak_rss_gb",
    "salmon_targets",
    "folded_units",
    "draw_units_stored",
    "U_size",
    "identifiability_minutes",
    "prior_overdispersion_folded_pe8",
    "peak_rss_gb_import",
    "productivity_minutes",
    "productivity_ok",
    "salmon_duplicate_aliases"
  ),
  value = c(
    t_total,
    t_import,
    ue_peak_rss_gb(),
    n_targets,
    length(units),
    length(keep_ids),
    length(u_ids),
    t_id,
    attr(od_pe, "prior"),
    rss_import,
    t_prod,
    as.numeric(prod_ok),
    sum(is_alias)
  )
)
ue_write(summary_tab, file.path(rep_dir, "run_summary.parquet"), ue)
ue_write_session(rep_dir, ue)
say(
  "done: ",
  round(t_total, 1),
  " min; peak RSS ",
  round(ue_peak_rss_gb(), 2),
  " GB"
)
