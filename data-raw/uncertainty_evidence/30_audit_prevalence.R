# WP-B (placenta audit), part 1: prevalence of fragile conclusions in the real
# preeclampsia (PE) against control comparison.
#
# What is measured: the sensitivity of conclusions to quantification
# uncertainty (self-consistency), not an error rate. Every naive conclusion
# of the 8 PE against the 15 independent control placentas (on the universe
# U of the S0 store, haplotype and patch copies folded) is refit on each of
# the 30 Gibbs draws (full refit: normalisation and dispersions re-estimated
# on the same features), and compared with the uncertainty-aware
# alternatives (divided counts, swish, scaled usage test). Each figure is set
# next to the same statistic for gene-level tests matched on the naive
# p-value, so that borderline p-values are not read as fragility.
#
# Families
#   F1 DTE        edgeR QL on transcript counts (filterByExpr); alternatives:
#                 counts divided by the RTA overdispersion (Baldoni 2024,
#                 package estimator on all folded units, equal to
#                 edgeR::catchSalmon), fishpond swish (seed 20261113)
#   F2 switches   isoform_switches() rule: a gene with a transcript up and one
#                 down at q < 0.05 and |log2FC| >= 1
#   F3 DTU        diffSplice backend of test_usage(), scale FALSE / TRUE
#   F4 NMD        test_productive_fraction() (min_gene_count 10) on counts,
#                 raw / divided, and delta_z of productive_fraction() over
#                 the draws (assay counts; repository code, amendment A18)
#   F5 enrichment run_enrichment() ranked by log2FC against lfc_z
#                 (secondary; Reactome C2 GMT of the preeclampsia article
#                 data, vignettes/articles/data/preeclampsia/)
#   G  gene level edgeR QL on gene counts (sum of the isoforms in U): the
#                 p-matched reference
# Confidence classes: add_uncertainty() + de_confidence() of the package on
# the whole of U (stage `package`), checked against the matrix code of
# 38_audit_common.R that 31_halfsplit_replication.R uses (stage `validate`).
#
# Stages (ISOFORMIC_UE_STAGE, default "all"; each writes its outputs and is
# skipped when they exist):
#   naive     naive and alternative fits (F1, F3, F4, G, divided, sex covariate)
#   swish     fishpond swish on transcripts and gene sums; F2 switch genes
#   package   package objects on U: add_uncertainty/de_confidence,
#             isoform_switches, test_usage, F5
#   productivity  F4 package objects on U: test_productive_fraction and the
#             delta_z of productive_fraction(). Amendment A18: run with
#             ISOFORMIC_UE_PKG=. (the repository code, not the snapshot),
#             because productive_fraction_sd and delta_z changed definition
#             after the snapshot
#   validate  matrix code against the package (validation.parquet)
#   perdraw   full refits of F1, F3, F4 and G on each of the 30 draws
#   null      20 null splits of the control units (7 v 8, sex-balanced)
#   summarise 30b_prevalence_summary.R (MUST REPORT (1), strata, null link)
#
# Inputs: prereg.yaml and the package snapshot (01_preregister.R); the S0
#   store (10_placenta_prepare.R); the GENCODE 50 annotation directory.
# Env vars:
#   ISOFORMIC_UE_OUT, ISOFORMIC_UE_PKG, ISOFORMIC_UE_THREADS  (ue_common.R)
#   ISOFORMIC_UE_STAGE    naive | swish | package | productivity | validate |
#                         perdraw | null | summarise | all
#   ISOFORMIC_UE_WORKERS  forked workers for the per-draw refits (default 3)
#   ISOFORMIC_UE_ANNOT    annotation directory, default
#                         ~/.cache/R/isoformic/gencode_v50_chr_patch_hapl_scaff_annotation_gff3
#   ISOFORMIC_UE_PILOT    "true": 1,500 random genes, output to
#                         placenta_audit_pilot/ (code check)
# Outputs: $ISOFORMIC_UE_OUT/placenta_audit/prevalence/*.parquet
# Memory (estimate): draws of U 57,294 x 23 x 30 x 8 B = 0.32 GB, moments of
#   650,973 units x 23 x 2 x 8 B = 0.24 GB, the package object on U with its
#   draws ~2-3 GB; per-draw workers ~0.5-1 GB each: < 8 GB with 3 workers.
# Measured per stage (peak RSS of the main process): naive 1 min, 2.2 GB;
#   swish 1 min, 3.7 GB; package 3 min, 4.4 GB (then with productivity);
#   productivity (repository code) 1 min, 3.9 GB; validate 1 min, 2.6 GB;
#   perdraw 3 min, 1.5 GB; null 2 min, 2.6 GB; summarise 2 min, 1.7 GB.
# Run from the repository root, inside a memory cap of about 8 GB (see
# README.md):
#   OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
#     pixi run --frozen Rscript data-raw/uncertainty_evidence/30_audit_prevalence.R
# (one-thread BLAS: forked workers x multithreaded OpenBLAS would use more
# than 4 threads; the evaluation ran each stage with ISOFORMIC_UE_STAGE=<stage>;
# the productivity stage with ISOFORMIC_UE_PKG=. as well, then validate and
# summarise again with the snapshot)

source("data-raw/uncertainty_evidence/ue_common.R")
ue <- ue_setup("30_audit_prevalence")
source("data-raw/uncertainty_evidence/38_audit_common.R")
stage <- Sys.getenv("ISOFORMIC_UE_STAGE", "all")
workers <- as.integer(Sys.getenv("ISOFORMIC_UE_WORKERS", "3"))
annot_dir <- Sys.getenv(
  "ISOFORMIC_UE_ANNOT",
  path.expand(
    "~/.cache/R/isoformic/gencode_v50_chr_patch_hapl_scaff_annotation_gff3"
  )
)
prev_dir <- wpb_dir("prevalence")
dir.create(prev_dir, recursive = TRUE, showWarnings = FALSE)
say <- wpb_logger(wpb_dir("logs", "30_audit_prevalence.log"))
run_stage <- function(s) stage %in% c(s, "all")
pq <- function(name) file.path(prev_dir, paste0(name, ".parquet"))
done <- function(name) file.exists(pq(name))
say(
  "stage %s, snapshot %s (sha256 %s), pilot %s",
  stage,
  ue$git_sha,
  ue$snapshot_sha256,
  wpb_pilot()
)

need_mom <- run_stage("naive") && !done("overdispersion_U")
st <- wpb_load_store(draws = TRUE, moments = need_mom, say = say)
cohort <- st$runs
groups <- st$groups
sex <- st$sex
U <- st$U
gene_of <- st$gene_of
counts_U <- st$counts
lib_full <- st$lib_size
arr <- st$arr
n_draw <- dim(arr)[3L]
tx <- st$tx
is_prod <- tx$productivity %in% c("productive", "altered_productive")
coding_genes <- unique(tx$gene_id[
  !is.na(tx$productivity) & tx$productivity != "noncoding"
])
in_b <- groups == "preeclampsia"
say(
  "cohort: %s; sex x group: %s",
  paste(table(groups), collapse = " control / PE "),
  paste(apply(table(sex, groups), 1, paste, collapse = "/"), collapse = "; ")
)

# Stage naive --------------------------------------------------------------------
if (run_stage("naive") && !done("naive_f4_div")) {
  seed <- wpb_seed(ue, 1L)
  set.seed(seed)
  if (!done("overdispersion_U")) {
    od <- wpb_overdispersion(st$mom, cohort)
    say(
      "overdispersion prior (23 cohort samples, all folded units): %.4f",
      attr(od, "prior")
    )
    ue_write(
      data.frame(
        transcript_id = U,
        overdispersion = unname(od[U]),
        prior = attr(od, "prior")
      ),
      pq("overdispersion_U"),
      ue,
      seed
    )
  }
  st$mom <- NULL
  invisible(gc())
  odt <- wpb_read(pq("overdispersion_U"))
  od <- stats::setNames(odt$overdispersion, odt$transcript_id)
  # F1 and G
  f1_raw <- wpb_edger(counts_U, groups)
  div_U <- wpb_divide(counts_U, od)
  f1_div <- wpb_edger(div_U, groups)
  f1_sex <- wpb_edger(counts_U, groups, covariates = data.frame(sex = sex))
  gcounts <- wpb_gene_sums(counts_U, gene_of[U])
  g_raw <- wpb_edger(gcounts, groups)
  # jarl-ignore unused_object: written below through get(nm)
  g_sex <- wpb_edger(gcounts, groups, covariates = data.frame(sex = sex))
  say(
    "F1: raw %d / %d tested, divided %d / %d, raw + sex %d calls at q < 0.05; G: %d / %d genes",
    sum(f1_raw$qvalue < 0.05),
    nrow(f1_raw),
    sum(f1_div$qvalue < 0.05),
    nrow(f1_div),
    sum(f1_sex$qvalue < 0.05),
    sum(g_raw$qvalue < 0.05),
    nrow(g_raw)
  )
  # F3 usage
  f3_raw <- wpb_usage(counts_U, gene_of, groups)
  f3_scaled <- wpb_usage(counts_U, gene_of, groups, od = od)
  say(
    "F3: raw %d / %d, scaled %d / %d transcripts at q < 0.05",
    sum(f3_raw$qvalue < 0.05),
    nrow(f3_raw),
    sum(f3_scaled$qvalue < 0.05),
    nrow(f3_scaled)
  )
  # F4 productive fraction (coding genes with mean gene count >= 10 per group)
  f4_genes <- wpb_prod_genes(counts_U, gene_of, groups, coding_genes, 10)
  f4_raw <- wpb_prod_test(counts_U, gene_of, is_prod, groups, f4_genes)
  f4_genes_div <- wpb_prod_genes(div_U, gene_of, groups, coding_genes, 10)
  f4_div <- wpb_prod_test(div_U, gene_of, is_prod, groups, f4_genes_div)
  say(
    "F4: raw %d / %d, divided %d / %d genes at q < 0.05",
    sum(f4_raw$qvalue < 0.05, na.rm = TRUE),
    nrow(f4_raw),
    sum(f4_div$qvalue < 0.05, na.rm = TRUE),
    nrow(f4_div)
  )
  for (nm in c(
    "f1_raw",
    "f1_div",
    "f1_sex",
    "g_raw",
    "g_sex",
    "f3_raw",
    "f3_scaled",
    "f4_raw"
  )) {
    ue_write(get(nm), pq(paste0("naive_", nm)), ue, seed)
  }
  ue_write(f4_div, pq("naive_f4_div"), ue, seed)
  say("naive done, peak RSS %.2f GB", ue_peak_rss_gb())
}

# Stage swish (seed 20261113, as specified): transcripts and gene sums ------------
if (run_stage("swish") && !done("naive_f2_switches")) {
  run_swish <- function(cnt, slice) {
    reps <- lapply(seq_len(n_draw), slice)
    names(reps) <- paste0("infRep", seq_len(n_draw))
    se <- SummarizedExperiment::SummarizedExperiment(
      assays = c(list(counts = cnt), reps),
      colData = S4Vectors::DataFrame(condition = groups, row.names = cohort)
    )
    rm(reps)
    se <- suppressMessages(fishpond::scaleInfReps(se, quiet = TRUE))
    se <- fishpond::labelKeep(se)
    se <- se[S4Vectors::mcols(se)$keep, ]
    set.seed(20261113)
    se <- suppressMessages(fishpond::swish(se, x = "condition", quiet = TRUE))
    mc <- as.data.frame(S4Vectors::mcols(se))
    out <- data.frame(
      feature_id = rownames(se),
      logFC = mc$log2FC,
      pvalue = mc$pvalue,
      qvalue = mc$qvalue,
      stat = mc$stat
    )
    rm(se)
    invisible(gc())
    out
  }
  if (!done("naive_f1_swish")) {
    f1_swish <- run_swish(counts_U, function(k) arr[,, k])
    ue_write(f1_swish, pq("naive_f1_swish"), ue, 20261113L)
    say(
      "swish transcripts: %d calls at q < 0.05 of %d kept, peak RSS %.2f GB",
      sum(f1_swish$qvalue < 0.05),
      nrow(f1_swish),
      ue_peak_rss_gb()
    )
  }
  if (!done("naive_g_swish")) {
    g_swish <- run_swish(wpb_gene_sums(counts_U, gene_of[U]), function(k) {
      wpb_gene_sums(arr[,, k], gene_of[U])
    })
    ue_write(g_swish, pq("naive_g_swish"), ue, 20261113L)
    say(
      "swish genes: %d calls at q < 0.05 of %d kept",
      sum(g_swish$qvalue < 0.05),
      nrow(g_swish)
    )
  }
  # F2 switches of raw, divided and swish DTE
  det_of <- function(f) {
    data.frame(
      transcript_id = f$feature_id,
      gene_id = gene_of[f$feature_id],
      logFC = f$logFC,
      qvalue = f$qvalue
    )
  }
  f2 <- rbind(
    data.frame(
      gene_id = wpb_switch_genes(det_of(wpb_read(pq("naive_f1_raw")))),
      method = "raw"
    ),
    data.frame(
      gene_id = wpb_switch_genes(det_of(wpb_read(pq("naive_f1_div")))),
      method = "divided"
    ),
    data.frame(
      gene_id = wpb_switch_genes(det_of(wpb_read(pq("naive_f1_swish")))),
      method = "swish"
    )
  )
  say(
    "F2 switches: %s",
    paste(names(table(f2$method)), table(f2$method), collapse = ", ")
  )
  ue_write(f2, pq("naive_f2_switches"), ue, wpb_seed(ue, 1L))
}

# The package object on U: naive F1 results and the 30 draws -----------------------
package_object <- function() {
  iso <- IsoformicExperiment(
    experiment_name = "placenta_U",
    annot_path = annot_dir,
    assays = list(counts = counts_U),
    col_data = data.frame(
      sample_id = cohort,
      condition = as.character(groups),
      library_size = unname(lib_full[cohort])
    )
  )
  meta <- iso@metadata
  meta$inf_reps <- lapply(seq_len(n_draw), function(k) arr[,, k])
  iso@metadata <- meta
  f1_raw <- wpb_read(pq("naive_f1_raw"))
  de_tx(iso) <- data.frame(
    transcript_id = f1_raw$feature_id,
    log2FC = f1_raw$logFC,
    pvalue = f1_raw$pvalue,
    qvalue = f1_raw$qvalue
  )
  iso
}

# Stage package --------------------------------------------------------------------
if (run_stage("package") && !done("package_confidence")) {
  seed <- wpb_seed(ue, 2L)
  set.seed(seed)
  iso <- package_object()
  if (stage == "package") {
    # the array is not needed again in this process
    rm(arr)
    st$arr <- NULL
    invisible(gc())
  }
  t0 <- Sys.time()
  iso <- add_uncertainty(
    iso,
    group = "condition",
    levels = c("control", "preeclampsia")
  )
  say(
    "add_uncertainty on U: %.1f min, peak RSS %.2f GB",
    as.numeric(difftime(Sys.time(), t0, units = "mins")),
    ue_peak_rss_gb()
  )
  conf <- as.data.frame(de_confidence(iso))
  conf <- conf[, !vapply(conf, is.list, logical(1)), drop = FALSE]
  ue_write(conf, pq("package_confidence"), ue, seed)
  amb <- as.data.frame(iso@dea[["ambiguity"]])
  ue_write(amb, pq("package_ambiguity_pairs"), ue, seed)
  unc <- as.data.frame(iso@dea[["uncertainty"]])
  ue_write(unc, pq("package_uncertainty"), ue, seed)
  say(
    "package classes on U: %s",
    paste(
      names(table(conf$confidence, useNA = "ifany")),
      table(conf$confidence, useNA = "ifany"),
      collapse = ", "
    )
  )
  # isoform_switches on the naive F1 results
  sw <- as.data.frame(isoform_switches(iso))
  ue_write(
    sw[, !vapply(sw, is.list, logical(1)), drop = FALSE],
    pq("package_isoform_switches"),
    ue,
    seed
  )
  # test_usage raw and scaled (scale uses the draws of the object, i.e. U)
  tu0 <- as.data.frame(test_usage(
    iso,
    group = "condition",
    levels = c("control", "preeclampsia"),
    scale = FALSE
  ))
  tu1 <- as.data.frame(test_usage(
    iso,
    group = "condition",
    levels = c("control", "preeclampsia"),
    scale = TRUE
  ))
  ue_write(tu0, pq("package_test_usage_raw"), ue, seed)
  ue_write(tu1, pq("package_test_usage_scaled"), ue, seed)
  # F5 enrichment (secondary): Reactome C2, log2FC against lfc_z
  # the GMT moved out of inst/extdata after the snapshot; the same file (sha256
  # c605faea...) is read from the repository
  gmt <- file.path(
    "vignettes",
    "articles",
    "data",
    "preeclampsia",
    "c2.cp.reactome.v2023.1.Hs.symbols.gmt.txt"
  )
  sets <- read_gmt(gmt)
  flat <- function(e, rank) {
    e <- as.data.frame(e)
    e$leading_edge <- vapply(e$leadingEdge, paste, character(1), collapse = ";")
    e <- e[, !vapply(e, is.list, logical(1)), drop = FALSE]
    e$rank_by <- rank
    e
  }
  set.seed(wpb_seed(ue, 200L))
  e_lfc <- run_enrichment(iso, sets, rank_by = "log2FC")
  set.seed(wpb_seed(ue, 200L))
  e_z <- run_enrichment(iso, sets, rank_by = "lfc_z")
  e_both <- rbind(flat(e_lfc, "log2FC"), flat(e_z, "lfc_z"))
  ue_write(e_both, pq("f5_enrichment"), ue, wpb_seed(ue, 200L))
  say("package stage done, peak RSS %.2f GB", ue_peak_rss_gb())
  rm(iso)
  invisible(gc())
}

# Stage productivity ---------------------------------------------------------------
# F4 package outputs. Amendment A18: rerun with the repository code
# (ISOFORMIC_UE_PKG=.) after productive_fraction_sd and delta_z changed
# definition; the other package outputs stay those of the snapshot.
if (run_stage("productivity") && !done("package_productive_fraction")) {
  seed <- wpb_seed(ue, 2L)
  set.seed(seed)
  iso <- package_object()
  # test_productive_fraction and productive_fraction delta_z (assay counts)
  tpf <- as.data.frame(test_productive_fraction(
    iso,
    group = "condition",
    levels = c("control", "preeclampsia"),
    min_gene_count = 10
  ))
  ue_write(tpf, pq("package_test_productive_fraction"), ue, seed)
  pf <- as.data.frame(productive_fraction(
    iso,
    group = "condition",
    levels = c("control", "preeclampsia"),
    assay = "counts",
    min_gene_count = 10
  ))
  ue_write(pf, pq("package_productive_fraction"), ue, seed)
  say("productivity stage done, peak RSS %.2f GB", ue_peak_rss_gb())
  rm(iso)
  invisible(gc())
}

# Stage validate ---------------------------------------------------------------------
if (run_stage("validate") && !done("validation")) {
  seed <- wpb_seed(ue, 3L)
  pk <- wpb_read(pq("package_confidence"))
  pk <- pk[pk$feature_type == "transcript", ]
  t0 <- Sys.time()
  mine <- wpb_confidence(
    arr,
    s = seq_along(cohort),
    in_b = in_b,
    lib = lib_full[cohort],
    ids = U,
    gene_of = gene_of
  )
  say(
    "matrix confidence on U: %.1f min",
    as.numeric(difftime(Sys.time(), t0, units = "mins"))
  )
  ue_write(mine$table, pq("confidence_matrix_U"), ue, seed)
  m <- mine$table[match(pk$feature_id, mine$table$transcript_id), ]
  both <- !is.na(pk$confidence) | !is.na(m$confidence)
  checks <- list()
  checks$confidence <- data.frame(
    check = "confidence class: matrix code (38_audit_common.R) vs package de_confidence() on all of U",
    n = sum(both),
    n_disagree = sum(
      xor(is.na(pk$confidence[both]), is.na(m$confidence[both])) |
        (pk$confidence[both] != m$confidence[both]) %in% TRUE
    ),
    max_abs_diff = max(abs(pk$sign_agreement - m$sign_agreement), na.rm = TRUE),
    max_abs_diff_2 = max(abs(pk$mean_inf_rv - m$mean_inf_rv), na.rm = TRUE),
    note = "max_abs_diff: sign_agreement; max_abs_diff_2: mean_inf_rv"
  )
  sw <- wpb_read(pq("package_isoform_switches"))
  pk_sw <- sw$gene_id[sw$switch_type == "opposite_direction"]
  f2 <- wpb_read(pq("naive_f2_switches"))
  my_sw <- f2$gene_id[f2$method == "raw"]
  checks$switches <- data.frame(
    check = "switch genes: wpb_switch_genes() vs isoform_switches() (naive F1)",
    n = length(union(pk_sw, my_sw)),
    n_disagree = length(setdiff(union(pk_sw, my_sw), intersect(pk_sw, my_sw))),
    max_abs_diff = NA_real_,
    max_abs_diff_2 = NA_real_,
    note = ""
  )
  tpf <- wpb_read(pq("package_test_productive_fraction"))
  f4 <- wpb_read(pq("naive_f4_raw"))
  j <- match(tpf$gene_id, f4$gene_id)
  dd <- abs(tpf$delta_productive_fraction - f4$delta[j]) > 1e-9
  checks$prodfrac <- data.frame(
    check = paste(
      "productive fraction test: wpb_prod_test() (S0 productivity classes) vs",
      "test_productive_fraction(min_gene_count = 10) (classes recomputed on U)"
    ),
    n = nrow(tpf),
    n_disagree = sum(is.na(j)) +
      sum(!f4$gene_id %in% tpf$gene_id) +
      sum(dd, na.rm = TRUE),
    max_abs_diff = max(abs(tpf$pvalue - f4$pvalue[j]), na.rm = TRUE),
    max_abs_diff_2 = max(
      abs(tpf$delta_productive_fraction - f4$delta[j]),
      na.rm = TRUE
    ),
    note = sprintf(
      "max_abs_diff: p-value; max_abs_diff_2: delta; %d genes differ in delta, %d of them with a p-value in either table, %d with q < 0.05 in either",
      sum(dd, na.rm = TRUE),
      sum(dd & (!is.na(tpf$pvalue) | !is.na(f4$pvalue[j])), na.rm = TRUE),
      sum(
        dd & ((tpf$padj < 0.05) %in% TRUE | (f4$qvalue[j] < 0.05) %in% TRUE),
        na.rm = TRUE
      )
    )
  )
  tu0 <- wpb_read(pq("package_test_usage_raw"))
  f3 <- wpb_read(pq("naive_f3_raw"))
  j <- match(tu0$transcript_id, f3$transcript_id)
  checks$usage_raw <- data.frame(
    check = "usage test raw: wpb_usage() vs test_usage(scale = FALSE)",
    n = nrow(tu0),
    n_disagree = sum(is.na(j)) + sum(!f3$transcript_id %in% tu0$transcript_id),
    max_abs_diff = max(abs(tu0$pvalue - f3$pvalue[j]), na.rm = TRUE),
    max_abs_diff_2 = max(abs(tu0$qvalue - f3$qvalue[j]), na.rm = TRUE),
    note = "max_abs_diff: p-value; max_abs_diff_2: q-value"
  )
  tu1 <- wpb_read(pq("package_test_usage_scaled"))
  f3s <- wpb_read(pq("naive_f3_scaled"))
  j <- match(tu1$transcript_id, f3s$transcript_id)
  checks$usage_scaled <- data.frame(
    check = paste(
      "usage test scaled: wpb_usage(od = all folded units) vs test_usage(scale = TRUE),",
      "whose overdispersion uses the draws of the object (U only): expected to differ slightly"
    ),
    n = nrow(tu1),
    n_disagree = sum(
      xor(tu1$qvalue < 0.05, f3s$qvalue[j] < 0.05),
      na.rm = TRUE
    ) +
      sum(is.na(j)),
    max_abs_diff = stats::median(
      abs(log10(tu1$pvalue) - log10(f3s$pvalue[j])),
      na.rm = TRUE
    ),
    max_abs_diff_2 = sum(tu1$qvalue < 0.05) - sum(f3s$qvalue < 0.05),
    note = "n_disagree: calls at q < 0.05 differing; max_abs_diff: median |log10 p| difference; max_abs_diff_2: package minus matrix call count"
  )
  val <- do.call(rbind, checks)
  print(val)
  ue_write(val, pq("validation"), ue, seed)
}

# Stage perdraw ------------------------------------------------------------------------
if (run_stage("perdraw") && !done("perdraw_f4")) {
  f1_raw <- wpb_read(pq("naive_f1_raw"))
  g_raw <- wpb_read(pq("naive_g_raw"))
  f3_raw <- wpb_read(pq("naive_f3_raw"))
  f4_raw <- wpb_read(pq("naive_f4_raw"))
  k_f1 <- f1_raw$feature_id
  k_g <- g_raw$feature_id
  k_f3 <- f3_raw$transcript_id
  k_f4 <- f4_raw$gene_id
  one_draw <- function(k) {
    m <- arr[,, k]
    f1 <- wpb_edger(m, groups, keep = k_f1)
    g <- wpb_edger(wpb_gene_sums(m, gene_of[U]), groups, keep = k_g)
    f3 <- wpb_usage(m, gene_of, groups, kept = k_f3, gene_p = FALSE)
    f4 <- wpb_prod_test(m, gene_of, is_prod, groups, k_f4)
    list(
      f1 = data.frame(
        draw = k,
        f1[, c("feature_id", "logFC", "pvalue", "qvalue")]
      ),
      g = data.frame(
        draw = k,
        g[, c("feature_id", "logFC", "pvalue", "qvalue")]
      ),
      f3 = data.frame(
        draw = k,
        f3[, c("transcript_id", "usage_log2FC", "pvalue", "qvalue")]
      ),
      f4 = data.frame(draw = k, f4[, c("gene_id", "delta", "pvalue", "qvalue")])
    )
  }
  t0 <- Sys.time()
  res <- parallel::mclapply(
    seq_len(n_draw),
    function(k) {
      out <- tryCatch(one_draw(k), error = function(e) e)
      say(
        "draw %d %s",
        k,
        if (inherits(out, "error")) conditionMessage(out) else "done"
      )
      out
    },
    mc.cores = workers,
    mc.preschedule = FALSE
  )
  bad <- vapply(res, function(r) inherits(r, "error") || is.null(r), logical(1))
  if (any(bad)) {
    stop(
      "per-draw refit failed for draws ",
      paste(which(bad), collapse = ","),
      call. = FALSE
    )
  }
  for (fam in c("f1", "g", "f3", "f4")) {
    ue_write(
      do.call(rbind, lapply(res, `[[`, fam)),
      pq(paste0("perdraw_", fam)),
      ue
    )
  }
  say(
    "perdraw done in %.1f min, peak RSS %.2f GB",
    as.numeric(difftime(Sys.time(), t0, units = "mins")),
    ue_peak_rss_gb()
  )
  rm(res)
  invisible(gc())
}

# Stage null: 20 null splits of the control units (non-circular link) --------------
if (run_stage("null") && !done("null_link_splits")) {
  ctrl <- cohort[groups == "control"]
  csex <- as.character(sex[groups == "control"])
  rows <- list()
  ids <- list()
  for (i in seq_len(20L)) {
    seed <- wpb_seed(ue, 100L + i)
    set.seed(seed)
    side <- unlist(lapply(split(ctrl, csex), function(r) {
      sample(r, floor(length(r) / 2))
    }))
    ib <- ctrl %in% side
    gr <- factor(ifelse(ib, "b", "a"), c("a", "b"))
    s_idx <- match(ctrl, cohort)
    fit <- wpb_edger(counts_U[, ctrl], gr)
    calls <- fit$feature_id[fit$pvalue < 0.001]
    cc <- if (length(calls)) {
      wpb_confidence(
        arr,
        s = s_idx,
        in_b = ib,
        lib = lib_full[ctrl],
        ids = calls,
        gene_of = gene_of
      )$table
    } else {
      NULL
    }
    q05 <- fit$feature_id[fit$qvalue < 0.05]
    if (length(calls)) {
      ids[[i]] <- data.frame(
        split = i,
        transcript_id = calls,
        pvalue = fit$pvalue[match(calls, fit$feature_id)],
        q_lt_0.05 = calls %in% q05,
        confidence = cc$confidence[match(calls, cc$transcript_id)]
      )
    }
    rows[[i]] <- data.frame(
      split = i,
      side_b = paste(sort(side), collapse = ","),
      n_tested = nrow(fit),
      n_calls_q05 = length(q05),
      n_calls_p001 = length(calls),
      n_flag_q05 = sum(!cc$confidence[cc$transcript_id %in% q05] %in% "high"),
      n_flag_p001 = sum(!cc$confidence %in% "high"),
      n_na_p001 = sum(is.na(cc$confidence)),
      seed = seed
    )
    say(
      "null split %d: %d calls at p < 0.001, %d at q < 0.05",
      i,
      length(calls),
      length(q05)
    )
  }
  ue_write(do.call(rbind, ids), pq("null_link_calls"), ue, wpb_seed(ue, 100L))
  ue_write(do.call(rbind, rows), pq("null_link_splits"), ue, wpb_seed(ue, 100L))
}

# Stage summarise ------------------------------------------------------------------------
if (run_stage("summarise")) {
  source("data-raw/uncertainty_evidence/30b_prevalence_summary.R")
}
ue_write_session(prev_dir, ue)
say(
  "30_audit_prevalence stage %s finished, peak RSS %.2f GB",
  stage,
  ue_peak_rss_gb()
)
