# WP-B (placenta audit), part 3: the pre-specified story genes.
#
# Purpose: show on named genes what the audit and the replication test mean,
# with rules fixed in prereg.yaml (`story_genes`) and in the WP-B plan before
# any result was seen, and reported whatever the outcome:
#   * FLT1 per isoform: naive (edgeR QL) log2FC with its 95% Wald interval,
#     swish log2FC with the 2.5-97.5% range of the 30 per-draw refits,
#     sig_frac and sign_agree over the per-draw refits, draw correlation
#     with each sibling (inferential_correlation() of the package),
#     de_confidence() class, held-out support across the 60 directed
#     half-splits, appearances among null calls (WP-A per-split tables when
#     present, else the 20 control-only splits of 30_audit_prevalence.R),
#     and whether the soluble / membrane attribution keeps its sign in every
#     draw.
#   * ENG and LEP at gene level, with their isoforms.
#   * Rule-selected genes: fragile (WP-B plan rule: smallest naive p among
#     F2/F4 conclusions of genes with gene-level q < 0.05 and sig_frac <
#     0.5), fragile (prereg rule), robust (prereg rule).
#   * plot_de_confidence() and plot_isoform_correlation() of the package.
# Definitions (prereg leaves them open; fixed here before reading results):
#   cohort_sign_support_in_v_halves (called "held-out support" in the
#     prereg and used by its rules) = share of the 60 validation halves in
#     which the test has the cohort sign and one-sided p < 0.05 (all 60,
#     whether or not the discovery half called it). It is NOT held out for
#     a cohort call: the sign comes from all 23 samples, and every
#     validation half is part of them.
#   replicated_when_called (the genuinely held-out figure) = among the
#     directed pairs whose discovery half called the feature (q < 0.10), the
#     share in which the validation half has the discovery sign and
#     one-sided p < 0.05 (n_called_discovery pairs).
#   FLT1 soluble = protein-coding isoforms with a CDS < 2,500 nt (the sFLT1
#     proteins are 665-733 aa); membrane = protein-coding isoforms with a CDS
#     >= 3,600 nt (>= 1,199 aa; full-length FLT1 is 1,338 aa); others
#     "other". The attribution is the log2 ratio PE / control of the summed
#     soluble and summed membrane CPM, per draw and on the point counts.
#
# Inputs: outputs of 30_audit_prevalence.R and 31_halfsplit_replication.R;
#   the S0 store; WP-A null tables ($ISOFORMIC_UE_OUT/null_placenta/dte/main/tx)
#   when present; the GENCODE 50 annotation directory.
# Env vars: ISOFORMIC_UE_OUT, ISOFORMIC_UE_PKG, ISOFORMIC_UE_THREADS,
#   ISOFORMIC_UE_ANNOT, ISOFORMIC_UE_PILOT.
# Outputs: $ISOFORMIC_UE_OUT/placenta_audit/story/*.parquet and *.png
# Memory (estimate): draws of U 0.32 GB plus a package object of the story
#   genes: < 3 GB. Measured: about 2 min, peak 1.6 GB.
# Run from the repository root, inside a memory cap of about 6 GB (see
# README.md):
#   pixi run --frozen Rscript data-raw/uncertainty_evidence/32_story.R

source("data-raw/uncertainty_evidence/ue_common.R")
ue <- ue_setup("32_story")
source("data-raw/uncertainty_evidence/38_audit_common.R")
seed <- wpb_seed(ue, 600L)
set.seed(seed)
annot_dir <- Sys.getenv(
  "ISOFORMIC_UE_ANNOT",
  path.expand(
    "~/.cache/R/isoformic/gencode_v50_chr_patch_hapl_scaff_annotation_gff3"
  )
)
story_dir <- wpb_dir("story")
dir.create(story_dir, recursive = TRUE, showWarnings = FALSE)
prev_dir <- wpb_dir("prevalence")
rep_dir <- wpb_dir("replication")
say <- wpb_logger(wpb_dir("logs", "32_story.log"))
rd <- function(name) wpb_read(file.path(prev_dir, paste0(name, ".parquet")))

st <- wpb_load_store(draws = TRUE, moments = FALSE, say = say)
U <- st$U
txu <- st$tx
gene_of <- st$gene_of
cohort <- st$runs
groups <- st$groups
arr <- st$arr
counts_U <- st$counts
lib <- st$lib_size
n_draw <- dim(arr)[3L]
f1_raw <- rd("naive_f1_raw")
f1_div <- rd("naive_f1_div")
f1_swish <- rd("naive_f1_swish")
g_raw <- rd("naive_g_raw")
g_swish <- rd("naive_g_swish")
conf <- rd("package_confidence")
conf <- conf[conf$feature_type == "transcript", ]
amb <- rd("package_ambiguity_pairs")
pd_f1 <- rd("perdraw_f1")
pd_g <- rd("perdraw_g")

perdraw_stats <- function(tab, ids, naive_lfc) {
  l <- matrix(NA_real_, length(ids), n_draw)
  q <- l
  i <- match(tab$feature_id, ids)
  ok <- !is.na(i)
  l[cbind(i[ok], tab$draw[ok])] <- tab$logFC[ok]
  q[cbind(i[ok], tab$draw[ok])] <- tab$qvalue[ok]
  agree <- sign(l) == sign(naive_lfc)
  agree[is.na(agree)] <- FALSE
  data.frame(
    sign_agree = ifelse(is.na(naive_lfc), NA, rowMeans(agree)),
    sig_frac = ifelse(
      is.na(naive_lfc),
      NA,
      rowMeans(agree & !is.na(q) & q < 0.05)
    ),
    draw_refit_lfc_lo = apply(l, 1, wpb_q, 0.025),
    draw_refit_lfc_hi = apply(l, 1, wpb_q, 0.975)
  )
}

# Held-out support across the directed half-splits
tracked <- do.call(
  rbind,
  lapply(
    list.files(
      file.path(rep_dir, "pairs"),
      pattern = "_tracked\\.parquet$",
      full.names = TRUE
    ),
    wpb_read
  )
)
support_of <- function(ids, level, cohort_lfc) {
  t <- tracked[tracked$level == level & tracked$feature_id %in% ids, ]
  t$cohort_sign <- sign(cohort_lfc[match(t$feature_id, ids)])
  t$v_p1 <- wpb_one_sided(t$v_logFC, t$v_p, t$cohort_sign)
  t$v_support <- !is.na(t$v_p1) & t$v_p1 < 0.05
  t$v_p1_d <- wpb_one_sided(t$v_logFC, t$v_p, sign(t$d_logFC))
  t$v_rep_d <- !is.na(t$v_p1_d) & t$v_p1_d < 0.05
  s <- do.call(
    rbind,
    lapply(split(t, t$feature_id), function(x) {
      dc <- x$d_called %in% TRUE
      data.frame(
        feature_id = x$feature_id[1],
        n_pairs = nrow(x),
        cohort_sign_support_in_v_halves = mean(x$v_support),
        n_called_discovery = sum(dc),
        replicated_when_called = if (any(dc)) mean(x$v_rep_d[dc]) else NA_real_
      )
    })
  )
  s[match(ids, s$feature_id), ]
}

# Null appearances: WP-A per-split tables when present, else the 20 splits of 30
wpa_dir <- file.path(ue_out_root(), "null_placenta", "dte", "main", "tx")
wpa_files <- list.files(wpa_dir, pattern = "\\.parquet$", full.names = TRUE)
null_of <- function(ids) {
  if (length(wpa_files)) {
    x <- wpb_query(
      sprintf(
        "SELECT transcript_id, split_id, M1_p FROM read_parquet([%s]) JOIN ids USING (transcript_id)",
        paste(sprintf("'%s'", wpa_files), collapse = ",")
      ),
      list(ids = data.frame(transcript_id = ids))
    )
    data.frame(
      transcript_id = ids,
      null_source = sprintf(
        "WP-A M1 (edgeR QL) per-split tables, %d splits",
        length(wpa_files)
      ),
      null_n_splits_tested = as.integer(table(factor(x$transcript_id, ids))),
      null_n_p_lt_0.001 = as.integer(table(factor(
        x$transcript_id[x$M1_p < 0.001],
        ids
      ))),
      null_n_p_lt_0.05 = as.integer(table(factor(
        x$transcript_id[x$M1_p < 0.05],
        ids
      )))
    )
  } else {
    nc <- rd("null_link_calls")
    data.frame(
      transcript_id = ids,
      null_source = "20 control-only splits of 30_audit_prevalence.R (calls at p < 0.001 only)",
      null_n_splits_tested = NA_integer_,
      null_n_p_lt_0.001 = vapply(
        ids,
        function(i) sum(nc$transcript_id == i),
        integer(1)
      ),
      null_n_p_lt_0.05 = NA_integer_
    )
  }
}

tx_story <- function(gene_names) {
  ids <- U[txu$gene_name %in% gene_names]
  f1 <- f1_raw[match(ids, f1_raw$feature_id), ]
  sw <- f1_swish[match(ids, f1_swish$feature_id), ]
  se <- abs(f1$logFC) / pmax(abs(f1$z), 1e-8)
  cf <- conf[match(ids, conf$feature_id), ]
  pdst <- perdraw_stats(pd_f1, ids, f1$logFC)
  sup <- support_of(ids, "transcript", f1$logFC)
  sib <- vapply(
    ids,
    function(i) {
      p <- amb[amb$transcript_id_1 == i | amb$transcript_id_2 == i, ]
      if (!nrow(p)) {
        return(NA_character_)
      }
      other <- ifelse(
        p$transcript_id_1 == i,
        p$transcript_id_2,
        p$transcript_id_1
      )
      nm <- txu$transcript_name[match(other, U)]
      o <- order(p$mean_cor)
      paste(sprintf("%s:%.2f", nm[o], p$mean_cor[o]), collapse = "; ")
    },
    character(1)
  )
  data.frame(
    gene_name = txu$gene_name[match(ids, U)],
    transcript_id = ids,
    transcript_name = txu$transcript_name[match(ids, U)],
    transcript_type = txu$transcript_type[match(ids, U)],
    mean_count = rowMeans(counts_U[ids, , drop = FALSE]),
    tested_edgeR = !is.na(f1$logFC),
    naive_log2FC = f1$logFC,
    naive_ci_lo = f1$logFC - 1.96 * se,
    naive_ci_hi = f1$logFC + 1.96 * se,
    naive_p = f1$pvalue,
    naive_q = f1$qvalue,
    divided_q = f1_div$qvalue[match(ids, f1_div$feature_id)],
    swish_log2FC = sw$logFC,
    swish_q = sw$qvalue,
    draw_refit_lfc_lo_2.5 = pdst$draw_refit_lfc_lo,
    draw_refit_lfc_hi_97.5 = pdst$draw_refit_lfc_hi,
    sign_agree = pdst$sign_agree,
    sig_frac = pdst$sig_frac,
    mean_inf_rv = cf$mean_inf_rv,
    de_confidence_sign_agreement = cf$sign_agreement,
    confidence = cf$confidence,
    confidence_reason = cf$reason,
    sibling_draw_correlations = sib,
    cohort_sign_support_in_v_halves = sup$cohort_sign_support_in_v_halves,
    n_called_discovery = sup$n_called_discovery,
    replicated_when_called = sup$replicated_when_called
  ) |>
    merge(null_of(ids), by = "transcript_id", sort = FALSE)
}

gene_story <- function(gene_names) {
  gids <- unique(txu$gene_id[txu$gene_name %in% gene_names])
  g <- g_raw[match(gids, g_raw$feature_id), ]
  gs <- g_swish[match(gids, g_swish$feature_id), ]
  se <- abs(g$logFC) / pmax(abs(g$z), 1e-8)
  pdst <- perdraw_stats(pd_g, gids, g$logFC)
  sup <- support_of(gids, "gene", g$logFC)
  data.frame(
    gene_name = txu$gene_name[match(gids, txu$gene_id)],
    gene_id = gids,
    n_isoforms_U = as.integer(table(gene_of)[gids]),
    naive_log2FC = g$logFC,
    naive_ci_lo = g$logFC - 1.96 * se,
    naive_ci_hi = g$logFC + 1.96 * se,
    naive_p = g$pvalue,
    naive_q = g$qvalue,
    swish_log2FC = gs$logFC,
    swish_q = gs$qvalue,
    draw_refit_lfc_lo_2.5 = pdst$draw_refit_lfc_lo,
    draw_refit_lfc_hi_97.5 = pdst$draw_refit_lfc_hi,
    sign_agree = pdst$sign_agree,
    sig_frac = pdst$sig_frac,
    cohort_sign_support_in_v_halves = sup$cohort_sign_support_in_v_halves,
    n_called_discovery = sup$n_called_discovery,
    replicated_when_called = sup$replicated_when_called
  )
}

# FLT1 ------------------------------------------------------------------------
annot_tx <- wpb_read(file.path(annot_dir, "annot_data_transcripts.parquet"))
flt1 <- tx_story("FLT1")
flt1$cds_length <- annot_tx$cds_length[match(
  flt1$transcript_id,
  annot_tx$transcript_id
)]
pc <- flt1$transcript_type == "protein_coding" & !is.na(flt1$cds_length)
flt1$form <- ifelse(
  pc & flt1$cds_length >= 3600,
  "membrane",
  ifelse(pc & flt1$cds_length < 2500, "soluble", "other")
)
ue_write(flt1, file.path(story_dir, "flt1_isoforms.parquet"), ue, seed)
say(
  "FLT1: %d isoforms in U (%s)",
  nrow(flt1),
  paste(flt1$transcript_name, flt1$form, flt1$confidence, collapse = "; ")
)

sol <- flt1$transcript_id[flt1$form == "soluble"]
mem <- flt1$transcript_id[flt1$form == "membrane"]
pe <- groups == "preeclampsia"
lfc_of <- function(m) {
  s <- colSums(m[sol, , drop = FALSE]) / lib[cohort] * 1e6
  b <- colSums(m[mem, , drop = FALSE]) / lib[cohort] * 1e6
  c(
    soluble_share_pe = mean((s / (s + b))[pe]),
    soluble_share_control = mean((s / (s + b))[!pe]),
    lfc_soluble = log2(mean(s[pe])) - log2(mean(s[!pe])),
    lfc_membrane = log2(mean(b[pe])) - log2(mean(b[!pe]))
  )
}
point <- lfc_of(counts_U)
dr <- t(vapply(seq_len(n_draw), function(k) lfc_of(arr[,, k]), numeric(4)))
diff_point <- point[["lfc_soluble"]] - point[["lfc_membrane"]]
diff_draw <- dr[, "lfc_soluble"] - dr[, "lfc_membrane"]
attr_tab <- data.frame(
  quantity = c(colnames(dr), "lfc_soluble_minus_membrane"),
  point = c(point, diff_point),
  draw_min = c(apply(dr, 2, min), min(diff_draw)),
  draw_max = c(apply(dr, 2, max), max(diff_draw)),
  same_sign_in_all_draws = c(
    apply(sweep(sign(dr), 2, sign(point), "=="), 2, all),
    all(sign(diff_draw) == sign(diff_point))
  ),
  n_soluble_isoforms = length(sol),
  soluble_isoforms = paste(
    flt1$transcript_name[flt1$form == "soluble"],
    collapse = ","
  ),
  n_membrane_isoforms = length(mem),
  membrane_isoforms = paste(
    flt1$transcript_name[flt1$form == "membrane"],
    collapse = ","
  )
)
ue_write(
  attr_tab,
  file.path(story_dir, "flt1_soluble_membrane.parquet"),
  ue,
  seed
)

# ENG, LEP ----------------------------------------------------------------------
genes_tab <- gene_story(c("FLT1", "ENG", "LEP"))
ue_write(
  genes_tab,
  file.path(story_dir, "story_genes_gene_level.parquet"),
  ue,
  seed
)
iso_tab <- tx_story(c("ENG", "LEP"))
ue_write(
  iso_tab,
  file.path(story_dir, "story_genes_isoforms.parquet"),
  ue,
  seed
)

# Rule-selected genes ------------------------------------------------------------
# Fragile, WP-B plan rule: smallest naive p among F2/F4 conclusions of genes
# with gene-level q < 0.05 and sig_frac < 0.5
c2 <- rd("conclusions_F2")
c4 <- rd("conclusions_F4")
cand <- rbind(
  c2[, c("family", "gene_id", "naive_p", "sig_frac", "gene_dge")],
  c4[, c("family", "gene_id", "naive_p", "sig_frac", "gene_dge")]
)
cand <- cand[cand$gene_dge & cand$sig_frac < 0.5, ]
cand <- cand[order(cand$naive_p), ]
fragile_task <- if (nrow(cand)) cand$gene_id[1] else NA_character_
# Fragile, prereg rule: among genes with an F1 q < 0.05 isoform, the gene with
# the most such isoforms whose calls all have the primary flag and are
# supported in < 50% of the validation halves (cohort_sign_support_in_v_halves, not held out) (ties: smallest min q)
f1c <- f1_raw[f1_raw$qvalue < 0.05, ]
f1c$gene_id <- gene_of[f1c$feature_id]
f1c$flag <- !conf$confidence[match(f1c$feature_id, conf$feature_id)] %in% "high"
f1c$cohort_sign_support_in_v_halves <- support_of(
  f1c$feature_id,
  "transcript",
  f1c$logFC
)$cohort_sign_support_in_v_halves
fr <- do.call(
  rbind,
  lapply(split(f1c, f1c$gene_id), function(x) {
    data.frame(
      gene_id = x$gene_id[1],
      n_calls = nrow(x),
      all_flag_and_low_support = all(
        x$flag & x$cohort_sign_support_in_v_halves < 0.5
      ),
      min_q = min(x$qvalue)
    )
  })
)
fr <- fr[fr$all_flag_and_low_support %in% TRUE, ]
fr <- fr[order(-fr$n_calls, fr$min_q), ]
fragile_prereg <- if (nrow(fr)) fr$gene_id[1] else NA_character_
# Robust, prereg rule: among genes with gene-level q < 0.05, the gene whose
# top isoform (smallest F1 p) has class high and support in the validation halves >= 0.9 (not held out),
# smallest gene q
gcalls <- g_raw[g_raw$qvalue < 0.05, ]
f1o <- f1_raw[order(f1_raw$pvalue), ]
f1o$gene_id <- gene_of[f1o$feature_id]
top_iso <- f1o[!duplicated(f1o$gene_id) & f1o$gene_id %in% gcalls$feature_id, ]
top_iso$confidence <- conf$confidence[match(
  top_iso$feature_id,
  conf$feature_id
)]
top_iso$cohort_sign_support_in_v_halves <- support_of(
  top_iso$feature_id,
  "transcript",
  top_iso$logFC
)$cohort_sign_support_in_v_halves
top_iso$gene_q <- gcalls$qvalue[match(top_iso$gene_id, gcalls$feature_id)]
rb <- top_iso[
  top_iso$confidence %in%
    "high" &
    top_iso$cohort_sign_support_in_v_halves >= 0.9,
]
robust_prereg <- if (nrow(rb)) {
  rb$gene_id[which.min(rb$gene_q)]
} else {
  NA_character_
}
sel <- data.frame(
  rule = c(
    "fragile (WP-B plan rule, F2/F4)",
    "fragile (prereg rule)",
    "robust (prereg rule)"
  ),
  gene_id = c(fragile_task, fragile_prereg, robust_prereg),
  n_candidates = c(nrow(cand), nrow(fr), nrow(rb))
)
sel$gene_name <- txu$gene_name[match(sel$gene_id, txu$gene_id)]
sel$family <- c(
  if (nrow(cand)) cand$family[1] else NA_character_,
  "F1_DTE",
  "G_gene_DGE"
)
print(sel)
ue_write(sel, file.path(story_dir, "rule_selected_genes.parquet"), ue, seed)
ue_write(
  cand,
  file.path(story_dir, "fragile_task_rule_candidates.parquet"),
  ue,
  seed
)
sel_names <- unique(stats::na.omit(sel$gene_name))
if (length(sel_names)) {
  ue_write(
    tx_story(sel_names),
    file.path(story_dir, "rule_selected_isoforms.parquet"),
    ue,
    seed
  )
  ue_write(
    gene_story(sel_names),
    file.path(story_dir, "rule_selected_gene_level.parquet"),
    ue,
    seed
  )
}

# Plots of the package -------------------------------------------------------------
plot_genes <- unique(c("FLT1", "ENG", "LEP", sel_names))
pids <- U[txu$gene_name %in% plot_genes]
iso <- IsoformicExperiment(
  experiment_name = "placenta_story",
  annot_path = annot_dir,
  assays = list(counts = counts_U[pids, ]),
  col_data = data.frame(
    sample_id = cohort,
    condition = as.character(groups),
    library_size = unname(lib[cohort])
  )
)
meta <- iso@metadata
meta$inf_reps <- lapply(seq_len(n_draw), function(k) arr[pids, , k])
iso@metadata <- meta
de_tx(iso) <- data.frame(
  transcript_id = pids,
  log2FC = f1_raw$logFC[match(pids, f1_raw$feature_id)],
  pvalue = f1_raw$pvalue[match(pids, f1_raw$feature_id)],
  qvalue = f1_raw$qvalue[match(pids, f1_raw$feature_id)]
)[pids %in% f1_raw$feature_id, ]
iso <- add_uncertainty(
  iso,
  group = "condition",
  levels = c("control", "preeclampsia")
)
pk <- as.data.frame(de_confidence(iso))
check <- data.frame(
  transcript_id = pk$feature_id,
  story_object_class = pk$confidence,
  whole_U_class = conf$confidence[match(pk$feature_id, conf$feature_id)]
)
check$same <- check$story_object_class %in%
  check$whole_U_class &
  (is.na(check$story_object_class) == is.na(check$whole_U_class))
say(
  "story object classes equal to the whole-U classes: %d of %d",
  sum(check$same),
  nrow(check)
)
ue_write(
  check,
  file.path(story_dir, "plot_object_class_check.parquet"),
  ue,
  seed
)
for (g in plot_genes) {
  p1 <- tryCatch(plot_de_confidence(iso, gene_name = g), error = function(e) e)
  if (inherits(p1, "error")) {
    say("plot_de_confidence %s: %s", g, conditionMessage(p1))
  } else {
    ggplot2::ggsave(
      file.path(story_dir, paste0("de_confidence_", g, ".png")),
      p1,
      width = 8,
      height = 5,
      dpi = 120
    )
  }
  p2 <- tryCatch(
    plot_isoform_correlation(iso, gene_name = g),
    error = function(e) e
  )
  if (inherits(p2, "error")) {
    say("plot_isoform_correlation %s: %s", g, conditionMessage(p2))
  } else {
    ggplot2::ggsave(
      file.path(story_dir, paste0("isoform_correlation_", g, ".png")),
      p2,
      width = 8,
      height = 7,
      dpi = 120
    )
  }
}
ue_write_session(story_dir, ue)
say("story done, peak RSS %.2f GB", ue_peak_rss_gb())
