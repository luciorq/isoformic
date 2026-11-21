# WP-B (placenta audit), part 2: held-out replication of real preeclampsia
# calls, the non-circular test of the uncertainty flags (decision rule R_B).
#
# Design (prereg.yaml `primary`, `splits`, `decision_rules$R_B`):
#   * 30 seeded half-splits of the 23 units of the S0 store (8 PE + the 15
#     primary control units; the shallower member of each A/B pair is not a
#     unit), sex-balanced: each side has 4 PE (1 female, 3 male) and 7 or 8
#     controls (3 male and 4 or 5 female). Each split is used in both
#     directions: 60 directed (discovery D, validation V) pairs.
#   * Discovery side only: edgeR QL DTE (filterByExpr on D) on transcript
#     counts of U, calls at q < 0.10 with the pre-specified fallback to the
#     top 500 by p when there are fewer than 100 calls; diffSplice DTU (same
#     rule); gene-level edgeR (control C3); confidence classes of the calls
#     (de_confidence() rules, add_uncertainty() defaults) from the discovery
#     draws only, with all 30 draws and with draws 1-15 (ablation);
#     ambiguity groups of the discovery draws (H3).
#   * Validation side: the same test on the held-out point counts (edgeR QL
#     on V, rows kept by filterByExpr on V plus every discovery call, so that
#     every call has a held-out p). A call replicates when the held-out
#     log2FC has the discovery sign and one-sided p < 0.05.
#   * Primary flag: class != high. Secondary flags: sign agreement < 0.95,
#     top InfRV tertile (tertiles of U from the discovery draws), ambiguous
#     partner, unique_fraction == 0 (S0 stratum, PE samples), not
#     identifiable (S0 per-sample tables of the PE samples of D).
#   * C1: label permutation of the discovery side within fetal-sex blocks
#     (100 permutations cycling over the 60 directed pairs); calls, classes
#     and held-out replication as above, with the fixed S0 strata.
#   * Every per-call table carries the fixed S0 strata (InfRV tertile of the
#     PE samples, expression decile) for the per-stratum comparisons.
#   * Tracked features (cohort calls of 30_audit_prevalence.R and the story
#     genes FLT1, ENG, LEP): discovery and held-out statistics in every pair,
#     for 32_story.R.
#   * The statistics (models, bootstrap, controls) are in
#     31b_replication_summary.R (stage summarise).
#
# Stages (ISOFORMIC_UE_STAGE: pairs | perm | summarise | all)
# Inputs: prereg.yaml and the snapshot; the S0 store; the naive tables of
#   30_audit_prevalence.R (tracked features).
# Env vars:
#   ISOFORMIC_UE_OUT, ISOFORMIC_UE_PKG, ISOFORMIC_UE_THREADS  (ue_common.R)
#   ISOFORMIC_UE_SPLITS   number of half-splits (default 30)
#   ISOFORMIC_UE_WORKERS  forked workers (default 3)
#   ISOFORMIC_UE_PERMS    C1 permutations (default 100)
#   ISOFORMIC_UE_BOOT     two-level bootstrap resamples (default 1000)
#   ISOFORMIC_UE_STAGE
#   OPENBLAS_NUM_THREADS, OMP_NUM_THREADS  set to 1: the glm bootstrap of the
#                         summary is otherwise multithreaded in every worker
#   ISOFORMIC_UE_PILOT    "true": pilot gene subset, placenta_audit_pilot/
# Outputs: $ISOFORMIC_UE_OUT/placenta_audit/replication/
#   splits.parquet, pairs/<pair>_{dte,dtu,gene,groups,tracked,info}.parquet,
#   perm/perm<i>.parquet, summary_*.parquet, pooled_dte.parquet
# Memory (estimate): draws of U 0.32 GB shared by forked workers; per worker
#   the draws of the discovery samples (57,294 x 12 x 30 x 8 B = 0.16 GB, two
#   copies) and the edgeR fits: < 1.5 GB per worker, < 6 GB with 3 workers.
# Measured with 3 workers (peak RSS of the main process): pairs 43 min,
#   1.5 GB; perm 15 min, 1.5 GB; summarise 33 min, 2.1 GB.
# Run from the repository root, inside a memory cap of about 8 GB (see
# README.md):
#   OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
#   ISOFORMIC_UE_SPLITS=30 ISOFORMIC_UE_WORKERS=3 \
#     pixi run --frozen Rscript data-raw/uncertainty_evidence/31_halfsplit_replication.R

source("data-raw/uncertainty_evidence/ue_common.R")
ue <- ue_setup("31_halfsplit_replication")
source("data-raw/uncertainty_evidence/38_audit_common.R")
stage <- Sys.getenv("ISOFORMIC_UE_STAGE", "all")
n_splits <- as.integer(Sys.getenv("ISOFORMIC_UE_SPLITS", "30"))
workers <- as.integer(Sys.getenv("ISOFORMIC_UE_WORKERS", "3"))
n_perm <- as.integer(Sys.getenv("ISOFORMIC_UE_PERMS", "100"))
n_boot <- as.integer(Sys.getenv("ISOFORMIC_UE_BOOT", "1000"))
run_stage <- function(s) stage %in% c(s, "all")
rep_dir <- wpb_dir("replication")
pair_dir <- file.path(rep_dir, "pairs")
perm_dir <- file.path(rep_dir, "perm")
dir.create(pair_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(perm_dir, recursive = TRUE, showWarnings = FALSE)
log_path <- wpb_dir("logs", "31_halfsplit_replication.log")
say <- wpb_logger(log_path)
say(
  "stage %s, splits %d, workers %d, snapshot %s",
  stage,
  n_splits,
  workers,
  ue$git_sha
)

st <- wpb_load_store(draws = stage != "summarise", moments = FALSE, say = say)
units <- st$units
runs <- st$runs
U <- st$U
gene_of <- st$gene_of
txu <- st$tx
counts_U <- st$counts
lib_full <- st$lib_size
arr <- st$arr
# the full strata table (650,973 units) is not needed after this point; freed
# so that forked workers do not copy it when their garbage collector runs
st$strata <- NULL
st$arr <- NULL
invisible(gc())
type_class <- ifelse(
  txu$transcript_type %in%
    c(
      "protein_coding",
      "retained_intron",
      "nonsense_mediated_decay",
      "lncRNA",
      "protein_coding_CDS_not_defined"
    ),
  txu$transcript_type,
  "other"
)
names(type_class) <- U
fixed_tertile <- stats::setNames(txu$infrv_tertile, U)
fixed_decile <- stats::setNames(txu$expr_decile, U)
zero_unique <- stats::setNames(txu$zero_unique, U)
pe_runs <- runs[units$condition == "preeclampsia"]
ident_by_sample <- lapply(stats::setNames(pe_runs, pe_runs), function(s) {
  x <- wpb_identifiability(s)
  stats::setNames(x$non_identifiable, x$transcript_id)
})

# Tracked features: cohort calls (30_audit_prevalence.R) and story genes
prev_dir <- wpb_dir("prevalence")
f1_path <- file.path(prev_dir, "naive_f1_raw.parquet")
if (!file.exists(f1_path)) {
  stop("Run 30_audit_prevalence.R (stage naive) first", call. = FALSE)
}
f1c <- wpb_read(f1_path)
gc0 <- wpb_read(file.path(prev_dir, "naive_g_raw.parquet"))
story_genes <- unique(txu$gene_id[txu$gene_name %in% c("FLT1", "ENG", "LEP")])
tracked_tx <- unique(c(
  f1c$feature_id[f1c$qvalue < 0.05],
  U[gene_of[U] %in% story_genes]
))
tracked_gene <- unique(c(gc0$feature_id[gc0$qvalue < 0.05], story_genes))
rm(f1c, gc0)

# Splits -----------------------------------------------------------------------------
splits_path <- file.path(rep_dir, "splits.parquet")
make_splits <- function(n) {
  seed <- wpb_seed(ue, 10L)
  set.seed(seed)
  blocks <- split(runs, paste(units$condition, units$fetal_sex))
  seen <- character()
  out <- list()
  tries <- 0L
  while (length(out) < n && tries < 100000L) {
    tries <- tries + 1L
    a <- sort(unlist(
      lapply(blocks, function(r) {
        k <- length(r) / 2
        if (k != floor(k)) {
          k <- floor(k) + sample(0:1, 1)
        }
        sample(r, k)
      }),
      use.names = FALSE
    ))
    b <- sort(setdiff(runs, a))
    key <- paste(
      sort(c(paste(a, collapse = ","), paste(b, collapse = ","))),
      collapse = "|"
    )
    if (key %in% seen) {
      next
    }
    seen <- c(seen, key)
    out[[length(out) + 1L]] <- a
  }
  do.call(
    rbind,
    lapply(seq_along(out), function(i) {
      data.frame(
        split = i,
        sample_id = runs,
        condition = units$condition,
        fetal_sex = units$fetal_sex,
        side = ifelse(runs %in% out[[i]], "A", "B")
      )
    })
  )
}
if (!file.exists(splits_path)) {
  sp <- make_splits(n_splits)
  ue_write(sp, splits_path, ue, wpb_seed(ue, 10L))
}
splits <- wpb_read(splits_path)
chk <- stats::aggregate(
  sample_id ~ split + side + condition,
  data = splits,
  FUN = length
)
stopifnot(chk$sample_id[chk$condition == "preeclampsia"] == 4L)
directed <- expand.grid(
  split = seq_len(n_splits),
  direction = c("AB", "BA"),
  stringsAsFactors = FALSE
)
directed$pair_id <- sprintf("s%02d_%s", directed$split, directed$direction)
sides_of <- function(sp, direction) {
  s <- splits[splits$split == sp, ]
  a <- s$sample_id[s$side == "A"]
  b <- s$sample_id[s$side == "B"]
  if (direction == "AB") list(D = a, V = b) else list(D = b, V = a)
}
grp_of <- function(r) {
  factor(units$condition[match(r, runs)], c("control", "preeclampsia"))
}

calls_rule <- function(ids, p, q, thr = 0.10, min_calls = 100L, top = 500L) {
  ok <- !is.na(p)
  hits <- ids[ok & q < thr]
  if (length(hits) >= min_calls) {
    return(list(ids = hits, fallback = FALSE))
  }
  o <- order(p)
  o <- o[ok[o]]
  list(ids = ids[o[seq_len(min(top, length(o)))]], fallback = TRUE)
}

# One directed pair (or a permutation of its discovery labels) --------------------
analyse_pair <- function(D, V, labels_D = NULL, perm = FALSE) {
  sD <- match(D, runs)
  gD <- labels_D %||% grp_of(D)
  gV <- grp_of(V)
  mD <- counts_U[, D]
  mV <- counts_U[, V]
  fD <- wpb_edger(mD, gD)
  cr <- calls_rule(fD$feature_id, fD$pvalue, fD$qvalue)
  ids <- cr$ids
  keepV <- edgeR::filterByExpr(mV, group = gV)
  keepV <- keepV | rownames(mV) %in% c(ids, if (!perm) tracked_tx)
  fV <- wpb_edger(mV, gV, keep = keepV)
  fd <- fD[match(ids, fD$feature_id), ]
  fv <- fV[match(ids, fV$feature_id), ]
  in_b <- gD == "preeclampsia"
  cf <- wpb_confidence(
    arr,
    s = sD,
    in_b = in_b,
    lib = lib_full[D],
    ids = ids,
    gene_of = gene_of
  )
  cf15 <- wpb_confidence(
    arr,
    s = sD,
    in_b = in_b,
    lib = lib_full[D],
    ids = ids,
    gene_of = gene_of,
    draws = 1:15
  )
  rvD <- cf$mean_inf_rv_all
  assessedD <- cf$mean_count_all >= 1
  terc <- stats::quantile(rvD[assessedD], c(1 / 3, 2 / 3), names = FALSE)
  pe_D <- intersect(D, pe_runs)
  ni <- vapply(
    ids,
    function(i) {
      v <- unlist(lapply(ident_by_sample[pe_D], function(x) x[i]))
      v <- v[!is.na(v)]
      if (length(v)) any(v) else NA
    },
    logical(1)
  )
  tab <- data.frame(
    transcript_id = ids,
    gene_id = unname(gene_of[ids]),
    fallback = cr$fallback,
    d_logFC = fd$logFC,
    d_p = fd$pvalue,
    d_q = fd$qvalue,
    d_absz = abs(fd$z),
    d_mean_count = rowMeans(mD[ids, , drop = FALSE]),
    n_isoforms = txu$n_isoforms_gene[match(ids, U)],
    type_class = unname(type_class[ids]),
    mean_inf_rv = cf$table$mean_inf_rv,
    sign_agreement = cf$table$sign_agreement,
    uncertain = cf$table$uncertain,
    ambiguous = cf$table$ambiguous,
    confidence = cf$table$confidence,
    sign_agreement15 = cf15$table$sign_agreement,
    confidence15 = cf15$table$confidence,
    infrv_top_tertile = unname(rvD[ids] >= terc[2]),
    infrv_top_tertile15 = NA,
    zero_unique = unname(zero_unique[ids]),
    non_identifiable = unname(ni),
    fixed_infrv_tertile = unname(fixed_tertile[ids]),
    fixed_expr_decile = unname(fixed_decile[ids]),
    v_logFC = fv$logFC,
    v_p = fv$pvalue
  )
  rv15 <- cf15$mean_inf_rv_all
  terc15 <- stats::quantile(
    rv15[cf15$mean_count_all >= 1],
    c(1 / 3, 2 / 3),
    names = FALSE
  )
  tab$infrv_top_tertile15 <- unname(rv15[ids] >= terc15[2])
  tab$v_p1 <- wpb_one_sided(tab$v_logFC, tab$v_p, sign(tab$d_logFC))
  tab$replicated <- !is.na(tab$v_p1) & tab$v_p1 < 0.05
  tab$sign_concordant <- !is.na(tab$v_logFC) &
    sign(tab$v_logFC) == sign(tab$d_logFC)
  tab$group_id <- NA_character_
  out <- list(dte = tab)
  if (perm) {
    return(out)
  }
  # H3: discovery ambiguity groups of the calls, summed and re-tested
  grp <- wpb_groups(cf$pairs)
  if (nrow(grp)) {
    members <- split(grp$transcript_id, grp$ambiguity_group)
    gsD <- t(vapply(
      members,
      function(m) colSums(mD[m, , drop = FALSE]),
      numeric(ncol(mD))
    ))
    gsV <- t(vapply(
      members,
      function(m) colSums(mV[m, , drop = FALSE]),
      numeric(ncol(mV))
    ))
    kD <- edgeR::filterByExpr(mD, group = gD)
    fgD <- wpb_edger(
      rbind(mD[kD, , drop = FALSE], gsD),
      gD,
      keep = rep(TRUE, sum(kD) + nrow(gsD))
    )
    fgV <- wpb_edger(
      rbind(mV[keepV, , drop = FALSE], gsV),
      gV,
      keep = rep(TRUE, sum(keepV) + nrow(gsV))
    )
    gd <- fgD[match(names(members), fgD$feature_id), ]
    gv <- fgV[match(names(members), fgV$feature_id), ]
    out$groups <- data.frame(
      group_id = names(members),
      members = vapply(members, paste, character(1), collapse = ";"),
      n_members = lengths(members),
      n_calls_in_group = vapply(
        members,
        function(m) sum(m %in% ids),
        integer(1)
      ),
      d_logFC = gd$logFC,
      d_p = gd$pvalue,
      v_logFC = gv$logFC,
      v_p = gv$pvalue
    )
    out$groups$v_p1 <- wpb_one_sided(gv$logFC, gv$pvalue, sign(gd$logFC))
    out$groups$replicated <- !is.na(out$groups$v_p1) & out$groups$v_p1 < 0.05
    out$dte$group_id <- grp$ambiguity_group[match(ids, grp$transcript_id)]
  }
  # Tracked transcripts (cohort calls and story genes)
  tr <- tracked_tx
  out$tracked <- data.frame(
    feature_id = tr,
    level = "transcript",
    d_logFC = fD$logFC[match(tr, fD$feature_id)],
    d_p = fD$pvalue[match(tr, fD$feature_id)],
    d_called = tr %in% ids,
    v_logFC = fV$logFC[match(tr, fV$feature_id)],
    v_p = fV$pvalue[match(tr, fV$feature_id)]
  )
  # Gene level (C3)
  gDc <- wpb_gene_sums(mD, gene_of[U])
  gVc <- wpb_gene_sums(mV, gene_of[U])
  ggD <- wpb_edger(gDc, gD)
  gcr <- calls_rule(ggD$feature_id, ggD$pvalue, ggD$qvalue)
  keepGV <- edgeR::filterByExpr(gVc, group = gV) |
    rownames(gVc) %in% c(gcr$ids, tracked_gene)
  ggV <- wpb_edger(gVc, gV, keep = keepGV)
  gmax <- tapply(rvD, gene_of[U], max)
  gd <- ggD[match(gcr$ids, ggD$feature_id), ]
  gv <- ggV[match(gcr$ids, ggV$feature_id), ]
  out$gene <- data.frame(
    gene_id = gcr$ids,
    fallback = gcr$fallback,
    d_logFC = gd$logFC,
    d_p = gd$pvalue,
    d_q = gd$qvalue,
    d_absz = abs(gd$z),
    d_mean_count = rowMeans(gDc[gcr$ids, , drop = FALSE]),
    max_inf_rv = as.vector(gmax[gcr$ids]),
    n_isoforms_U = as.integer(table(gene_of)[gcr$ids]),
    v_logFC = gv$logFC,
    v_p = gv$pvalue
  )
  out$gene$v_p1 <- wpb_one_sided(gv$logFC, gv$pvalue, sign(gd$logFC))
  out$gene$replicated <- !is.na(out$gene$v_p1) & out$gene$v_p1 < 0.05
  trg <- tracked_gene
  out$tracked <- rbind(
    out$tracked,
    data.frame(
      feature_id = trg,
      level = "gene",
      d_logFC = ggD$logFC[match(trg, ggD$feature_id)],
      d_p = ggD$pvalue[match(trg, ggD$feature_id)],
      d_called = trg %in% gcr$ids,
      v_logFC = ggV$logFC[match(trg, ggV$feature_id)],
      v_p = ggV$pvalue[match(trg, ggV$feature_id)]
    )
  )
  # DTU (diffSplice, raw counts; V keeps its own filter)
  uD <- wpb_usage(mD, gene_of, gD, gene_p = FALSE)
  uV <- wpb_usage(mV, gene_of, gV, gene_p = FALSE)
  ucr <- calls_rule(uD$transcript_id, uD$pvalue, uD$qvalue)
  ud <- uD[match(ucr$ids, uD$transcript_id), ]
  uv <- uV[match(ucr$ids, uV$transcript_id), ]
  cfu <- wpb_confidence(
    arr,
    s = sD,
    in_b = in_b,
    lib = lib_full[D],
    ids = ucr$ids,
    gene_of = gene_of
  )$table
  out$dtu <- data.frame(
    transcript_id = ucr$ids,
    gene_id = unname(gene_of[ucr$ids]),
    fallback = ucr$fallback,
    d_usage = ud$usage_log2FC,
    d_p = ud$pvalue,
    d_q = ud$qvalue,
    d_mean_count = rowMeans(mD[ucr$ids, , drop = FALSE]),
    confidence = cfu$confidence,
    mean_inf_rv = cfu$mean_inf_rv,
    sign_agreement = cfu$sign_agreement,
    v_usage = uv$usage_log2FC,
    v_p = uv$pvalue,
    v_tested = ucr$ids %in% uV$transcript_id
  )
  out$dtu$v_p1 <- wpb_one_sided(
    uv$usage_log2FC,
    uv$pvalue,
    sign(ud$usage_log2FC)
  )
  out$dtu$replicated <- !is.na(out$dtu$v_p1) & out$dtu$v_p1 < 0.05
  out$info <- data.frame(
    n_D = length(D),
    n_V = length(V),
    n_dte_tested = nrow(fD),
    n_dte_q10 = sum(fD$qvalue < 0.10),
    n_dte_calls = length(ids),
    dte_fallback = cr$fallback,
    n_dtu_tested = nrow(uD),
    n_dtu_q10 = sum(uD$qvalue < 0.10),
    dtu_fallback = ucr$fallback,
    n_gene_tested = nrow(ggD),
    n_gene_q10 = sum(ggD$qvalue < 0.10),
    gene_fallback = gcr$fallback,
    infrv_tertile_cut_low_mid = terc[1],
    infrv_tertile_cut_mid_high = terc[2],
    D = paste(D, collapse = ","),
    V = paste(V, collapse = ",")
  )
  out
}

# Stage pairs ----------------------------------------------------------------------------
if (run_stage("pairs")) {
  todo <- directed[
    !file.exists(file.path(pair_dir, paste0(directed$pair_id, "_dte.parquet"))),
  ]
  say("%d directed pairs to run", nrow(todo))
  invisible(parallel::mclapply(
    seq_len(nrow(todo)),
    function(i) {
      row <- todo[i, ]
      t0 <- Sys.time()
      sd <- sides_of(row$split, row$direction)
      res <- tryCatch(analyse_pair(sd$D, sd$V), error = function(e) e)
      if (inherits(res, "error")) {
        say("%s ERROR %s", row$pair_id, conditionMessage(res))
        return(NULL)
      }
      for (nm in c("gene", "dtu", "groups", "info", "tracked")) {
        if (!is.null(res[[nm]])) {
          x <- res[[nm]]
          x$pair_id <- row$pair_id
          ue_write(
            x,
            file.path(pair_dir, paste0(row$pair_id, "_", nm, ".parquet")),
            ue
          )
        }
      }
      x <- res$dte
      x$pair_id <- row$pair_id
      # written last: its presence marks the pair as done
      ue_write(x, file.path(pair_dir, paste0(row$pair_id, "_dte.parquet")), ue)
      say(
        "%s done in %.1f min (%d DTE calls, fallback %s; %d q < 0.10), RSS %.2f GB",
        row$pair_id,
        as.numeric(difftime(Sys.time(), t0, units = "mins")),
        nrow(x),
        x$fallback[1],
        res$info$n_dte_q10,
        ue_peak_rss_gb()
      )
      NULL
    },
    mc.cores = workers,
    mc.preschedule = FALSE
  ))
}

# Stage perm (C1) ----------------------------------------------------------------------------
if (run_stage("perm")) {
  todo <- seq_len(n_perm)
  todo <- todo[
    !file.exists(file.path(perm_dir, sprintf("perm%03d.parquet", todo)))
  ]
  say("%d permutations to run", length(todo))
  invisible(parallel::mclapply(
    todo,
    function(i) {
      seed <- wpb_seed(ue, 500L + i)
      set.seed(seed)
      row <- directed[(i - 1L) %% nrow(directed) + 1L, ]
      sd <- sides_of(row$split, row$direction)
      g <- grp_of(sd$D)
      sx <- units$fetal_sex[match(sd$D, runs)]
      gp <- g
      for (b in unique(sx)) {
        ix <- which(sx == b)
        gp[ix] <- g[ix][sample.int(length(ix))]
      }
      res <- tryCatch(
        analyse_pair(sd$D, sd$V, labels_D = gp, perm = TRUE),
        error = function(e) e
      )
      if (inherits(res, "error")) {
        say("perm %d ERROR %s", i, conditionMessage(res))
        return(NULL)
      }
      x <- res$dte
      x$perm <- i
      x$pair_id <- row$pair_id
      x$n_relabelled <- sum(gp != g)
      ue_write(x, file.path(perm_dir, sprintf("perm%03d.parquet", i)), ue, seed)
      say(
        "perm %d (%s) done: %d calls, %d relabelled",
        i,
        row$pair_id,
        nrow(x),
        sum(gp != g)
      )
      NULL
    },
    mc.cores = workers,
    mc.preschedule = FALSE
  ))
}

# Stage summarise --------------------------------------------------------------------------------
if (run_stage("summarise")) {
  source("data-raw/uncertainty_evidence/31b_replication_summary.R")
}
ue_write_session(rep_dir, ue)
say(
  "31_halfsplit_replication stage %s finished, peak RSS %.2f GB",
  stage,
  ue_peak_rss_gb()
)
