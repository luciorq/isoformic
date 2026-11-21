# Pre-registration and frozen code snapshot of the uncertainty evidence study.
#
# Purpose: fix, before any test is run, (a) the package code every later
# script loads and (b) the seeds, split counts, methods, thresholds, the
# universe U, the primary flag and metric, the decision rules, the story
# genes and the stop rules. The maintainer commits prereg.yaml (copied to
# data-raw/uncertainty_evidence/prereg.yaml); scripts never change it after a
# test has run. Later scripts read it through ue_common.R.
#
# (a) Snapshot: `git archive <HEAD sha>` is written as a tarball, hashed
#     (sha256) and extracted to $ISOFORMIC_UE_OUT/pkg_snapshot/<sha>/. This
#     reads committed objects only (no checkout, the working tree is not
#     touched), so uncommitted edits of R/ by other workflows are excluded.
# (b) $ISOFORMIC_UE_OUT/prereg.yaml and a copy in this directory.
#
# Inputs: the git repository (HEAD). Nothing else.
# Env vars:
#   ISOFORMIC_UE_OUT       output root (required)
#   ISOFORMIC_V50_DIR      work directory of the Salmon runs (required; only
#                          recorded in prereg.yaml `data`)
#   ISOFORMIC_UE_THREADS   not used here
#   ISOFORMIC_UE_PKG       not used here (this script creates the snapshot)
#   ISOFORMIC_UE_REF       commit to snapshot, default HEAD
#   ISOFORMIC_UE_FORCE     "true" to overwrite an existing prereg.yaml
#                          (refused once any later step has written results)
#
# Resources: seconds, well under 1 GB (git archive of the package).
# Run from the repository root:
#   pixi run --frozen Rscript data-raw/uncertainty_evidence/01_preregister.R
# Then the maintainer commits data-raw/uncertainty_evidence/prereg.yaml,
# with the absolute paths it records replaced by ${ISOFORMIC_UE_OUT} and
# ${ISOFORMIC_V50_DIR} (the copy committed for this evaluation is redacted
# that way: amendment A10 of prereg_amendments.yaml). The prereg.yaml in
# $ISOFORMIC_UE_OUT keeps the real paths, which the later scripts read.

need_env <- function(name) {
  value <- Sys.getenv(name)
  if (!nzchar(value)) {
    stop(
      "Set the environment variable ",
      name,
      "; see data-raw/uncertainty_evidence/README.md",
      call. = FALSE
    )
  }
  value
}
out <- need_env("ISOFORMIC_UE_OUT")
v50 <- need_env("ISOFORMIC_V50_DIR")
ref <- Sys.getenv("ISOFORMIC_UE_REF", "HEAD")
force <- identical(tolower(Sys.getenv("ISOFORMIC_UE_FORCE")), "true")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
prereg_path <- file.path(out, "prereg.yaml")
if (file.exists(prereg_path) && !force) {
  stop(
    "prereg.yaml exists; set ISOFORMIC_UE_FORCE=true to rewrite it: ",
    prereg_path,
    call. = FALSE
  )
}
later <- setdiff(
  list.dirs(out, recursive = FALSE, full.names = FALSE),
  "pkg_snapshot"
)
if (force && length(later)) {
  stop(
    "Results already exist (",
    paste(later, collapse = ", "),
    "): the pre-registration is frozen.",
    call. = FALSE
  )
}

# (a) Code snapshot ------------------------------------------------------------
sha <- system2("git", c("rev-parse", ref), stdout = TRUE)
snap_root <- file.path(out, "pkg_snapshot")
dir.create(snap_root, recursive = TRUE, showWarnings = FALSE)
tar_path <- file.path(snap_root, paste0(sha, ".tar"))
snap_dir <- file.path(snap_root, sha)
if (!file.exists(tar_path)) {
  status <- system2("git", c("archive", "--format=tar", "-o", tar_path, sha))
  stopifnot(status == 0L)
}
tar_sha256 <- strsplit(system2("sha256sum", tar_path, stdout = TRUE), " ")[[1]][
  1
]
if (!dir.exists(snap_dir)) {
  dir.create(snap_dir)
  stopifnot(utils::untar(tar_path, exdir = snap_dir) == 0L)
  Sys.chmod(list.files(snap_dir, recursive = TRUE, full.names = TRUE), "0444")
}
dirty <- system2("git", c("status", "--porcelain", "--", "R"), stdout = TRUE)

# (b) Pre-registration ---------------------------------------------------------
pe <- sprintf("SRR%d", 11498039:11498046)
ctrl <- sprintf("SRR%d", 11498063:11498079)

prereg <- list(
  study = "isoformic uncertainty evidence (placenta null, replication, knockdown)",
  written = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  snapshot = list(
    git_ref = ref,
    git_sha = sha,
    tar_path = tar_path,
    tar_sha256 = tar_sha256,
    path = snap_dir,
    uncommitted_changes_in_R_at_snapshot = length(dirty),
    note = "git archive of the commit: uncommitted edits are not in the snapshot"
  ),
  seeds = list(
    base = 20261117L,
    rule = "seed of a unit = base + 1000 * package_index + unit_index",
    package_index = list(S0 = 0L, A = 1L, B = 2L, C = 3L, D = 4L, prior = 5L),
    placenta_split_sampling = 20261117L,
    sex_kmeans = 20261118L,
    validation_subsets = 20261119L
  ),
  data = list(
    placenta_quants = file.path(v50, "human", "placenta", "quants"),
    gtf = file.path(
      v50,
      "human",
      "ref",
      "gencode.v50.chr_patch_hapl_scaff.annotation.gtf.gz"
    ),
    preeclampsia = pe,
    control = ctrl,
    technical_pairs = list(
      p27 = c("SRR11498064", "SRR11498065"),
      p38 = c("SRR11498076", "SRR11498077")
    ),
    control_units = "17 controls minus the shallower member (fewer mapped fragments) of each A/B pair = 15 units; sensitivity: the other member",
    null_samples = "control units only; the 8 preeclampsia samples define U and the strata and are never tested in the null"
  ),
  haplotype_folding = list(
    primary = "sum counts and every draw of patch/haplotype copies onto the primary-chromosome transcript of the same gene_name and transcript_name; Salmon duplicate clusters retained on a non-primary copy are renamed to the primary member",
    sensitivity = "unfolded store; non-primary transcripts dropped"
  ),
  universe = list(
    name = "U",
    defined_from = "preeclampsia samples only (8)",
    contigs = "primary chromosomes chr1-22, chrX, chrY, chrM",
    cpm_floor = 1,
    cpm_library_size = "sum of folded NumReads of the sample",
    min_samples = 3L,
    rule = "folded CPM >= 1 in >= 3 of the 8 preeclampsia samples"
  ),
  strata = list(
    infrv = "mean over the 8 PE samples of per-sample InfRV (fishpond form, pc 5, shift 0.01) of the folded 30 draws; tertiles on U (low, mid, high) plus uncertain = InfRV >= 0.5",
    infrv_uncertain_cutoff = 0.5,
    identifiability = "equivalence_identifiability() on the 8 PE samples; non_identifiable = !identifiable (not identifiable in at least one PE sample with fragments)",
    unique_fraction = "equivalence_support() on the 8 PE samples; zero_unique = unique_fraction == 0",
    isoforms_per_gene = "annotated transcripts of the gene on primary chromosomes after folding: 1, 2-4, 5-9, 10+",
    expression_decile = "deciles of mean PE CPM on U",
    gene_level = c(
      "max isoform mean InfRV in the gene (on U)",
      "any non-identifiable isoform in U"
    ),
    robustness = "same strata from draws 1-15 only"
  ),
  methods = list(
    M1 = "edgeR QL (glmQLFit/glmQLFTest), transcript counts, no uncertainty",
    M2 = "edgeR QL on counts divided by the RTA overdispersion of the split samples (Baldoni 2024)",
    M3 = "limma-voom, transcript counts",
    M4 = "limma-voom, divided counts",
    M5 = "DESeq2 Wald, transcript counts (pilot splits only, stop rule)",
    M6 = "fishpond swish on the 30 draws (not at 3 vs 3: 20 permutations)",
    M7 = "M1 calls kept only when de_confidence class == high",
    M8 = "diffSplice (edgeR) usage, raw counts",
    M9 = "diffSplice (edgeR) usage, divided counts",
    M10 = "satuRn and DEXSeq usage (pilot splits only, stop rule)",
    M11 = "edgeR QL, gene counts (sum of isoforms), the gene-level comparator"
  ),
  thresholds = list(
    q = c(0.05, 0.10),
    p = c(0.01, 0.001),
    sign_agreement_unstable = 0.95,
    de_confidence = "isoformic snapshot defaults"
  ),
  splits = list(
    null_sizes = c(3L, 5L, 7L),
    null_splits_per_size = 100L,
    pilot_splits = 5L,
    pilot_fallback_splits_per_size = 50L,
    pilot_fallback_rule = "use 50 per size if the pilot projects > 24 h total under the lock",
    swish_sizes = c(5L, 7L),
    deseq2_dexseq_satuRn_splits = 10L,
    replication_half_splits = 30L,
    replication_directions = 2L,
    replication_directed_total = 60L,
    jackknife = "leave one control unit out; Monte Carlo correction of the split-sampling variance; 28-per-unit check"
  ),
  primary = list(
    flag = "de_confidence class != high",
    replication_metric = "a discovery replicates when the held-out half has the same sign and one-sided p < 0.05",
    replication_bootstrap = "two-level: half-splits, then transcripts within splits",
    controls = c(
      "placebo flag (random flag of the same size per split)",
      "count-matched flag (same mean-count distribution)",
      "permutation baseline of the labels"
    )
  ),
  decision_rules = list(
    R_A = "Ignoring uncertainty inflates false isoform calls: the M1 top-InfRV-tertile adjusted rate ratio (quasi-Poisson with expression decile and isoform class; CMH as a check) has a jackknife 95% CI excluding 1, AND M2 or M6 reduces it with a CI excluding 0, AND M11 (gene level) shows no comparable excess.",
    R_B = "Flags predict non-replication: the adjusted OR of the primary flag has a 95% CI excluding 1 in the placenta half-splits, AND the flag beats the count-matched flag.",
    negative_results = "A failed rule is reported as a headline, including the Baldoni 2024 pattern (raw counts calibrated, gain only in power), an adjusted OR near 1, and an InfRV >= 0.5 flag that is too broad."
  ),
  story_genes = list(
    isoform_level = "FLT1",
    gene_level = c("ENG", "LEP"),
    fragile_rule = "among genes with an M1 q < 0.05 isoform in the full PE vs control comparison, the gene with the most such isoforms whose calls all have the primary flag and replicate in < 50% of directed half-splits (ties: smallest min q)",
    robust_rule = "among genes with M11 q < 0.05, the gene whose top isoform has class high and replicates in >= 90% of directed half-splits, with the smallest M11 q"
  ),
  stop_rules = list(
    dexseq_minutes_per_split = 20,
    deseq2_minutes_per_split = 5,
    measured_by = "pilot"
  ),
  compute = list(
    lock = "shared lock file of the heavy jobs (machine specific)",
    lock_rule = "> 8 GB RAM, or > 4 cores for > 10 min",
    threads_outside_lock = 4L
  )
)
yaml::write_yaml(prereg, prereg_path)
file.copy(
  prereg_path,
  "data-raw/uncertainty_evidence/prereg.yaml",
  overwrite = TRUE
)
cat(
  "snapshot:",
  snap_dir,
  "\nsha256:",
  tar_sha256,
  "\nprereg:",
  prereg_path,
  "\n"
)
