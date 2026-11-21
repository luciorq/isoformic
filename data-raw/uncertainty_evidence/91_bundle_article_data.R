# 91_bundle_article_data.R: copy the article tables into the package
# repository, with evidence_numbers.parquet pruned to the keys the article
# uses.
#
# Purpose: vignettes/articles/isoformic-uncertainty-evidence.qmd reads four
# small tables from vignettes/articles/data/uncertainty_evidence/. This
# script copies them from $ISOFORMIC_UE_OUT/article_tables/ (written by
# 90_article_tables.R) without changing any value. evidence_numbers.parquet
# keeps only the keys that the two articles name literally, plus the keys that
# its figures and tables build in loops (the patterns below); the article
# stops on a missing key, so rendering it checks that nothing was pruned
# by mistake.
#
# Inputs: $ISOFORMIC_UE_OUT/article_tables/{evidence_numbers,
#   flt1_isoforms, story_genes, replication_flags}.parquet;
#   vignettes/articles/isoformic-uncertainty-evidence.qmd and
#   isoformic-uncertainty.qmd (the keys they name)
# Outputs: vignettes/articles/data/uncertainty_evidence/ (the same four
#   files; Parquet through DuckDB, zstd)
# Env vars: ISOFORMIC_UE_OUT (required).
# Resources: seconds, < 1 GB. It reads the tables with the installed or
# loaded isoformic (read_parquet(), write_parquet()), not the snapshot.
# Run from the repository root:
#   pixi run --frozen Rscript data-raw/uncertainty_evidence/91_bundle_article_data.R

out_root <- Sys.getenv("ISOFORMIC_UE_OUT")
if (!nzchar(out_root)) {
  stop("set ISOFORMIC_UE_OUT", call. = FALSE)
}
src <- file.path(out_root, "article_tables")
dest <- file.path("vignettes", "articles", "data", "uncertainty_evidence")
# the evidence article, and the uncertainty tutorial, which quotes a few of
# the same numbers in its pointers to it
articles <- file.path(
  "vignettes",
  "articles",
  c("isoformic-uncertainty-evidence.qmd", "isoformic-uncertainty.qmd")
)
dir.create(dest, recursive = TRUE, showWarnings = FALSE)

ev <- isoformic::read_parquet(file.path(src, "evidence_numbers.parquet"))
txt <- paste(unlist(lapply(articles, readLines)), collapse = "\n")
lit <- regmatches(txt, gregexpr('"(design|S0|A|B|C|D|R)\\.[^"]+"', txt))[[1]]
lit <- unique(gsub('"', "", lit))
# keys built in loops: the null rate-ratio figure (raw jackknife, placenta;
# mouse methods and subsets), the thinning tables, and the supplementary
# tables of every null method and of the robustness variants
patterns <- c(
  "^A\\.adj_rr_rawjk\\.[357]v[357]\\.M[126]\\.prereg\\.top$",
  "^C\\.M[126]\\.adjusted_rate_ratio_high_low(_batch_balanced)?$",
  "^A\\.thin\\.(7v7_louo|5v5)\\.(M1|M2|M6|M9a|M10a)\\.(all|top_tertile)\\.q_0\\.05\\.(fdr_pooled|total_calls|total_false)$",
  "^A\\.(mean_calls_q05|inflation_p01)\\.[357]v[357]\\.M[0-9]+[a-z]*\\.own\\.all$",
  "^A\\.adj_rr\\.[357]v[357]\\.M[0-9]+\\.prereg\\.top$",
  "^A\\.robust\\.[^.]+\\.M1\\.(all\\.all|chr_xy\\.FALSE|infrv_tertile\\.(high|low))$",
  # the od-numbers chunk: raw against divided counts and the prior swaps
  paste0(
    "^D\\.impact\\.(",
    "placenta_8_PE_vs_15_control_units_DTE_edgeR_QL\\.raw_counts_-_divided_counts_genome-wide_prior\\.(separate|common)",
    "|GSE(129540|75489)_DTE_edgeR_QL_divided_counts\\.panel_prior_-_genome_prior_swapped\\.common",
    "|placenta_8_PE_vs_15_control_units_DTE_edgeR_QL_divided_counts\\.genome_prior_-_panel_median\\.common",
    ")\\.[a-z_]+$"
  )
)
from_patterns <- ev$key[Reduce(`|`, lapply(patterns, grepl, x = ev$key))]
# helpers in the article that build keys from a short name
helper_keys <- c(
  paste0("A.adj_rr.", c("3v3.M1", "5v5.M1", "7v7.M1", "7v7.M6"), ".prereg.top"),
  paste0("A.adj_rr_rawjk.", c("5v5", "7v7"), ".M1.prereg.top"),
  paste0(
    "B.cluster.cluster.",
    c(
      "gene.unstable_or",
      "split_x_gene.unstable_or",
      "gene.unstable_or_signconc",
      "gene.primary_or",
      "gene.primary_minus_c4"
    )
  )
)
keep <- unique(c(intersect(lit, ev$key), from_patterns, helper_keys))
missing <- setdiff(helper_keys, ev$key)
if (length(missing)) {
  stop(
    "keys missing from the evidence table: ",
    paste(missing, collapse = ", "),
    call. = FALSE
  )
}
pruned <- ev[ev$key %in% keep, ]
pruned <- pruned[order(match(pruned$key, ev$key)), ]
isoformic::write_parquet(pruned, file.path(dest, "evidence_numbers.parquet"))
for (f in c(
  "flt1_isoforms.parquet",
  "story_genes.parquet",
  "replication_flags.parquet"
)) {
  isoformic::write_parquet(
    isoformic::read_parquet(file.path(src, f)),
    file.path(dest, f)
  )
}
message(sprintf(
  "91: kept %d of %d evidence keys; %s KB in %s",
  nrow(pruned),
  nrow(ev),
  format(round(sum(file.size(list.files(dest, full.names = TRUE))) / 1024, 1)),
  dest
))
