# Simulate barcoded long reads for Oarfish single-cell mode
#
# Sixty cells in three samples (two controls and one treated, 20 cells each)
# of 150 reads each, drawn from the 141 transcripts of
# inst/extdata/example_transcripts.fa.gz like simulate_long_reads.R does. The
# cell barcode is the first field of the read name ("<barcode>_<n>"); the run
# script moves it to the CB:Z tag that Oarfish reads. The treated cells raise
# three isoforms. A table of cells and samples is written next to the reads.
#
# Usage: Rscript data-raw/oarfish/simulate_single_cell.R OUTDIR
args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args)) args[[1L]] else "."
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

fa <- readLines(gzfile("inst/extdata/example_transcripts.fa.gz"))
grp <- cumsum(startsWith(fa, ">"))
tx <- split(fa, grp)
seqs <- vapply(tx, function(x) paste(x[-1L], collapse = ""), "")

set.seed(2027)
n_tx <- length(seqs)
base_weight <- stats::rlnorm(n_tx, meanlog = 0, sdlog = 1.2)
boosted <- c(2L, 5L, 9L)
samples <- c("control_1", "control_2", "treated_1")
bases <- c("A", "C", "G", "T")
mutate_read <- function(s) {
  b <- strsplit(s, "")[[1L]]
  u <- stats::runif(length(b))
  sub <- u < 0.03
  b[sub] <- sample(bases, sum(sub), replace = TRUE)
  del <- u >= 0.03 & u < 0.045
  ins <- u >= 0.045 & u < 0.06
  b[ins] <- paste0(b[ins], sample(bases, sum(ins), replace = TRUE))
  paste(b[!del], collapse = "")
}

cells <- data.frame(
  barcode = vapply(
    seq_len(60L),
    function(i) paste(sample(bases, 16L, replace = TRUE), collapse = ""),
    ""
  ),
  sample_id = rep(samples, each = 20L)
)
fq <- character()
for (i in seq_len(nrow(cells))) {
  w <- base_weight * stats::rlnorm(n_tx, 0, 0.3)
  if (startsWith(cells$sample_id[[i]], "treated")) {
    w[boosted] <- w[boosted] * 4
  }
  pick <- sample(n_tx, 150L, replace = TRUE, prob = w)
  for (k in seq_along(pick)) {
    s <- seqs[[pick[k]]]
    len <- nchar(s)
    read <- mutate_read(substr(
      s,
      1L + floor(stats::runif(1L, 0, 0.3) * len),
      len
    ))
    fq <- c(
      fq,
      paste0("@", cells$barcode[[i]], "_", k),
      read,
      "+",
      strrep("5", nchar(read))
    )
  }
}
writeLines(fq, gzfile(file.path(out_dir, "cells.fastq.gz")))
utils::write.table(
  cells,
  file.path(out_dir, "cells.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
message(nrow(cells), " cells, ", length(fq) / 4L, " reads")
