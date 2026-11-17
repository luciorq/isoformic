# Simulate Oarfish input: ONT-like cDNA reads from the bundled transcripts
#
# Three samples (two controls, one treated) of 3,000 reads each, drawn from
# the 141 transcripts of inst/extdata/example_transcripts.fa.gz with random 5'
# truncation and a 6% error rate (substitutions, insertions, deletions). The
# treated sample raises three isoforms. The reads only exist to give Oarfish
# something to quantify; they are not a model of any real library.
#
# Usage: Rscript data-raw/oarfish/simulate_long_reads.R OUTDIR
args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args)) args[[1L]] else "."
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

fa <- readLines(gzfile("inst/extdata/example_transcripts.fa.gz"))
is_header <- startsWith(fa, ">")
grp <- cumsum(is_header)
tx <- split(fa, grp)
names(tx) <- vapply(tx, function(x) sub("^>", "", x[[1L]]), "")
seqs <- vapply(tx, function(x) paste(x[-1L], collapse = ""), "")

set.seed(2026)
n_tx <- length(seqs)
base_weight <- stats::rlnorm(n_tx, meanlog = 0, sdlog = 1.2)
boosted <- c(2L, 5L, 9L)
samples <- list(control_1 = 1, control_2 = 1, treated_1 = 1)

mutate_read <- function(s) {
  b <- strsplit(s, "")[[1L]]
  u <- stats::runif(length(b))
  bases <- c("A", "C", "G", "T")
  sub <- u < 0.03
  b[sub] <- sample(bases, sum(sub), replace = TRUE)
  del <- u >= 0.03 & u < 0.045
  ins <- u >= 0.045 & u < 0.06
  b[ins] <- paste0(b[ins], sample(bases, sum(ins), replace = TRUE))
  paste(b[!del], collapse = "")
}

for (nm in names(samples)) {
  w <- base_weight * stats::rlnorm(n_tx, 0, 0.15)
  if (startsWith(nm, "treated")) {
    w[boosted] <- w[boosted] * 4
  }
  # Reads per transcript are proportional to abundance and length-independent
  n_reads <- 3000L
  pick <- sample(n_tx, n_reads, replace = TRUE, prob = w)
  lines <- character(2L * n_reads)
  for (i in seq_len(n_reads)) {
    s <- seqs[[pick[i]]]
    len <- nchar(s)
    from <- 1L + floor(stats::runif(1L, 0, 0.3) * len)
    read <- mutate_read(substr(s, from, len))
    lines[2L * i - 1L] <- paste0("@", nm, "_read", i)
    lines[2L * i] <- read
  }
  fq <- character(4L * n_reads)
  fq[seq(1L, by = 4L, length.out = n_reads)] <- lines[seq(
    1L,
    by = 2L,
    length.out = n_reads
  )]
  fq[seq(2L, by = 4L, length.out = n_reads)] <- lines[seq(
    2L,
    by = 2L,
    length.out = n_reads
  )]
  fq[seq(3L, by = 4L, length.out = n_reads)] <- "+"
  fq[seq(4L, by = 4L, length.out = n_reads)] <- vapply(
    nchar(lines[seq(2L, by = 2L, length.out = n_reads)]),
    function(n) strrep("5", n),
    ""
  )
  writeLines(fq, gzfile(file.path(out_dir, paste0(nm, ".fastq.gz"))))
  message(nm, ": ", n_reads, " reads")
}
