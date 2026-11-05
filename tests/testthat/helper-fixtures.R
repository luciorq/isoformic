# Shared fixtures for the test suite --------------------------------------------

mini_annotation <- function(file = "mini.gff", envir = parent.frame()) {
  out <- withr::local_tempdir(.local_envir = envir)
  suppressMessages(prepare_isoformic_annotation(
    input_path = testthat::test_path("testdata", file),
    output_path = out
  ))
}

mini_assay <- function() {
  matrix(
    c(10, 20, 30, 40, 50, 60, 1, 2, 3, 4, 6, 8),
    nrow = 3,
    dimnames = list(c("tx1", "tx2", "tx3"), c("s1", "s2", "s3", "s4"))
  )
}

mini_col_data <- function() {
  data.frame(
    sample_id = c("s1", "s2", "s3", "s4"),
    condition = c("a", "a", "b", "b"),
    stringsAsFactors = FALSE
  )
}

mini_experiment <- function(envir = parent.frame(), with_dea = TRUE) {
  iso <- IsoformicExperiment(
    experiment_name = "mini",
    annot_path = mini_annotation(envir = envir),
    assays = list(abundance = mini_assay(), counts = mini_assay() * 2),
    col_data = mini_col_data()
  )
  if (isTRUE(with_dea)) {
    de_tx(iso) <- data.frame(
      transcript_id = c("tx1", "tx2", "tx3"),
      log2FoldChange = c(2, -1.5, 0.1),
      pvalue = c(0.001, 0.01, 0.5),
      padj = c(0.003, 0.02, 0.6)
    )
    de_gene(iso) <- data.frame(
      gene_id = c("gene1", "gene2"),
      log2FC = c(1.2, 0),
      pvalue = c(0.01, 0.9),
      qvalue = c(0.02, 0.95)
    )
  }
  iso
}

# Deterministic synthetic genome matching testdata/mini_cds.gff3 ---------------
#
# chr1 (25 kb). The canonical isoform tA1 carries an ORF from genomic 1100
# (ATG) to 8300 (TAA) built from T-free codons so that no stop codon can occur
# in any frame inside its exons. A uORF (ATG GCA GCA TAA) sits at 1010-1021.
# The tA3 cassette exon (4000-4100) and the tA4 retained intron (3201-4999)
# are T-free up to a planted in-frame TAA at 4046-4048 and 3301-3303.
mini_genome_seq <- function() {
  set.seed(42)
  acgt <- c("A", "C", "G", "T")
  acg <- c("A", "C", "G")
  g <- sample(acgt, 25000, replace = TRUE)
  set_seq <- function(g, start, chars) {
    g[start:(start + length(chars) - 1L)] <- chars
    g
  }
  # 5' UTR of exon 1 (1000-1099) without T except the planted uORF
  g <- set_seq(g, 1000, sample(acg, 100, replace = TRUE))
  g <- set_seq(g, 1010, strsplit("ATGGCAGCATAA", "")[[1L]])
  # Canonical CDS exons: T-free, ATG at 1100, TAA at 8298-8300
  for (r in list(c(1100, 1200), c(3000, 3200), c(5000, 5200), c(8000, 8300))) {
    g <- set_seq(
      g,
      r[[1L]],
      sample(acg, r[[2L]] - r[[1L]] + 1L, replace = TRUE)
    )
  }
  g <- set_seq(g, 1100, c("A", "T", "G"))
  g <- set_seq(g, 8298, c("T", "A", "A"))
  # tA3 cassette exon: T-free, in-frame TAA at 4046
  g <- set_seq(g, 4000, sample(acg, 101, replace = TRUE))
  g <- set_seq(g, 4046, c("T", "A", "A"))
  # tA4 retained intron: T-free until the planted in-frame TAA at 3301
  g <- set_seq(g, 3201, sample(acg, 100, replace = TRUE))
  g <- set_seq(g, 3301, c("T", "A", "A"))
  paste(g, collapse = "")
}

mini_genome_fasta <- function(
  envir = parent.frame(),
  seq_name = "chr1",
  gzip = FALSE
) {
  seq <- mini_genome_seq()
  lines <- substring(
    seq,
    seq(1, nchar(seq), by = 60),
    pmin(seq(60, nchar(seq) + 59, by = 60), nchar(seq))
  )
  path <- withr::local_tempfile(
    fileext = if (gzip) ".fa.gz" else ".fa",
    .local_envir = envir
  )
  con <- if (gzip) gzfile(path, "w") else file(path, "w")
  writeLines(
    c(paste0(">", seq_name, " synthetic"), lines, ">chrUn_decoy", "ACGT"),
    con
  )
  close(con)
  path
}

# Synthetic Salmon output with Gibbs samples and equivalence classes -----------
#
# Transcripts of testdata/mini_apa.gff3 in six samples (s1-s3 group "a", s4-s6
# group "b"). Gene gA: tA1 and tA2 share their reads (draws sum to a constant,
# so their correlation is -1) and swap expression between the groups; tA3 is
# well determined and four times higher in b; tA4 is strongly overdispersed;
# tA5 is not expressed; tA6 is stable. Gene gB: tB1 and tB2 are independent.
salmon_uncertainty_fixture <- function(envir = parent.frame(), n_draws = 30L) {
  withr::local_seed(11L, .local_envir = envir)
  root <- withr::local_tempdir(.local_envir = envir)
  tx <- c("tA1", "tA2", "tA3", "tA4", "tA5", "tA6", "tB1", "tB2")
  samples <- paste0("s", 1:6)
  in_b <- samples %in% c("s4", "s5", "s6")
  draws <- array(
    0,
    dim = c(length(tx), length(samples), n_draws),
    dimnames = list(tx, samples, NULL)
  )
  for (s in seq_along(samples)) {
    noise <- exp(stats::rnorm(1L, 0, 0.1))
    p0 <- if (in_b[s]) 0.7 else 0.3
    p <- stats::rbeta(n_draws, 20 * p0, 20 * (1 - p0))
    draws["tA1", s, ] <- 200 * p
    draws["tA2", s, ] <- 200 * (1 - p)
    draws["tA3", s, ] <- stats::rpois(
      n_draws,
      (if (in_b[s]) 80 else 20) * noise
    )
    draws["tA4", s, ] <- stats::rnbinom(n_draws, size = 3, mu = 30 * noise)
    draws["tA5", s, ] <- stats::rpois(n_draws, 0.2)
    draws["tA6", s, ] <- stats::rpois(n_draws, 50 * noise)
    draws["tB1", s, ] <- stats::rpois(n_draws, 40 * noise)
    draws["tB2", s, ] <- stats::rpois(n_draws, 40 * noise)
  }
  mapped <- c(90, 88, 30, 91, 89, 92)
  dirs <- character()
  for (s in seq_along(samples)) {
    dir <- fs::path(root, paste0(samples[s], "_quant"))
    aux <- fs::path(dir, "aux_info")
    boot <- fs::path(aux, "bootstrap")
    fs::dir_create(boot)
    counts <- rowMeans(draws[, s, ])
    utils::write.table(
      data.frame(
        Name = tx,
        Length = 1000L,
        EffectiveLength = 800,
        TPM = counts / sum(counts) * 1e6,
        NumReads = counts
      ),
      fs::path(dir, "quant.sf"),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
    names_con <- gzfile(fs::path(boot, "names.tsv.gz"), "wb")
    writeBin(charToRaw(paste0(paste(tx, collapse = "\t"), "\n")), names_con)
    close(names_con)
    boot_con <- gzfile(fs::path(boot, "bootstraps.gz"), "wb")
    writeBin(as.vector(draws[, s, ]), boot_con, size = 8L, endian = "little")
    close(boot_con)
    writeLines(
      c(
        "{",
        '    "salmon_version": "2.8.0",',
        '    "samp_type": "gibbs",',
        '    "opt_type": "vb",',
        '    "library_types": [',
        '        "IU"',
        "    ],",
        '    "detected_library_type": "IU",',
        '    "seq_bias_correct": true,',
        '    "gc_bias_correct": true,',
        '    "pos_bias_correct": true,',
        '    "frag_length_mean": 250.5,',
        '    "frag_length_sd": 50.1,',
        paste0('    "num_valid_targets": ', length(tx), ","),
        '    "em_converged": true,',
        '    "num_processed": 1000,',
        paste0('    "num_mapped": ', mapped[s] * 10, ","),
        paste0('    "percent_mapped": ', mapped[s], ","),
        paste0('    "num_bootstraps": ', n_draws, ","),
        '    "call": "quant"',
        "}"
      ),
      fs::path(aux, "meta_info.json")
    )
    eq <- c(
      as.character(length(tx)),
      "8",
      tx,
      "2\t0\t1\t120",
      "1\t0\t40",
      "1\t1\t40",
      "1\t2\t20",
      "2\t2\t3\t0.5\t0.5\t10",
      "1\t3\t30",
      "1\t5\t50",
      "1\t6\t40",
      "1\t7\t40"
    )
    eq_con <- gzfile(fs::path(aux, "eq_classes.txt.gz"), "wt")
    writeLines(eq, eq_con)
    close(eq_con)
    dirs[samples[s]] <- as.character(dir)
  }
  list(
    dirs = dirs,
    draws = draws,
    samples = samples,
    col_data = data.frame(
      sample_id = samples,
      condition = ifelse(in_b, "b", "a")
    )
  )
}

# An IsoformicExperiment on the synthetic Salmon output, with DE results
uncertainty_experiment <- function(
  envir = parent.frame(),
  inferential = FALSE
) {
  fx <- salmon_uncertainty_fixture(envir = envir)
  iso <- read_salmon(
    fx$dirs,
    annot_path = mini_annotation("mini_apa.gff3", envir = envir),
    col_data = fx$col_data,
    inferential = inferential
  )
  de_tx(iso) <- data.frame(
    transcript_id = c("tA1", "tA2", "tA3", "tB1"),
    log2FC = c(1.2, -1.2, 1.8, 0.1),
    pvalue = c(0.001, 0.001, 0.0001, 0.6),
    qvalue = c(0.01, 0.01, 0.001, 0.7)
  )
  de_gene(iso) <- data.frame(
    gene_id = c("gA", "gB"),
    log2FC = c(0.1, 0),
    pvalue = c(0.6, 0.9),
    qvalue = c(0.7, 0.95)
  )
  list(iso = iso, fixture = fx)
}

# Many genes with three isoforms each, for tests of usage; genes g1 to g8 switch
# between two isoforms in group b. With three groups, genes g9 to g16 switch
# between two other isoforms in group c, and group c is like a for g1 to g8.
usage_experiment <- function(
  n_genes = 60L,
  n_switch = 8L,
  n_groups = 2L,
  envir = parent.frame()
) {
  gff <- withr::local_tempfile(fileext = ".gff3", .local_envir = envir)
  lines <- "##gff-version 3"
  for (g in seq_len(n_genes)) {
    start <- 1000L + (g - 1L) * 20000L
    lines <- c(
      lines,
      sprintf(
        "chr1\tT\tgene\t%d\t%d\t.\t+\t.\tID=g%d;gene_id=g%d;gene_name=GENE%d;gene_type=protein_coding;",
        start,
        start + 9000L,
        g,
        g,
        g
      )
    )
    for (k in 1:3) {
      id <- sprintf("g%dt%d", g, k)
      lines <- c(
        lines,
        sprintf(
          paste0(
            "chr1\tT\ttranscript\t%d\t%d\t.\t+\t.\tID=%s;Parent=g%d;gene_id=g%d;",
            "transcript_id=%s;transcript_type=protein_coding;",
            "transcript_name=GENE%d-20%d;"
          ),
          start,
          start + 9000L,
          id,
          g,
          g,
          id,
          g,
          k
        ),
        sprintf(
          "chr1\tT\texon\t%d\t%d\t.\t+\t.\tID=%se1;Parent=%s;transcript_id=%s;exon_number=1;",
          start,
          start + 200L,
          id,
          id,
          id
        ),
        sprintf(
          "chr1\tT\texon\t%d\t%d\t.\t+\t.\tID=%se2;Parent=%s;transcript_id=%s;exon_number=2;",
          start + 1000L * k,
          start + 1000L * k + 200L,
          id,
          id,
          id
        )
      )
    }
  }
  writeLines(lines, gff)
  annot <- withr::local_tempdir(.local_envir = envir)
  annot <- suppressMessages(
    prepare_isoformic_annotation(input_path = gff, output_path = annot)
  )
  ids <- as.vector(outer(seq_len(n_genes), 1:3, function(g, k) {
    sprintf("g%dt%d", g, k)
  }))
  n_samples <- 3L * n_groups
  samples <- paste0("s", seq_len(n_samples))
  set.seed(42)
  base <- rep(c(200, 100, 50), each = n_genes)
  mu <- matrix(base, nrow = length(ids), ncol = n_samples)
  # In the switching genes, the second isoform loses its share to the third
  switching <- ids %in%
    c(
      sprintf("g%dt2", seq_len(n_switch)),
      sprintf("g%dt3", seq_len(n_switch))
    )
  mu[ids %in% sprintf("g%dt2", seq_len(n_switch)), 4:6] <- 25
  mu[ids %in% sprintf("g%dt3", seq_len(n_switch)), 4:6] <- 125
  c_genes <- n_switch + seq_len(n_switch)
  if (n_groups >= 3L) {
    # In group c, the first isoform loses its share to the second
    mu[ids %in% sprintf("g%dt1", c_genes), 7:9] <- 60
    mu[ids %in% sprintf("g%dt2", c_genes), 7:9] <- 240
  }
  counts <- matrix(
    stats::rnbinom(length(mu), mu = mu, size = 30),
    nrow = length(ids),
    dimnames = list(ids, samples)
  )
  iso <- IsoformicExperiment(
    experiment_name = "usage",
    annot_path = annot,
    assays = list(
      counts = tibble::as_tibble(counts, rownames = "transcript_id")
    ),
    col_data = data.frame(
      sample_id = samples,
      condition = rep(letters[seq_len(n_groups)], each = 3)
    )
  )
  list(
    iso = iso,
    switching = switching,
    switch_genes = sprintf("g%d", seq_len(n_switch)),
    c_genes = sprintf("g%d", c_genes)
  )
}
