# Run Salmon Quantification

Perform transcript quantification using Salmon's
selective-alignment-based mode from raw RNA-seq reads.

## Usage

``` r
salmon_quant(
  input_r1,
  input_r2 = NULL,
  index_path = "salmon_index",
  output_dir = "quant_output",
  num_threads = 8,
  num_gibbs = 30,
  min_score_fraction = "0.65",
  env_name = "salmon-env"
)
```

## Arguments

- input_r1:

  Path to the FASTQ file for read 1 (or single-end reads).

- input_r2:

  Optional path to the FASTQ file for read 2 (paired-end reads).

- index_path:

  Path to the Salmon index directory (default is "salmon_index").

- output_dir:

  Directory to save the quantification output (default is
  "quant_output").

- num_threads:

  Number of threads to use (default is 8).

- num_gibbs:

  Number of Gibbs samples (inferential replicates) drawn per sample
  (default is 30). They quantify how well each transcript abundance is
  determined by the reads and are used by
  [`add_uncertainty()`](https://luciorq.github.io/isoformic/reference/uncertainty.md)
  and by `swish`; 30 is enough for those summaries, and more only costs
  time and disk space. Set it to `0` to skip them.

- min_score_fraction:

  Minimum score fraction for alignments (default is "0.65").

- env_name:

  Name of the conda environment with Salmon installed (default is
  "salmon-env").

## Value

processx style output list.

## Details

The command enables the sequence, GC and positional bias corrections,
allows dovetailing mates and soft-clipped overhangs, writes the
equivalence classes with their weights (`--dumpEq --dumpEqWeights`) and
draws `num_gibbs` Gibbs samples (`aux_info/bootstrap/`) so that
inferential replicates are available for `swish`. Selective alignment
and the variational Bayes optimiser are Salmon defaults and are not
passed explicitly; Salmon 2 quantifies deterministically (identical
output for any thread count). Single-end reads are passed with
`--unmatedReads`.
