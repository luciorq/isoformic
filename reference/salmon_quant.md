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
  num_gibbs = 100,
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

  Number of Gibbs samples for uncertainty estimation (default is 100).

- min_score_fraction:

  Minimum score fraction for alignments (default is "0.65").

- env_name:

  Name of the conda environment with Salmon installed (default is
  "salmon-env").

## Value

processx style output list.
