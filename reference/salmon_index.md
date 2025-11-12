# Build Salmon Index

Create a Salmon index from a reference transcriptome FASTA file.

## Usage

``` r
salmon_index(
  fasta_path,
  index_path = "salmon_index",
  kmer_len = 31,
  num_threads = 2,
  env_name = "salmon-env",
  is_gencode = FALSE,
  decoy_fasta = NULL,
  clip_poly_a = TRUE
)
```

## Arguments

- fasta_path:

  Path to the reference transcriptome FASTA file.

- index_path:

  Directory path to save the Salmon index (default is "salmon_index").

- kmer_len:

  K-mer length for the index (default is 31).

- num_threads:

  Number of threads to use (default is 2).

- env_name:

  Name of the conda environment with Salmon installed (default is
  "salmon-env").

- is_gencode:

  Logical indicating if the FASTA is from GENCODE (default is `FALSE`).

- decoy_fasta:

  Optional path to a FASTA file containing decoy sequences (default is
  `NULL`).

- clip_poly_a:

  Logical indicating whether to clip poly-A tails (default is `TRUE`).

## Value

processx style output list.
