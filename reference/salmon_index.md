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

  Logical indicating if the FASTA is from GENCODE (default is `FALSE`);
  the transcript names are then truncated at the first `|` so that they
  match the annotation.

- decoy_fasta:

  Optional path to a FASTA file containing decoy sequences (default is
  `NULL`).

- clip_poly_a:

  Logical indicating whether to clip poly-A tails (default is `TRUE`).

## Value

processx style output list.

## Details

The index must be built with the same Salmon major version used for
[`salmon_quant()`](https://luciorq.github.io/isoformic/reference/salmon_quant.md):
Salmon 2 (the Rust rewrite) cannot read indices built by Salmon 1.x and
rejects them with an error. Rebuild the index after upgrading.
