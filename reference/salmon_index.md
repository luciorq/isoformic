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

  Optional path to a FASTA file, plain or compressed with gzip, xz or
  bzip2, with the decoy sequences (default is `NULL`, an index without
  decoys). Usually the genome assembly the transcripts come from (for
  GENCODE, the `GRCh38.primary_assembly.genome.fa.gz` of the same
  release). A full human genome makes the index several times larger and
  slower to build (tens of GB of memory and disk instead of a few), and
  needs free space in
  [`tempdir()`](https://rdrr.io/r/base/tempfile.html) for a copy of the
  genome while the index is built.

- clip_poly_a:

  Logical indicating whether to clip poly-A tails (default is `TRUE`).

## Value

processx style output list.

## Details

The index must be built with the same Salmon major version used for
[`salmon_quant()`](https://luciorq.github.io/isoformic/reference/salmon_quant.md):
Salmon 2 (the Rust rewrite) cannot read indices built by Salmon 1.x and
rejects them with an error. Rebuild the index after upgrading.

Transcripts with identical sequences are all kept (`--keepDuplicates`),
so that every annotated transcript is quantified.

With `decoy_fasta`, the index is decoy-aware, as recommended by the
Salmon authors: the function writes, in a temporary directory, a
`gentrome` FASTA (the transcripts followed by the decoy sequences) and a
`decoys.txt` file with the names of the decoy sequences (the first word
of each FASTA header), passes them to `salmon index` as `--transcripts`
and `--decoys`, and removes them afterwards. Reads that align better to
a decoy than to a transcript are then discarded instead of being
assigned to a transcript. The `gentrome` is gzipped when `decoy_fasta`
is (the decoy file is appended byte by byte, and concatenated gzip
streams are valid gzip), so a gzipped genome is not decompressed to
disk; the transcript FASTA, plain or gzipped, is re-encoded to match. A
decoy compressed with xz or bzip2 is decompressed into a plain
`gentrome`. Decoy names that are also transcript names are an error.
