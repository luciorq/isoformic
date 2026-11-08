# Extract Transcript Sequences into an Annotation Directory

Adds the spliced sequence of every transcript of an annotation directory
as `annot_data_sequences.parquet`, read by
[`transcript_orfs()`](https://luciorq.github.io/isoformic/reference/transcript_orfs.md),
[`coding_potential()`](https://luciorq.github.io/isoformic/reference/transcript_orfs.md)
and the sequence-aware rules of
[`transcript_productivity()`](https://luciorq.github.io/isoformic/reference/productivity.md).
Sequences come either from a **genome FASTA** (`genome_fasta`; exons are
cut out of each chromosome, reverse complemented on the minus strand and
concatenated in transcript order, so novel transcripts from long-read or
assembly annotations are covered) or from a **transcript FASTA**
(`transcript_fasta`, for example the GENCODE `transcripts.fa.gz`;
sequences are matched by the first token of the header, with and without
the identifier version).

The genome is streamed in chunks (plain or gzip-compressed FASTA, no
index needed) and only one chromosome is held in memory at a time; a
human genome takes a few minutes. Chromosome names are matched to the
annotation `seqid` directly, or after adding or removing a `chr` prefix.
Soft-masked (lower-case) bases are upper-cased.

## Usage

``` r
extract_transcript_sequences(
  x,
  genome_fasta = NULL,
  transcript_fasta = NULL,
  overwrite = FALSE,
  chunk_bytes = 6.4e+07
)

has_sequences(x)

transcript_sequences(x, transcript_id = NULL)
```

## Arguments

- x:

  An annotation directory path or an
  [IsoformicExperiment](https://luciorq.github.io/isoformic/reference/IsoformicExperiment.md).

- genome_fasta:

  Path to a genome FASTA file (`.fa`, `.fa.gz`), see
  `download_reference(file_type = "genome_fasta")`.

- transcript_fasta:

  Path to a transcript FASTA file, see
  `download_reference(file_type = "fasta")`.

- overwrite:

  Logical, replace existing sequences (default `FALSE`).

- chunk_bytes:

  Number of bytes read per chunk while streaming (default 64 MB).

- transcript_id:

  Optional transcript identifiers to read; `NULL` reads all sequences.

## Value

The path of the written Parquet file, invisibly. The table has columns
`transcript_id`, `sequence` and `seq_length`.

## Details

The sequences table is optional: annotation directories without it keep
working, and `has_sequences()` tells whether it is present. Transcripts
whose chromosome is absent from the FASTA (or whose identifier is absent
from a transcript FASTA) are reported and left out.

## See also

[`transcript_orfs()`](https://luciorq.github.io/isoformic/reference/transcript_orfs.md),
[`coding_potential()`](https://luciorq.github.io/isoformic/reference/transcript_orfs.md)

## Examples

``` r
gff <- system.file("extdata", "example.gff3", package = "isoformic")
fasta <- system.file("extdata", "example_transcripts.fa.gz", package = "isoformic")
annot <- prepare_isoformic_annotation(gff, output_path = tempdir(), annot_name = "seq_demo")
#> ℹ Parsing example.gff3 (GFF) into /tmp/RtmpV2WK1u/seq_demo.
extract_transcript_sequences(annot, transcript_fasta = fasta, overwrite = TRUE)
#> ✔ Wrote 72 transcript sequences to
#>   /tmp/RtmpV2WK1u/seq_demo/annot_data_sequences.parquet.
has_sequences(annot)
#> [1] TRUE
head(transcript_sequences(annot)[, c("transcript_id", "seq_length")])
#> # A tibble: 6 × 2
#>   transcript_id     seq_length
#>   <chr>                  <int>
#> 1 ENST00000282397.9       7123
#> 2 ENST00000615611.4       1476
#> 3 ENST00000617835.4        434
#> 4 ENST00000543394.2       1275
#> 5 ENST00000540678.2       1927
#> 6 ENST00000541932.5       2969

if (FALSE) { # \dontrun{
# From a genome FASTA (covers novel transcripts of any annotation)
genome <- download_reference(version = "49", file_type = "genome_fasta", region = "primary")
extract_transcript_sequences(annot, genome_fasta = genome)
} # }
```
