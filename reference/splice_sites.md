# Splice-Site Strength from the Annotation and the Genome

Scores the donor (5') and acceptor (3') splice site of every intron with
position weight matrices built from the annotation itself, and reports
the weakest site of every transcript. A weak site is one way for an
isoform to be rare or to be a splicing error; the score is evidence to
read next to
[`transcript_productivity()`](https://luciorq.github.io/isoformic/reference/productivity.md),
not a classifier.

- `extract_splice_sites(x, genome_fasta)` cuts the sequences around the
  splice sites out of a genome FASTA and stores them in the annotation
  directory.
  [`extract_transcript_sequences()`](https://luciorq.github.io/isoformic/reference/transcript_sequences.md)
  does the same when it is given a genome FASTA, so this is only needed
  for directories whose sequences came from a transcript FASTA.

- `splice_site_strength(x)` returns one row per transcript, or one row
  per intron with `by = "intron"`.

## Usage

``` r
extract_splice_sites(x, genome_fasta, overwrite = FALSE, chunk_bytes = 6.4e+07)

has_splice_sites(x)

splice_site_strength(
  self,
  gene_name = NULL,
  by = c("transcript", "intron"),
  ...
)
```

## Arguments

- x:

  An annotation directory path or an
  [IsoformicExperiment](https://luciorq.github.io/isoformic/reference/IsoformicExperiment.md).

- genome_fasta:

  Path to a genome FASTA file (`.fa`, `.fa.gz`).

- overwrite:

  Logical, replace an existing table (default `FALSE`).

- chunk_bytes:

  Number of bytes read per chunk while streaming.

- self:

  An annotation directory path or an
  [IsoformicExperiment](https://luciorq.github.io/isoformic/reference/IsoformicExperiment.md).
  An experiment reports its own transcripts.

- gene_name:

  Optional gene names or IDs to keep.

- by:

  `"transcript"` (default) or `"intron"`.

- ...:

  Not used; reserved for methods.

## Value

- `extract_splice_sites()`: the path of the written Parquet file,
  invisibly. The table has `seqid`, `strand`, `intron_start`,
  `intron_end` (1-based, inclusive), `donor_seq` and `acceptor_seq`. An
  existing table is kept, with a message, unless `overwrite = TRUE`.

- `has_splice_sites()`: `TRUE` when the table is present.

- `splice_site_strength()`: a tibble. By transcript: `transcript_id`,
  `gene_id`, `gene_name`, `n_introns`, `n_scored`, `weakest_donor`,
  `weakest_acceptor`, `mean_donor`, `mean_acceptor` and
  `n_noncanonical`; transcripts without introns have `NA` scores. By
  intron: one row per transcript and intron with the coordinates, the
  two window sequences, `donor_score`, `acceptor_score` and
  `canonical_dinucleotides`.

## Details

The donor window holds the last 3 exonic and the first 6 intronic bases,
the acceptor window the last 20 intronic and the first 3 exonic bases,
in the direction of transcription (the windows of MaxEntScan, Yeo and
Burge 2004). The matrices are trained on the unique introns of the
canonical transcripts (`is_canonical`, see
[`prepare_isoformic_annotation()`](https://luciorq.github.io/isoformic/reference/prepare_isoformic_annotation.md)),
or on all introns when the annotation marks fewer than 100 canonical
introns. A score is the sum over the window of
`log2(frequency / background)`, with the background taken from the base
composition of the training windows and half a pseudocount per base.
Positive means closer to the consensus.

The matrices are trained once per annotation directory and kept for the
session. With `gene_name`, only the introns of those genes are read and
scored, with the same matrices, so a transcript gets the same scores as
in the genome-wide table. On an
[IsoformicExperiment](https://luciorq.github.io/isoformic/reference/IsoformicExperiment.md),
the table by transcript is cached in the object, and later calls with
`gene_name` take their rows from it.

Limits. The matrices assume that positions are independent, so they are
weaker than MaxEntScan or a sequence model, and they describe the major
(GT-AG) spliceosome: U12 introns and other non-canonical sites get low
scores whatever their real efficiency (`n_noncanonical` counts the
introns that are not GT-AG). Canonical transcripts train the matrices,
so their own scores are slightly optimistic. Windows that run off the
sequence or hold a base other than A, C, G or T are not scored.

## References

Yeo, G., Burge, C.B. (2004). Maximum entropy modeling of short sequence
motifs with applications to RNA splicing signals. *Journal of
Computational Biology* 11(2-3), 377-394.

## See also

[`extract_transcript_sequences()`](https://luciorq.github.io/isoformic/reference/transcript_sequences.md),
[`transcript_productivity()`](https://luciorq.github.io/isoformic/reference/productivity.md)

## Examples

``` r
gff <- system.file("extdata", "example.gff3.gz", package = "isoformic")
annot <- prepare_isoformic_annotation(
  gff,
  output_path = tempdir(),
  annot_name = "splice_demo"
)
#> ℹ Parsing example.gff3.gz (GFF) into /tmp/RtmpPCTCTT/splice_demo.
has_splice_sites(annot)
#> [1] FALSE
if (FALSE) { # \dontrun{
genome <- download_reference(file_type = "genome_fasta", region = "primary")
extract_splice_sites(annot, genome_fasta = genome)
splice_site_strength(annot)
} # }
```
