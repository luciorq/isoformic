# Open Reading Frames of Transcripts from Their Sequences

Computes, for every assayed transcript (every transcript of the
annotation when `self` is a path), an open reading frame (ORF) and the
sequence-level features that decide whether it is productive. Requires
transcript sequences in the annotation directory (see
[`extract_transcript_sequences()`](https://luciorq.github.io/isoformic/reference/transcript_sequences.md)).
The ORF of each transcript comes from the first applicable source
(`orf_source`):

1.  `"annotation"`: a complete annotated CDS. Its stop codon is added
    when the annotation excludes it (as GENCODE does) and the ORF is
    re-read from the annotated start codon to the first in-frame stop,
    so a CDS that does not match the sequence is corrected
    (`annotation_cds_agrees = FALSE`). Transcripts tagged `seleno`,
    `ribosomal_slippage` or `non_ATG_start` keep the annotated CDS.

2.  `"projected"`: the start codon of the gene's canonical isoform lies
    in an exon of the transcript and reads `ATG` there; translation
    proceeds from that position to the first in-frame stop (the ORFanage
    approach). This is what rescues novel isoforms without CDS.

3.  `"annotation_partial"`: an annotated CDS flagged incomplete
    (`cds_start_NF` / `cds_end_NF`) whose start could not be projected;
    the annotated span is kept, and when its 3' end is known the NMD
    rule still applies (GENCODE annotates many NMD transcripts this
    way).

4.  `"longest"`: none of the above; the longest `ATG`-to-stop ORF of at
    least `min_orf_length` nt is used and `coding_probability` tells
    whether it is believable.

5.  `"none"`: no ORF of at least `min_orf_length` nt.

From the ORF the function derives the stop codon (`NA` when translation
runs into the transcript end, a non-stop decay substrate), the UTR
lengths, uORFs in the 5' UTR (`n_uorf`, complete upstream ORFs;
`n_overlapping_uorf`, out-of-frame upstream ORFs overlapping the main
start), the genomic coordinates of start and stop, whether start and
stop are shared with the canonical isoform and whether the reading frame
at the stop codon is in frame with the canonical ORF
(`frame_vs_canonical`), and the NMD class:

- `"ptc"`: the stop lies more than `nmd_rule` nt upstream of the last
  exon-exon junction (the 50-nt rule), NMD predicted.

- `"ptc_start_proximal"`: same, but the ORF is shorter than
  `start_proximal` nt; still called unproductive, with the note that
  such transcripts can escape NMD by translation re-initiation.

- `"non_stop"`: no in-frame stop codon.

- `"long_3utr"`: no premature stop, but a 3' UTR longer than `long_utr3`
  nt (default `Inf`, so the flag is off unless asked for).

- `"none"`: no NMD feature.

- `NA`: no stop codon in a transcript whose 3' end is annotated as
  incomplete (`cds_end_NF` / `mRNA_end_NF`, column `end_incomplete`), so
  nothing can be said.

`coding_potential()` scores every transcript's longest ORF with the
Fickett TESTCODE statistic and an in-frame hexamer usage score, combined
with the ORF length and the fraction of the transcript the ORF covers by
a logistic model trained on GENCODE 50 human coding versus lncRNA
transcripts (see `data-raw/build_coding_model.R` and
[`build_coding_model()`](https://luciorq.github.io/isoformic/reference/build_coding_model.md));
`coding_probability` above `coding_cutoff` (0.364, the CPAT human
cutoff) is called coding.

## Usage

``` r
transcript_orfs(
  self,
  gene_name = NULL,
  min_orf_length = 90L,
  nmd_rule = 50,
  start_proximal = 150,
  long_utr3 = Inf,
  ...
)

coding_potential(
  self,
  gene_name = NULL,
  min_orf_length = 90L,
  coding_cutoff = 0.364,
  ...
)
```

## Arguments

- self:

  An
  [IsoformicExperiment](https://luciorq.github.io/isoformic/reference/IsoformicExperiment.md)
  or an annotation directory path.

- gene_name:

  Optional genes (names or IDs) to restrict the result to.

- min_orf_length:

  Minimum ORF length in nt, stop codon included (default 90, i.e. 29
  amino acids).

- nmd_rule:

  Distance (nt) from the stop codon to the last junction above which NMD
  is predicted; default 50.

- start_proximal:

  ORF length (nt) below which a premature stop is classed
  `ptc_start_proximal` instead of `ptc`; default 150.

- long_utr3:

  3' UTR length (nt) above which `long_3utr` is flagged; default `Inf`.

- ...:

  Not used; reserved for methods.

- coding_cutoff:

  Coding probability above which a transcript is called coding (default
  0.364).

## Value

`transcript_orfs()` returns a tibble with one row per transcript:
`transcript_id`, `transcript_name`, `gene_id`, `gene_name`,
`has_sequence`, `orf_source`, `annotation_cds_agrees`, `orf_tx_start`
and `orf_tx_end` (ORF span in transcript coordinates, stop codon
included), `orf_length`, `protein_length`, `start_codon`, `stop_codon`,
`orf_genomic_start`, `orf_genomic_end`, `utr5_length`, `utr3_length`,
`stop_to_last_junction`, `nmd_class`, `end_incomplete`, `n_uorf`,
`n_overlapping_uorf`, `canonical_transcript_id`,
`start_same_as_canonical`, `stop_same_as_canonical`,
`frame_vs_canonical` (`"in_frame"`, `"frameshift"` or `"novel_sequence"`
when the stop lies outside the canonical isoform), and the coding
potential columns `longest_orf_length`, `orf_coverage`, `fickett_score`,
`hexamer_score` and `coding_probability`. `coding_potential()` returns
`transcript_id`, `transcript_name`, `gene_id`, `gene_name`, the five
coding potential columns and `is_coding`.

## Stored ORF table

The analysis takes minutes for a whole GENCODE release, so the rows
computed are written to `annot_data_orfs.parquet` in the annotation
directory and read back by later calls, sessions and objects that use
the same annotation. A stored row is used when the canonical isoform of
the gene, among the transcripts of the call, is the one it was computed
with, so the result is the same with and without the table. A gene with
a transcript that has no such row is computed whole, and its rows are
stored when its canonical isoform is the one of the whole annotation.
The table is computed again when the transcripts, exons or sequences of
the annotation, `min_orf_length`, the coding potential model
([`build_coding_model()`](https://luciorq.github.io/isoformic/reference/build_coding_model.md))
or the package version change. It is not written into a read-only
directory or the example annotation of the installed package, and
`options(isoformic.orf_cache = FALSE)` turns it off. The warning about
sequences that differ in length from their exons is given when the rows
are computed, not when they are read back.

## See also

[`transcript_productivity()`](https://luciorq.github.io/isoformic/reference/productivity.md),
[`transcript_features()`](https://luciorq.github.io/isoformic/reference/transcript_features.md)

## Examples

``` r
iso <- example_isoformic()
orfs <- transcript_orfs(iso, gene_name = "RBPJ")
orfs[, c("transcript_name", "orf_source", "orf_length", "nmd_class", "coding_probability")]
#> # A tibble: 65 × 5
#>    transcript_name orf_source         orf_length nmd_class coding_probability
#>    <chr>           <chr>                   <dbl> <chr>                  <dbl>
#>  1 RBPJ-222        annotation_partial        283 NA                    0.0809
#>  2 RBPJ-247        annotation                177 ptc                   0.855 
#>  3 RBPJ-243        longest                  1098 none                  0.885 
#>  4 RBPJ-201        annotation               1503 none                  0.984 
#>  5 RBPJ-213        annotation_partial        291 NA                    0.101 
#>  6 RBPJ-223        annotation               1503 none                  0.984 
#>  7 RBPJ-203        annotation               1398 none                  0.980 
#>  8 RBPJ-206        annotation               1503 none                  0.959 
#>  9 RBPJ-225        longest                  1140 none                  0.909 
#> 10 RBPJ-229        annotation_partial        145 NA                    0.0200
#> # ℹ 55 more rows
```
