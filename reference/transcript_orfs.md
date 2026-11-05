# Open Reading Frames of Transcripts from Their Sequences

Computes, for every assayed transcript, an open reading frame (ORF) and
the sequence-level features that decide whether it is productive.
Requires transcript sequences in the annotation directory (see
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
by a logistic model trained on GENCODE coding versus lncRNA transcripts
(see `data-raw/build_coding_model.R`); `coding_probability` above
`coding_cutoff` (0.364, the CPAT human cutoff) is called coding.

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

  Passed on to methods.

- coding_cutoff:

  Coding probability above which a transcript is called coding (default
  0.364).

## Value

A tibble with one row per transcript.

## See also

[`transcript_productivity()`](https://luciorq.github.io/isoformic/reference/productivity.md),
[`transcript_features()`](https://luciorq.github.io/isoformic/reference/transcript_features.md)

## Examples

``` r
iso <- example_isoformic()
orfs <- transcript_orfs(iso, gene_name = "RBPJ")
#> Warning: ! 1 transcript has a sequence length differing from the exon length; its ORF is
#>   skipped.
orfs[, c("transcript_name", "orf_source", "orf_length", "nmd_class", "coding_probability")]
#> # A tibble: 10 × 5
#>    transcript_name orf_source         orf_length nmd_class coding_probability
#>    <chr>           <chr>                   <dbl> <chr>                  <dbl>
#>  1 RBPJ-222        annotation_partial        283 NA                    0.0723
#>  2 RBPJ-201        annotation               1503 none                  0.975 
#>  3 RBPJ-206        annotation               1503 none                  0.984 
#>  4 RBPJ-225        projected                1140 none                  0.950 
#>  5 RBPJ-204        annotation               1464 none                  0.980 
#>  6 RBPJ-205        annotation               1461 none                  0.968 
#>  7 RBPJ-214        annotation               1398 none                  0.976 
#>  8 RBPJ-208        annotation               1149 none                  0.931 
#>  9 RBPJ-202        annotation               1461 none                  0.968 
#> 10 RBPJ-210        longest                   462 none                  0.548 
```
