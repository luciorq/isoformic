# Productivity of Transcripts and Genes

The biotype split used by the workflow (protein coding versus retained
intron, nonsense-mediated decay, ...) is a proxy. These functions make
productivity a computed, graded property of every transcript and
summarise it per gene:

- `transcript_productivity()` classifies each assayed transcript as
  `"productive"` (complete ORF, no predicted NMD, same ORF as the
  canonical isoform), `"altered_productive"` (complete ORF, no NMD, but
  a different ORF than the canonical isoform: alternative start,
  in-frame change, frameshift or stop in novel sequence),
  `"unproductive"` (no ORF, premature stop codon by the 50-nt rule, no
  stop codon, an unproductive biotype such as `nonsense_mediated_decay`,
  `retained_intron`, `protein_coding_CDS_not_defined`,
  `protein_coding_LoF`, or a longest ORF without coding potential),
  `"noncoding"` (transcripts of genes without any coding isoform) or
  `"uncertain"` (incomplete CDS annotation without a premature stop, a
  credible ORF unrelated to the canonical one, FANTOM
  `uncertain_coding`). The `reason` column explains each call and
  `source` tells whether it came from the sequence (`"sequence"`, when
  the annotation has transcript sequences, see
  [`extract_transcript_sequences()`](https://luciorq.github.io/isoformic/reference/transcript_sequences.md)
  and
  [`transcript_orfs()`](https://luciorq.github.io/isoformic/reference/transcript_orfs.md)),
  the annotation (`"annotation"`) or external evidence (`"evidence"`).
  Evidence columns take precedence, in this order: a
  `productivity_evidence` column added with
  [`add_annotation_columns()`](https://luciorq.github.io/isoformic/reference/add_annotation_columns.md),
  a validated NMD exon
  ([`add_nmd_exons()`](https://luciorq.github.io/isoformic/reference/productivity_evidence.md)),
  a Ribo-seq ORF matching the transcript's ORF
  ([`add_riboseq_orfs()`](https://luciorq.github.io/isoformic/reference/productivity_evidence.md)),
  an APPRIS principal isoform
  ([`add_appris()`](https://luciorq.github.io/isoformic/reference/productivity_evidence.md)).

- `productive_fraction()` is the share of a gene's expression carried by
  productive (and altered productive) isoforms, per sample or per group,
  with the difference between the last and the first group level. A gene
  whose expression is flat while its productive fraction drops is a
  post-transcriptional regulation call that gene-level differential
  expression misses. With inferential replicates (`inf_reps`, Salmon
  bootstraps or Gibbs samples) it reports the uncertainty of each
  fraction, and with a `reference`
  ([`gtex_reference()`](https://luciorq.github.io/isoformic/reference/gtex_reference.md))
  the difference to the healthy-tissue baseline.

- `productivity_fractions()` splits the gene's expression over all
  productivity classes and, for unproductive isoforms, over the
  mechanism (`nmd`, `retained_intron`, `no_orf`, `non_stop`, `other`).

- `test_productive_fraction()` tests the change in productive fraction
  between two groups on counts, with a quasi-binomial model per gene
  (pooled proportions per group, Pearson over-dispersion, t test on the
  log odds ratio, BH adjustment). A gene with fewer than 5 reads in its
  smaller class (productive or not), pooled over both groups, carries
  too little information on usage: it gets a `log_odds_ratio` of 0 and
  no statistic or p-value, and does not count in the BH adjustment. A
  gene without reads in one group gets `NA` there too.

- `run_productivity_enrichment()` ranks genes by the change in
  productive fraction between the last and the first group level
  (`rank_by = "delta"`), by that change over its standard deviation
  (`rank_by = "delta_z"`) or by the statistic of
  `test_productive_fraction()` (`rank_by = "statistic"`), and runs gene
  set enrichment (`fgsea`) on gene symbols, which captures regulated
  unproductive splicing (AS-NMD) directly. It takes an
  [IsoformicExperiment](https://luciorq.github.io/isoformic/reference/IsoformicExperiment.md),
  a `productive_fraction()` table or a `test_productive_fraction()`
  table, and one or several gene set collections.

- `poison_exons()` lists exons that occur only in unproductive isoforms
  of genes that also have productive isoforms and that do not overlap
  any exon of a productive isoform (cassette poison / NMD exons, as
  opposed to retained introns or alternative splice sites);
  [`plot_genomic_context()`](https://luciorq.github.io/isoformic/reference/plot_genomic_context.md)
  outlines them with `highlight_poison_exons`.

- `benchmark_productivity()` measures the rules against the annotation
  biotypes: recall of `nonsense_mediated_decay` and `retained_intron`
  transcripts by the sequence/structure rules alone, and the share of
  `protein_coding` and canonical transcripts called productive.

- `productivity_palette()` gives colours for the classes.

## Usage

``` r
productivity_palette()

transcript_productivity(
  self,
  gene_name = NULL,
  nmd_rule = 50,
  use_biotype = TRUE,
  use_sequence = TRUE,
  use_evidence = TRUE,
  coding_cutoff = 0.364,
  ...
)

productive_fraction(
  self,
  group = NULL,
  levels = NULL,
  assay = "abundance",
  include_altered = TRUE,
  as_matrix = FALSE,
  gene_name = NULL,
  min_expr = 0,
  min_gene_count = 0,
  inf_reps = NULL,
  reference = NULL,
  ...
)

productivity_fractions(
  self,
  group = NULL,
  levels = NULL,
  assay = "abundance",
  min_expr = 0,
  gene_name = NULL,
  ...
)

test_productive_fraction(
  self,
  group,
  levels = NULL,
  assay = "counts",
  include_altered = TRUE,
  min_expr = 0,
  min_gene_count = 0,
  gene_name = NULL,
  ...
)

run_productivity_enrichment(
  self,
  genesets_list,
  group = NULL,
  levels = NULL,
  assay = "abundance",
  min_size = 15,
  max_size = 500,
  rank_by = c("delta", "delta_z", "statistic"),
  include_altered = TRUE,
  min_expr = 0,
  min_gene_count = 0,
  inf_reps = NULL,
  workers = getOption("isoformic.workers", 1L),
  ...
)

poison_exons(self, gene_name = NULL, ...)

benchmark_productivity(self, nmd_rule = 50, use_sequence = TRUE, ...)
```

## Arguments

- self:

  An
  [IsoformicExperiment](https://luciorq.github.io/isoformic/reference/IsoformicExperiment.md).
  `run_productivity_enrichment()` also takes the table of
  `productive_fraction()` for two or more groups, or of
  `test_productive_fraction()` for `rank_by = "statistic"`, so that one
  table serves several gene set collections; the arguments that build
  the table (`group`, `levels`, `assay`, `include_altered`, `min_expr`,
  `min_gene_count`, `inf_reps`) are then ignored.

- gene_name:

  Optional genes (names or IDs) to restrict the result to.

- nmd_rule:

  Distance (nt) from the stop codon to the last junction above which NMD
  is predicted; default 50.

- use_biotype:

  Logical, also use unproductive biotypes and FANTOM coding status
  (default `TRUE`).

- use_sequence:

  Logical, use the sequence-derived ORFs of
  [`transcript_orfs()`](https://luciorq.github.io/isoformic/reference/transcript_orfs.md)
  when the annotation has sequences (default `TRUE`).

- use_evidence:

  Logical, apply the evidence columns (default `TRUE`).

- coding_cutoff:

  Coding probability below which a longest ORF (one that could not be
  projected from the canonical start) is not believed (default 0.364).

- ...:

  Used only by `run_productivity_enrichment()`, which passes them to
  [`fgsea::fgseaMultilevel()`](https://rdrr.io/pkg/fgsea/man/fgseaMultilevel.html);
  the other functions do not use them (reserved for methods). Each
  `fgsea` call runs with
  [`BiocParallel::SerialParam()`](https://rdrr.io/pkg/BiocParallel/man/SerialParam-class.html)
  unless `nproc` or `BPPARAM` is given; prefer `workers` (see *Workers*
  in
  [`run_enrichment()`](https://luciorq.github.io/isoformic/reference/enrichment.md)).

- group:

  Name of a `col_data` column defining groups; `NULL` for per-sample
  values. `run_productivity_enrichment()` requires it for an
  [IsoformicExperiment](https://luciorq.github.io/isoformic/reference/IsoformicExperiment.md)
  and ignores it for a `productive_fraction()` table.

- levels:

  Levels of `group` to use, in this order; the change is `last - first`.
  Samples of other groups are ignored, as if the object had been subset
  to these samples, so every returned value comes from the samples kept.
  `NULL` takes all the levels of `group`. The tests
  (`test_productive_fraction()`,
  [`test_apa()`](https://luciorq.github.io/isoformic/reference/polyadenylation.md))
  need exactly two levels.

- assay:

  Transcript assay to use (default `"abundance"`; `"counts"` for the
  test).

- include_altered:

  Logical, count `altered_productive` isoforms as productive in
  `productive_fraction()` (default `TRUE`).

- as_matrix:

  Logical, return the genes x samples matrix of productive fractions
  instead of the summary table.

- min_expr:

  Minimum assay value for an isoform to contribute (default 0): isoforms
  below it count as zero.

- min_gene_count:

  Minimum mean count of a gene, the sum of its isoforms in the
  `"counts"` assay, in every group (over all the samples kept without
  `group`); genes below it are dropped. Needs a `"counts"` assay
  whatever `assay` is. The default 0 keeps every gene. See Details.

- inf_reps:

  Optional list of transcript matrices (inferential replicates with the
  transcripts in rows and the samples of the assay in columns, as
  counts, as
  [`read_inferential_replicates()`](https://luciorq.github.io/isoformic/reference/inferential-io.md)
  returns them); `NULL` uses `x@metadata$inf_reps` when present and
  `FALSE` skips the replicates, which leaves out the standard deviations
  and `delta_z` and is much faster. For an `assay` other than `"counts"`
  the draws are put on its scale with the `"counts"` assay (see
  Details); an object without a `"counts"` assay takes the draws to be
  on the scale of `assay`. `run_productivity_enrichment()` reads the
  replicates only for `rank_by = "delta_z"`.

- reference:

  Optional tibble from
  [`gtex_reference()`](https://luciorq.github.io/isoformic/reference/gtex_reference.md),
  or a GTEx tissue name, giving the healthy-tissue productive fraction
  per gene.

- genesets_list:

  Named list of gene symbol vectors (see
  [`read_gmt()`](https://luciorq.github.io/isoformic/reference/enrichment.md)),
  or a named list of such lists: each collection is then tested on its
  own and the result has a `collection` column (see
  [`run_enrichment()`](https://luciorq.github.io/isoformic/reference/enrichment.md)).

- min_size, max_size:

  Gene set size limits passed to `fgsea`.

- rank_by:

  Statistic that ranks the genes in `run_productivity_enrichment()`:
  `"delta"` (default), the change in productive fraction, `"delta_z"`,
  that change divided by its standard deviation over the inferential
  replicates, which moves the genes with an uncertain change towards the
  middle of the ranking, or `"statistic"`, the test of
  `test_productive_fraction()` on the `"counts"` assay between the first
  and the last level of `levels`, as a signed normal quantile of its
  p-value (see Details); `assay` and `inf_reps` are then not used.
  `"delta_z"` needs inferential replicates; a change that is the same in
  every replicate (`delta_z` is `NA` in the `productive_fraction()`
  table) is ranked with the largest absolute `delta_z` of the table (at
  least 1), with the sign of the change. It reflects the quantification
  uncertainty only, not the variation between samples or the number of
  reads; `"statistic"` reflects both (see Details).

- workers:

  Number of worker processes that run the `fgsea` calls of
  `run_productivity_enrichment()` (one per gene set collection) in
  parallel: a whole number, or `"auto"` for the CPUs available to the
  session. Defaults to the `isoformic.workers` option, or 1 (serial).
  Needs the `mirai` package above 1; see *Workers* in
  [`run_enrichment()`](https://luciorq.github.io/isoformic/reference/enrichment.md),
  also for the random seeds.

## Value

- `transcript_productivity()`: one row per transcript with
  `transcript_id`, `transcript_name`, `transcript_type`, `gene_id`,
  `gene_name`, `productivity` (a factor with the levels of
  `productivity_palette()`), `reason`, `source`, `has_cds`,
  `cds_complete`, `orf_length`, `nmd_predicted`,
  `stop_to_last_junction`, `orf_fraction_of_canonical`, `is_canonical`
  and, with sequences, the ORF columns behind the call (`orf_source`,
  `nmd_class`, `frame_vs_canonical`, `coding_probability`, `n_uorf`).

- `productive_fraction()`: one row per coding gene with `gene_id`,
  `gene_name`, `mean_gene_expr`, `productive_fraction` (mean over the
  samples kept), with `group` one `productive_fraction_<level>` column
  per level and `delta_productive_fraction`, with inferential replicates
  `productive_fraction_sd` (the standard deviation over the draws, see
  Details) and, with `group`, `delta_productive_fraction_sd` and
  `delta_z`, and with `reference` the `reference_*` columns and
  `delta_vs_reference`. A genes x samples matrix with
  `as_matrix = TRUE`.

- `productivity_fractions()`: one row per gene and class present
  (`productivity`, a factor), with `fraction` and, with `group`, the
  `fraction_<level>` and `delta_fraction` columns.

- `test_productive_fraction()`: one row per coding gene with
  `mean_gene_expr`, the two `productive_fraction_<level>` columns,
  `delta_productive_fraction`, `log_odds_ratio`, `se`, `statistic`,
  `dispersion`, `df`, `pvalue` and `padj`, ordered by `pvalue`.

- `run_productivity_enrichment()`: the `fgsea` result with the columns
  `experiment` (`"productive_fraction"`), `rank_by` and, for a list of
  collections, `collection`.

- `poison_exons()`: one row per exon with `gene_id`, `gene_name`,
  `seqid`, `start_pos`, `end_pos`, `strand`, `n_unproductive`,
  `transcript_ids` and `exon_length`.

- `benchmark_productivity()`: a list with the `summary` (`metric`, `n`,
  `value`, `source`) and `confusion` (`transcript_type`, `productivity`,
  `n`) tibbles.

- `productivity_palette()`: a named character vector of colours.

## Details

**Low-count genes.** A productive fraction computed from a few reads is
mostly sampling noise: a gene with 4 reads moves by 0.25 when one read
changes isoform. `delta_z` does not catch this, because the inferential
replicates measure how the reads are split between isoforms, not how
many reads were sampled. `min_gene_count` removes genes whose mean count
(the sum of the counts of their isoforms in the `"counts"` assay) is
below it in any group, before the statistics and the BH adjustment; with
`rank_by = "statistic"`, before the ranking. The filter is off by
default (0); 10 reads is a reasonable start. The rule asks for enough
reads in every group, because a fraction is compared between groups: it
is stricter than the transcript filter of
[`test_usage()`](https://luciorq.github.io/isoformic/reference/test_usage.md)
(`min_count` reads in as many samples as the smallest group), and
`min_expr` sets single isoform values to zero rather than removing
genes. `rank_by = "statistic"` accounts for the reads of every gene and
for the variation between its samples. Its dispersion is estimated for
each gene alone, on `n - 2` degrees of freedom and not moderated, so
with three samples per group the ranking is still noisy. The t statistic
is put on the normal scale through its p-value
(`sign(statistic) * qnorm(pvalue / 2, lower.tail = FALSE)`), because the
degrees of freedom differ between genes when samples without reads drop
out. Its sign is that of the pooled change (`delta_productive_fraction`
of `test_productive_fraction()`), which can differ from the mean of the
per-sample fractions in `productive_fraction()` when the samples differ
much in depth. Genes with fewer than 5 reads in one class get no
statistic and are left out of the ranking, so that a change in the
expression of the whole gene does not rank it.

**Uncertainty of the fractions.** With inferential replicates,
`productive_fraction()` computes the productive fraction of every gene
and sample in each draw. The draws are counts; for another `assay` each
draw is first put on the scale of that assay, multiplied by
`assay / counts` of the same transcript and sample (for TPM, the inverse
of the effective length times the scaling of the sample), so that the
fractions of the draws and of the assay, and `min_expr`, are on one
scale. Where the point count or value is zero (Gibbs draws often give
reads to such cells) the factor is estimated: from the `length` assay
when the assay follows it (TPM of
[`read_salmon()`](https://luciorq.github.io/isoformic/reference/read_salmon.md)
and
[`read_oarfish()`](https://luciorq.github.io/isoformic/reference/read_oarfish.md)),
otherwise from a transcript effect plus a sample effect on the log scale
fitted on the other cells, and for a transcript without any count from
its annotated length. `productive_fraction_sd` is the sample standard
deviation of a gene's fraction over the draws, around the mean of the
draws, averaged over the samples kept; `delta_productive_fraction_sd` is
the standard deviation over the draws of the change between the group
means, and `delta_z` is `delta_productive_fraction` divided by it. In
every draw, only the samples in which the assay gives the gene a
fraction are used. Both measure the spread of the draws, not their
distance to the point estimate: the mean of Gibbs draws can differ from
the EM estimate of the assay.

## Examples

``` r
iso <- example_isoformic()
prod <- transcript_productivity(iso, gene_name = "RBPJ")
prod[, c("transcript_name", "transcript_type", "productivity", "reason", "source")]
#> # A tibble: 65 × 5
#>    transcript_name transcript_type         productivity       reason      source
#>    <chr>           <chr>                   <fct>              <chr>       <chr> 
#>  1 RBPJ-222        protein_coding          uncertain          incomplete… seque…
#>  2 RBPJ-247        nonsense_mediated_decay unproductive       unproducti… annot…
#>  3 RBPJ-243        retained_intron         unproductive       unproducti… annot…
#>  4 RBPJ-201        protein_coding          altered_productive alternativ… seque…
#>  5 RBPJ-213        protein_coding          uncertain          incomplete… seque…
#>  6 RBPJ-223        protein_coding          altered_productive alternativ… seque…
#>  7 RBPJ-203        protein_coding          altered_productive alternativ… seque…
#>  8 RBPJ-206        protein_coding          altered_productive alternativ… seque…
#>  9 RBPJ-225        retained_intron         unproductive       unproducti… annot…
#> 10 RBPJ-229        protein_coding          uncertain          incomplete… seque…
#> # ℹ 55 more rows
head(productive_fraction(iso, group = "condition"))
#> # A tibble: 6 × 7
#>   gene_id    gene_name mean_gene_expr productive_fraction productive_fraction_…¹
#>   <chr>      <chr>              <dbl>               <dbl>                  <dbl>
#> 1 ENSG00000… SNX6               47.8                0.337                  0.154
#> 2 ENSG00000… ERO1A              17.9                0.830                  0.655
#> 3 ENSG00000… DNAJB6             52.5                0.659                  0.496
#> 4 ENSG00000… GPT2                8.34               0.511                  0.377
#> 5 ENSG00000… PIP5K1A           108.                 0.408                  0.541
#> 6 ENSG00000… MYO7B              10.2                0.291                  0.414
#> # ℹ abbreviated name: ¹​productive_fraction_preeclampsia
#> # ℹ 2 more variables: productive_fraction_control <dbl>,
#> #   delta_productive_fraction <dbl>
head(productivity_fractions(iso, group = "condition"))
#> # A tibble: 6 × 7
#>   gene_id gene_name productivity fraction fraction_preeclampsia fraction_control
#>   <chr>   <chr>     <fct>           <dbl>                 <dbl>            <dbl>
#> 1 ENSG00… FSTL3     productive    0.889               0.926               0.871 
#> 2 ENSG00… FSTL3     altered_pro…  0.0452              0.0404              0.0474
#> 3 ENSG00… FSTL3     nmd           0.0210              0.0294              0.0171
#> 4 ENSG00… FSTL3     no_orf        0.00700             0.0000164           0.0103
#> 5 ENSG00… FSTL3     uncertain     0.0382              0.00450             0.0541
#> 6 ENSG00… CDC42     productive    0.0212              0.0120              0.0255
#> # ℹ 1 more variable: delta_fraction <dbl>
head(test_productive_fraction(iso, group = "condition"))
#> # A tibble: 6 × 13
#>   gene_id gene_name mean_gene_expr productive_fraction_…¹ productive_fraction_…²
#>   <chr>   <chr>              <dbl>                  <dbl>                  <dbl>
#> 1 ENSG00… MYO7B               593.                  0.647                  0.354
#> 2 ENSG00… CCSAP              5483.                  0.961                  0.977
#> 3 ENSG00… TPT1              24353.                  0.923                  0.813
#> 4 ENSG00… DNAJB6             1408.                  0.599                  0.711
#> 5 ENSG00… PIP5K1A            3240.                  0.815                  0.682
#> 6 ENSG00… PSMA3               889.                  0.618                  0.504
#> # ℹ abbreviated names: ¹​productive_fraction_preeclampsia,
#> #   ²​productive_fraction_control
#> # ℹ 8 more variables: delta_productive_fraction <dbl>, log_odds_ratio <dbl>,
#> #   se <dbl>, statistic <dbl>, dispersion <dbl>, df <dbl>, pvalue <dbl>,
#> #   padj <dbl>
poison_exons(iso, gene_name = "RBPJ")
#> # A tibble: 11 × 9
#>    gene_id            gene_name seqid start_pos  end_pos strand n_unproductive
#>    <chr>              <chr>     <chr>     <int>    <int> <chr>           <int>
#>  1 ENSG00000168214.24 RBPJ      chr4   26163458 26163614 +                   1
#>  2 ENSG00000168214.24 RBPJ      chr4   26163530 26163614 +                   1
#>  3 ENSG00000168214.24 RBPJ      chr4   26263617 26263714 +                   1
#>  4 ENSG00000168214.24 RBPJ      chr4   26321228 26321295 +                   1
#>  5 ENSG00000168214.24 RBPJ      chr4   26367874 26368145 +                   1
#>  6 ENSG00000168214.24 RBPJ      chr4   26385917 26386072 +                   1
#>  7 ENSG00000168214.24 RBPJ      chr4   26408060 26408174 +                   3
#>  8 ENSG00000168214.24 RBPJ      chr4   26410042 26410292 +                   2
#>  9 ENSG00000168214.24 RBPJ      chr4   26412930 26412987 +                   4
#> 10 ENSG00000168214.24 RBPJ      chr4   26414049 26414155 +                   3
#> 11 ENSG00000168214.24 RBPJ      chr4   26422132 26422156 +                   2
#> # ℹ 2 more variables: transcript_ids <chr>, exon_length <int>
benchmark_productivity(iso)$summary
#> # A tibble: 5 × 4
#>   metric                                  n value source  
#>   <chr>                               <int> <dbl> <chr>   
#> 1 nmd_biotype_called_unproductive       353 1     sequence
#> 2 retained_intron_called_unproductive   125 0.68  sequence
#> 3 protein_coding_called_productive     1226 0.900 sequence
#> 4 canonical_called_productive            57 0.982 sequence
#> 5 lncrna_called_noncoding               107 1     sequence
```
