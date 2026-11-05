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
  `retained_intron`, `protein_coding_CDS_not_defined`, or a longest ORF
  without coding potential), `"noncoding"` (transcripts of genes without
  any coding isoform) or `"uncertain"` (incomplete CDS annotation
  without a premature stop, a credible ORF unrelated to the canonical
  one, FANTOM `uncertain_coding`). The `reason` column explains each
  call and `source` tells whether it came from the sequence
  (`"sequence"`, when the annotation has transcript sequences, see
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
  log odds ratio, BH adjustment).

- `run_productivity_enrichment()` ranks genes by the change in
  productive fraction between two groups and runs gene set enrichment
  (`fgsea`) on gene symbols, which captures regulated unproductive
  splicing (AS-NMD) directly.

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
  gene_name = NULL,
  ...
)

run_productivity_enrichment(
  self,
  genesets_list,
  group,
  levels = NULL,
  assay = "abundance",
  min_size = 15,
  max_size = 500,
  ...
)

poison_exons(self, gene_name = NULL, ...)

benchmark_productivity(self, nmd_rule = 50, use_sequence = TRUE, ...)
```

## Arguments

- self:

  An
  [IsoformicExperiment](https://luciorq.github.io/isoformic/reference/IsoformicExperiment.md).

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

  Passed on to methods or to
  [`fgsea::fgseaMultilevel()`](https://rdrr.io/pkg/fgsea/man/fgseaMultilevel.html).

- group:

  Name of a `col_data` column defining groups; `NULL` for per-sample
  values.

- levels:

  Order of the group levels; the change is `last - first`. For the tests
  (`test_productive_fraction()`,
  [`test_apa()`](https://luciorq.github.io/isoformic/reference/polyadenylation.md)),
  the two groups to compare; samples of other groups are ignored.

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

- inf_reps:

  Optional list of transcript matrices (inferential replicates with the
  same dimnames as the assay); `NULL` uses `metadata(x)$inf_reps` when
  present.

- reference:

  Optional tibble from
  [`gtex_reference()`](https://luciorq.github.io/isoformic/reference/gtex_reference.md),
  or a GTEx tissue name, giving the healthy-tissue productive fraction
  per gene.

- genesets_list:

  Named list of gene symbol vectors (see
  [`read_gmt()`](https://luciorq.github.io/isoformic/reference/enrichment.md)).

- min_size, max_size:

  Gene set size limits passed to `fgsea`.

## Value

Tibbles, a matrix for `productive_fraction(as_matrix = TRUE)`, a named
character vector for `productivity_palette()`, a list with `summary` and
`confusion` tibbles for `benchmark_productivity()`.

## Examples

``` r
iso <- example_isoformic()
prod <- transcript_productivity(iso, gene_name = "RBPJ")
#> Warning: ! 1 transcript has a sequence length differing from the exon length; its ORF is
#>   skipped.
prod[, c("transcript_name", "transcript_type", "productivity", "reason", "source")]
#> # A tibble: 10 × 5
#>    transcript_name transcript_type productivity       reason              source
#>    <chr>           <chr>           <fct>              <chr>               <chr> 
#>  1 RBPJ-222        protein_coding  uncertain          incomplete 3' end,… seque…
#>  2 RBPJ-201        protein_coding  productive         complete ORF ident… seque…
#>  3 RBPJ-206        protein_coding  productive         complete ORF ident… seque…
#>  4 RBPJ-225        retained_intron unproductive       unproductive bioty… annot…
#>  5 RBPJ-204        protein_coding  altered_productive alternative start … seque…
#>  6 RBPJ-205        protein_coding  altered_productive alternative start … seque…
#>  7 RBPJ-214        protein_coding  altered_productive alternative start … seque…
#>  8 RBPJ-208        protein_coding  altered_productive alternative start … seque…
#>  9 RBPJ-202        protein_coding  altered_productive alternative start … seque…
#> 10 RBPJ-210        retained_intron unproductive       unproductive bioty… annot…
head(productive_fraction(iso, group = "condition"))
#> # A tibble: 6 × 7
#>   gene_id    gene_name mean_gene_expr productive_fraction productive_fraction_…¹
#>   <chr>      <chr>              <dbl>               <dbl>                  <dbl>
#> 1 ENSG00000… PNCK                50.6               0.251                  0.450
#> 2 ENSG00000… REV3L              831.                0.265                  0.458
#> 3 ENSG00000… PSMA3              860.                0.527                  0.614
#> 4 ENSG00000… ANO8               208.                0.668                  0.747
#> 5 ENSG00000… GDPD3              152.                0.792                  0.862
#> 6 ENSG00000… TSG101            1002.                0.847                  0.912
#> # ℹ abbreviated name: ¹​productive_fraction_preeclampsia
#> # ℹ 2 more variables: productive_fraction_control <dbl>,
#> #   delta_productive_fraction <dbl>
head(productivity_fractions(iso, group = "condition"))
#> # A tibble: 6 × 7
#>   gene_id gene_name productivity fraction fraction_preeclampsia fraction_control
#>   <chr>   <chr>     <fct>           <dbl>                 <dbl>            <dbl>
#> 1 ENSG00… UPF1      productive     0.132                 0.141           0.128  
#> 2 ENSG00… UPF1      altered_pro…   0.829                 0.823           0.832  
#> 3 ENSG00… UPF1      retained_in…   0.0386                0.0364          0.0397 
#> 4 ENSG00… SPPL2B    productive     0.399                 0.350           0.422  
#> 5 ENSG00… SPPL2B    altered_pro…   0.307                 0.377           0.274  
#> 6 ENSG00… SPPL2B    nmd            0.0111                0.0138          0.00984
#> # ℹ 1 more variable: delta_fraction <dbl>
head(test_productive_fraction(iso, group = "condition"))
#> # A tibble: 6 × 13
#>   gene_id gene_name mean_gene_expr productive_fraction_…¹ productive_fraction_…²
#>   <chr>   <chr>              <dbl>                  <dbl>                  <dbl>
#> 1 ENSG00… TSG101            20031.                  0.910                  0.814
#> 2 ENSG00… FAM120A           73031.                  0.933                  0.970
#> 3 ENSG00… POLR2F             7000.                  0.729                  0.634
#> 4 ENSG00… REV3L             16616.                  0.465                  0.183
#> 5 ENSG00… MED24             23507.                  0.822                  0.869
#> 6 ENSG00… U2AF2             40683.                  0.874                  0.840
#> # ℹ abbreviated names: ¹​productive_fraction_preeclampsia,
#> #   ²​productive_fraction_control
#> # ℹ 8 more variables: delta_productive_fraction <dbl>, log_odds_ratio <dbl>,
#> #   se <dbl>, statistic <dbl>, dispersion <dbl>, df <dbl>, pvalue <dbl>,
#> #   padj <dbl>
poison_exons(iso, gene_name = "RBPJ")
#> # A tibble: 0 × 9
#> # ℹ 9 variables: gene_id <chr>, gene_name <chr>, seqid <chr>, start_pos <int>,
#> #   end_pos <int>, strand <chr>, n_unproductive <int>, transcript_ids <chr>,
#> #   exon_length <int>
benchmark_productivity(iso)$summary
#> # A tibble: 5 × 4
#>   metric                                  n value source  
#>   <chr>                               <int> <dbl> <chr>   
#> 1 nmd_biotype_called_unproductive        38 0.974 sequence
#> 2 retained_intron_called_unproductive   145 0.607 sequence
#> 3 protein_coding_called_productive      252 0.782 sequence
#> 4 canonical_called_productive            67 1     sequence
#> 5 lncrna_called_noncoding                19 1     sequence
```
