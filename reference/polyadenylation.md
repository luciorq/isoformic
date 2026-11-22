# Poly(A) Sites and Alternative Polyadenylation of Transcripts

`polya_sites()` describes the 3' end of every assayed transcript:

- `tes_pos`, `tes_group` (alternative 3' ends of the gene, from
  [`transcript_features()`](https://luciorq.github.io/isoformic/reference/transcript_features.md)),
  `tes_rank` (1 = most proximal 3' end along the direction of
  transcription), `n_tes_groups`, `is_proximal`, `is_distal`.

- `apa_type` relative to the canonical isoform: `"same_end"`,
  `"tandem_utr"` (different 3' end, same stop codon: 3' UTR shortening
  or lengthening), `"alternative_last_exon"` (different last exon),
  `"intronic"` (3' end inside an intron of the canonical isoform,
  typically a composite terminal exon) or `"other"`; `NA` for a gene
  without a canonical isoform.

- `utr3_length` (from
  [`transcript_orfs()`](https://luciorq.github.io/isoformic/reference/transcript_orfs.md)
  when sequences are available, else from the annotated CDS).

- `polya_signal` and `polya_signal_distance`: the poly(A) signal hexamer
  (`AATAAA`, `ATTAAA` or one of the weaker variants) found in the last
  `signal_window` nt of the transcript sequence, and its distance to the
  3' end (canonical signals sit 10-30 nt upstream); `NA` without
  sequences, and when no hexamer lies in that window (common, because
  annotated 3' ends are often not the exact cleavage site).

- With `polya_atlas`, the nearest atlas site within `window` bp of the
  3' end: `pas_site_id`, `pas_class`, `pas_tpm`, `pas_distance`,
  `pas_supported`. The atlas is the frozen PolyASite 2.0 table
  (`polya_atlas = "polyasite"`, see
  [frozen_data](https://luciorq.github.io/isoformic/reference/frozen_data.md)),
  a data frame with `seqid`, `cleavage_pos`, `strand` (and optionally
  `site_id`, `site_class`, `tpm`), or a GENCODE `polyAs.gff3.gz` file
  (`download_reference(file_type = "polya")`), which `read_polya_gff()`
  reads into such a table.

`apa_usage()` quantifies per gene, for genes with at least two 3' end
groups: the share of expression ending at the distal-most 3' end
(`distal_fraction`) and the abundance-weighted 3' UTR length
(`weighted_utr3_length`), per sample or per group with the change
between the last and the first group level. `test_apa()` tests the
change in distal usage between two groups with the quasi-binomial test
used by
[`test_productive_fraction()`](https://luciorq.github.io/isoformic/reference/productivity.md),
on counts. A gene with fewer than 5 reads at its distal or at its other
ends, pooled over both groups, is not tested (no statistic or p-value).

`run_apa_enrichment()` runs a gene set enrichment analysis
([`fgsea::fgseaMultilevel()`](https://rdrr.io/pkg/fgsea/man/fgseaMultilevel.html))
on the genes ranked by their change in distal usage
(`rank_by = "delta"`, the `delta_distal_fraction` of `apa_usage()`) or
in weighted 3' UTR length (`rank_by = "utr3_length"`,
`delta_weighted_utr3_length`): a positive score is a gene set whose
genes move to distal 3' ends, or to longer 3' UTRs, in the last group
level. It also takes the table of `apa_usage()` for two or more groups,
so that one table serves several calls. As in
[`run_enrichment()`](https://luciorq.github.io/isoformic/reference/enrichment.md),
`genesets_list` can hold several collections, each tested on its own. A
gene symbol that occurs more than once keeps its largest absolute
change, and all the gene sets are returned (there is no p-value cutoff).

## Usage

``` r
polya_sites(
  self,
  gene_name = NULL,
  polya_atlas = NULL,
  window = 50L,
  signal_window = 50L,
  ...
)

read_polya_gff(path)

apa_usage(
  self,
  group = NULL,
  levels = NULL,
  assay = "abundance",
  min_expr = 0,
  as_matrix = FALSE,
  gene_name = NULL,
  ...
)

test_apa(
  self,
  group,
  levels = NULL,
  assay = "counts",
  min_expr = 0,
  gene_name = NULL,
  ...
)

run_apa_enrichment(
  self,
  genesets_list,
  group = NULL,
  levels = NULL,
  assay = "abundance",
  min_size = 15,
  max_size = 500,
  rank_by = c("delta", "utr3_length"),
  min_expr = 0,
  workers = getOption("isoformic.workers", 1L),
  ...
)
```

## Arguments

- self:

  An
  [IsoformicExperiment](https://luciorq.github.io/isoformic/reference/IsoformicExperiment.md);
  for `run_apa_enrichment()` also the table of `apa_usage()` for two or
  more groups.

- gene_name:

  Optional genes (names or IDs) to restrict the result to.

- polya_atlas:

  `NULL` (no atlas), `"polyasite"`, a data frame or a file path (Parquet
  table or GENCODE poly(A) GFF3).

- window:

  Distance in bp within which an atlas site supports a 3' end (default
  50).

- signal_window:

  Length of the 3' terminal sequence searched for a poly(A) signal
  (default 50 nt).

- ...:

  Not used by `polya_sites()`, `apa_usage()` and `test_apa()`; reserved
  for methods. `run_apa_enrichment()` passes them to
  [`fgsea::fgseaMultilevel()`](https://rdrr.io/pkg/fgsea/man/fgseaMultilevel.html),
  which runs with
  [`BiocParallel::SerialParam()`](https://rdrr.io/pkg/BiocParallel/man/SerialParam-class.html)
  unless `nproc` or `BPPARAM` is given; prefer `workers` (see *Workers*
  in
  [`run_enrichment()`](https://luciorq.github.io/isoformic/reference/enrichment.md)).

- path:

  Path to a GENCODE `polyAs.gff3.gz` file.

- group:

  Name of a `col_data` column defining groups; `NULL` for per-sample
  values. `test_apa()` and `run_apa_enrichment()` on an object need it.
  With an `apa_usage()` table, `run_apa_enrichment()` ignores `group`,
  `levels`, `assay` and `min_expr`.

- levels:

  Levels of `group` to use, in this order; the change is `last - first`.
  Samples of other groups are ignored, as if the object had been subset
  to these samples, so every returned value comes from the samples kept.
  `NULL` takes all the levels of `group`. `test_apa()` needs exactly two
  levels.

- assay:

  Transcript assay to use (default `"abundance"`; `"counts"` for
  `test_apa()`).

- min_expr:

  Minimum assay value for an isoform to count (default 0): values below
  it count as zero.

- as_matrix:

  Logical, return the genes x samples matrix of distal fractions instead
  of the summary table.

- genesets_list:

  Named list of character vectors (gene symbols), as returned by
  [`read_gmt()`](https://luciorq.github.io/isoformic/reference/enrichment.md),
  or a named list of such lists: each collection is then tested on its
  own and the result has a `collection` column (see
  [`run_enrichment()`](https://luciorq.github.io/isoformic/reference/enrichment.md)).

- min_size, max_size:

  Minimum and maximum gene set size passed to `fgsea`.

- rank_by:

  Statistic that ranks the genes in `run_apa_enrichment()`: `"delta"`
  (default), the change in distal fraction, or `"utr3_length"`, the
  change in weighted 3' UTR length.

- workers:

  Number of worker processes that run the `fgsea` calls of
  `run_apa_enrichment()` (one per gene set collection) in parallel: a
  whole number, or `"auto"` for the CPUs available to the session.
  Defaults to the `isoformic.workers` option, or 1 (serial). Needs the
  `mirai` package above 1; see *Workers* in
  [`run_enrichment()`](https://luciorq.github.io/isoformic/reference/enrichment.md),
  also for the random seeds.

## Value

Tibbles with one row per transcript (`polya_sites()`) or per gene
(`apa_usage()`, `test_apa()`); a matrix for `as_matrix = TRUE`.
`run_apa_enrichment()` returns the `fgsea` table with the columns
`experiment` (`"distal_fraction"` or `"weighted_utr3_length"`),
`rank_by` and, for several collections, `collection`. `read_polya_gff()`
returns a tibble with one row per poly(A) feature (`site_id`, `seqid`,
`start_pos`, `end_pos`, `cleavage_pos`, `strand`, `tpm`, `site_class`,
`signals`).

## Examples

``` r
iso <- example_isoformic()
pa <- polya_sites(iso, gene_name = "RBPJ")
pa[, c("transcript_name", "tes_rank", "apa_type", "utr3_length", "polya_signal")]
#> # A tibble: 65 × 5
#>    transcript_name tes_rank apa_type              utr3_length polya_signal
#>    <chr>              <int> <chr>                       <dbl> <chr>       
#>  1 RBPJ-222               3 alternative_last_exon          NA NA          
#>  2 RBPJ-247              12 alternative_last_exon         936 NA          
#>  3 RBPJ-243              12 alternative_last_exon          NA NA          
#>  4 RBPJ-201              19 same_end                     4124 NA          
#>  5 RBPJ-213               3 alternative_last_exon          NA NA          
#>  6 RBPJ-223              19 same_end                     4124 NA          
#>  7 RBPJ-203              19 same_end                     4124 NA          
#>  8 RBPJ-206              15 tandem_utr                      0 NA          
#>  9 RBPJ-225              12 alternative_last_exon          NA NA          
#> 10 RBPJ-229               2 alternative_last_exon          NA NA          
#> # ℹ 55 more rows
usage <- apa_usage(iso, group = "condition")
head(usage)
#> # A tibble: 6 × 13
#>   gene_id        gene_name n_tes_groups apa_types mean_gene_expr distal_fraction
#>   <chr>          <chr>            <int> <chr>              <dbl>           <dbl>
#> 1 ENSG000002352… ENSG0000…            3 intronic            9.62           0.152
#> 2 ENSG000001979… ERO1A               15 alternat…          17.9            0.795
#> 3 ENSG000001295… SNX6                12 alternat…          47.8            0.313
#> 4 ENSG000001433… PIP5K1A             11 alternat…         108.             0.498
#> 5 ENSG000001661… GPT2                 5 alternat…           8.34           0.578
#> 6 ENSG000001078… SFXN3                5 alternat…          54.2            0.598
#> # ℹ 7 more variables: weighted_utr3_length <dbl>,
#> #   distal_fraction_preeclampsia <dbl>, distal_fraction_control <dbl>,
#> #   delta_distal_fraction <dbl>, weighted_utr3_length_preeclampsia <dbl>,
#> #   weighted_utr3_length_control <dbl>, delta_weighted_utr3_length <dbl>
head(test_apa(iso, group = "condition"))
#> # A tibble: 6 × 12
#>   gene_id            gene_name distal_fraction_preeclam…¹ distal_fraction_cont…²
#>   <chr>              <chr>                          <dbl>                  <dbl>
#> 1 ENSG00000143398.22 PIP5K1A                      0.942                  0.862  
#> 2 ENSG00000070831.20 CDC42                        0.131                  0.233  
#> 3 ENSG00000154429.12 CCSAP                        0.961                  0.977  
#> 4 ENSG00000124942.15 AHNAK                        0.00112                0.00275
#> 5 ENSG00000113916.19 BCL6                         0.885                  0.908  
#> 6 ENSG00000197930.15 ERO1A                        0.824                  0.969  
#> # ℹ abbreviated names: ¹​distal_fraction_preeclampsia, ²​distal_fraction_control
#> # ℹ 8 more variables: delta_distal_fraction <dbl>, log_odds_ratio <dbl>,
#> #   se <dbl>, statistic <dbl>, dispersion <dbl>, df <dbl>, pvalue <dbl>,
#> #   padj <dbl>

# Rank genes from the usage table computed above
if (requireNamespace("fgsea", quietly = TRUE)) {
  gmt <- system.file("extdata", "reactome_subset.gmt", package = "isoformic")
  head(run_apa_enrichment(usage, read_gmt(gmt), min_size = 2))
}
#> # A tibble: 6 × 10
#>   pathway          pval  padj log2err     ES    NES  size leadingEdge experiment
#>   <chr>           <dbl> <dbl>   <dbl>  <dbl>  <dbl> <int> <list>      <chr>     
#> 1 REACTOME_ADAP… 0.789  0.866  0.0496  0.616  0.824     2 <chr [1]>   distal_fr…
#> 2 REACTOME_CELL… 0.0461 0.866  0.257   0.864  1.35      8 <chr [1]>   distal_fr…
#> 3 REACTOME_CELL… 0.511  0.866  0.0975 -0.726 -1.04      4 <chr [2]>   distal_fr…
#> 4 REACTOME_CELL… 0.538  0.866  0.0843 -0.755 -1.02      2 <chr [2]>   distal_fr…
#> 5 REACTOME_CELL… 0.657  0.866  0.0557  0.700  0.951     3 <chr [1]>   distal_fr…
#> 6 REACTOME_CELL… 0.657  0.866  0.0557  0.700  0.951     3 <chr [1]>   distal_fr…
#> # ℹ 1 more variable: rank_by <chr>
```
