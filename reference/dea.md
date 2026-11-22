# Differential Expression Results in an IsoformicExperiment

An
[IsoformicExperiment](https://luciorq.github.io/isoformic/reference/IsoformicExperiment.md)
stores differential expression results for transcripts (`det`) and genes
(`deg`) in its `dea` property, together with the significance cutoffs
used to flag features. Results can come from any tool: on assignment the
tables are standardised to the column names `log2FC`, `pvalue` and
`qvalue` (see Details), checked against the annotation, and stored as
tibbles.

- `de_tx(x)` / `de_gene(x)`: results with an added `det_sig` / `deg_sig`
  column (`"yes"`/`"no"`) computed from the cutoffs, sorted by absolute
  `log2FC`. Cutoffs can be overridden per call.

- `de_tx(x) <- value` / `de_gene(x) <- value`: set results.

- `de_cutoffs(x)` / `de_cutoffs(x) <- value`: the stored cutoffs, a list
  with `log2fc` (default 1), `pvalue` (default 0.05) and `use_fdr`
  (default `TRUE`, i.e. compare `qvalue` rather than `pvalue`).

- `de_table(x)`: a long table combining gene and transcript results, one
  row per feature, with `feature_id`, `feature_name`, `feature_type`
  (`"gene"` or the transcript type), `gene_id`, `gene_name`, `log2FC`,
  `pvalue`, `qvalue` and `is_de`. This is the input of
  [`plot_log2fc()`](https://luciorq.github.io/isoformic/reference/plot_log2fc.md).

## Usage

``` r
de_tx(
  self,
  log2fc_cutoff = NULL,
  pvalue_cutoff = NULL,
  use_fdr = NULL,
  sig_only = FALSE,
  ...
)

de_gene(
  self,
  log2fc_cutoff = NULL,
  pvalue_cutoff = NULL,
  use_fdr = NULL,
  sig_only = FALSE,
  ...
)

de_tx(self) <- value

de_gene(self) <- value

de_cutoffs(self, ...)

de_cutoffs(self) <- value

de_table(self, ...)
```

## Arguments

- self:

  An
  [IsoformicExperiment](https://luciorq.github.io/isoformic/reference/IsoformicExperiment.md).

- log2fc_cutoff, pvalue_cutoff, use_fdr:

  Optional per-call overrides of the stored cutoffs.

- sig_only:

  Logical, return only significant features.

- ...:

  Not used; reserved for methods.

- value:

  A data frame of results (`de_tx<-`, `de_gene<-`) or a list of cutoffs
  (`de_cutoffs<-`).

## Value

Tibbles, or the modified object for replacement functions.

## Details

Recognised column synonyms: `log2FoldChange`/`logFC`/`log2fc`/`log2_fc`
for `log2FC`; `padj`/`FDR`/`adj.P.Val`/`svalue`/`fdr`/`q_value` for
`qvalue`; `PValue`/`P.Value`/`pval`/`p_value` for `pvalue`. The ID
column (`transcript_id` or `gene_id`) is also found as `feature_id`,
`id` or `ID` (and `tx_id`/`transcript` or `gene`). Row names are used as
IDs when no ID column is present, and GENCODE-style `|`-suffixed
transcript IDs are truncated. When `qvalue` is absent it is computed
with Benjamini-Hochberg adjustment of `pvalue`.

`de_cutoffs<-` takes a named list with any of `log2fc`, `pvalue` and
`use_fdr` (also as `log2FC`/`lfc`, `padj` and `fdr`); missing entries
take the defaults and other names are an error.

The `dea` property also holds the uncertainty tables of
[`add_uncertainty()`](https://luciorq.github.io/isoformic/reference/uncertainty.md)
(`uncertainty`, `ambiguity`, and `confidence`, the fold-change table of
`add_uncertainty(group = )`).

## Examples

``` r
iso <- example_isoformic()
de_cutoffs(iso)
#> $log2fc
#> [1] 1
#> 
#> $pvalue
#> [1] 0.05
#> 
#> $use_fdr
#> [1] TRUE
#> 
head(de_tx(iso))
#> # A tibble: 6 × 5
#>   transcript_id     log2FC       pvalue    qvalue det_sig
#>   <chr>              <dbl>        <dbl>     <dbl> <chr>  
#> 1 ENST00000308868.5   6.73 0.0000000789 0.0000198 yes    
#> 2 ENST00000965599.2   6.42 0.0000000789 0.0000198 yes    
#> 3 ENST00000450046.3  -6.05 0.0000000789 0.0000198 yes    
#> 4 ENST00001134423.1   5.95 0.0460       0.193     no     
#> 5 ENST00001002978.1   5.64 0.000384     0.00772   yes    
#> 6 ENST00001035611.1  -5.52 0.0000000789 0.0000198 yes    
head(de_gene(iso, sig_only = TRUE))
#> # A tibble: 6 × 8
#>   gene_id            log10mean  stat log2FC      pvalue   locfdr  qvalue deg_sig
#>   <chr>                  <dbl> <dbl>  <dbl>       <dbl>    <dbl>   <dbl> <chr>  
#> 1 ENSG00000174697.7       3.97  68     6.58 0.000000451 0.000450 3.57e-4 yes    
#> 2 ENSG00000070404.12      4.13  68     4.85 0.000000451 0.000450 3.57e-4 yes    
#> 3 ENSG00000169495.6       3.73  68     4.40 0.000000451 0.000450 3.57e-4 yes    
#> 4 ENSG00000102755.15      4.67  64.6   3.47 0.0000135   0.00285  1.38e-3 yes    
#> 5 ENSG00000169994.22      2.79  68     3.16 0.000000451 0.000450 3.57e-4 yes    
#> 6 ENSG00000104827.13      4.24  68     3.05 0.000000451 0.000450 3.57e-4 yes    
head(de_table(iso))
#> # A tibble: 6 × 9
#>   feature_id  feature_name feature_type gene_id gene_name log2FC  pvalue  qvalue
#>   <chr>       <chr>        <chr>        <chr>   <chr>      <dbl>   <dbl>   <dbl>
#> 1 ENSG000001… LEP          gene         ENSG00… LEP         6.58 4.51e-7 3.57e-4
#> 2 ENSG000000… FSTL3        gene         ENSG00… FSTL3       4.85 4.51e-7 3.57e-4
#> 3 ENSG000001… HTRA4        gene         ENSG00… HTRA4       4.40 4.51e-7 3.57e-4
#> 4 ENSG000001… FLT1         gene         ENSG00… FLT1        3.47 1.35e-5 1.38e-3
#> 5 ENSG000001… MYO7B        gene         ENSG00… MYO7B       3.16 4.51e-7 3.57e-4
#> 6 ENSG000001… CGB3         gene         ENSG00… CGB3        3.05 4.51e-7 3.57e-4
#> # ℹ 1 more variable: is_de <chr>

de_cutoffs(iso) <- list(log2fc = 0.5)
sum(de_tx(iso)$det_sig == "yes")
#> [1] 307
```
