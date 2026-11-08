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

  Passed on to methods.

- value:

  A data frame of results (`de_tx<-`, `de_gene<-`) or a list of cutoffs
  (`de_cutoffs<-`).

## Value

Tibbles, or the modified object for replacement functions.

## Details

Recognised column synonyms: `log2FoldChange`/`logFC` for `log2FC`;
`padj`/`FDR`/`adj.P.Val`/`svalue` for `qvalue`; `PValue`/`P.Value` for
`pvalue`. Row names are used as IDs when no ID column is present, and
GENCODE-style `|`-suffixed transcript IDs are truncated. When `qvalue`
is absent it is computed with Benjamini-Hochberg adjustment of `pvalue`.

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
#>   transcript_id     log2FC      pvalue    qvalue det_sig
#>   <chr>              <dbl>       <dbl>     <dbl> <chr>  
#> 1 ENST00000638463.1   8.57 0.0226      0.102     no     
#> 2 ENST00000557508.5  -5.47 0.000000138 0.0000236 yes    
#> 3 ENST00000166139.9   5.01 0.000000138 0.0000236 yes    
#> 4 ENST00000592947.5   4.85 0.00000332  0.000211  yes    
#> 5 ENST00000508871.5  -4.61 0.0116      0.0658    no     
#> 6 ENST00000615840.4   3.82 0.000000138 0.0000236 yes    
head(de_gene(iso, sig_only = TRUE))
#> # A tibble: 6 × 8
#>   gene_id            log10mean  stat log2FC      pvalue   locfdr  qvalue deg_sig
#>   <chr>                  <dbl> <dbl>  <dbl>       <dbl>    <dbl>   <dbl> <chr>  
#> 1 ENSG00000070404.10      4.20  68     5.23 0.000000462  5.48e-5 2.58e-5 yes    
#> 2 ENSG00000130822.15      1.78  64.5   3.39 0.00000139   1.33e-4 4.69e-5 yes    
#> 3 ENSG00000102755.12      4.66  63.4   3.37 0.00000277   2.56e-4 8.56e-5 yes    
#> 4 ENSG00000070388.11      1.89  68.0   2.56 0.000000462  5.48e-5 2.58e-5 yes    
#> 5 ENSG00000061656.10      2.15  67.1   2.54 0.000000462  5.48e-5 2.58e-5 yes    
#> 6 ENSG00000106991.13      4.30  67     2.49 0.000000462  5.48e-5 2.58e-5 yes    
head(de_table(iso))
#> # A tibble: 6 × 9
#>   feature_id  feature_name feature_type gene_id gene_name log2FC  pvalue  qvalue
#>   <chr>       <chr>        <chr>        <chr>   <chr>      <dbl>   <dbl>   <dbl>
#> 1 ENSG000000… FSTL3        gene         ENSG00… FSTL3       5.23 4.62e-7 2.58e-5
#> 2 ENSG000001… PNCK         gene         ENSG00… PNCK        3.39 1.39e-6 4.69e-5
#> 3 ENSG000001… FLT1         gene         ENSG00… FLT1        3.37 2.77e-6 8.56e-5
#> 4 ENSG000000… FGF22        gene         ENSG00… FGF22       2.56 4.62e-7 2.58e-5
#> 5 ENSG000000… SPAG4        gene         ENSG00… SPAG4       2.54 4.62e-7 2.58e-5
#> 6 ENSG000001… ENG          gene         ENSG00… ENG         2.49 4.62e-7 2.58e-5
#> # ℹ 1 more variable: is_de <chr>

de_cutoffs(iso) <- list(log2fc = 0.5)
sum(de_tx(iso)$det_sig == "yes")
#> [1] 135
```
