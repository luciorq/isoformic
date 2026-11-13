# Alternative polyadenylation

A change of 3’ end does not change the protein but changes the 3’ UTR,
its miRNA and RNA-binding protein sites, and sometimes the last exon.
This article ranks the 3’ ends of every gene, classifies each isoform by
the kind of alternative polyadenylation (APA) it represents, quantifies
the usage of the distal site and tests its shift between conditions.

It starts from the object built in the
[introduction](https://luciorq.github.io/isoformic/articles/isoformic-introduction.md),
with the transcript sequences added in the [productivity
article](https://luciorq.github.io/isoformic/articles/isoformic-productivity.md)
(they are optional here; without them the poly(A) signal columns are
`NA`):

``` r

library(isoformic)
library(dplyr)

extdata <- fs::path_package("isoformic", "extdata")
annot_dir <- prepare_isoformic_annotation(download_reference(version = "33", file_type = "gff"))
extract_transcript_sequences(
  annot_dir,
  transcript_fasta = download_reference(version = "33", file_type = "fasta")
)
iso <- IsoformicExperiment(
  experiment_name = "preeclampsia_placenta",
  annot_path = annot_dir,
  assays = list(abundance = read_parquet(fs::path(extdata, "pe_abundance.parquet"))),
  col_data = read_parquet(fs::path(extdata, "pe_col_data.parquet"))
)
```

## Poly(A) sites per isoform

[`polya_sites()`](https://luciorq.github.io/isoformic/reference/polyadenylation.md)
ranks the 3’ ends of each gene along the direction of transcription
(`tes_rank`, 1 for the most proximal), classifies every isoform relative
to the canonical one and reads the poly(A) signal hexamer closest to the
end from the sequence:

- `same_end`: ends at the canonical site;
- `tandem_utr`: same last exon, different cleavage site;
- `alternative_last_exon`: ends in another last exon;
- `intronic`: ends inside an intron of the canonical isoform;
- `other`.

``` r

pa <- polya_sites(iso, gene_name = "RBPJ")
#> Warning: ! 91 transcripts have a sequence length differing from the exon length; their
#>   ORFs are skipped.
pa[, c("transcript_name", "tes_rank", "n_tes_groups", "apa_type", "utr3_length", "polya_signal", "polya_signal_distance")]
#> # A tibble: 10 × 7
#>    transcript_name tes_rank n_tes_groups apa_type       utr3_length polya_signal
#>    <chr>              <int>        <int> <chr>                <dbl> <chr>       
#>  1 RBPJ-222               1            6 alternative_l…          NA <NA>        
#>  2 RBPJ-201               6            6 same_end              4121 <NA>        
#>  3 RBPJ-206               3            6 tandem_utr               0 <NA>        
#>  4 RBPJ-225               2            6 alternative_l…         347 <NA>        
#>  5 RBPJ-204               4            6 tandem_utr              59 AAAAAG      
#>  6 RBPJ-205               6            6 same_end              4121 <NA>        
#>  7 RBPJ-214               5            6 tandem_utr             176 <NA>        
#>  8 RBPJ-208               5            6 other                  357 <NA>        
#>  9 RBPJ-202               6            6 same_end              4124 <NA>        
#> 10 RBPJ-210               5            6 tandem_utr              NA <NA>        
#> # ℹ 1 more variable: polya_signal_distance <int>
```

With `polya_atlas = "polyasite"` the ends are checked against the
PolyASite 2.0 atlas (one of the frozen datasets, downloaded once); the
GENCODE poly(A) feature file (`download_reference(file_type = "polya")`)
or any table of cleavage positions work the same way:

``` r

polya_sites(iso, polya_atlas = "polyasite")[, c("transcript_name", "pas_supported", "pas_class", "pas_distance")]
```

## Distal usage per gene

[`apa_usage()`](https://luciorq.github.io/isoformic/reference/polyadenylation.md)
quantifies, per gene and group, the share of expression ending at the
distal-most site and the abundance-weighted 3’ UTR length;
[`test_apa()`](https://luciorq.github.io/isoformic/reference/polyadenylation.md)
tests the shift in distal usage between two groups with the same
quasi-binomial model as
[`test_productive_fraction()`](https://luciorq.github.io/isoformic/reference/productivity.md).
The test runs on `counts` when the object has them; this object only
carries TPM.

``` r

head(apa_usage(iso, group = "condition"))
#> # A tibble: 6 × 13
#>   gene_id        gene_name n_tes_groups apa_types mean_gene_expr distal_fraction
#>   <chr>          <chr>            <int> <chr>              <dbl>           <dbl>
#> 1 ENSG000001302… RPL36                3 alternat…        1693.             0.423
#> 2 ENSG000001334… MYO18B               2 intronic            6.75           0.364
#> 3 ENSG000001203… SEC16B               2 alternat…           9.63           0.439
#> 4 ENSG000002550… ZFP91-CN…            2 other               9.66           0.276
#> 5 ENSG000002127… LINC02693            2 other              11.2            0.381
#> 6 ENSG000001830… RGPD6                2 alternat…         741.             0.820
#> # ℹ 7 more variables: weighted_utr3_length <dbl>,
#> #   distal_fraction_preeclampsia <dbl>, distal_fraction_control <dbl>,
#> #   delta_distal_fraction <dbl>, weighted_utr3_length_preeclampsia <dbl>,
#> #   weighted_utr3_length_control <dbl>, delta_weighted_utr3_length <dbl>
head(test_apa(iso, group = "condition", assay = "abundance"))
#> # A tibble: 6 × 12
#>   gene_id            gene_name distal_fraction_preeclam…¹ distal_fraction_cont…²
#>   <chr>              <chr>                          <dbl>                  <dbl>
#> 1 ENSG00000164587.13 RPS14                          0.485                 0.116 
#> 2 ENSG00000128609.15 NDUFA5                         0.376                 0.0833
#> 3 ENSG00000156508.18 EEF1A1                         0.178                 0.0630
#> 4 ENSG00000132475.10 H3-3B                          0.504                 0.669 
#> 5 ENSG00000170144.20 HNRNPA3                        0.518                 0.667 
#> 6 ENSG00000152926.15 ZNF117                         0.517                 0.343 
#> # ℹ abbreviated names: ¹​distal_fraction_preeclampsia, ²​distal_fraction_control
#> # ℹ 8 more variables: delta_distal_fraction <dbl>, log_odds_ratio <dbl>,
#> #   se <dbl>, statistic <dbl>, dispersion <dbl>, df <dbl>, pvalue <dbl>,
#> #   padj <dbl>
```

## Session information

``` r

sessionInfo()
#> R version 4.6.2 (2026-10-29)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 26.04.1 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.32.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] dplyr_1.2.1     isoformic_0.1.3
#> 
#> loaded via a namespace (and not attached):
#>  [1] vctrs_0.7.3      cli_3.6.6        knitr_1.52       rlang_1.3.0     
#>  [5] xfun_0.62        DBI_1.3.0        otel_0.3.0       generics_0.1.4  
#>  [9] S7_0.2.2         jsonlite_2.0.0   glue_1.8.1       htmltools_0.5.9 
#> [13] rmarkdown_2.32   evaluate_1.0.5   tibble_3.3.1     fastmap_1.2.0   
#> [17] yaml_2.3.12      lifecycle_1.0.5  duckdb_2.0.0     compiler_4.6.2  
#> [21] fs_2.1.0         pkgconfig_2.0.3  digest_0.6.39    R6_2.6.1        
#> [25] utf8_1.2.6       tidyselect_1.2.1 pillar_1.11.1    curl_8.1.0      
#> [29] magrittr_2.0.5   tools_4.6.2      withr_3.0.3
```
