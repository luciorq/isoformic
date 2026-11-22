# Alternative polyadenylation

A change of 3’ end does not change the protein but changes the 3’ UTR,
its miRNA and RNA-binding protein sites, and sometimes the last exon.
This article ranks the 3’ ends of every gene, classifies each isoform by
the kind of alternative polyadenylation (APA) it represents, quantifies
the usage of the distal site, tests its shift between conditions and
tests gene sets for a shift of their genes.

It starts from the object built in the
[introduction](https://luciorq.github.io/isoformic/articles/isoformic-introduction.md)
and adds the transcript sequences, as in the [productivity
article](https://luciorq.github.io/isoformic/articles/isoformic-productivity.md)
(they are optional here; without them the poly(A) signal columns are
`NA`, and with them a signal is `NA` when no signal hexamer lies within
`signal_window` nucleotides, 50 by default, of the annotated 3’ end;
this is common, because annotated ends are often not the exact cleavage
site):

``` r

library(isoformic)
library(dplyr)

data_dir <- "data/preeclampsia"
annot_dir <- prepare_isoformic_annotation(download_reference(version = "50", file_type = "gff"))
extract_transcript_sequences(
  annot_dir,
  transcript_fasta = download_reference(version = "50", file_type = "fasta")
)
iso <- IsoformicExperiment(
  experiment_name = "preeclampsia_placenta",
  annot_path = annot_dir,
  assays = list(abundance = read_parquet(fs::path(data_dir, "pe_abundance.parquet"))),
  col_data = read_parquet(fs::path(data_dir, "pe_col_data.parquet"))
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
pa[, c("transcript_name", "tes_rank", "n_tes_groups", "apa_type", "utr3_length", "polya_signal", "polya_signal_distance")]
#> # A tibble: 24 × 7
#>    transcript_name tes_rank n_tes_groups apa_type       utr3_length polya_signal
#>    <chr>              <int>        <int> <chr>                <dbl> <chr>       
#>  1 RBPJ-222               1            8 alternative_l…          NA <NA>        
#>  2 RBPJ-243               2            8 alternative_l…          NA <NA>        
#>  3 RBPJ-201               8            8 tandem_utr            4124 <NA>        
#>  4 RBPJ-203               8            8 tandem_utr            4124 <NA>        
#>  5 RBPJ-206               3            8 tandem_utr               0 <NA>        
#>  6 RBPJ-264               2            8 intronic               277 <NA>        
#>  7 RBPJ-258               8            8 other                 4332 <NA>        
#>  8 RBPJ-204               8            8 tandem_utr            4124 <NA>        
#>  9 RBPJ-262               5            8 other                  895 AGTAAA      
#> 10 RBPJ-253               8            8 tandem_utr            4124 <NA>        
#> # ℹ 14 more rows
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
#> 1 ENSG000002358… HLA-DPA1             2 other             0.357            0.362
#> 2 ENSG000002835… PRRT1B               2 other             0.302            0.583
#> 3 ENSG000002764… MYH11                3 intronic…         0.721            0.632
#> 4 ENSG000001640… MST1R                2 other             0.0641           0.5  
#> 5 ENSG000001482… SURF2                2 tandem_u…         6.89             0.520
#> 6 ENSG000002795… ZNF43                2 tandem_u…         0.736            0.510
#> # ℹ 7 more variables: weighted_utr3_length <dbl>,
#> #   distal_fraction_preeclampsia <dbl>, distal_fraction_control <dbl>,
#> #   delta_distal_fraction <dbl>, weighted_utr3_length_preeclampsia <dbl>,
#> #   weighted_utr3_length_control <dbl>, delta_weighted_utr3_length <dbl>
head(test_apa(iso, group = "condition", assay = "abundance"))
#> # A tibble: 6 × 12
#>   gene_id            gene_name distal_fraction_preeclam…¹ distal_fraction_cont…²
#>   <chr>              <chr>                          <dbl>                  <dbl>
#> 1 ENSG00000101166.18 PRELID3B                      0.851                  0.324 
#> 2 ENSG00000186652.12 PRG2                          0.0110                 0.0346
#> 3 ENSG00000163283.7  ALPP                          0.605                  0.787 
#> 4 ENSG00000037280.18 FLT4                          0.413                  0.618 
#> 5 ENSG00000129596.7  CDO1                          0.646                  0.897 
#> 6 ENSG00000137818.14 RPLP1                         0.504                  0.883 
#> # ℹ abbreviated names: ¹​distal_fraction_preeclampsia, ²​distal_fraction_control
#> # ℹ 8 more variables: delta_distal_fraction <dbl>, log_odds_ratio <dbl>,
#> #   se <dbl>, statistic <dbl>, dispersion <dbl>, df <dbl>, pvalue <dbl>,
#> #   padj <dbl>
```

`levels` chooses the groups and their order, as in the productivity
functions: the change is the last level minus the first, and the samples
of other groups are left out.

## Enrichment on 3’ end usage

[`run_apa_enrichment()`](https://luciorq.github.io/isoformic/reference/polyadenylation.md)
ranks the genes by their change in distal usage (`rank_by = "delta"`,
the `delta_distal_fraction` of
[`apa_usage()`](https://luciorq.github.io/isoformic/reference/polyadenylation.md))
or in weighted 3’ UTR length (`rank_by = "utr3_length"`) and runs gene
set enrichment (`fgsea`) on gene symbols. A positive score is a gene set
whose genes move to distal 3’ ends, or to longer 3’ UTRs, in the last
group level. It also takes the table of
[`apa_usage()`](https://luciorq.github.io/isoformic/reference/polyadenylation.md)
for two or more groups, so that one table serves several calls, and one
collection of gene sets or a named list of collections (see the
[enrichment
article](https://luciorq.github.io/isoformic/articles/isoformic-enrichment.md)).

``` r

genesets <- read_gmt(fs::path(data_dir, "c2.cp.reactome.v2023.1.Hs.symbols.gmt.txt"))
usage <- apa_usage(iso, group = "condition")
set.seed(42)
apa_enrichment <- run_apa_enrichment(usage, genesets)
head(apa_enrichment[order(apa_enrichment$padj), c("pathway", "NES", "padj")])
#> # A tibble: 6 × 3
#>   pathway                                                             NES   padj
#>   <chr>                                                             <dbl>  <dbl>
#> 1 REACTOME_SELENOAMINO_ACID_METABOLISM                              -1.80 0.0157
#> 2 REACTOME_EUKARYOTIC_TRANSLATION_ELONGATION                        -1.90 0.0164
#> 3 REACTOME_INTRA_GOLGI_TRAFFIC                                       1.78 0.0276
#> 4 REACTOME_NONSENSE_MEDIATED_DECAY_NMD                              -1.72 0.0370
#> 5 REACTOME_SARS_COV_2_MODULATES_HOST_TRANSLATION_MACHINERY          -1.89 0.0642
#> 6 REACTOME_SRP_DEPENDENT_COTRANSLATIONAL_PROTEIN_TARGETING_TO_MEMB… -1.61 0.0754
set.seed(42)
utr_enrichment <- run_apa_enrichment(usage, genesets, rank_by = "utr3_length")
head(utr_enrichment[order(utr_enrichment$padj), c("pathway", "NES", "padj")])
#> # A tibble: 6 × 3
#>   pathway                                                             NES  padj
#>   <chr>                                                             <dbl> <dbl>
#> 1 REACTOME_ANTIGEN_PROCESSING_UBIQUITINATION_PROTEASOME_DEGRADATION  1.35 0.202
#> 2 REACTOME_EPHA_MEDIATED_GROWTH_CONE_COLLAPSE                        1.70 0.202
#> 3 REACTOME_MEMBRANE_TRAFFICKING                                      1.26 0.202
#> 4 REACTOME_VESICLE_MEDIATED_TRANSPORT                                1.25 0.202
#> 5 REACTOME_INTRA_GOLGI_AND_RETROGRADE_GOLGI_TO_ER_TRAFFIC            1.41 0.287
#> 6 REACTOME_CLASS_I_MHC_MEDIATED_ANTIGEN_PROCESSING_PRESENTATION      1.23 0.418
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
#>  [1] Matrix_1.7-6        gtable_0.3.6        jsonlite_2.0.0     
#>  [4] compiler_4.6.2      tidyselect_1.2.1    Rcpp_1.1.2         
#>  [7] parallel_4.6.2      scales_1.4.0        BiocParallel_1.46.0
#> [10] yaml_2.3.12         fastmap_1.2.0       lattice_0.23-1     
#> [13] ggplot2_4.0.4       R6_2.6.1            generics_0.1.4     
#> [16] curl_8.1.0          knitr_1.52          tibble_3.3.1       
#> [19] DBI_1.3.0           RColorBrewer_1.1-3  pillar_1.11.1      
#> [22] rlang_1.3.0         fastmatch_1.1-8     utf8_1.2.6         
#> [25] xfun_0.62           fs_2.1.0            S7_1.0.0           
#> [28] otel_0.3.0          cli_3.6.6           withr_3.0.3        
#> [31] magrittr_2.0.5      grid_4.6.2          digest_0.6.39      
#> [34] fgsea_1.38.0        cowplot_1.2.0       lifecycle_1.0.5    
#> [37] vctrs_0.7.3         data.table_1.18.6.1 evaluate_1.0.5     
#> [40] glue_1.8.1          farver_2.1.2        duckdb_2.0.0       
#> [43] codetools_0.2-20    rmarkdown_2.32      tools_4.6.2        
#> [46] pkgconfig_2.0.3     htmltools_0.5.9
```
