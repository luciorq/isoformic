# IsoformicExperiment: Transcript, Gene and Exon Level Data in One Object

`IsoformicExperiment` is an S7 class modelled after
`SummarizedExperiment`, but organised around *transcriptomic diversity*:
the same samples are described at three levels that are kept consistent
by construction.

- **Transcript level**: assays are numeric matrices with transcript IDs
  as row names and sample IDs as column names (for example `abundance`,
  `counts`, `length`). All assays share the same dimnames.

- **Gene level**: gene assays are derived from the transcript assays by
  summing the isoforms of each gene (abundance-weighted means for
  `length`), cached, and can be overridden with gene_assays\<-.

- **Exon level**: exon structures come from the annotation and are
  available through
  [`annot_data_exons()`](https://luciorq.github.io/isoformic/reference/IsoformicExperiment-accessors.md),
  [`exon_table()`](https://luciorq.github.io/isoformic/reference/IsoformicExperiment-accessors.md)
  and
  [`exon_diversity()`](https://luciorq.github.io/isoformic/reference/diversity.md).

Annotation is never copied into the object. The `annot_path` property
points to an *annotation directory* of Parquet files created by
[`prepare_isoformic_annotation()`](https://luciorq.github.io/isoformic/reference/prepare_isoformic_annotation.md)
(a GFF3/GTF file is accepted and converted on the fly). Tables are read
with DuckDB on first access and cached in memory keyed by file path,
size and modification time, so a re-created annotation is picked up
automatically and several objects sharing one annotation share one copy
of the tables.

Derived results (gene assays, row data, isoform fractions) are cached
per object in the `.cache` environment, keyed by a hash of their inputs,
and therefore survive copies but never go stale. See
[`clear_cache()`](https://luciorq.github.io/isoformic/reference/IsoformicExperiment-accessors.md).

Objects persist to disk as a directory of Parquet files with
[`write_isoformic()`](https://luciorq.github.io/isoformic/reference/isoformic-io.md)
and are restored with
[`read_isoformic()`](https://luciorq.github.io/isoformic/reference/isoformic-io.md).

## Usage

``` r
IsoformicExperiment(
  experiment_name = NULL,
  annot_path = NULL,
  assays = list(),
  col_data = NULL,
  gene_assays = list(),
  dea = list(),
  gsea = list(),
  metadata = list(),
  data_path = NULL,
  tx_type_palette = NULL
)
```

## Arguments

- experiment_name:

  Single string naming the experiment.

- annot_path:

  Annotation directory or GFF3/GTF file.

- assays:

  Named list of transcript-level matrices (or data frames with a
  `transcript_id` column).

- col_data:

  Data frame of sample metadata with a `sample_id` column (or sample IDs
  as row names / first column).

- gene_assays:

  Optional named list of gene-level matrices.

- dea:

  Optional list of differential expression results.

- gsea:

  Optional named list of enrichment results.

- metadata:

  Optional named list of free-form metadata.

- data_path:

  Optional on-disk location of the experiment.

- tx_type_palette:

  Optional named colour vector per transcript type.

## Value

An `IsoformicExperiment` object.

## Properties

- `experiment_name`:

  Single string. Generated when not supplied.

- `annot_path`:

  Path to the annotation directory, or `NULL` for an empty object.

- `data_path`:

  Directory the object was written to or read from (see
  [`write_isoformic()`](https://luciorq.github.io/isoformic/reference/isoformic-io.md)),
  or `NULL`.

- `assays`:

  Named list of transcript-level matrices.

- `gene_assays`:

  Named list of gene-level matrices. Derived from `assays` when not set
  explicitly.

- `col_data`:

  Tibble of sample metadata; its `sample_id` column matches the assay
  column names.

- `dea`:

  List with differential expression results: `det` (transcripts), `deg`
  (genes) and `cutoffs`. Use de_tx\<- and de_gene\<- to set them.

- `gsea`:

  Named list of enrichment result tables.

- `metadata`:

  Free-form named list.

- `annot_metadata`:

  Named list describing the annotation source, read from the annotation
  directory unless set explicitly.

- `tx_type_palette`:

  Named character vector of colours per transcript type, used by all
  plots. Defaults to
  [`tx_type_palette()`](https://luciorq.github.io/isoformic/reference/tx_type_palette.md).

- `annot_data_transcripts`, `annot_data_genes`, `annot_data_exons`:

  Read-only, cached annotation tables.

## See also

[`example_isoformic()`](https://luciorq.github.io/isoformic/reference/example_isoformic.md)
for a ready-made object,
[`read_salmon()`](https://luciorq.github.io/isoformic/reference/read_salmon.md)
and
[`as_isoformic()`](https://luciorq.github.io/isoformic/reference/isoformic-conversion.md)
to create one from quantifications,
[IsoformicExperiment-accessors](https://luciorq.github.io/isoformic/reference/IsoformicExperiment-accessors.md)
for the accessor generics.

## Examples

``` r
iso <- example_isoformic()
iso
#> <IsoformicExperiment> pe_placenta_example
#>   annotation : gencode_v33_example_subset (69 genes, 892 transcripts, 7,032 exons)
#>   samples    : 25 [condition]
#>   transcripts: 499 assayed from 69 genes
#>   assays     : abundance, counts
#>   gene assays: derived from transcript assays on demand
#>   DE results : det (499 transcripts), deg (69 genes)
dim(iso)
#> [1] 499  25
assay_names(iso)
#> [1] "abundance" "counts"   
head(row_data(iso))
#> # A tibble: 6 × 19
#>   transcript_id     transcript_name transcript_type  gene_id gene_name gene_type
#>   <chr>             <chr>           <chr>            <chr>   <chr>     <chr>    
#> 1 ENST00000374514.8 LYPLA2-205      protein_coding   ENSG00… LYPLA2    protein_…
#> 2 ENST00000420982.5 LYPLA2-206      protein_coding   ENSG00… LYPLA2    protein_…
#> 3 ENST00000374502.7 LYPLA2-202      protein_coding   ENSG00… LYPLA2    protein_…
#> 4 ENST00000374501.1 LYPLA2-201      protein_coding   ENSG00… LYPLA2    protein_…
#> 5 ENST00000436063.7 DNTTIP2-202     protein_coding   ENSG00… DNTTIP2   protein_…
#> 6 ENST00000496535.6 DNTTIP2-205     processed_trans… ENSG00… DNTTIP2   protein_…
#> # ℹ 13 more variables: source_gene_id <chr>, source_transcript_id <chr>,
#> #   tags <chr>, coding_status <chr>, transcript_support_level <chr>,
#> #   tie_score <dbl>, cds_start <int>, cds_end <int>, cds_length <int>,
#> #   is_canonical <lgl>, haplotype <chr>, n_exons <int>, tx_length <int>
head(gene_assay(iso, "abundance"))
#>                    SRR11498039 SRR11498040 SRR11498041 SRR11498042 SRR11498043
#> ENSG00000005007.12    3253.512    2563.848    2656.175    3295.886    4497.997
#> ENSG00000005206.17    1727.018    1515.008    1807.224    1551.435     986.158
#> ENSG00000007384.15    1970.727    1612.095    2188.516    2179.135    1910.157
#> ENSG00000007541.16     550.154     526.270     745.990     723.957     498.251
#> ENSG00000008838.20    1058.643    1018.991     978.694    1061.083    1152.704
#> ENSG00000009413.16     724.057     599.019     675.888     696.338     419.116
#>                    SRR11498044 SRR11498045 SRR11498046 SRR11498063 SRR11498064
#> ENSG00000005007.12    4732.194    2686.639    3904.654    3335.450    3699.352
#> ENSG00000005206.17    1850.509     609.902    1492.657    2069.531    2048.620
#> ENSG00000007384.15    2708.052    2264.347    3182.521    2966.487    3382.728
#> ENSG00000007541.16     962.738     245.692     969.690     958.960    1220.914
#> ENSG00000008838.20    1700.362    1002.685    1392.596    1419.984    1543.717
#> ENSG00000009413.16     824.227    1010.559    1235.516     837.198    1042.775
#>                    SRR11498065 SRR11498066 SRR11498067 SRR11498068 SRR11498069
#> ENSG00000005007.12    2405.434    2521.298    2951.488    2412.240    3125.850
#> ENSG00000005206.17    1502.081    1117.840    1625.626    1202.581    1971.119
#> ENSG00000007384.15    2154.739    1830.455    2151.931    2016.009    2502.582
#> ENSG00000007541.16     758.946     604.632     712.162     563.255     840.969
#> ENSG00000008838.20    1133.716     996.395    1158.077     867.589    1375.700
#> ENSG00000009413.16     981.165     736.211     844.466     923.871     996.592
#>                    SRR11498070 SRR11498071 SRR11498072 SRR11498073 SRR11498074
#> ENSG00000005007.12    4497.310    2265.091    3065.700    2692.038    3329.921
#> ENSG00000005206.17    2303.939    1076.992    1372.618    1172.162    1693.740
#> ENSG00000007384.15    3482.066    1743.023    2339.072    2031.263    2347.569
#> ENSG00000007541.16    1170.131     528.531     626.839     608.772     803.845
#> ENSG00000008838.20    1783.556     927.626     953.281    1095.061    1114.967
#> ENSG00000009413.16    1375.634     715.627     462.376     620.009     733.193
#>                    SRR11498075 SRR11498076 SRR11498077 SRR11498078 SRR11498079
#> ENSG00000005007.12    3200.173    3211.790    2886.842    2522.062    3415.826
#> ENSG00000005206.17    1728.208    1544.569    1448.556    1375.868    1866.532
#> ENSG00000007384.15    2321.575    2703.154    2317.905    2011.830    3311.084
#> ENSG00000007541.16     701.155     760.047     646.552     600.129     922.636
#> ENSG00000008838.20    1117.793    1124.630    1119.471     952.758    1333.459
#> ENSG00000009413.16     924.570     742.140     629.304     927.810    1091.680

# Subset like a SummarizedExperiment
iso[1:10, c("SRR11498039", "SRR11498063")]
#> <IsoformicExperiment> pe_placenta_example
#>   annotation : gencode_v33_example_subset (69 genes, 892 transcripts, 7,032 exons)
#>   samples    : 2 [condition]
#>   transcripts: 10 assayed from 3 genes
#>   assays     : abundance, counts
#>   gene assays: derived from transcript assays on demand
#>   DE results : det (499 transcripts), deg (69 genes)

# Build one from matrices
annot <- system.file("extdata", "example_experiment", "annotation", package = "isoformic")
tpm <- assay(iso, "abundance")[1:5, 1:3]
small <- IsoformicExperiment(
  experiment_name = "tiny",
  annot_path = annot,
  assays = list(abundance = tpm),
  col_data = data.frame(sample_id = colnames(tpm), condition = c("a", "a", "b"))
)
small
#> <IsoformicExperiment> tiny
#>   annotation : gencode_v33_example_subset (69 genes, 892 transcripts, 7,032 exons)
#>   samples    : 3 [condition]
#>   transcripts: 5 assayed from 2 genes
#>   assays     : abundance
#>   gene assays: derived from transcript assays on demand
#>   DE results : <none>
```
