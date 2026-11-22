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

Derived results (gene assays, row data, isoform fractions, ORFs,
transcript features and productivity classes) are cached per object in
the `.cache` environment, keyed by a hash of their inputs, and therefore
survive copies but never go stale. See
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

  Single string naming the experiment; a random name is generated when
  `NULL`.

- annot_path:

  Annotation directory or GFF3/GTF file.

- assays:

  Named list of transcript-level matrices (or data frames with a
  `transcript_id` column).

- col_data:

  Data frame of sample metadata with a `sample_id` column (or sample IDs
  as row names / first column). `NULL` builds it from the column names
  of the assays.

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
  (genes) and `cutoffs`, and the uncertainty tables of
  [`add_uncertainty()`](https://luciorq.github.io/isoformic/reference/uncertainty.md).
  Use de_tx\<- and de_gene\<- to set the results.

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

## Subsetting

`x[i, j]` keeps the transcripts `i` and the samples `j`, given as IDs,
positions or logical vectors; a missing `i` or `j` keeps all of them,
and the subset must hold at least one transcript and one sample.
Everything stored by transcript or by sample follows: the assays, the
sample table, stored gene assays (the genes of the kept transcripts) and
the inferential replicates in `x@metadata$inf_reps`. The uncertainty
tables of
[`add_uncertainty()`](https://luciorq.github.io/isoformic/reference/uncertainty.md)
are restricted to the kept transcripts (a collapsed ambiguity group
stays only with all its members); when samples are dropped they are not
recomputed, and a message says so. Differential expression and
enrichment results are kept as they are.

The inferential replicates are not copied: the `metadata$inf_reps` of a
subset holds the draws of `x` and a `selection` attribute, a list with
the IDs of the kept `transcripts` and `samples`. The functions of the
package (and
[`inferential_replicates_table()`](https://luciorq.github.io/isoformic/reference/inferential-io.md))
read the draws through that selection, and
[`write_isoformic()`](https://luciorq.github.io/isoformic/reference/isoformic-io.md)
writes only the selected rows and columns, so the results are those of a
copy of the kept draws, while the memory of a subset does not grow with
the draws. The trade-off: the draws of `x` (every transcript and sample)
stay in memory as long as a subset made from it exists, even after `x`
itself is removed, and they go with every serialised copy of the subset:
[`saveRDS()`](https://rdrr.io/r/base/readRDS.html),
[`serialize()`](https://rdrr.io/r/base/serialize.html) and the objects
sent to parallel workers that do not share memory (such as `future` or
`mirai` workers) carry every transcript and sample of the draws of `x`,
while
[`write_isoformic()`](https://luciorq.github.io/isoformic/reference/isoformic-io.md)
writes only the selected ones. Code that indexes the matrices of
`metadata$inf_reps` directly sees the draws of `x`. A list operation on
`metadata$inf_reps` (such as
[`lapply()`](https://rdrr.io/r/base/lapply.html) or `[`) drops the
`selection` attribute and returns the draws of `x`: the functions that
take the object still read only its kept transcripts and samples, while
code that gets only the list, such as
[`inferential_replicates_table()`](https://luciorq.github.io/isoformic/reference/inferential-io.md)
without its `x` argument, reads every transcript and sample of `x`. When
no kept transcript has draws, the subset keeps only the replicate
moments (with zero-row draws), so that
[`scale_counts()`](https://luciorq.github.io/isoformic/reference/inferential-overdispersion.md)
and
[`inferential_overdispersion()`](https://luciorq.github.io/isoformic/reference/inferential-overdispersion.md)
still work, but the functions that need the draws treat it as an object
without inferential replicates; without moments, `metadata$inf_reps` is
`NULL`. To replace the draws by plain matrices of the kept draws (and
free those of `x`), write the subset with
[`write_isoformic()`](https://luciorq.github.io/isoformic/reference/isoformic-io.md)
and read it back; to drop them, set `metadata$inf_reps <- NULL`. To
compare two groups of a larger object, pass `levels` to the function
instead of subsetting: the samples of the other groups are ignored. A
subset with other transcripts starts with an empty cache of derived
tables (the transcript annotation join is shared).

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
#>   annotation : gencode_v50_example_subset (69 genes, 1,914 transcripts, 18,656 exons)
#>   samples    : 25 [condition]
#>   transcripts: 1,911 assayed from 69 genes
#>   assays     : abundance, counts
#>   gene assays: derived from transcript assays on demand
#>   DE results : det (710 transcripts), deg (69 genes)
dim(iso)
#> [1] 1911   25
assay_names(iso)
#> [1] "abundance" "counts"   
head(row_data(iso))
#> # A tibble: 6 × 19
#>   transcript_id     transcript_name transcript_type gene_id  gene_name gene_type
#>   <chr>             <chr>           <chr>           <chr>    <chr>     <chr>    
#> 1 ENST00001127150.1 CDC42-254       protein_coding  ENSG000… CDC42     protein_…
#> 2 ENST00001127151.1 CDC42-255       protein_coding  ENSG000… CDC42     protein_…
#> 3 ENST00001091268.1 CDC42-249       protein_coding  ENSG000… CDC42     protein_…
#> 4 ENST00001127149.1 CDC42-253       protein_coding  ENSG000… CDC42     protein_…
#> 5 ENST00000695857.2 CDC42-218       protein_coding  ENSG000… CDC42     protein_…
#> 6 ENST00000695860.2 CDC42-221       protein_coding  ENSG000… CDC42     protein_…
#> # ℹ 13 more variables: source_gene_id <chr>, source_transcript_id <chr>,
#> #   tags <chr>, coding_status <chr>, transcript_support_level <chr>,
#> #   tie_score <dbl>, cds_start <int>, cds_end <int>, cds_length <int>,
#> #   is_canonical <lgl>, haplotype <chr>, n_exons <int>, tx_length <int>
head(gene_assay(iso, "abundance"))
#>                    SRR11498039 SRR11498040 SRR11498041 SRR11498042 SRR11498043
#> ENSG00000070404.12  2569.68384   299.07857  1417.10732   603.17133   471.63517
#> ENSG00000070831.20   122.70026    99.10083    80.20022    86.70524   124.04831
#> ENSG00000092199.21   266.25891   326.93814   276.23063   273.64675   514.47320
#> ENSG00000100567.15    83.58850    89.24681    83.38840    83.66714    89.96554
#> ENSG00000101181.20    19.77858    22.22232    26.51377    22.33130    22.33756
#> ENSG00000102755.15  6674.69963   294.43483  2383.73606   559.48987   632.91902
#>                    SRR11498044 SRR11498045 SRR11498046 SRR11498063 SRR11498064
#> ENSG00000070404.12   310.76213   896.73093  3187.15377    29.18200    31.38722
#> ENSG00000070831.20    78.78730   204.01240   107.61872    52.03164    35.00855
#> ENSG00000092199.21   289.21006   525.55022   327.02047   267.16888   228.18136
#> ENSG00000100567.15    84.23372   116.61350    96.74124   166.90329   113.95220
#> ENSG00000101181.20    21.93102    30.98572    19.77071    29.64593    26.51956
#> ENSG00000102755.15  1899.57688  3828.13284  3580.22704   168.79459   124.46575
#>                    SRR11498065 SRR11498066 SRR11498067 SRR11498068 SRR11498069
#> ENSG00000070404.12    48.99216    14.16329    66.08116    53.78511    13.59521
#> ENSG00000070831.20    38.44176    48.57682    51.06573    43.95347    43.87882
#> ENSG00000092199.21   239.72877   240.21498   277.28130   254.05305   257.29341
#> ENSG00000100567.15   118.80365   121.72739   161.33785   131.05597   111.49771
#> ENSG00000101181.20    21.75279    20.98340    28.30428    28.21587    23.79867
#> ENSG00000102755.15   142.90353   145.31683   220.13165   341.21090   115.23941
#>                    SRR11498070 SRR11498071 SRR11498072 SRR11498073 SRR11498074
#> ENSG00000070404.12    44.79135   136.64611    40.09664    17.48312    20.45126
#> ENSG00000070831.20    46.05098    43.95159    53.90587    55.45721    47.51016
#> ENSG00000092199.21   214.11104   227.47375   246.52902   220.25141   260.88590
#> ENSG00000100567.15   106.60685   124.09180   115.31332   105.69139   105.46941
#> ENSG00000101181.20    17.86774    23.33720    21.65134    21.76331    23.31649
#> ENSG00000102755.15   177.58189   447.94442   270.29311   169.66820   496.04255
#>                    SRR11498075 SRR11498076 SRR11498077 SRR11498078 SRR11498079
#> ENSG00000070404.12    18.72177    49.55750    19.41653   125.57957    17.47104
#> ENSG00000070831.20    36.86573    49.52355    47.49247    42.47217    35.72663
#> ENSG00000092199.21   217.86640   225.80300   225.53977   243.13445   238.39553
#> ENSG00000100567.15   102.48072   114.81176   114.87302   136.99088   106.77130
#> ENSG00000101181.20    23.07231    24.28106    24.41587    25.10409    22.52178
#> ENSG00000102755.15   192.84398   202.08497   209.51895   298.36576   134.90583

# Subset like a SummarizedExperiment. Stored gene assays, inferential
# replicates and uncertainty tables follow the transcripts and samples kept
iso[1:10, c("SRR11498039", "SRR11498063")]
#> <IsoformicExperiment> pe_placenta_example
#>   annotation : gencode_v50_example_subset (69 genes, 1,914 transcripts, 18,656 exons)
#>   samples    : 2 [condition]
#>   transcripts: 10 assayed from 1 genes
#>   assays     : abundance, counts
#>   gene assays: derived from transcript assays on demand
#>   DE results : det (710 transcripts), deg (69 genes)

# Build one from matrices
annot <- system.file("extdata", "example_experiment", "annotation", package = "isoformic")
tpm <- tx_assay(iso, "abundance")[1:5, 1:3]
small <- IsoformicExperiment(
  experiment_name = "tiny",
  annot_path = annot,
  assays = list(abundance = tpm),
  col_data = data.frame(sample_id = colnames(tpm), condition = c("a", "a", "b"))
)
small
#> <IsoformicExperiment> tiny
#>   annotation : gencode_v50_example_subset (69 genes, 1,914 transcripts, 18,656 exons)
#>   samples    : 3 [condition]
#>   transcripts: 5 assayed from 1 genes
#>   assays     : abundance
#>   gene assays: derived from transcript assays on demand
#>   DE results : <none>
```
