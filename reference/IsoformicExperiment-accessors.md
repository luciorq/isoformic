# Accessors for IsoformicExperiment Objects

S7 generics giving `SummarizedExperiment`-style access to the three
levels of an
[IsoformicExperiment](https://luciorq.github.io/isoformic/reference/IsoformicExperiment.md).
All of them return in-memory tibbles or matrices; the annotation tables
are read from Parquet on first use and cached (see the class
documentation).

**Samples**

- `col_data()`, `col_data<-`: sample metadata (`sample_id` first).

- `col_names()`: sample IDs of the assays.

**Transcript level**

- `assays()`, `assay(x, name)`, `assay_names()`, `assay<-`: transcript
  matrices. `assay()` without a name returns the first assay.

- `row_names()`: transcript IDs of the assays.

- `row_data()`: transcript annotation (with gene columns) for the
  assayed transcripts, in assay row order.

- `tx_annot()`: the same join for *all* annotated transcripts.

- `tx_to_gene()`: two-column transcript/gene mapping.

**Gene level**

- `gene_assays()`, `gene_assay(x, name)`, `gene_assays<-`: gene
  matrices, derived by summing isoforms unless set explicitly.

- `gene_data()`: gene annotation for the genes represented in the
  assays.

- `gene_names()`: gene IDs of the gene assays.

**Annotation**

- `annot_data_transcripts()`, `annot_data_genes()`,
  `annot_data_exons()`: the raw annotation tables.

- `annot_data()`: transcript-centric join of transcripts and genes
  (exons optional).

- `annot_metadata()`: named list describing the annotation source.

- `exon_table(x, gene_name)`: exons of the transcripts of selected
  genes.

**Housekeeping**

- `clear_cache()`: drop cached derived data for the object (and, with
  `annotation = TRUE`, the shared annotation tables).

## Usage

``` r
col_data(self, ...)

col_data(self) <- value

col_names(self, ...)

assays(self, ...)

assay_names(self, ...)

assay(self, name = NULL, ...)

assay(self, name) <- value

row_names(self, ...)

row_data(self, ...)

tx_annot(self, ...)

tx_to_gene(self, ...)

gene_assays(self, ...)

gene_assays(self) <- value

gene_assay(self, name = NULL, ...)

gene_names(self, ...)

gene_data(self, ...)

annot_data_transcripts(self, ...)

annot_data_genes(self, ...)

annot_data_exons(self, ...)

annot_data(self, with_exons = FALSE, ...)

annot_metadata(self, ...)

exon_table(self, gene_name, ...)

clear_cache(self, annotation = FALSE, ...)
```

## Arguments

- self:

  An
  [IsoformicExperiment](https://luciorq.github.io/isoformic/reference/IsoformicExperiment.md).

- ...:

  Passed on to methods.

- value:

  Replacement value.

- name:

  Name of the assay. `NULL` selects the first assay.

- with_exons:

  Logical, also join the exon table (one row per exon).

- gene_name:

  Character vector of gene names or gene IDs.

- annotation:

  Logical, also clear the shared annotation cache.

## Value

See the individual descriptions above.

## Examples

``` r
iso <- example_isoformic()
col_data(iso)
#> # A tibble: 25 × 2
#>    sample_id   condition   
#>    <chr>       <chr>       
#>  1 SRR11498039 preeclampsia
#>  2 SRR11498040 preeclampsia
#>  3 SRR11498041 preeclampsia
#>  4 SRR11498042 preeclampsia
#>  5 SRR11498043 preeclampsia
#>  6 SRR11498044 preeclampsia
#>  7 SRR11498045 preeclampsia
#>  8 SRR11498046 preeclampsia
#>  9 SRR11498063 control     
#> 10 SRR11498064 control     
#> # ℹ 15 more rows
assay_names(iso)
#> [1] "abundance" "counts"   
assay(iso, "abundance")[1:3, 1:3]
#>                   SRR11498039 SRR11498040 SRR11498041
#> ENST00000374514.8     895.871     850.783     956.829
#> ENST00000420982.5       0.000       0.000       0.000
#> ENST00000374502.7      20.998      16.324      17.151
gene_assay(iso, "abundance")[1:3, 1:3]
#>                    SRR11498039 SRR11498040 SRR11498041
#> ENSG00000005007.12    3253.512    2563.848    2656.175
#> ENSG00000005206.17    1727.018    1515.008    1807.224
#> ENSG00000007384.15    1970.727    1612.095    2188.516
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
head(gene_data(iso))
#> # A tibble: 6 × 14
#>   gene_id   gene_name gene_type source_gene_id gene_class gene_category dhs_type
#>   <chr>     <chr>     <chr>     <chr>          <chr>      <chr>         <chr>   
#> 1 ENSG0000… UPF1      protein_… NA             NA         NA            NA      
#> 2 ENSG0000… SPPL2B    protein_… NA             NA         NA            NA      
#> 3 ENSG0000… RHBDF1    protein_… NA             NA         NA            NA      
#> 4 ENSG0000… PIGQ      protein_… NA             NA         NA            NA      
#> 5 ENSG0000… MED24     protein_… NA             NA         NA            NA      
#> 6 ENSG0000… REV3L     protein_… NA             NA         NA            NA      
#> # ℹ 7 more variables: coding_status <chr>, cumulative_support <chr>,
#> #   haplotype <chr>, seqid <chr>, start_pos <int>, end_pos <int>, strand <chr>
head(exon_table(iso, "RBPJ"))
#> # A tibble: 6 × 12
#>   gene_id      gene_name transcript_id transcript_name transcript_type haplotype
#>   <chr>        <chr>     <chr>         <chr>           <chr>           <chr>    
#> 1 ENSG0000016… RBPJ      ENST00000342… RBPJ-201        protein_coding  NA       
#> 2 ENSG0000016… RBPJ      ENST00000342… RBPJ-201        protein_coding  NA       
#> 3 ENSG0000016… RBPJ      ENST00000342… RBPJ-201        protein_coding  NA       
#> 4 ENSG0000016… RBPJ      ENST00000342… RBPJ-201        protein_coding  NA       
#> 5 ENSG0000016… RBPJ      ENST00000342… RBPJ-201        protein_coding  NA       
#> 6 ENSG0000016… RBPJ      ENST00000342… RBPJ-201        protein_coding  NA       
#> # ℹ 6 more variables: exon_id <chr>, exon_number <int>, seqid <chr>,
#> #   start_pos <int>, end_pos <int>, strand <chr>
```
