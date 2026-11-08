# Other annotations: T2T-CHM13, Ensembl, RefSeq, haplotypes, FANTOM CAT and long reads

The
[introduction](https://luciorq.github.io/isoformic/articles/isoformic-introduction.md)
uses GENCODE, but
[`prepare_isoformic_annotation()`](https://luciorq.github.io/isoformic/reference/prepare_isoformic_annotation.md)
does not depend on the GENCODE vocabulary. Genes are features of a
gene-like type or parents of transcripts, transcripts are features of a
transcript-like type (`transcript`, `mRNA`, `lnc_RNA`, …) or parents of
exons, and identifiers are taken from `gene_id`/`transcript_id`
attributes or from `ID`/`Parent`. This article shows how to obtain and
prepare the other references the downloader knows about. Most chunks are
not evaluated because they download gigabytes of data; every file lands
in the package cache (`?isoformic-options`).

``` r

library(isoformic)
```

## Telomere-to-telomere assembly (T2T-CHM13)

The same object works on the T2T-CHM13v2.0 assembly, with the NCBI
RefSeq annotation or the Liftoff annotations of the T2T consortium:

``` r

# NCBI RefSeq annotation of T2T-CHM13v2.0: GFF3 + matching transcript FASTA
t2t_gff <- download_reference(reference = "t2t", file_type = "gff")
t2t_fasta <- download_reference(reference = "t2t", file_type = "fasta")
annot_t2t <- prepare_isoformic_annotation(t2t_gff)

# T2T consortium Liftoff annotations keep the GRCh38 identifiers of each
# lifted gene and transcript in `source_gene_id` / `source_transcript_id`
annot_lift <- prepare_isoformic_annotation(
  "chm13v2.0_GENCODEv35_CAT_Liftoff.vep.gff3.gz"
)
read_parquet(fs::path(annot_lift, "annot_data_transcripts.parquet"), n_max = 3)
```

For lifted annotations, build the transcript FASTA for Salmon from the
assembly with
`gffread -w transcripts.fa -g chm13v2.0.fa annotation.gff3`, and use
`source_transcript_id` to join results obtained on GRCh38 (for example
gene sets defined on GENCODE identifiers).
`extract_transcript_sequences(annot, genome_fasta = )` derives the
transcript sequences from the assembly for the sequence-aware analyses.

## Ensembl, NCBI and MANE

The same downloader covers Ensembl (any division), any NCBI assembly
accession and MANE, and arbitrary URLs:

``` r

ens_gff <- download_reference(reference = "ensembl")                       # current human release
ath_gff <- download_reference(reference = "ensembl", division = "plants",
                              species = "arabidopsis_thaliana")
grch38_refseq <- download_reference(reference = "ncbi", accession = "GCF_000001405.40")
mane_gff <- download_reference(reference = "mane")                         # canonical transcripts
```

Transcripts carry their `tags` (`MANE_Select`, `Ensembl_canonical`,
`basic`, …) and an `is_canonical` flag whatever the source, so
[`dominant_isoform()`](https://luciorq.github.io/isoformic/reference/diversity.md)
and
[`transcript_features()`](https://luciorq.github.io/isoformic/reference/transcript_features.md)
know the canonical isoform of each gene.

## Haplotype-resolved assemblies (trio binning, HPRC)

Phased assemblies come with one annotation per haplotype.
[`merge_annotations()`](https://luciorq.github.io/isoformic/reference/merge_annotations.md)
combines them into one annotation directory whose identifiers carry a
`_<haplotype>` suffix (rename the transcript FASTA headers the same way
before indexing), and
[`haplotype_summary()`](https://luciorq.github.io/isoformic/reference/diversity.md)
reports the expression of each gene per haplotype, the haplotype
fraction and whether the dominant isoform differs between haplotypes:

``` r

files <- download_hprc_annotation("HG00438")           # HG00438.1.chm13.gff3.gz, HG00438.2.chm13.gff3.gz
merged <- merge_annotations(files, annot_name = "HG00438_chm13")
iso_hap <- read_salmon(quant_dirs, annot_path = merged)
haplotype_summary(iso_hap, group = "condition")
```

## FANTOM CAT and FANTOM6

The FANTOM CAT transcriptome assembly (the catalog used by FANTOM6) is
read with its own variables: `geneClass` becomes the gene type and the
type of transcripts without a biotype, and `geneCategory`, `DHS_type`,
`coding_status`, `cumulative_support` and `TIEScore` are kept; the
companion tables add the CAGE-defined promoter of every transcript
(`cage_cluster_id`) and the gene information columns:

``` r

gtf <- download_reference(reference = "fantom", version = "lv3_robust", file_type = "gtf")
tables <- download_fantom_cat_tables("lv3_robust")
annot_fantom <- prepare_isoformic_annotation(gtf)
add_fantom_cat_tables(annot_fantom, tables[["id_mapping"]], tables[["gene_info"]])
```

## Long-read transcript catalogs

StringTie, Bambu, FLAIR or IsoQuant GTF files are read like any
annotation; novel transcripts have no biotype until SQANTI3 structural
categories are added, after which they behave as transcript types in
every plot and summary:

``` r

annot_lr <- prepare_isoformic_annotation("flair_isoforms.gtf")
add_sqanti_classification(annot_lr, "flair_isoforms_classification.txt")
```

Novel transcripts have no annotated CDS, so the [sequence-aware
ORFs](https://luciorq.github.io/isoformic/articles/isoformic-productivity.md)
derived from the genome are what classify them as productive or not.

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
#> [1] isoformic_0.1.2.9022
#> 
#> loaded via a namespace (and not attached):
#>  [1] compiler_4.6.2  fastmap_1.2.0   cli_3.6.6       tools_4.6.2    
#>  [5] htmltools_0.5.9 otel_0.3.0      yaml_2.3.12     rmarkdown_2.32 
#>  [9] knitr_1.52      S7_0.2.2        jsonlite_2.0.0  xfun_0.62      
#> [13] digest_0.6.39   rlang_1.3.0     evaluate_1.0.5
```
