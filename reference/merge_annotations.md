# Merge Per-Haplotype Annotations Into One Annotation Directory

Phased (trio-binned or pangenome) assemblies come with one annotation
per haplotype, for example the HPRC CAT files `HG00438.1.chm13.gff3.gz`
and `HG00438.2.chm13.gff3.gz`. `merge_annotations()` converts each file
(or reuses an existing annotation directory), tags every table with the
haplotype label, makes identifiers unique by appending `_<haplotype>`,
and writes one combined annotation directory. Transcript, gene and exon
IDs therefore become `<id>_<haplotype>`, the convention expected for the
row names of a haplotype-resolved
[IsoformicExperiment](https://luciorq.github.io/isoformic/reference/IsoformicExperiment.md)
(rename the sequences of the per-haplotype transcript FASTA files the
same way before building the Salmon index).

The `source_gene_id` / `source_transcript_id` columns (GRCh38
identifiers kept by Liftoff/CAT) are left untouched, which is what
[`haplotype_summary()`](https://luciorq.github.io/isoformic/reference/diversity.md)
uses to compare haplotypes of the same gene.

## Usage

``` r
merge_annotations(
  annot_paths,
  output_path = NULL,
  annot_name = "merged_haplotypes",
  id_suffix = TRUE,
  overwrite = FALSE
)
```

## Arguments

- annot_paths:

  Named character vector: names are haplotype labels, values are
  GFF3/GTF files or annotation directories.

- output_path:

  Directory under which the merged annotation directory is created.
  `NULL` uses
  [`get_isoformic_cache()`](https://luciorq.github.io/isoformic/reference/get_isoformic_cache.md).

- annot_name:

  Name of the merged annotation directory.

- id_suffix:

  Logical. Append `_<haplotype>` to gene, transcript and exon
  identifiers (default `TRUE`). Set to `FALSE` when the identifiers are
  already unique across haplotypes.

- overwrite:

  Logical. Replace an existing merged directory.

## Value

The path to the merged annotation directory.

## Examples

``` r
gff <- system.file("extdata", "example.gff3", package = "isoformic")
merged <- merge_annotations(
  c(hap1 = gff, hap2 = gff),
  output_path = tempdir(),
  annot_name = "example_diploid",
  overwrite = TRUE
)
#> ℹ Using cached annotation directory /tmp/RtmpV2WK1u/example.
#> ℹ Using cached annotation directory /tmp/RtmpV2WK1u/example.
tx <- read_parquet(file.path(merged, "annot_data_transcripts.parquet"))
table(tx$haplotype)
#> 
#> hap1 hap2 
#>   72   72 
head(tx$transcript_id)
#> [1] "ENST00000282397.9_hap1" "ENST00000615611.4_hap1" "ENST00000617835.4_hap1"
#> [4] "ENST00000543394.2_hap1" "ENST00000540678.2_hap1" "ENST00000541932.5_hap1"
```
