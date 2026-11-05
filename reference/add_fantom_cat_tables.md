# Attach FANTOM CAT Tables to an Annotation

FANTOM CAT distributes, next to the GTF, an ID mapping table (`geneID`,
`transcriptID`, `CAGEClusterID`) and a gene information table
(`geneClass`, `geneCategory`, `TIEScore_tier`, `splicingType`,
`directionality`, transposable element associations, ...). This function
joins them into the annotation directory: the CAGE cluster (the
CAGE-defined promoter of each transcript) becomes `cage_cluster_id`,
with its coordinates parsed into `cage_cluster_seqid`,
`cage_cluster_start`, `cage_cluster_end`, and the gene information
columns are added to the gene table in snake case. See
[`download_fantom_cat_tables()`](https://luciorq.github.io/isoformic/reference/download_fantom_cat_tables.md).

## Usage

``` r
add_fantom_cat_tables(x, id_mapping, gene_info = NULL)
```

## Arguments

- x:

  An annotation directory path or an
  [IsoformicExperiment](https://luciorq.github.io/isoformic/reference/IsoformicExperiment.md).

- id_mapping:

  Path (or data frame) of the ID mapping table.

- gene_info:

  Optional path (or data frame) of the gene info table.

## Value

The annotation directory path, invisibly.

## Examples

``` r
if (FALSE) { # \dontrun{
gtf <- download_reference(reference = "fantom", version = "lv3_robust")
tables <- download_fantom_cat_tables("lv3_robust")
annot <- prepare_isoformic_annotation(gtf)
add_fantom_cat_tables(annot, tables[["id_mapping"]], tables[["gene_info"]])
} # }
```
