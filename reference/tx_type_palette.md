# Transcript Type Color Palette

Default colour palette for transcript types: the GENCODE/Ensembl
biotypes used by the workflow, the `gene` level, generic labels for
annotations without biotypes (`unknown`, `novel`) and the SQANTI3
structural categories assigned by
[`add_sqanti_classification()`](https://luciorq.github.io/isoformic/reference/add_sqanti_classification.md)
to long-read transcript catalogs.

## Usage

``` r
tx_type_palette()
```

## Value

A named character vector of hexadecimal colours.

## Examples

``` r
head(tx_type_palette())
#>                           gene                 protein_coding 
#>                      "#fb8072"                      "#a6d854" 
#>                retained_intron protein_coding_CDS_not_defined 
#>                      "#8da0cb"                      "#fc8d62" 
#>        nonsense_mediated_decay                         lncRNA 
#>                      "#66c2a5"                      "#e78ac3" 
```
