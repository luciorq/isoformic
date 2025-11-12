# Create ContextData Object

This function instantiates a `ContextData` object containing information
about the genomic context of transcripts. It requires a GFF file for
gene annotation and constructs a `TxDb` object from it. The function
also prepares an annotation table and updates the transcript names in
the `TxDb` object. The `ContextData` object can then be used in
conjunction with `IsoformicExperiment` for transcriptomic analyses.

## Usage

``` r
create_context_data(
  gff_file,
  ...,
  organism,
  orgdb_package,
  bsgenome_package,
  tx_type_palette = NULL
)
```

## Arguments

- gff_file:

  Character string specifying the path to a GFF file containing the gene
  annotation.

- ...:

  These dots are for future extensions and must be empty.

- organism:

  Character string specifying the organism name (e.g., "Homo sapiens").

- orgdb_package:

  Character string specifying the name of the organism database package
  (e.g., "org.Hs.eg.db").

- bsgenome_package:

  Character string specifying the name of the BSgenome package (e.g.,
  "BSgenome.Hsapiens.UCSC.hg38").

- tx_type_palette:

  Named character vector specifying the color palette for transcript
  types.
