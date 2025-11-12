# ContextData Class

The `ContextData` class holds information on the genomic context of
transcripts and their annotation in a `IsoformicExperiment` analysis.

## Usage

``` r
ContextData(
  gene_name = character(0),
  gff_file = character(0),
  txdb = NULL,
  annotation_table = NULL,
  annotation_name = character(0),
  assembly_name = character(0),
  ideogram_assembly = character(0),
  organism = character(0),
  orgdb_package = character(0),
  bsgenome_package = character(0),
  tx_type_palette = character(0)
)
```

## Arguments

- gene_name:

  Character vector of gene names to print.

- gff_file:

  Character string specifying the path to a GFF file containing the gene
  annotation.

- txdb:

  A `TxDb` object containing the transcript database.

- annotation_table:

  A data frame containing the transcript annotation.

- annotation_name:

  Character string specifying the name of the annotation (e.g.,
  "Ensembl_v104").

- assembly_name:

  Character string specifying the name of the genome assembly (e.g.,
  "GRCh38").

- ideogram_assembly:

  Character string specifying the name of the genome assembly for
  ideogram plotting (e.g., "hg38").

- organism:

  Character string specifying the organism name (e.g., "Homo sapiens").

- orgdb_package:

  Character string specifying the name of the organism database package
  (e.g., "org.Hs.eg.db").

- bsgenome_package:

  Character string specifying the name of the BSgenome package (e.g.,
  "BSgenome.Hsapiens.UCSC.hg38").

- tx_type_palette:

  Character vector specifying the color palette for transcript types.

## Details

The preferred way to construct an object of this class is through the
[`create_context_data()`](https://luciorq.github.io/isoformic/reference/create_context_data.md)
function.
