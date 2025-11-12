# IsoformicExperiment Class

The `IsoformicExperiment` class encapsulates the core data structure for
transcriptomic analyses in the `isoformic` package. It holds the path to
the dataset, sample metadata, and provides access to transcript, gene,
and exon annotations through properties. The preferred way to construct
an object of this class is through the `IsoformicExperiment()`
constructor.

## Usage

``` r
IsoformicExperiment(
  experiment_name = NA_character_,
  data_path = NULL,
  annot_path = NULL,
  assay = NULL,
  col_data = NULL,
  annot_metadata = NULL,
  dea = NULL,
  gsea = NULL,
  tx_type_palette = character(0)
)

col_data(self, ...)

annot_data_transcripts(self, ...)

annot_data_genes(self, ...)

annot_data_exons(self, ...)

annot_data(self, ...)

annot_row_names(self, ...)

col_names(self, ...)

row_names(self, ...)

tx_to_gene(self, ...)

row_data(self, ...)

tx_annot(self, ...)

de_tx(self, ...)

de_gene(self, ...)
```

## Arguments

- experiment_name:

  Character string specifying the name of the experiment. This name is
  used for caching the assays experiment. If a name is not provided a
  random identifier is used.

- data_path:

  Character string specifying the path to the data directory.

- annot_path:

  Character string specifying the path to the annotation database
  directory.

- assay:

  A list of matrices or data frames containing assay data, with
  transcript IDs as row names and sample IDs as column names. Each
  element of the list represents a different assay (e.g., TPM, counts).

- col_data:

  A data frame containing sample metadata. First column must be
  `sample_id` matching the column names of the assays.

- annot_metadata:

  A list containing metadata about the annotation, such as source,
  version, and date.

- dea:

  A list containing differential expression analysis results for
  transcripts and genes.

- gsea:

  A list containing gene set enrichment analysis results.

- tx_type_palette:

  A named character vector specifying the color palette for different
  transcript types.

- self:

  An `IsoformicExperiment` object.

- ...:

  Additional arguments passed to methods.

- annot_data_transcripts:

  A property that retrieves transcript annotation data.

- annot_data_genes:

  A property that retrieves gene annotation data.

- annot_data_exons:

  A property that retrieves exon annotation data.

- annot_data:

  A property that aggregates transcript, gene, and exon annotation data.
