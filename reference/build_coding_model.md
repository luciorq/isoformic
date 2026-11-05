# Train the Coding Potential Model

Builds the hexamer log-odds table and the logistic regression used by
[`coding_potential()`](https://luciorq.github.io/isoformic/reference/transcript_orfs.md)
from a set of coding ORF sequences and non-coding transcript sequences,
and writes both as Parquet files. The shipped model was trained on
GENCODE human protein-coding CDS versus lncRNA transcripts by
`data-raw/build_coding_model.R`; rerun that script to refresh it for a
new GENCODE release or another organism.

## Usage

``` r
build_coding_model(
  coding_cds,
  coding,
  noncoding,
  output_dir,
  min_orf_length = 90L,
  pseudocount = 1
)
```

## Arguments

- coding_cds:

  Character vector of CDS sequences (start to stop codon), used to build
  the hexamer table.

- coding:

  Character vector of full-length coding transcript sequences (mRNAs),
  the positive class of the logistic model.

- noncoding:

  Character vector of non-coding transcript sequences, the negative
  class (and the hexamer background).

- output_dir:

  Directory receiving `hexamer_logodds.parquet` and
  `coding_model_coefficients.parquet`.

- min_orf_length:

  Minimum ORF length used when scoring the training sequences.

- pseudocount:

  Pseudocount added to hexamer counts.

## Value

A list with the hexamer table, the fitted `glm` and the training
accuracy, invisibly.

## Examples

``` r
if (FALSE) { # \dontrun{
build_coding_model(cds_sequences, mrna_sequences, lncrna_sequences, output_dir = "coding_model")
} # }
```
