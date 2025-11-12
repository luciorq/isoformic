# Create Transcript-to-Gene Relationship Table

Extracts a transcript-to-gene mapping table from GENCODE annotation
files, such as the transcriptome FASTA file. Currently, only FASTA files
are supported.

## Usage

``` r
make_tx_to_gene(file_path, file_type = c("fasta", "gff", "gtf"))
```

## Arguments

- file_path:

  A character string specifying the path to the reference file (e.g.,
  GENCODE FASTA file).

- file_type:

  A character string specifying the type of the reference file.
  Currently, only `"fasta"` is supported. Default is `"fasta"`.

## Value

A `tibble` containing the transcript-to-gene mapping information,
including transcript IDs, gene IDs, transcript names, gene names, and
transcript types.

## Details

The function reads the headers of the FASTA file and extracts relevant
information to create a mapping table. For GTF or GFF3 files, support is
not yet implemented.

## Examples

``` r
if (FALSE) { # \dontrun{
# Assuming you have downloaded the GENCODE transcriptome FASTA file:
fasta_file <- download_reference(
  version = "43",
  organism = "human",
  file_type = "fasta",
  output_path = "data-raw"
)

# Create the transcript-to-gene mapping table
tx_to_gene <- make_tx_to_gene(file_path = fasta_file, file_type = "fasta")

# View the first few rows
utils::head(tx_to_gene)
} # }
```
