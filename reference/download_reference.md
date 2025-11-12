# Download Reference Files from GENCODE

Downloads reference annotation files from the GENCODE database for human
or mouse genomes. Supports downloading GFF, GTF, transcriptome FASTA,
and genome FASTA files. The function handles directory creation and
checks for existing files to avoid redundant downloads.

## Usage

``` r
download_reference(
  version = "49",
  reference = "gencode",
  organism = c("human", "mouse"),
  file_type = c("gff", "gtf", "fasta", "genome_fasta"),
  output_path = ":cache:",
  timeout_limit = 3600,
  method = "auto"
)
```

## Arguments

- version:

  A character string specifying the GENCODE release version. For mouse
  references, include the letter 'M' in the version string (e.g.,
  `"M38"`). Default is `"49"`.

- reference:

  A character string specifying the source of the reference file.
  Currently, only `"gencode"` is supported. Default is `"gencode"`.

- organism:

  A character string specifying the organism. Valid options are
  `"human"` or `"mouse"`.

- file_type:

  A character string specifying the type of file to download. Valid
  options are `"gff"`, `"gtf"`, `"fasta"` or `"genome_fasta"`. Defaults
  to `"gff"`. **Note:** `"fasta"` refers to the transcriptome FASTA
  file. `"genome_fasta"` refers to the whole genome sequence FASTA file.

- output_path:

  A character string specifying the directory where the downloaded file
  will be saved. Defaults to `":cache:"`. Cache path is defined by the
  `[isoformic::get_isoformic_cache()]` function.

- timeout_limit:

  A numeric value specifying the maximum time in seconds for the
  download to complete. This argument takes precedence over
  `options("timeout")`. Defaults to `3600` seconds (1 hour).

- method:

  A character string specifying the method used by
  [`utils::download.file()`](https://rdrr.io/r/utils/download.file.html).
  Defaults to `"auto"`.

## Value

A character string with the full path to the downloaded file.

## Details

The function constructs the appropriate download URL based on the
specified organism, version, and file type, and downloads the file to
the specified output path, being a user cache by default. If the file
already exists in the output directory, the function will not download
it again and will return the existing file path. The function requires
an internet connection and handles timeout settings to prevent download
interruptions.

## Note

Currently, only `"gencode"` reference files are supported. The `"mane"`
reference is not implemented yet.

## Examples

``` r
if (FALSE) { # \dontrun{
# Download human GFF file for GENCODE release 49
gff_file <- download_reference(
  version = "49",
  organism = "human",
  file_type = "gff",
  output_path = ":cache:"
)

# Download mouse GFF file for GENCODE release M38
gff_file_mouse <- download_reference(
  version = "M38",
  organism = "mouse",
  file_type = "gff",
  output_path = ":cache:"
)

# Download human transcriptome FASTA file for GENCODE release 49
fasta_file <- download_reference(
  version = "49",
  organism = "human",
  file_type = "fasta",
  output_path = ":cache:"
)
} # }
```
