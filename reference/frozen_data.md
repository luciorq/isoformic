# Frozen Reference Datasets

`isoformic` ships evidence tables derived from public resources, frozen
at a given version so that results are reproducible:

- `appris`: APPRIS principal/alternative isoform labels and method
  scores (functional residues, structure, conservation, domains) per
  transcript, from the APPRIS release files.

- `appris_domains`: Pfam domain genomic coordinates per transcript
  (APPRIS SPADE).

- `riboseq_orfs`: GENCODE Ribo-seq ORF catalog, phase 1 and phase 2,
  with exon blocks and compatible transcripts.

- `nmd_exons`: NMD-inducing cassette exons validated by NMD inhibition
  and nanopore sequencing (Karousis et al. 2021, Genome Biology,
  Supplementary Table 1).

- `polyasite`: PolyASite 2.0 poly(A) site atlas (GRCh38) with site class
  and poly(A) signals.

- `gtex_transcript_tpm`: GTEx v8 median transcript TPM per tissue
  (RSEM).

`frozen_datasets()` lists them with version, source URLs and the MD5 of
the frozen file. `fetch_frozen_data()` returns the local path of a
dataset, downloading the frozen file from the first mirror that has it
(`frozen_mirrors()`: any URL in `options(isoformic.frozen_mirrors = )`,
the package's GitHub release, then the Zenodo snapshot of the processed
files) or, when that is not possible (`build = "auto"`) or asked for
(`build = "always"`), rebuilding it from the source files with the
matching `build_frozen_*()` function. `read_frozen_data()` fetches and
reads it. `bundle_frozen_data()` assembles the snapshot that is
deposited on Zenodo: every frozen file, the coding potential model, an
`MD5SUMS` file and a README stating versions and sources.

## Usage

``` r
frozen_datasets()

frozen_mirrors(
  release_tag = "frozen-data",
  zenodo_record = frozen_zenodo_record()
)

fetch_frozen_data(
  name,
  output_path = ":cache:",
  build = c("auto", "never", "always"),
  verify = TRUE,
  release_tag = "frozen-data",
  mirrors = frozen_mirrors(release_tag),
  ...
)

read_frozen_data(
  name,
  output_path = ":cache:",
  build = c("auto", "never", "always"),
  ...
)

bundle_frozen_data(
  output_dir,
  output_path = ":cache:",
  build = c("auto", "never", "always"),
  include_coding_model = TRUE,
  verify = TRUE
)
```

## Arguments

- release_tag:

  Tag of the GitHub release hosting the frozen files.

- zenodo_record:

  Zenodo record id (or `https://zenodo.org/records/<id>` URL) of a
  deposit holding the frozen files; defaults to the record published
  with the package, `NA` when there is none.

- name:

  Dataset name, one of `names(frozen_datasets())`.

- output_path:

  Directory holding the frozen files; `":cache:"` (default) uses
  `get_isoformic_cache("frozen")`.

- build:

  `"auto"` (download, else build from sources), `"never"` or `"always"`.

- verify:

  Logical, check the MD5 of the file against the registry (default
  `TRUE`; a mismatch is a warning, since a rebuilt file from a newer
  source is still usable).

- mirrors:

  Character vector of base URLs tried in order for the frozen file;
  default `frozen_mirrors()`.

- ...:

  Reserved for future use; must be empty.

- output_dir:

  Directory receiving the bundle (created if needed).

- include_coding_model:

  Logical, also copy the coding potential model files (default `TRUE`).

## Value

`frozen_datasets()` returns a tibble; `fetch_frozen_data()` the file
path, invisibly; `read_frozen_data()` a tibble; `frozen_mirrors()` a
character vector; `bundle_frozen_data()` the bundle directory,
invisibly.

## Examples

``` r
frozen_datasets()[, c("name", "version", "file")]
#> # A tibble: 6 × 3
#>   name                version       file                                        
#>   <chr>               <chr>         <chr>                                       
#> 1 appris              e113v49       appris_e113v49.parquet                      
#> 2 appris_domains      e113v49       appris_domains_e113v49.parquet              
#> 3 riboseq_orfs        phase1_phase2 gencode_riboseq_orfs_phase1_phase2.parquet  
#> 4 nmd_exons           karousis2021  nmd_exons_karousis2021.parquet              
#> 5 polyasite           2.0_GRCh38.96 polyasite_2.0_GRCh38.96.parquet             
#> 6 gtex_transcript_tpm v8            gtex_v8_transcript_tpm_median_by_tissue.par…
if (FALSE) { # \dontrun{
appris <- read_frozen_data("appris")
head(appris)
} # }
```
