# Build the Frozen Datasets from Their Sources

Each function downloads the source files of one frozen dataset (see
[`frozen_datasets()`](https://luciorq.github.io/isoformic/reference/frozen_data.md))
into the sources cache, converts them to the Parquet table used by the
package and returns the file path. They are what
[`fetch_frozen_data()`](https://luciorq.github.io/isoformic/reference/frozen_data.md)
runs when a frozen file cannot be downloaded, and what
`data-raw/build_frozen_data.R` runs to publish a new version. The
resulting tables:

- `build_frozen_appris()`: one row per transcript with
  `appris_annotation` (`PRINCIPAL:1` to `ALTERNATIVE:2`, `MINOR`),
  `appris_score`, the method scores (`firestar_residues`,
  `matador_structure`, `corsair_conservation`, `spade_domain_score`,
  `thump_tm_helices`, `crash_signal`, `trifid_score`), `protein_length`,
  `ccds_id` and `mane` label.

- `build_frozen_appris_domains()`: one row per Pfam domain per
  transcript with its genomic span (`seqid`, `start_pos`, `end_pos`,
  `strand`), `pfam_acc`, `pfam_name`, `evalue`, `pep_start`, `pep_end`.

- `build_frozen_riboseq_orfs()`: one row per Ribo-seq ORF with `phase`,
  genomic span, `block_starts`/`block_ends` (1-based, `;`-separated),
  `orf_type`, `gene_id`, `gene_name`, `transcript_ids` (compatible
  transcripts, phase 1 only), `replicated` and `studies`.

- `build_frozen_nmd_exons()`: one row per validated NMD-inducing exon
  with `seqid`, `start_pos`, `end_pos`, `strand`, the flanking exons,
  `nmd_transcript_ids`, `non_nmd_transcript_ids`, `gene_id`,
  `gene_name`. Reading the supplementary spreadsheet needs the `readxl`
  package.

- `build_frozen_polyasite()`: one row per poly(A) site cluster with
  `seqid`, `start_pos`, `end_pos`, `cleavage_pos`, `strand`, `tpm`,
  `frac_samples`, `n_protocols`, `site_class` (TE terminal exon, EX
  exon, IN intron, DS downstream, AE/AI/AU antisense, IG intergenic) and
  `signals`.

- `build_frozen_gtex()`: one row per transcript and tissue with
  `median_tpm` and `n_samples` (GTEx v8, RSEM transcript TPM;
  identifiers are GENCODE v26, stored without version).

## Usage

``` r
build_frozen_appris(output_path = ":cache:", sources_dir = NULL)

build_frozen_appris_domains(output_path = ":cache:", sources_dir = NULL)

build_frozen_riboseq_orfs(output_path = ":cache:", sources_dir = NULL)

build_frozen_nmd_exons(output_path = ":cache:", sources_dir = NULL)

build_frozen_polyasite(output_path = ":cache:", sources_dir = NULL)

build_frozen_gtex(output_path = ":cache:", sources_dir = NULL)
```

## Arguments

- output_path:

  Directory receiving the Parquet file (default: the frozen-data cache).

- sources_dir:

  Directory holding (or receiving) the source files; default:
  `get_isoformic_cache("sources", <name>)`.

## Value

The path of the written Parquet file, invisibly.

## Examples

``` r
if (FALSE) { # \dontrun{
build_frozen_riboseq_orfs()
} # }
```
