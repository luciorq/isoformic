# Download Reference Annotations and Sequences

Downloads annotation (GFF3/GTF) and sequence (transcript or genome
FASTA) files from the main public providers:

- `"gencode"` (default): any GENCODE release for human (GRCh38) or mouse
  (GRCm39), including future releases.

- `"ensembl"`: any Ensembl release, for any species of any division
  (`vertebrates`, `plants`, `metazoa`, `fungi`, `protists`).

- `"ncbi"`: any NCBI assembly by accession (for example
  `GCF_009914755.1` for T2T-CHM13v2.0, `GCF_000001405.40` for
  GRCh38.p14).

- `"t2t"`: shorthand for
  `reference = "ncbi", accession = "GCF_009914755.1"`
  (telomere-to-telomere T2T-CHM13v2.0).

- `"mane"`: MANE (Matched Annotation from NCBI and EMBL-EBI) releases,
  the canonical transcript per protein-coding gene, in Ensembl or RefSeq
  coordinates and identifiers.

- `"fantom"`: the FANTOM CAT transcriptome assembly (Hon et al. 2017;
  GRCh37/hg19) at one of its levels (`version = "lv3_robust"` by
  default; `"lv1_raw"`, `"lv2_permissive"`, `"lv4_stringent"`), GTF
  only. FANTOM6 analyses use the lv3 robust catalog. The companion ID
  mapping, gene information and CAGE cluster tables come from
  [`download_fantom_cat_tables()`](https://luciorq.github.io/isoformic/reference/download_fantom_cat_tables.md).

Instead of assuming file names, the function lists the release directory
on the server and picks the matching file, so new releases and new
genome patch versions (for example `GRCh38.p15`) work without a package
update. Files already present in `output_path` are not downloaded again.
Any other file reachable by URL (for example the T2T consortium Liftoff
GFF3 files) can be cached the same way with
[`download_reference_url()`](https://luciorq.github.io/isoformic/reference/download_reference_url.md),
and HPRC per-haplotype annotations with
[`download_hprc_annotation()`](https://luciorq.github.io/isoformic/reference/download_hprc_annotation.md).

## Usage

``` r
download_reference(
  version = NULL,
  reference = c("gencode", "ensembl", "ncbi", "t2t", "mane", "fantom"),
  organism = c("human", "mouse"),
  file_type = c("gff", "gtf", "fasta", "ncrna_fasta", "genome_fasta", "polya"),
  region = c("all", "chr", "primary"),
  species = NULL,
  division = c("vertebrates", "plants", "metazoa", "fungi", "protists"),
  accession = NULL,
  mane_source = c("ensembl", "refseq"),
  output_path = ":cache:",
  timeout_limit = 3600,
  method = "auto"
)
```

## Arguments

- version:

  Release version. GENCODE: the release number (`"49"`, `49`, or `"M38"`
  / `"38"` for mouse; default `"49"`). Ensembl: the release number
  (`"116"`); `NULL` (default) uses the current release. MANE: `"1.5"`;
  `NULL` uses the current release. NCBI/T2T: ignored (the accession
  identifies the assembly version).

- reference:

  One of `"gencode"`, `"ensembl"`, `"ncbi"`, `"t2t"`, `"mane"`,
  `"fantom"`.

- organism:

  `"human"` (default) or `"mouse"`. For Ensembl it sets the default
  `species` (`homo_sapiens` / `mus_musculus`).

- file_type:

  `"gff"` (default), `"gtf"`, `"fasta"` (transcript sequences; for
  Ensembl the `cdna.all` file), `"ncrna_fasta"` (Ensembl non-coding
  transcript sequences), `"genome_fasta"` (genome sequence; not provided
  by MANE) or `"polya"` (GENCODE manually annotated poly(A) sites and
  signals, `polyAs.gff3.gz`; see
  [`polya_sites()`](https://luciorq.github.io/isoformic/reference/polyadenylation.md)).

- region:

  Annotation flavour: `"all"` (default; reference chromosomes, patches,
  haplotypes and scaffolds), `"chr"` (reference chromosomes only) or
  `"primary"` (primary assembly). For `file_type = "genome_fasta"`,
  `"primary"` selects the primary assembly genome and the other values
  the full patched / toplevel genome. Ignored for NCBI and MANE.

- species:

  Ensembl species name in lower snake case (`"homo_sapiens"`,
  `"arabidopsis_thaliana"`). Defaults from `organism`.

- division:

  Ensembl division: `"vertebrates"` (default, the main Ensembl site),
  `"plants"`, `"metazoa"`, `"fungi"` or `"protists"`.

- accession:

  NCBI assembly accession (`GCF_...` or `GCA_...`), required for
  `reference = "ncbi"`.

- mane_source:

  Identifier space of MANE files: `"ensembl"` (default) or `"refseq"`.

- output_path:

  Directory where the file is saved. Defaults to `":cache:"`, the
  package cache (see
  [`get_isoformic_cache()`](https://luciorq.github.io/isoformic/reference/get_isoformic_cache.md)).

- timeout_limit:

  Maximum time in seconds for the download. Takes precedence over
  `options("timeout")`. Defaults to `3600`.

- method:

  Method used by
  [`utils::download.file()`](https://rdrr.io/r/utils/download.file.html).
  Defaults to `"auto"`.

## Value

The full path to the downloaded (or already present) file.

## Errors

The function fails early with a specific condition class:

- `isoformic_download_reference_no_internet_error`: no connection.

- `isoformic_download_reference_release_not_found_error`: the release
  (or species, assembly) directory does not exist; for GENCODE the
  message lists the latest available releases.

- `isoformic_download_reference_file_not_found_error`: the release
  exists but has no file of the requested type; the message lists the
  candidate files found.

- `isoformic_download_reference_error`: the transfer failed or produced
  an empty file (partial files are removed).

## Matching annotation and sequences

Salmon indices must be built from the transcript FASTA of the same
provider and release as the annotation directory, so that transcript IDs
match: GENCODE FASTA headers are `|`-separated (handled by
`ignore_after_bar`), Ensembl `cdna.all` headers carry the versioned
transcript ID first, NCBI `rna.fna` headers carry the RefSeq accession,
and MANE provides both identifier spaces. For lifted annotations (T2T
consortium Liftoff, HPRC CAT) extract the FASTA from the assembly with
`gffread -w transcripts.fa -g assembly.fa annotation.gff3`.

## Examples

``` r
if (FALSE) { # \dontrun{
# GENCODE
gff_file <- download_reference(version = "49", file_type = "gff")
gtf_file <- download_reference(version = "50", file_type = "gtf", region = "chr")
gff_file_mouse <- download_reference(version = "M38", organism = "mouse")
genome_file <- download_reference(version = "49", file_type = "genome_fasta")

# Ensembl: current human release, Arabidopsis from Ensembl Plants
ens_gff <- download_reference(reference = "ensembl")
ens_cdna <- download_reference(reference = "ensembl", file_type = "fasta")
ath_gff <- download_reference(
  reference = "ensembl", division = "plants", species = "arabidopsis_thaliana"
)

# NCBI assemblies by accession, and the T2T-CHM13v2.0 shorthand
t2t_gff <- download_reference(reference = "t2t", file_type = "gff")
t2t_fasta <- download_reference(reference = "t2t", file_type = "fasta")
grch38_gff <- download_reference(reference = "ncbi", accession = "GCF_000001405.40")

# MANE canonical transcripts (Ensembl identifiers)
mane_gff <- download_reference(reference = "mane", file_type = "gff")
c(gff_file, gtf_file, gff_file_mouse, genome_file, ens_gff, ens_cdna,
  ath_gff, t2t_gff, t2t_fasta, grch38_gff, mane_gff)
} # }
```
