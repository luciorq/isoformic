# Attach External Productivity Evidence to an Annotation

These functions add evidence columns to the transcript table of the
annotation directory of `x` from the frozen datasets
([frozen_data](https://luciorq.github.io/isoformic/reference/frozen_data.md)),
so that
[`transcript_productivity()`](https://luciorq.github.io/isoformic/reference/productivity.md)
can use observed translation and decay instead of rules alone.
Identifiers are matched without version.

- `add_appris()` adds `appris_annotation` (`PRINCIPAL:1` ...
  `ALTERNATIVE:2`, `MINOR`), `appris_principal`, `appris_score` and the
  method scores (`firestar_residues`, `matador_structure`,
  `corsair_conservation`, `spade_domain_score`, `thump_tm_helices`,
  `trifid_score`). Principal isoforms are productive evidence.

- `add_riboseq_orfs()` adds `riboseq_n_orfs` (Ribo-seq ORFs compatible
  with the transcript's exon structure), `riboseq_orf_types`,
  `riboseq_main_orf` (`TRUE` when a translated ORF coincides with the
  transcript's own ORF, from
  [`transcript_orfs()`](https://luciorq.github.io/isoformic/reference/transcript_orfs.md)
  when sequences are available or else from the annotated CDS) and
  `riboseq_uorf`. A translated main ORF is productive evidence. Phase 1
  ORFs list their compatible transcripts; all ORFs are also matched by
  exon-block compatibility (every block inside an exon, consecutive
  blocks in consecutive exons, internal block boundaries on exon
  boundaries).

- `add_nmd_exons()` adds `nmd_validated_exon`, `TRUE` when the
  transcript contains an exon with the exact coordinates of an
  NMD-inducing cassette exon validated by NMD inhibition (Karousis et
  al. 2021), or is one of the listed NMD isoforms. This is unproductive
  evidence.

Each function returns the annotation directory path invisibly and
invalidates the object caches of `x` (the next call to a method rereads
the annotation).

## Usage

``` r
add_appris(x, table = NULL)

add_riboseq_orfs(x, table = NULL)

add_nmd_exons(x, table = NULL)
```

## Arguments

- x:

  An annotation directory path or an
  [IsoformicExperiment](https://luciorq.github.io/isoformic/reference/IsoformicExperiment.md).

- table:

  A data frame with the dataset, or a path to its Parquet file; `NULL`
  (default) fetches the frozen dataset with
  [`fetch_frozen_data()`](https://luciorq.github.io/isoformic/reference/frozen_data.md).

## Value

The annotation directory path, invisibly.

## Examples

``` r
if (FALSE) { # \dontrun{
iso <- example_isoformic()
add_appris(iso)
add_riboseq_orfs(iso)
add_nmd_exons(iso)
prod <- transcript_productivity(iso, gene_name = "RBPJ")
prod[, c("transcript_name", "productivity", "reason", "source")]
} # }
```
