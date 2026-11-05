# Protein Domain Retention of Isoforms

Compares the Pfam domains of every isoform with those of the canonical
isoform of its gene, using the APPRIS SPADE domain coordinates (frozen
dataset `appris_domains`,
[frozen_data](https://luciorq.github.io/isoformic/reference/frozen_data.md)).
For isoforms present in APPRIS the domain sets are compared directly.
For other isoforms (novel transcripts) the genomic span of each
canonical domain is projected onto the isoform: a domain is `retained`
when the isoform covers the same exonic bases over the domain span and
its ORF (see
[`transcript_orfs()`](https://luciorq.github.io/isoformic/reference/transcript_orfs.md),
or the annotated CDS) extends over it without a frameshift, `partial`
when only part of the span is covered, and `lost` otherwise.

An `altered_productive` isoform that keeps all domains of the canonical
protein is a milder event than one that loses a domain; the
`domain_retention_fraction` column grades this and is included in
[`transcript_productivity()`](https://luciorq.github.io/isoformic/reference/productivity.md)
output when the dataset is available.

## Usage

``` r
domain_retention(self, gene_name = NULL, domains = NULL, ...)
```

## Arguments

- self:

  An
  [IsoformicExperiment](https://luciorq.github.io/isoformic/reference/IsoformicExperiment.md).

- gene_name:

  Optional genes (names or IDs) to restrict the result to.

- domains:

  Data frame or Parquet path of the domain table; `NULL` fetches the
  frozen `appris_domains` dataset.

- ...:

  Passed on to methods.

## Value

A tibble with one row per transcript: `n_domains_canonical`,
`n_retained`, `n_partial`, `n_lost`, `domains_lost` (Pfam names) and
`domain_retention_fraction` (`NA` when the canonical isoform has no
domain).

## Examples

``` r
if (FALSE) { # \dontrun{
iso <- example_isoformic()
domain_retention(iso, gene_name = "RBPJ")
} # }
```
