# isoformic

`isoformic` is a workflow for isoform-level biological interpretation of
transcriptomic data. It uses annotated transcripts and their types
(protein coding, retained intron, nonsense-mediated decay, …) to
describe how expression is distributed across the isoforms of each gene,
and how that distribution changes between conditions.

The package is organised around the `IsoformicExperiment` class: a
`SummarizedExperiment`-like container that keeps transcript, gene and
exon level data for the same samples consistent, backed by Parquet files
on disk and read through DuckDB.

## Installation

``` r

install.packages("isoformic", repos = c("https://luciorq.r-universe.dev", getOption("repos")))
```

The development version can be installed from
[GitHub](https://github.com/luciorq/isoformic):

``` r

# install.packages("pak")
pak::pak("luciorq/isoformic")
```

## Quick start

``` r

library(isoformic)

# A bundled experiment: isoforms of 69 genes, 25 placenta samples, GENCODE v33
iso <- example_isoformic()
iso
#> <IsoformicExperiment> pe_placenta_example
#>   annotation : gencode_v33_example_subset (69 genes, 892 transcripts, 7,032 exons)
#>   samples    : 25 [condition]
#>   transcripts: 499 assayed from 69 genes
#>   assays     : abundance, counts
#>   gene assays: derived from transcript assays on demand
#>   DE results : det (499 transcripts), deg (69 genes)

# Three levels, one object
assay(iso, "abundance")[1:3, 1:3]   # transcripts x samples
gene_assay(iso, "abundance")[1:3, 1:3] # genes x samples, summed from isoforms
head(exon_table(iso, "RBPJ"))        # exon structures from the annotation

# Transcriptomic diversity
head(tx_diversity(iso))               # expressed isoforms, dominant isoform, effective number
head(transcript_features(iso))        # ORF, UTRs, NMD rule, promoters, splicing events vs canonical
head(transcript_productivity(iso))    # productive / altered / unproductive / noncoding, with reasons
head(productive_fraction(iso, group = "condition")) # productive share of each gene and its change
dominant_isoform(iso, group = "condition") # dominance switches between groups
isoform_switches(iso)                 # transcript vs gene level differential expression

# Plots
plot_log2fc(iso, feature = "RBPJ")
plot_tx_expr(iso, gene_name = "RBPJ", group = "condition")
plot_genomic_context(iso, gene_name = "RBPJ")
```

Any GFF3/GTF annotation works: GENCODE (GRCh38, GRCm39), Ensembl (all
divisions), NCBI RefSeq (any assembly accession), MANE, FANTOM CAT /
FANTOM6 (`download_reference(reference = "fantom")`,
[`add_fantom_cat_tables()`](https://luciorq.github.io/isoformic/reference/add_fantom_cat_tables.md)),
the telomere-to-telomere T2T-CHM13 annotations
(`download_reference(reference = "t2t")`, or the consortium’s
Liftoff/CAT files whose GRCh38 identifiers are kept), HPRC per-haplotype
CAT annotations
([`download_hprc_annotation()`](https://luciorq.github.io/isoformic/reference/download_hprc_annotation.md) +
[`merge_annotations()`](https://luciorq.github.io/isoformic/reference/merge_annotations.md) +
[`haplotype_summary()`](https://luciorq.github.io/isoformic/reference/diversity.md)),
and long-read catalogs (StringTie, Bambu, FLAIR, with SQANTI3 categories
through
[`add_sqanti_classification()`](https://luciorq.github.io/isoformic/reference/add_sqanti_classification.md)).

Your own data enters through one of three doors:

``` r

annot <- prepare_isoformic_annotation("gencode.v49.annotation.gff3.gz") # once, cached as Parquet
iso <- read_salmon(c(sample1 = "quants/sample1", sample2 = "quants/sample2"), annot_path = annot)
iso <- as_isoformic(se_from_tximeta, annot_path = annot)
iso <- IsoformicExperiment(annot_path = annot, assays = list(abundance = tpm_matrix), col_data = samples)

de_tx(iso) <- swish_results    # any tool's columns are standardised
de_gene(iso) <- deseq2_results
iso <- write_isoformic(iso, "my_experiment") # Parquet directory, restored with read_isoformic()
```

See the articles on the [package
website](https://luciorq.github.io/isoformic/) for the full workflow,
from Salmon quantification to enrichment of transcript types.

## How confident are we in an isoform-level result?

Salmon draws replicates from the posterior of the abundances
([`salmon_quant()`](https://luciorq.github.io/isoformic/reference/salmon_quant.md)
asks for 30 Gibbs samples) and writes the equivalence classes.
`read_salmon(inferential = TRUE)` or
[`add_uncertainty()`](https://luciorq.github.io/isoformic/reference/uncertainty.md)
turn them into the uncertainty of each isoform
([`tx_uncertainty()`](https://luciorq.github.io/isoformic/reference/uncertainty.md)),
the pairs of isoforms that cannot be told apart
([`isoform_ambiguity()`](https://luciorq.github.io/isoformic/reference/uncertainty.md))
and the stability of every fold-change over the draws
([`de_confidence()`](https://luciorq.github.io/isoformic/reference/uncertainty.md),
with a `high`/`moderate`/`low` class), with plots
([`plot_de_confidence()`](https://luciorq.github.io/isoformic/reference/uncertainty-plots.md),
[`plot_isoform_correlation()`](https://luciorq.github.io/isoformic/reference/uncertainty-plots.md)),
summary tables
([`uncertainty_report()`](https://luciorq.github.io/isoformic/reference/uncertainty.md))
and a methods paragraph
([`uncertainty_methods()`](https://luciorq.github.io/isoformic/reference/uncertainty.md))
for reporting. `isoform_switches(confident_only = TRUE)` keeps only the
confident calls.
[`inferential_overdispersion()`](https://luciorq.github.io/isoformic/reference/inferential-overdispersion.md)
and
[`scale_counts()`](https://luciorq.github.io/isoformic/reference/inferential-overdispersion.md)
divide the read-ambiguity overdispersion out of the counts,
[`test_usage()`](https://luciorq.github.io/isoformic/reference/test_usage.md)
tests differential transcript usage with `edgeR` `diffSplice`, `satuRn`
or `DEXSeq` on those counts,
[`read_terminus_groups()`](https://luciorq.github.io/isoformic/reference/read_terminus_groups.md)
brings in the groups of Terminus, and
[`read_oarfish()`](https://luciorq.github.io/isoformic/reference/read_oarfish.md)
reads long-read quantifications from Oarfish with their bootstrap
replicates. The [uncertainty
article](https://luciorq.github.io/isoformic/articles/isoformic-uncertainty.html)
works through two real ESRP knockdown datasets (prostate cancer
GSE129540 and H358 cells GSE75489): `SCRIB-207` is a high confidence
call, while the `CD44` and `CTNND1` switches have a certain direction
but isoforms that the reads cannot tell apart.

## Productive versus unproductive transcription

With transcript sequences in the annotation
([`extract_transcript_sequences()`](https://luciorq.github.io/isoformic/reference/transcript_sequences.md)
from a genome or transcript FASTA), productivity is computed from the
sequence rather than read from the biotype: the canonical start codon is
projected onto every isoform, translated, and checked for premature
stops, frameshifts, non-stop decay and uORFs. External evidence (APPRIS
principal isoforms, GENCODE Ribo-seq ORFs, validated NMD exons, GTEx
tissue baselines, Pfam domain retention) comes from frozen, versioned
datasets:

``` r

iso <- example_isoformic()                     # ships with sequences
transcript_orfs(iso, gene_name = "RBPJ")       # ORF, NMD class, coding probability
transcript_productivity(iso, gene_name = "RBPJ")
test_productive_fraction(iso, group = "condition")   # quasi-binomial test per gene
polya_sites(iso, gene_name = "RBPJ")           # alternative polyadenylation events
apa_usage(iso, group = "condition")

add_appris(iso); add_riboseq_orfs(iso); add_nmd_exons(iso)   # frozen datasets
productive_fraction(iso, group = "condition", reference = "Whole Blood")
```

## Citation

If you use `isoformic` in publications, please cite:

Izabela Mamede, Lucio R Queiroz, Carlos Mata-Machado, Júlia Teixeira
Rodrigues, Thomaz Luscher-Dias, Nayara E de Toledo, Paulo P Amaral,
Luigi Marchionni, Gloria R Franco, **Isoformic**: a workflow for
transcript-level RNA-seq interpretation, *NAR Genomics and
Bioinformatics*, Volume 7, Issue 4, December 2025, lqaf176,
<https://doi.org/10.1093/nargab/lqaf176>

## Motivation

![Types of Transcripts](https://i.imgur.com/UWoAr0k.png)

Types of Transcripts

Differential expression at the isoform level, rather than at the gene
level, captures regulation that gene-level summaries hide: alternative
splicing and alternative transcription start sites produce isoforms with
distinct structures and functions, and a gene whose total expression
does not change can still switch from a protein coding isoform to an
unproductive one. `isoformic` makes those changes visible by keeping
transcript types, exon structures and gene-level summaries together and
by testing each transcript category on its own during functional
enrichment.
