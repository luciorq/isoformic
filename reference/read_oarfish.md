# Read Oarfish Quantifications Into an IsoformicExperiment

Oarfish quantifies long reads (Oxford Nanopore and PacBio) against a
transcriptome. For every sample it writes, next to a common output
prefix, `<prefix>.quant` (transcript name, length and estimated read
count), `<prefix>.ambig_info.tsv` (reads that map to the transcript
alone and reads shared with others), `<prefix>.meta_info.json` and, with
`--num-bootstraps`, `<prefix>.infreps.pq` (the bootstrap replicates as a
Parquet table). `read_oarfish()` reads them like
[`read_salmon()`](https://luciorq.github.io/isoformic/reference/read_salmon.md)
reads Salmon output: an
[IsoformicExperiment](https://luciorq.github.io/isoformic/reference/IsoformicExperiment.md)
with `counts`, `length` and `abundance` assays, and, with
`inferential = TRUE`, the uncertainty tables of
[`add_uncertainty()`](https://luciorq.github.io/isoformic/reference/uncertainty.md).

## Usage

``` r
read_oarfish(
  paths,
  annot_path,
  col_data = NULL,
  experiment_name = NULL,
  ignore_after_bar = TRUE,
  drop_unannotated = FALSE,
  inferential = FALSE
)
```

## Arguments

- paths:

  Output prefixes of Oarfish (the value of `-o`), paths to the `.quant`
  files, or directories holding a single `.quant` file. Names, when
  present, are the sample IDs; otherwise the prefix is used.

- annot_path:

  Annotation directory or GFF3/GTF file (see
  [`prepare_isoformic_annotation()`](https://luciorq.github.io/isoformic/reference/prepare_isoformic_annotation.md)).

- col_data:

  Optional data frame of sample metadata whose first column (or a
  `sample_id` column) matches the sample IDs. When `NULL`, a minimal
  table with just `sample_id` is created.

- experiment_name:

  Optional experiment name.

- ignore_after_bar:

  Logical. Keep only the part of transcript names before the first `|`
  (GENCODE FASTA headers used as Salmon references).

- drop_unannotated:

  Logical. Drop transcripts absent from the annotation (with a message)
  instead of failing.

- inferential:

  Logical. Read the bootstrap replicates and the ambiguity counts and
  compute the uncertainty tables with
  [`add_uncertainty()`](https://luciorq.github.io/isoformic/reference/uncertainty.md).
  The samples must have been quantified with `--num-bootstraps`.

## Value

An
[IsoformicExperiment](https://luciorq.github.io/isoformic/reference/IsoformicExperiment.md).

## Details

Oarfish does not write TPM; the `abundance` assay is computed from the
read counts and the transcript lengths (reads per base, scaled to one
million). The ambiguity counts of `.ambig_info.tsv` play the role of the
equivalence classes of Salmon in
[`equivalence_support()`](https://luciorq.github.io/isoformic/reference/inferential-io.md):
`unique_fraction` is the share of reads mapped to the transcript alone,
while `n_classes` is not available. The functions
[`read_inferential_replicates()`](https://luciorq.github.io/isoformic/reference/inferential-io.md),
[`equivalence_support()`](https://luciorq.github.io/isoformic/reference/inferential-io.md),
[`salmon_run_qc()`](https://luciorq.github.io/isoformic/reference/inferential-io.md)
and
[`add_uncertainty()`](https://luciorq.github.io/isoformic/reference/uncertainty.md)
accept Oarfish prefixes as `paths`, but not together with Salmon
directories.

## Examples

``` r
prefixes <- file.path(
  system.file("extdata", "oarfish", package = "isoformic"),
  c("control_1", "control_2", "treated_1")
)
names(prefixes) <- basename(prefixes)
annot <- system.file(
  "extdata", "example_experiment", "annotation",
  package = "isoformic"
)
iso <- read_oarfish(prefixes, annot_path = annot, inferential = TRUE)
iso
#> <IsoformicExperiment> isoformic_edcd44d335
#>   annotation : gencode_v33_example_subset (69 genes, 892 transcripts, 7,032 exons)
#>   samples    : 3 [oarfish_num_processed, oarfish_num_mapped, oarfish_percent_mapped, oarfish_replicate_type, oarfish_num_replicates, oarfish_optimiser, oarfish_alignment_source, oarfish_prob_model, oarfish_qc_flag]
#>   transcripts: 72 assayed from 3 genes
#>   assays     : abundance, counts, length
#>   gene assays: derived from transcript assays on demand
#>   DE results : <none>
head(tx_uncertainty(iso))
#> # A tibble: 6 × 18
#>   transcript_id     gene_id gene_name transcript_name transcript_type mean_count
#>   <chr>             <chr>   <chr>     <chr>           <chr>                <dbl>
#> 1 ENST00000512351.5 ENSG00… RBPJ      RBPJ-222        protein_coding       59.9 
#> 2 ENST00000510778.6 ENSG00… RBPJ      RBPJ-218        protein_coding       13.7 
#> 3 ENST00000506956.5 ENSG00… RBPJ      RBPJ-213        protein_coding       21.4 
#> 4 ENST00000512671.5 ENSG00… RBPJ      RBPJ-223        protein_coding        1.19
#> 5 ENST00000345843.7 ENSG00… RBPJ      RBPJ-203        protein_coding        8.92
#> 6 ENST00000342295.5 ENSG00… RBPJ      RBPJ-201        protein_coding       41.3 
#> # ℹ 12 more variables: assessed <lgl>, mean_inf_rv <dbl>, max_inf_rv <dbl>,
#> #   mean_inf_cv <dbl>, uncertain <lgl>, min_sibling_cor <dbl>,
#> #   n_ambiguous_partners <int>, ambiguity_group <chr>, n_classes <dbl>,
#> #   n_fragments <dbl>, n_unique_fragments <dbl>, unique_fraction <dbl>
```
