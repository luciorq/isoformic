# From Salmon quantifications to an IsoformicExperiment

This article goes from FASTQ files to an `IsoformicExperiment` with
transcript-level and gene-level differential expression results. The
[introduction](https://luciorq.github.io/isoformic/articles/isoformic-introduction.md)
starts one step later, from ready-made tables, and the other articles
analyse the object.

## Quantifying with Salmon

We recommend [Salmon](https://github.com/COMBINE-lab/salmon) ([Patro et
al. 2017](#ref-salmon)) in selective-alignment mode against the GENCODE
transcriptome, with inferential replicates so that `swish` ([Zhu et al.
2019](#ref-swish)) can be used for transcript-level differential
expression.

Salmon 2 is a from-scratch rewrite in Rust (2.8.0 at the time of
writing). The workflow (`salmon index`, `salmon quant`, `quant.sf`) and
the output files are the same, so
[`read_salmon()`](https://luciorq.github.io/isoformic/reference/read_salmon.md),
`tximeta`, `tximport` and `swish` work unchanged, but the index format
is new: an index built with Salmon 1.x is rejected and must be rebuilt
(the last C++ release, 1.12.0, remains available as the `salmon-cpp`
conda package). Other changes worth knowing: selective alignment is the
only mapping mode (`--validateMappings` is a no-op), quantification is
deterministic by default (identical output for any thread count), decoy
sequences no longer appear in `quant.sf`, and `--ignoreTxVersion`
matches transcript identifiers without their version suffix when a gene
map is given. The Salmon [migration
guide](https://github.com/COMBINE-lab/salmon/blob/master/MIGRATION.md)
lists every removed or ignored option.

A typical single-sample invocation from a shell is below: `--libType A`
detects the library type, the index is built from the GENCODE transcript
FASTA with `salmon index --gencode` (which truncates the transcript
names at the first `|`), `--posBias --seqBias --gcBias` correct the
biases, `--dumpEq --dumpEqWeights` keep the equivalence classes and
`--numGibbsSamples` writes the inferential replicates that `swish` needs
into `aux_info/bootstrap/`.

``` bash
salmon index \
  --transcripts gencode.v50.transcripts.fa.gz \
  --index salmon_index \
  --gencode \
  --threads 8

salmon quant \
  --libType A \
  --index salmon_index \
  --mates1 /path/SRR..._1.fastq.gz \
  --mates2 /path/SRR..._2.fastq.gz \
  --output /path/SRR..._quant \
  --threads 8 \
  --allowDovetail --softclipOverhangs \
  --dumpEq --dumpEqWeights \
  --posBias --seqBias --gcBias \
  --numGibbsSamples 30
```

The same can be driven from R through the `condathis` package, which
installs Salmon from Bioconda into a conda environment (Linux and macOS,
or WSL2 on Windows).
[`salmon_index()`](https://luciorq.github.io/isoformic/reference/salmon_index.md)
and
[`salmon_quant()`](https://luciorq.github.io/isoformic/reference/salmon_quant.md)
run the commands above; the first call creates the environment with
Salmon 2.8.0:

``` r

library(isoformic)
tx_fasta <- download_reference(version = "50", file_type = "fasta")
gff_path <- download_reference(version = "50", file_type = "gff")

salmon_index(fasta_path = tx_fasta, index_path = "salmon_index", is_gencode = TRUE, num_threads = 8)
salmon_quant(
  input_r1 = "sample_1.fastq.gz",
  input_r2 = "sample_2.fastq.gz",
  index_path = "salmon_index",
  output_dir = "quants/sample_quant",
  num_threads = 8
)
```

[`salmon_quant()`](https://luciorq.github.io/isoformic/reference/salmon_quant.md)
draws 30 Gibbs samples per sample by default (`num_gibbs = 0` skips
them). Besides `swish`, they quantify how well the reads determine each
isoform, which isoforms cannot be told apart and how stable the
direction of each differential expression call is over the replicates;
`read_salmon(inferential = TRUE)` reads them and the [uncertainty
article](https://luciorq.github.io/isoformic/articles/isoformic-uncertainty.md)
explains how to use and report them.

For a decoy-aware index (recommended by the Salmon authors), pass the
genome FASTA as `decoy_fasta`, for GENCODE the
`GRCh38.primary_assembly.genome.fa.gz` of the same release:
[`salmon_index()`](https://luciorq.github.io/isoformic/reference/salmon_index.md)
writes the `gentrome` FASTA (transcripts followed by the genome) and
`decoys.txt` (the genome sequence names) in a temporary directory, runs
`salmon index --decoys decoys.txt` and removes them. The genome makes
the index much larger and slower to build, and needs free disk space in
[`tempdir()`](https://rdrr.io/r/base/tempfile.html) for a copy of the
genome while it runs.

## Reading quantifications directly

[`read_salmon()`](https://luciorq.github.io/isoformic/reference/read_salmon.md)
reads any number of `quant.sf` files with DuckDB and returns an
`IsoformicExperiment` with `abundance` (TPM), `counts` and `length`
assays. The package ships the Salmon output of six placenta samples
(three with preeclampsia, three controls) quantified against the GENCODE
v50 transcripts of 15 genes.

``` r

library(isoformic)

quants_dir <- fs::path_package("isoformic", "extdata", "mini_quants")
quant_dirs <- fs::dir_ls(quants_dir, type = "directory")
names(quant_dirs) <- sub("_quant$", "", fs::path_file(quant_dirs))

col_data <- tibble::tibble(
  sample_id = names(quant_dirs),
  condition = factor(c(rep("preeclampsia", 3), rep("control", 3)), levels = c("control", "preeclampsia"))
)
col_data
#> # A tibble: 6 × 2
#>   sample_id   condition   
#>   <chr>       <fct>       
#> 1 SRR11498039 preeclampsia
#> 2 SRR11498040 preeclampsia
#> 3 SRR11498041 preeclampsia
#> 4 SRR11498077 control     
#> 5 SRR11498078 control     
#> 6 SRR11498079 control
```

The annotation must match the reference used for quantification, GENCODE
v50 here. For real data it is downloaded once into the package cache and
converted to Parquet once:

``` r

gff_path <- download_reference(version = "50", file_type = "gff")
annot_dir <- prepare_isoformic_annotation(gff_path)
```

The package bundles the subset of that annotation covering the
transcripts of the mini dataset.

``` r

annot_dir <- fs::path_package("isoformic", "extdata", "mini_annotation")

iso <- read_salmon(quant_dirs, annot_path = annot_dir, col_data = col_data, experiment_name = "mini_placenta")
iso
#> <IsoformicExperiment> mini_placenta
#>   annotation : gencode_v50_mini_subset (15 genes, 158 transcripts, 1,802 exons)
#>   samples    : 6 [condition]
#>   transcripts: 158 assayed from 15 genes
#>   assays     : abundance, counts, length
#>   gene assays: derived from transcript assays on demand
#>   DE results : <none>
assay_names(iso)
#> [1] "abundance" "counts"    "length"
head(row_data(iso))
#> # A tibble: 6 × 19
#>   transcript_id      transcript_name transcript_type gene_id gene_name gene_type
#>   <chr>              <chr>           <chr>           <chr>   <chr>     <chr>    
#> 1 ENST00000216797.10 NFKBIA-201      protein_coding  ENSG00… NFKBIA    protein_…
#> 2 ENST00000225831.4  CCL2-201        protein_coding  ENSG00… CCL2      protein_…
#> 3 ENST00000258743.10 IL6-201         protein_coding  ENSG00… IL6       protein_…
#> 4 ENST00000263341.7  IL1B-201        protein_coding  ENSG00… IL1B      protein_…
#> 5 ENST00000264832.8  ICAM1-201       protein_coding  ENSG00… ICAM1     protein_…
#> 6 ENST00000306602.3  CXCL10-201      protein_coding  ENSG00… CXCL10    protein_…
#> # ℹ 13 more variables: source_gene_id <chr>, source_transcript_id <chr>,
#> #   tags <chr>, coding_status <chr>, transcript_support_level <chr>,
#> #   tie_score <dbl>, cds_start <int>, cds_end <int>, cds_length <int>,
#> #   is_canonical <lgl>, haplotype <chr>, n_exons <int>, tx_length <int>
```

## Transcript-level differential expression with swish

`swish` needs the inferential replicates, which `tximeta` imports from
the Salmon output (the bundled quantifications carry their 20 Gibbs
draws separately, in `mini_draws.parquet`, to keep the package small, so
`tximeta` is told to skip them).

``` r

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(fishpond)
})

tximeta_coldata <- data.frame(
  files = fs::path(quant_dirs, "quant.sf"),
  names = names(quant_dirs),
  condition = col_data$condition
)
se_tx <- tximeta::tximeta(
  coldata = tximeta_coldata,
  type = "salmon",
  txOut = TRUE,
  skipMeta = TRUE,
  useHub = FALSE,
  ignoreAfterBar = TRUE,
  skipSeqinfo = TRUE,
  dropInfReps = TRUE # the draws are not in the bundled directories
)
mini_draws <- as_inferential_replicates(
  read_parquet(fs::path_package("isoformic", "extdata", "mini_draws.parquet"))
)
names(mini_draws) <- paste0("infRep", seq_along(mini_draws))
SummarizedExperiment::assays(se_tx) <- c(
  as.list(SummarizedExperiment::assays(se_tx)),
  mini_draws
)

y <- fishpond::scaleInfReps(se_tx)
y <- fishpond::labelKeep(y, minCount = 10, minN = 3)
set.seed(20261029) # swish uses permutations
y <- fishpond::swish(y, x = "condition")

swish_results <- mcols(y) |>
  as.data.frame() |>
  tibble::rownames_to_column("transcript_id")
head(swish_results)
#>        transcript_id log10mean  keep stat log2FC pvalue locfdr qvalue
#> 1 ENST00000367468.10 1.0366158 FALSE   NA     NA     NA     NA     NA
#> 2  ENST00000681605.1 0.7665459 FALSE   NA     NA     NA     NA     NA
#> 3  ENST00000680451.1 1.0361596 FALSE   NA     NA     NA     NA     NA
#> 4  ENST00000490885.6 1.0145195 FALSE   NA     NA     NA     NA     NA
#> 5  ENST00000559627.1 0.8826379 FALSE   NA     NA     NA     NA     NA
#> 6  ENST00000466691.1 0.3940242 FALSE   NA     NA     NA     NA     NA
```

The results are attached to the object; `tximeta` objects can also be
converted directly with `as_isoformic(se_tx, annot_path = annot_dir)`.

``` r

de_tx(iso) <- swish_results
head(de_tx(iso))
#> # A tibble: 6 × 9
#>   transcript_id     log10mean keep   stat log2FC  pvalue locfdr qvalue det_sig
#>   <chr>                 <dbl> <lgl> <dbl>  <dbl>   <dbl>  <dbl>  <dbl> <chr>  
#> 1 ENST00000907726.2      1.97 TRUE   4.3    2.69 0.00545  0.176  0.18  no     
#> 2 ENST00000306602.3      1.45 TRUE   0.5    1.60 0.779    1      0.951 no     
#> 3 ENST00000580907.6      1.13 TRUE  -4.1   -1.46 0.0155   0.276  0.255 no     
#> 4 ENST00000361099.8      1.65 TRUE  -1.7   -1.43 0.343    1      0.566 no     
#> 5 ENST00000225831.4      1.69 TRUE  -2.9   -1.36 0.0830   0.599  0.404 no     
#> 6 ENST00000540176.6      1.51 TRUE   3.35   1.31 0.0448   0.446  0.296 no
```

## Gene-level differential expression with DESeq2

[`as_summarized_experiment()`](https://luciorq.github.io/isoformic/reference/isoformic-conversion.md)
exports the gene-level assays (isoforms summed, lengths averaged by
abundance) so that `DESeq2` can run on them.

``` r

se_gene <- as_summarized_experiment(iso, level = "gene")
SummarizedExperiment::assay(se_gene, "counts") <- round(
  SummarizedExperiment::assay(se_gene, "counts")
)
se_gene$condition <- factor(se_gene$condition, levels = c("control", "preeclampsia"))

dds <- DESeq2::DESeqDataSet(se_gene, design = ~condition)
dds <- dds[rowSums(DESeq2::counts(dds) >= 10) >= 3, ]
dds <- DESeq2::DESeq(dds, quiet = TRUE)
res <- DESeq2::results(dds, contrast = c("condition", "preeclampsia", "control"))
shrink_type <- if (rlang::is_installed("apeglm")) "apeglm" else "normal"
res_shr <- DESeq2::lfcShrink(dds, coef = "condition_preeclampsia_vs_control", res = res, type = shrink_type, quiet = TRUE)

de_gene(iso) <- as.data.frame(res_shr)
head(de_gene(iso))
#> # A tibble: 6 × 7
#>   gene_id            baseMean log2FC lfcSE   pvalue  qvalue deg_sig
#>   <chr>                 <dbl>  <dbl> <dbl>    <dbl>   <dbl> <chr>  
#> 1 ENSG00000184557.6     256.   1.16  0.430 0.000601 0.00722 yes    
#> 2 ENSG00000108691.11     97.1 -0.802 0.899 0.0193   0.116   no     
#> 3 ENSG00000090339.11    250.   0.435 0.392 0.0806   0.242   no     
#> 4 ENSG00000274233.4     106.  -0.297 0.422 0.170    0.407   no     
#> 5 ENSG00000073756.13     47.2  0.257 0.509 0.0405   0.162   no     
#> 6 ENSG00000115415.23    779.  -0.180 0.213 0.318    0.504   no
iso
#> <IsoformicExperiment> mini_placenta
#>   annotation : gencode_v50_mini_subset (15 genes, 158 transcripts, 1,802 exons)
#>   samples    : 6 [condition]
#>   transcripts: 158 assayed from 15 genes
#>   assays     : abundance, counts, length
#>   gene assays: derived from transcript assays on demand
#>   DE results : det (158 transcripts), deg (12 genes)
```

## Looking at the results

``` r

top_gene <- de_gene(iso) |>
  dplyr::arrange(qvalue) |>
  dplyr::slice_head(n = 1) |>
  dplyr::pull(gene_id)
gene_name <- gene_data(iso)$gene_name[gene_data(iso)$gene_id == top_gene]
gene_name
#> [1] "SOCS3"

plot_log2fc(iso, feature = gene_name)
```

![](isoformic-preprocess_files/figure-html/unnamed-chunk-8-1.png)

``` r

plot_tx_expr(iso, gene_name = gene_name, group = "condition")
```

![](isoformic-preprocess_files/figure-html/unnamed-chunk-8-2.png)

``` r

plot_genomic_context(iso, gene_name = gene_name)
```

![](isoformic-preprocess_files/figure-html/unnamed-chunk-8-3.png)

``` r

isoform_switches(iso) |> head()
#> # A tibble: 0 × 12
#> # ℹ 12 variables: gene_id <chr>, gene_name <chr>, gene_type <chr>,
#> #   n_tx_tested <int>, n_det_sig <int>, n_det_up <int>, n_det_down <int>,
#> #   det_sig_types <chr>, max_abs_log2FC <dbl>, gene_log2FC <dbl>,
#> #   deg_sig <chr>, switch_type <chr>
```

## Saving

The object, with both result tables, is written as a directory of
Parquet files and restored with
[`read_isoformic()`](https://luciorq.github.io/isoformic/reference/isoformic-io.md).

``` r

out_dir <- fs::path(tempdir(), "mini_placenta")
iso <- write_isoformic(iso, out_dir, overwrite = TRUE)
iso@data_path
#> [1] "/tmp/RtmpPPTSB0/mini_placenta"
fs::dir_ls(out_dir)
#> /tmp/RtmpPPTSB0/mini_placenta/assays
#> /tmp/RtmpPPTSB0/mini_placenta/col_data.parquet
#> /tmp/RtmpPPTSB0/mini_placenta/dea
#> /tmp/RtmpPPTSB0/mini_placenta/experiment.parquet
#> /tmp/RtmpPPTSB0/mini_placenta/tx_type_palette.parquet
```

## Long reads with Oarfish

Oxford Nanopore and PacBio reads are quantified with Oarfish, which
writes bootstrap replicates and a table of unique and shared reads next
to the counts.
[`read_oarfish()`](https://luciorq.github.io/isoformic/reference/read_oarfish.md)
reads them like
[`read_salmon()`](https://luciorq.github.io/isoformic/reference/read_salmon.md)
and the uncertainty functions of the [uncertainty
article](https://luciorq.github.io/isoformic/articles/isoformic-uncertainty.md)
work on them unchanged. The three samples below are simulated nanopore
reads of the 141 GENCODE v50 transcripts of RBPJ, PNCK and FLT1
(`example_transcripts.fa.gz`), quantified with Oarfish 0.10.3 and 30
bootstrap replicates (`data-raw/oarfish/`), so they only show the file
layout and the functions, not a real library.

``` r

oarfish_dir <- fs::path_package("isoformic", "extdata", "oarfish")
prefixes <- fs::path(oarfish_dir, c("control_1", "control_2", "treated_1"))
names(prefixes) <- fs::path_file(prefixes)
fs::dir_ls(oarfish_dir, regexp = "treated_1")
#> /home/runner/work/_temp/Library/isoformic/extdata/oarfish/treated_1.ambig_info.tsv
#> /home/runner/work/_temp/Library/isoformic/extdata/oarfish/treated_1.infreps.pq
#> /home/runner/work/_temp/Library/isoformic/extdata/oarfish/treated_1.meta_info.json
#> /home/runner/work/_temp/Library/isoformic/extdata/oarfish/treated_1.quant

iso_long <- read_oarfish(
  prefixes,
  annot_path = fs::path_package("isoformic", "extdata", "example_experiment", "annotation"),
  inferential = TRUE
)
iso_long
#> <IsoformicExperiment> isoformic_95e441dbd6
#>   annotation : gencode_v50_example_subset (69 genes, 1,914 transcripts, 18,656 exons)
#>   samples    : 3 [oarfish_num_processed, oarfish_num_mapped, oarfish_percent_mapped, oarfish_replicate_type, oarfish_num_replicates, oarfish_optimiser, oarfish_alignment_source, oarfish_prob_model, oarfish_qc_flag]
#>   transcripts: 141 assayed from 3 genes
#>   assays     : abundance, counts, length
#>   gene assays: derived from transcript assays on demand
#>   DE results : <none>
tx_uncertainty(iso_long) |>
  dplyr::arrange(dplyr::desc(mean_inf_rv)) |>
  dplyr::select(transcript_name, mean_count, mean_inf_rv, unique_fraction) |>
  head(6)
#> # A tibble: 6 × 4
#>   transcript_name mean_count mean_inf_rv unique_fraction
#>   <chr>                <dbl>       <dbl>           <dbl>
#> 1 RBPJ-202              62.3       28.7          0      
#> 2 RBPJ-204              40.9       23.3          0      
#> 3 RBPJ-238              60.7       15.3          0      
#> 4 PNCK-233              18.5        7.71         0      
#> 5 RBPJ-240              33.8        7.58         0.00127
#> 6 PNCK-239              20.2        6.65         0
```

The Oarfish output has no TPM, so the `abundance` assay is computed from
the read counts and the transcript lengths. `unique_fraction` comes from
the unique and shared read counts of `.ambig_info.tsv`; Oarfish has no
equivalence classes, so `n_classes` is not available. Samples quantified
by Oarfish and by Salmon cannot be read together.

In single-cell mode (`--single-cell`, on a transcriptome-aligned BAM
whose reads carry their cell barcode in the `CB:Z` tag) Oarfish writes a
cell by transcript count matrix and no replicates.
[`read_oarfish_sc()`](https://luciorq.github.io/isoformic/reference/read_oarfish_sc_counts.md)
sums the cells of each sample, given a table of barcodes and samples, so
the isoform analyses apply to samples or to cell types (give each
combination of sample and cell type its own `sample_id`). The example
has 60 simulated cells in three samples:

``` r

sc_dir <- fs::path_package("isoformic", "extdata", "oarfish_sc")
cells <- read.delim(fs::path(sc_dir, "cells.tsv"))
head(cells, 3)
#>            barcode sample_id
#> 1 CACATCGAACAGCTAA control_1
#> 2 ATACAAGACTGGCTCG control_1
#> 3 GTTCATGACTTAACCG control_1
iso_sc <- read_oarfish_sc(
  fs::path(sc_dir, "out"),
  cells = cells,
  annot_path = fs::path_package("isoformic", "extdata", "example_experiment", "annotation")
)
col_data(iso_sc)
#> # A tibble: 3 × 3
#>   sample_id n_cells n_reads
#>   <chr>       <int>   <dbl>
#> 1 control_1      20   3000.
#> 2 control_2      20   3000.
#> 3 treated_1      20   3000.
```

## Next steps

From here the object is analysed like the one of the introduction:
[transcriptomic
diversity](https://luciorq.github.io/isoformic/articles/isoformic-diversity.md),
[isoform structure and
productivity](https://luciorq.github.io/isoformic/articles/isoformic-productivity.md),
[alternative
polyadenylation](https://luciorq.github.io/isoformic/articles/isoformic-polyadenylation.md),
[enrichment by transcript
type](https://luciorq.github.io/isoformic/articles/isoformic-enrichment.md)
and
[plots](https://luciorq.github.io/isoformic/articles/isoformic-plots.md).
The inferential replicates imported above can be stored in
`iso@metadata$inf_reps` so that
[`productive_fraction()`](https://luciorq.github.io/isoformic/reference/productivity.md)
reports the uncertainty of each fraction.

## Citing isoformic

If you use `isoformic` in publications, please cite Mamede et al.
([2025](#ref-isoformic)); `citation("isoformic")` prints the reference.

## Session information

``` r

sessionInfo()
#> R version 4.6.2 (2026-10-29)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 26.04.1 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.32.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats4    stats     graphics  grDevices utils     datasets  methods  
#> [8] base     
#> 
#> other attached packages:
#>  [1] fishpond_2.18.0             SummarizedExperiment_1.42.0
#>  [3] Biobase_2.72.0              GenomicRanges_1.64.0       
#>  [5] Seqinfo_1.2.0               IRanges_2.46.0             
#>  [7] S4Vectors_0.50.3            BiocGenerics_0.58.1        
#>  [9] generics_0.1.4              MatrixGenerics_1.24.0      
#> [11] matrixStats_1.5.0           isoformic_0.1.3            
#> 
#> loaded via a namespace (and not attached):
#>   [1] RColorBrewer_1.1-3          jsonlite_2.0.0             
#>   [3] tximport_1.40.0             magrittr_2.0.5             
#>   [5] GenomicFeatures_1.64.0      farver_2.1.2               
#>   [7] rmarkdown_2.32              fs_2.1.0                   
#>   [9] BiocIO_1.22.0               vctrs_0.7.3                
#>  [11] memoise_2.0.1               Rsamtools_2.28.0           
#>  [13] RCurl_1.98-1.20             htmltools_0.5.9            
#>  [15] S4Arrays_1.12.1             progress_1.2.3             
#>  [17] AnnotationHub_4.2.2         curl_8.1.0                 
#>  [19] SparseArray_1.12.3          plyr_1.8.9                 
#>  [21] httr2_1.3.0                 cachem_1.1.0               
#>  [23] GenomicAlignments_1.48.0    lifecycle_1.0.5            
#>  [25] pkgconfig_2.0.3             Matrix_1.7-6               
#>  [27] R6_2.6.1                    fastmap_1.2.0              
#>  [29] digest_0.6.39               numDeriv_2016.8-1.1        
#>  [31] AnnotationDbi_1.74.0        DESeq2_1.52.0              
#>  [33] tximeta_1.30.1              RSQLite_3.53.3             
#>  [35] filelock_1.0.3              labeling_0.4.3             
#>  [37] httr_1.4.9                  abind_1.4-8                
#>  [39] compiler_4.6.2              bit64_4.8.6                
#>  [41] withr_3.0.3                 S7_1.0.0                   
#>  [43] BiocParallel_1.46.0         DBI_1.3.0                  
#>  [45] duckdb_2.0.0                biomaRt_2.68.6             
#>  [47] MASS_7.3-66                 rappdirs_0.3.4             
#>  [49] DelayedArray_0.38.2         rjson_0.2.23               
#>  [51] gtools_3.9.5                tools_4.6.2                
#>  [53] otel_0.3.0                  glue_1.8.1                 
#>  [55] restfulr_0.0.17             grid_4.6.2                 
#>  [57] reshape2_1.4.5              gtable_0.3.6               
#>  [59] tzdb_0.5.0                  ensembldb_2.36.1           
#>  [61] hms_1.1.4                   utf8_1.2.6                 
#>  [63] XVector_0.52.0              ggrepel_0.9.8              
#>  [65] BiocVersion_3.23.1          pillar_1.11.1              
#>  [67] stringr_1.6.0               emdbook_1.3.14             
#>  [69] vroom_1.7.1                 splines_4.6.2              
#>  [71] dplyr_1.2.1                 BiocFileCache_3.2.0        
#>  [73] lattice_0.23-1              rtracklayer_1.72.0         
#>  [75] bit_4.6.0                   tidyselect_1.2.1           
#>  [77] SingleCellExperiment_1.34.0 locfit_1.5-9.12            
#>  [79] Biostrings_2.80.2           knitr_1.52                 
#>  [81] ProtGenerics_1.44.0         xfun_0.62                  
#>  [83] stringi_1.8.9               UCSC.utils_1.8.0           
#>  [85] lazyeval_0.2.3              yaml_2.3.12                
#>  [87] evaluate_1.0.5              codetools_0.2-20           
#>  [89] cigarillo_1.2.1             bbmle_1.0.26               
#>  [91] tibble_3.3.1                qvalue_2.44.0              
#>  [93] BiocManager_1.30.27         cli_3.6.6                  
#>  [95] Rcpp_1.1.2                  GenomeInfoDb_1.48.0        
#>  [97] svMisc_1.4.3                dbplyr_2.6.0               
#>  [99] coda_0.19-4.1               png_0.1-9                  
#> [101] bdsmatrix_1.3-7             XML_3.99-0.25              
#> [103] parallel_4.6.2              ggplot2_4.0.4              
#> [105] readr_2.2.0                 blob_1.3.0                 
#> [107] prettyunits_1.2.0           AnnotationFilter_1.36.0    
#> [109] bitops_1.1-0                txdbmaker_1.8.0            
#> [111] mvtnorm_1.4-2               apeglm_1.34.0              
#> [113] scales_1.4.0                crayon_1.5.3               
#> [115] rlang_1.3.0                 KEGGREST_1.52.2
```

## References

Mamede, Izabela, Lucio R. Queiroz, Carlos Mata-Machado, et al. 2025.
“Isoformic: A Workflow for Transcript-Level RNA-seq Interpretation.”
*NAR Genomics and Bioinformatics* 7 (4): lqaf176.
<https://doi.org/10.1093/nargab/lqaf176>.

Patro, R, G Duggal, MI Love, RA Irizarry, and C Kingsford. 2017. “Salmon
Provides Fast and Bias-Aware Quantification of Transcript Expression.”
*Nature Methods* 14: 417–19. <https://doi.org/10.1038/nmeth.4197>.

Zhu, A, A Srivastava, JG Ibrahim, R Patro, and MI Love. 2019.
“Nonparametric expression analysis using inferential replicate counts.”
*Nucleic Acids Research* 47: e105. <https://doi.org/10.1093/nar/gkz622>.
