
# isoformic <img src="man/figures/logo.png" align="right" height="136" />

<!-- badges: start -->
<!-- badges: end -->

Welcome to `isoformic`, a workflow for isoform-level biological interpretation of transcriptomic data.
This workflow uses known annotated transcripts to produce biologically relevant results based on the different types of transcripts for any comparison of case versus control transcriptomic data.

## Installation

You can install the development version of `isoformic` from [GitHub][github-ref] with:

``` r
if (requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("luciorq/isoformic")
```

## Overview

![trancript_types](https://i.imgur.com/UWoAr0k.png)

## Motivation

Performing differential expression analysis at the isoform level, rather than the gene level, is important for gaining a more detailed understanding of gene regulation and functional differences in biological systems.
While gene-level analysis provides insights into the overall expression changes of genes, it fails to capture the complex dynamics occurring at the isoform level.
Isoforms, resulting from alternative splicing or alternative transcription start sites, can possess distinct structural features and functional properties.
By examining isoform-level expression changes, researchers can identify and characterize specific isoforms that may have unique roles in cellular processes, such as isoform-specific protein-protein interactions or protein functions. This granularity enables a more accurate interpretation of complex biological phenomena, such as tissue-specific expression patterns, cell differentiation, or disease progression.
Additionally, isoform-level analysis can help uncover regulatory mechanisms and identify potential biomarkers or therapeutic targets that may be missed when solely relying on gene-level analysis.
Therefore, considering the isoform-level expression changes provides a more comprehensive view of transcriptional dynamics and enhances our understanding of gene regulation in a given biological context.

---

[github-ref]: https://github.com/
