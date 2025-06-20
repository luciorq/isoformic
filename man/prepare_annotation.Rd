% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prepare_annotation.R
\name{prepare_annotation}
\alias{prepare_annotation}
\title{Prepare Annotation}
\usage{
prepare_annotation(file_path, file_type = c("gtf", "gff"))
}
\arguments{
\item{file_path}{Path to annotation file.}

\item{file_type}{Character indicating the type of file to download.
One of \code{"gtf"} or \code{"gff"}. Defaults to \code{"gtf"}.}
}
\description{
Prepare annotation to be imported as \code{rowRanges} and \code{rowData} for both
Genes, Transcripts and Exons based Position Annotation Table.
From a GTF or GFF3 annotation file.
}
