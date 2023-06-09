% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{make_tx_to_gene}
\alias{make_tx_to_gene}
\title{Create Transcript to Gene Relationship Table}
\usage{
make_tx_to_gene(file_path, file_type = "fasta")
}
\arguments{
\item{file_path}{Path to file containing the reference.}

\item{file_type}{One of gff, fasta, or gtf.}
}
\description{
This function can extract a transcript to gene relationship table
from GENCODE annotation files, such as the transcriptome FASTA and
GFF3 or GTF annotation files.
}
