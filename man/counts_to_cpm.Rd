% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/counts_to_TPM.R
\name{counts_to_cpm}
\alias{counts_to_cpm}
\title{Convert counts to counts per million (CPM).}
\usage{
counts_to_cpm(counts)
}
\arguments{
\item{counts}{A numeric matrix of raw feature counts i.e.
fragments assigned to each gene.}
}
\value{
cpm A numeric matrix normalized by library size
}
\description{
Convert a numeric matrix of features (rows) and conditions (columns) with
raw feature counts to counts per million.
}
