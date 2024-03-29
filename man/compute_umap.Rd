% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/projections.R
\name{compute_umap}
\alias{compute_umap}
\title{A function to perform *UMAP* projection of data}
\usage{
compute_umap(se, assay = "data")
}
\arguments{
\item{se}{A \code{\link{SummarizedExperiment}} object}

\item{assay}{A name of the slot data is stored}
}
\value{
A \code{\link{SummarizedExperiment}}  object with *UMAP* projection
  stored in \code{metadata(se[[umap]]} slot
}
\description{
A function to perform *UMAP* projection of data
}
