% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classifier_train.R
\name{filter_features_AUC}
\alias{filter_features_AUC}
\title{A function for filtering-based feature selection, using ROC/AUC as a metric}
\usage{
filter_features_AUC(
  dd,
  cls,
  threshold = 0.8,
  direction = c("two.sided", "greater", "less")
)
}
\arguments{
\item{dd}{A \code{matrix} or \code{data.frame}
where each row is an observation}

\item{cls}{A set of classes, either in \code{factor}}

\item{threshold}{A minimum value of AUC for filtering}

\item{direction}{Direction of AUC statistical test: "two.sided", "greater" or "less"}
}
\value{
A \code{data.frame} of \code{feature}, \code{diff}, \code{score}
        where \code{score} > \code{threshold}
}
\description{
A function for filtering-based feature selection, using ROC/AUC as a metric
}
