% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classifier_train.R
\name{cv_loop_train}
\alias{cv_loop_train}
\title{iterates training a model with CV (serial version)}
\usage{
cv_loop_train(
  data,
  cls,
  stratify = NA,
  fitControl,
  K = 25,
  resampling_rate = 0.8,
  n_features = NA,
  filter_method = c("ROC", "WILCOX"),
  filter_direction = c("two.sided", "greater", "less"),
  filter_threshold_diff = 1,
  filter_threshold_score = 0.8,
  observation_weights = NULL,
  feature_weights = c("uniform", "weighted")
)
}
\arguments{
\item{data}{input matrix, of dimension \code{nobs x nvars}; each row is an
observation vector. Since this is an input to \code{\link{glmnet}}, it
should be the format that can be used with \code{\link{glmnet}}}

\item{cls}{class labels}

\item{K}{...}

\item{resampling_rate}{...}

\item{n_features}{...}

\item{filter_method}{...}

\item{filter_direction}{...}

\item{filter_threshold_diff}{...}

\item{filter_threshold_score}{...}

\item{observation_weights}{...}

\item{feature_weights}{...}
}
\value{
a list of trained models -- see \code{\link{cv_loop_train_iter}}
}
\description{
iterates training a model with CV (serial version)
}
