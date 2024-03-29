% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classifier_model_summary.R
\name{predict.type.multiclass}
\alias{predict.type.multiclass}
\title{apply a list of models to a set of samples}
\usage{
\method{predict}{type.multiclass}(models, data, final.class = NA, min.prob = 0.25)
}
\arguments{
\item{models}{a list of models each which is an outcome of
\code{\link{cv_loop_train_iter}} or \code{\link{cv_final_train}}}

\item{data}{input matrix, of dimension \code{nobs x nvars}; each row is an
observation vector. Since this is an input to \code{\link{glmnet}}, it
should be the format that can be used with \code{\link{glmnet}}}

\item{final.class}{the name (column heading) of final class outcomes,
majority among multiple classes}

\item{min.prob}{minimum probability for a sample to be classified to a class
(default: 0.25)}
}
\value{
a table where each rows is predicted outcomes of applying predictive
  models to a sample
}
\description{
returns a table where each rows is predicted outcomes of a sample, by
applying predictive models
}
