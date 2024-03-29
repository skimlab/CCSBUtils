% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classifier_train.R
\name{train_cv_model_glmnet}
\alias{train_cv_model_glmnet}
\title{A function train a model with CV}
\usage{
train_cv_model_glmnet(
  data_train,
  cls_train,
  fitControl,
  observation_weights = NULL,
  penalty.factor = 1
)
}
\arguments{
\item{data_train}{input matrix, of dimension nobs x nvars; each row is an
observation vector. Since this is an input to \code{glmnet}, it should be
the format that can be used with \code{glmnet}}

\item{cls_train}{class labels}

\item{fitControl}{A list of training parameters.  See
\code{\link{caret::trainControl}} for detail}

\item{penalty.factor}{Separate penalty factors can be applied to each
coefficient.  See \code{\link{glmnet::glmnet}} for detail}
}
\value{
A list returned from \code{caret::train}
}
\description{
A function train a model with CV
}
