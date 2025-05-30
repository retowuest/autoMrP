% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ebma.R
\name{ebma_mc_draws}
\alias{ebma_mc_draws}
\title{EBMA multicore tuning - parallelises over draws.}
\usage{
ebma_mc_draws(
  train.preds,
  train.y,
  ebma.fold,
  y,
  L1.x,
  L2.x,
  L2.unit,
  L2.reg,
  pc.names,
  model.bs,
  model.pca,
  model.lasso,
  model.gb,
  model.svm,
  model.mrp,
  tol,
  n.draws,
  cores,
  preds_all,
  post.strat,
  dv_type,
  deep.mrp
)
}
\arguments{
\item{train.preds}{Predictions of classifiers on the classifier training
data. A tibble.}

\item{train.y}{Outcome variable of the classifier training data. A numeric
vector.}

\item{ebma.fold}{New data for EBMA tuning. A list containing the the data
that must not have been used in classifier training.}

\item{y}{Outcome variable. A character vector containing the column names of
the outcome variable. A character scalar containing the column name of
the outcome variable in \code{survey}.}

\item{L1.x}{Individual-level covariates. A character vector containing the
column names of the individual-level variables in \code{survey} and
\code{census} used to predict outcome \code{y}. Note that geographic unit
is specified in argument \code{L2.unit}.}

\item{L2.x}{Context-level covariates. A character vector containing the
column names of the context-level variables in \code{survey} and
\code{census} used to predict outcome \code{y}. To exclude context-level
variables, set \code{L2.x = NULL}.}

\item{L2.unit}{Geographic unit. A character scalar containing the column
name of the geographic unit in \code{survey} and \code{census} at which
outcomes should be aggregated.}

\item{L2.reg}{Geographic region. A character scalar containing the column
name of the geographic region in \code{survey} and \code{census} by which
geographic units are grouped (\code{L2.unit} must be nested within
\code{L2.reg}). Default is \code{NULL}.}

\item{pc.names}{Principal Component Variable names. A character vector
containing the names of the context-level principal components variables.}

\item{model.bs}{The tuned model from the multilevel regression with best
subset selection classifier. An \code{\link[lme4]{glmer}} object.}

\item{model.pca}{The tuned model from the multilevel regression with
principal components as context-level predictors classifier. An
\code{\link[lme4]{glmer}} object.}

\item{model.lasso}{The tuned model from the multilevel regression with L1
regularization classifier. A \code{\link[glmmLasso]{glmmLasso}} object.}

\item{model.gb}{The tuned model from the gradient boosting classifier. A
\code{\link[gbm]{gbm}} object.}

\item{model.svm}{The tuned model from the support vector machine classifier.
An \code{\link[e1071]{svm}} object.}

\item{model.mrp}{The standard MrP model. An \code{\link[lme4]{glmer}} object}

\item{tol}{EBMA tolerance. A numeric vector containing the tolerance values
for improvements in the log-likelihood before the EM algorithm stops
optimization. Values should range at least from \eqn{0.01} to \eqn{0.001}.
Default is \code{c(0.01, 0.005, 0.001, 0.0005, 0.0001, 0.00005, 0.00001)}.
Passed on from \code{ebma.tol}.}

\item{n.draws}{EBMA number of samples. An integer-valued scalar specifying
the number of bootstrapped samples to be drawn from the EBMA fold and used
for tuning EBMA. Default is \eqn{100}. Passed on from \code{ebma.n.draws}.}

\item{cores}{The number of cores to be used. An integer indicating the number
of processor cores used for parallel computing. Default is 1.}

\item{post.strat}{Post-stratification results. A list containing the best
models for each of the tuned classifiers, the individual level predictions
on the data classifier trainig data and the post-stratified context-level
predictions.}

\item{dv_type}{The type of the depenedent variable. A character string.
Either "binary" or "linear".}

\item{deep.mrp}{Deep MRP classifier. A logical argument indicating whether
the deep MRP classifier should be used for best subset prediction. Setting
\code{deep.mrp = TRUE} will include all interactions of L1.x in the best
subset classifier. Default is \code{FALSE}.}
}
\value{
The classifier weights. A numeric vector.
}
\description{
\code{ebma_mc_draws} is called from within \code{ebma}. It tunes using
multiple cores.
}
