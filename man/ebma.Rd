% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ebma.R
\name{ebma}
\alias{ebma}
\title{Bayesian Ensemble Model Averaging EBMA}
\usage{
ebma(
  ebma.fold,
  y,
  L1.x,
  L2.x,
  L2.unit,
  L2.reg,
  pc.names,
  post.strat,
  n.draws,
  tol,
  best.subset.opt,
  pca.opt,
  lasso.opt,
  gb.opt,
  svm.opt,
  deep.mrp,
  verbose,
  cores,
  preds_all
)
}
\arguments{
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

\item{post.strat}{Post-stratification results. A list containing the best
models for each of the tuned classifiers, the individual level predictions
on the data classifier trainig data and the post-stratified context-level
predictions.}

\item{n.draws}{EBMA number of samples. An integer-valued scalar specifying
the number of bootstrapped samples to be drawn from the EBMA fold and used
for tuning EBMA. Default is \eqn{100}. Passed on from \code{ebma.n.draws}.}

\item{tol}{EBMA tolerance. A numeric vector containing the tolerance values
for improvements in the log-likelihood before the EM algorithm stops
optimization. Values should range at least from \eqn{0.01} to \eqn{0.001}.
Default is \code{c(0.01, 0.005, 0.001, 0.0005, 0.0001, 0.00005, 0.00001)}.
Passed on from \code{ebma.tol}.}

\item{best.subset.opt}{Tuned best subset parameters. A list returned from
\code{run_best_subset()}.}

\item{pca.opt}{Tuned best subset with principal components parameters. A list
returned from \code{run_pca()}.}

\item{lasso.opt}{Tuned lasso parameters. A list returned from
\code{run_lasso()}.}

\item{gb.opt}{Tuned gradient tree boosting parameters. A list returned from
\code{run_gb()}.}

\item{svm.opt}{Tuned support vector machine parameters. A list returned from
\code{run_svm()}.}

\item{deep.mrp}{Deep MRP classifier. A logical argument indicating whether
the deep MRP classifier should be used for best subset prediction. Setting
\code{deep.mrp = TRUE} will include all interactions of L1.x in the best
subset classifier. Default is \code{FALSE}.}

\item{verbose}{Verbose output. A logical argument indicating whether or not
verbose output should be printed. Default is \code{FALSE}.}

\item{cores}{The number of cores to be used. An integer indicating the number
of processor cores used for parallel computing. Default is 1.}
}
\description{
\code{ebma} tunes EBMA and generates weights for classifier averaging.
}
