% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{loss_function}
\alias{loss_function}
\title{Estimates loss value.}
\usage{
loss_function(
  pred,
  data.valid,
  loss.unit = c("individuals", "L2 units"),
  loss.fun = c("MSE", "MAE", "cross-entropy"),
  y,
  L2.unit
)
}
\arguments{
\item{pred}{Predictions of outcome. A numeric vector of outcome predictions.}

\item{data.valid}{Test data set. A tibble of data that was not used for
prediction.}

\item{loss.unit}{Loss function unit. A character-valued scalar indicating
whether performance loss should be evaluated at the level of individual
respondents (\code{individuals}) or geographic units (\code{L2 units}).
Default is \code{individuals}.}

\item{loss.fun}{Loss function. A character-valued scalar indicating whether
prediction loss should be measured by the mean squared error (\code{MSE})
or the mean absolute error (\code{MAE}). Default is \code{MSE}.}

\item{y}{Outcome variable. A character vector containing the column names of
the outcome variable.}

\item{L2.unit}{Geographic unit. A character scalar containing the column name
of the geographic unit in \code{survey} and \code{census} at which outcomes
should be aggregated.}
}
\value{
Returns a tibble with number of rows equal to the number of loss
  functions tested (defaults to 4 for cross-entropy, f1, MSE, and msfe). The
  number of columns is 2 where the first is called measure and contains the
  names of the loss-functions and the second is called value and contains the
  loss-function scores.
}
\description{
\code{loss_function()} estimates the loss based on a loss function.
}
