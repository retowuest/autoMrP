% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/census_data.R
\docType{data}
\name{census}
\alias{census}
\title{Quasi census data.}
\format{
A data frame with 2934 rows and 13 variables:
\describe{
  \item{state}{U.S. state}
  \item{L2.unit}{U.S. state id}
  \item{region}{U.S. region (four categories: 1 = Northeast; 2 = Midwest; 3 = South; 4 = West)}
  \item{L1x1}{Age group (four categories)}
  \item{L1x2}{Education level (four categories)}
  \item{L1x3}{Gender-race combination (six categories)}
  \item{proportion}{State-level proportion of respondents of that ideal type in the population}
  \item{L2.x1}{State-level share of votes for the Republican candidate in the previous presidential election}
  \item{L2.x2}{State-level percentage of Evangelical Protestant or Mormon respondents}
  \item{L2.x3}{State-level percentage of the population living in urban areas}
  \item{L2.x4}{State-level unemployment rate}
  \item{L2.x5}{State-level share of Hispanics}
  \item{L2.x6}{State-level share of Whites}
}
}
\source{
The data set (excluding L2.x3, L2.x4, L2.x5, L2.x6) is taken from the
  article: Buttice, Matthew K, and Benjamin Highton. 2013. "How does
  multilevel regression and poststrat-stratification perform with
  conventional national surveys?" Political Analysis 21(4): 449-467. L2.x3,
  L2.x3, L2.x4, L2.x5 and L2.x6 are available at
  \url{https://www.census.gov}.
}
\usage{
census
}
\description{
The census file is generated from the full 2008 Cooperative Congressional Election Studies
item cc418_1 by dissaggregating the 64 ideal type combinations of the individual level variables
L1x1, L2x2 and L1x3. A row is an ideal type in a given state.
}
\keyword{datasets}
