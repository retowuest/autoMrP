#' A sample from the 2008 Cooperative Congressional Election Studies (item cc418_1) with additional state-level variables.
#'
#' @docType data
#'
#' @usage data(survey)
#'
#' The survey question asked: "Would you approve of the use of U.S. military troops in order to ensure the supply of oil?" This 2008 CCES item cc418_1 contains 36,832 respondents. This sample mimics a typical national survey. It contains at least 5 respondents from each state but is otherwise a random sample.
#'
#' @format A data frame with 1500 rows and 13 variables:
#' \describe{
#'   \item{YES}{1 if individual supports use of troops; 0 otherwise}
#'   \item{L1x1}{Age group (four categories)}
#'   \item{L1x2}{Education level (four categories)}
#'   \item{L1x3}{Gender-race combination (six categories)}
#'   \item{state}{The state that the respondent is from}
#'   \item{L2.unit}{An aditional  variable indicating the state that the respondent is from}
#'   \item{region}{Region in the U.S. that the respondent is from; 1 = something, 2 = something, 3 = something, 4 = something}
#'   \item{L2.x1}{The share of votes for the Republican candidate in the previous presidential election}
#'   \item{L2.x2}{The percentage of Evangelical Protestant or Mormon respondents}
#'   \item{L2.x3}{The percentage of the population living in urban areas}
#'   \item{L2.x4}{The unemployment rate}
#'   \item{L2.x5}{The share of Hispanics}
#'   \item{L2.x6}{The share of Whites}
#'   ...
#' }
#' @source \url{https://cces.gov.harvard.edu/pages/welcome-cooperative-congressional-election-study}
"survey"
