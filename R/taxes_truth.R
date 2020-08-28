#' Sample of tax rates item from the 2008 National Annenberg Election Studies.
#'
#' The 2008 National Annenberg Election Studies (NAES) item (CBb01) asked: "I'm
#' going to read you some options about federal income taxes. Please tell me
#' which one comes closest to your view on what we should be doing about federal
#' income taxes: (1) Cut taxes; (2) Keep taxes as they are; (3) Raise taxes if
#' necessary; (4) None of these; (998) Don't know; (999) No answer. Category (3)
#' was turned into a 'raise taxes response,' categories (1) and (2) were
#' combined into a 'do not raise taxes' response. The original item from the
#' phone and online surveys contains 50,483 respondents. This sample mimics a
#' typical national survey. It contains at least 5 respondents from each state
#' but is otherwise a random sample.
#'
#'
#' @format A data frame with 1500 rows and 13 variables:
#' \describe{
#'   \item{YES}{1 if individual supports raising taxes; 0 otherwise}
#'   \item{L1x1}{Age group (four categories: 1 = 18-29; 2 = 30-44; 3 = 45-64; 4 = 65+)}
#'   \item{L1x2}{Education level (four categories: 1 = < high school; 2 = high school graduate; 3 = some college; 4 = college graduate)}
#'   \item{L1x3}{Gender-race combination (six categories: 1 = white male; 2 = black male; 3 = hispanic male; 4 = white female; 5 = black female; 6 = hispanic female)}
#'   \item{state}{U.S. state}
#'   \item{L2.unit}{U.S. state id}
#'   \item{region}{U.S. region (four categories: 1 = Northeast; 2 = Midwest; 3 = South; 4 = West)}
#'   \item{L2.x1}{State-level share of votes for the Republican candidate in the previous presidential election}
#'   \item{L2.x2}{State-level percentage of Evangelical Protestant or Mormon respondents}
#'   \item{L2.x3}{State-level percentage of the population living in urban areas}
#'   \item{L2.x4}{State-level unemployment rate}
#'   \item{L2.x5}{State-level share of Hispanics}
#'   \item{L2.x6}{State-level share of Whites}
#' }
#' @usage data(taxes_survey)
#' @source The data set (excluding L2.x3, L2.x4, L2.x5, L2.x6) is taken from the
#'   article: Buttice, Matthew K, and Benjamin Highton. 2013. "How does
#'   multilevel regression and poststrat-stratification perform with
#'   conventional national surveys?" Political Analysis 21(4): 449-467. It is a
#'   random sample with at least 5 respondents per state. L2.x3, L2.x3, L2.x4,
#'   L2.x5 and L2.x6 are available at \url{https://www.census.gov}.
"taxes_survey"
