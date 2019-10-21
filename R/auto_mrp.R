#' Improve MrP through ensemble learning.
#'
#' \code{auto_MrP} improves the prediction performance of multilevel regression
#' with post-stratification (MrP) by combining a number of machine learning
#' methods through ensemble bayesian model averaging (EBMA).
#'
#' @param y Outcome variable. A character scalar containing the column name of
#'   the outcome variable.
#' @param x Predictor variables. A character vector of column names corresponding
#'   to the variables used to predict the outcome variable.
#' @param survey Survey data. A data.frame containing the y and x column names.
#' @param census Census data. A data.frame containing the x column names.
#' @param geo.unit Geographic unit. A character scalar indicating the column
#'   name of the geographic unit at which outcomes should be aggregated.
#' @param proportion Proportion of ideal types. A character vector indicating
#'   the column name of the variable containing the proportions of the ideal
#'   types in a geographic unit.
#' @return
#' @keywords MRP multilevel regression post-stratification machine learning
#'   EBMA ensemble bayesian model averaging
#' @examples
#' @export

auto_MrP <- function(y, x, survey, census, geo.unit,
                     proportion = "None", set.seed = NULL) {
  # Set seed
  set.seed(set.seed)

  # Error and warning checks
  if (!all(x %in% colnames(survey))) {
    stop(paste("Variable '", x[which(!(x %in% colnames(survey)))],
               "' is not in your survey data.", sep = ""))
  }
  if(!all(x %in% colnames(census))) {
    stop(paste("Variable '", x[which(!(x %in% colnames(census)))],
               "' is not in your census data.", sep = ""))
  }
  if(!(y %in% colnames(survey))) {
    stop(paste("Outcome '", y,
               "' is not in your survey data.", sep = ""))
  }
  if(!(geo.unit %in% x)) {
    stop(paste("The geographic unit '", geo.unit,
               "' is not among your predictor variables x.", sep = ""))
  }

  # Calculate proportions of ideal types if not provided
  if(proportion == "None") {
    census <- census %>%
      group_by_(.dots = x) %>%
      summarise(n = n())
  } else {
    census$n <- census[[proportion]]
  }


}
