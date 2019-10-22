#' Improve MrP through ensemble learning.
#'
#' \code{auto_MrP} improves the prediction performance of multilevel regression
#' with post-stratification (MrP) by combining a number of machine learning
#' methods through ensemble bayesian model averaging (EBMA).
#'
#' @param y Outcome variable. A character scalar containing the column name of
#'   the outcome variable.
#' @param L1.x Individual-level covariates. A character vector of column names
#'   corresponding to the individual-level variables used to predict the outcome
#'   variable.
#' @param L2.x Context-level covariates. A character vector of column names
#'   corresponding to the context-level variables used to predict the outcome
#'   variable.
#' @param survey Survey data. A data.frame containing the y and x column names.
#' @param census Census data. A data.frame containing the x column names.
#' @param geo.unit Geographic unit. A character scalar indicating the column
#'   name of the geographic unit at which outcomes should be aggregated.
#' @param n Bin size for ideal types. A character vector indicating the column
#'   name of the variable containing the bin size for ideal types in a geographic
#'   unit.
#' @return
#' @keywords MRP multilevel regression post-stratification machine learning
#'   EBMA ensemble bayesian model averaging
#' @examples
#' @export

auto_MrP <- function(y, L1.x, L2.x, survey, census, geo.unit,
                     proportion = "None", set.seed = NULL) {
  # Set seed
  set.seed(set.seed)

  # Error and warning checks
  if (!all(L1.x %in% colnames(survey))) {
    stop(paste("Individual-level variable(s) '",
               L1.x[which(!(L1.x %in% colnames(survey)))],
               "' is/are not in your survey data.", sep = ""))
  }
  if(!all(L1.x %in% colnames(census))) {
    stop(paste("Individual-level variable(s) '",
               L1.x[which(!(L1.x %in% colnames(census)))],
               "' is/are not in your census data.", sep = ""))
  }
  if (!all(L2.x %in% colnames(survey))) {
    stop(paste("Individual-level variable(s) '",
               L2.x[which(!(L2.x %in% colnames(survey)))],
               "' is/are not in your survey data.", sep = ""))
  }
  if(!all(L2.x %in% colnames(census))) {
    stop(paste("Individual-level variable(s) '",
               L2.x[which(!(L2.x %in% colnames(census)))],
               "' is/are not in your census data.", sep = ""))
  }
  if(!(y %in% colnames(survey))) {
    stop(paste("Outcome '", y,
               "' is not in your survey data.", sep = ""))
  }
  if(!(geo.unit %in% colnames(survey))) {
    stop(paste("The geographic unit '", geo.unit,
               "' is not in your survey data.", sep = ""))
  }
  if(!(geo.unit %in% colnames(census))) {
    stop(paste("The geographic unit '", geo.unit,
               "' is not in your census data.", sep = ""))
  }

  # Calculate bin size for ideal types if not provided in census data
  if(n == "None") {
    census <- census %>%
      dplyr::group_by(.dots = L1.x) %>%
      dplyr::summarise(n = dplyr::n())
  } else {
    census$n <- census[[n]]
  }

  # Prepare data ---------------------------------------------------------------

  # Create folds ---------------------------------------------------------------

  # Run individual classifiers -------------------------------------------------

  # Classifier: best subset
  best_subset_out <- best_subset(train.data = )

  # Classifier: lasso

}
