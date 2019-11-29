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
#'   variable. Must include geographic unit.
#' @param L2.x Context-level covariates. A character vector of column names
#'   corresponding to the context-level variables used to predict the outcome
#'   variable.
#' @param L2.unit Geographic unit. A character scalar indicating the column
#'   name of the geographic unit at which outcomes should be aggregated.
#' @param survey Survey data. A data.frame containing the y and x column names.
#' @param census Census data. A data.frame containing the x column names.
#' @param bin.size Bin size for ideal types. A character vector indicating the
#'   column name of the variable in census containing the bin size for ideal
#'   types in a geographic unit.
#' @param ebma.size Size of EBMA hold-out fold. A rational number in the open
#'   unit interval indicating the share of respondents to be contained in the
#'   EBMA hold-out fold. If left unspecified (NULL), then ebma.size is set to
#'   1/4 of the survey sample size.
#' @param k.folds Number of folds. An integer-valued scalar indicating the
#'   number of folds to be used for cross-validation. Defaults to the value of 5.
#' @param cv.sampling Sampling method. A character-valued scalar indicating
#'   whether sampling in the creation of cross-validation folds should be done
#'   by respondents or geographic units. Default is by units.
#' @param seed Seed. An integer-valued scalar to control random number
#'   generation. If left unspecified (NULL), then seed is set to 12345.
#' @return
#' @keywords MRP multilevel regression post-stratification machine learning
#'   EBMA ensemble bayesian model averaging
#' @examples
#' @export

auto_MrP <- function(y, L1.x, L2.x, L2.unit, survey, census,
                     bin.size = "None", ebma.size = NULL, k.folds = 5,
                     cv.sampling = "units", seed = NULL) {
  # Set seed
  if (is.null(seed)) {
    set.seed(12345)
  } else {
    set.seed(seed)
  }

  # Error and warning checks
  if (!all(L1.x %in% colnames(survey))) {
    stop(paste("Individual-level variable(s) '",
               L1.x[which(!(L1.x %in% colnames(survey)))],
               "' is/are not in your survey data.", sep = ""))
  }

  if (!all(L1.x %in% colnames(census))) {
    stop(paste("Individual-level variable(s) '",
               L1.x[which(!(L1.x %in% colnames(census)))],
               "' is/are not in your census data.", sep = ""))
  }

  if (!all(L2.x %in% colnames(survey))) {
    stop(paste("Context-level variable(s) '",
               L2.x[which(!(L2.x %in% colnames(survey)))],
               "' is/are not in your survey data.", sep = ""))
  }

  if (!all(L2.x %in% colnames(census))) {
    stop(paste("Context-level variable(s) '",
               L2.x[which(!(L2.x %in% colnames(census)))],
               "' is/are not in your census data.", sep = ""))
  }

  if (!(y %in% colnames(survey))) {
    stop(paste("Outcome '", y,
               "' is not in your survey data.", sep = ""))
  }

  if (!(L2.unit %in% colnames(survey))) {
    stop(paste("The geographic unit '", L2.unit,
               "' is not in your survey data.", sep = ""))
  }

  if (!(L2.unit %in% colnames(census))) {
    stop(paste("The geographic unit '", L2.unit,
               "' is not in your census data.", sep = ""))
  }

  if (is.null(ebma.size)) {
    ebma.size <- round(nrow(survey) / 4, digits = 0)
  } else if (is.numeric(ebma.size) & ebma.size > 0 & ebma.size < 1) {
    ebma.size <- round(nrow(survey) * ebma.size, digits = 0)
  } else {
    stop("ebma.size must be a rational number in the open unit interval.")
  }

  if (!(is.numeric(k.folds) & k.folds == round(k.folds, digits = 0))) {
    stop("k.folds must be an integer number.")
  }

  if (!cv.sampling %in% c("respondents", "units")) {
    stop("cv.sampling must take either the value 'respondents' or 'units'.")
  }

  # ------------------------------- Prepare data -------------------------------

  # If not provided in census data, calculate bin size for each ideal type
  if (bin.size == "None") {
    census <- census %>%
      dplyr::group_by(.dots = L1.x) %>%
      dplyr::summarise(n = dplyr::n())
  } else {
    census$n <- census[[bin.size]]
  }

  # Scale context-level variables in survey and census data
  survey[, L2.x] <- scale(survey[, L2.x], center = TRUE, scale = TRUE)
  census[, L2.x] <- scale(census[, L2.x], center = TRUE, scale = TRUE)

  # ------------------------------- Create folds -------------------------------

  # EBMA hold-out fold
  ebma_folding_out <- ebma_folding(data = survey,
                                   L2.unit = L2.unit,
                                   ebma.size = ebma.size)

  ebma_fold <- ebma_folding_out$ebma_fold
  cv_data <- ebma_folding_out$cv_data

  # K folds for cross-validation
  cv_folds <- cv_folding(data = cv_data,
                         L2.unit = L2.unit,
                         k.folds = k.folds,
                         cv.sampling = cv.sampling)

  # ------------------------ Run individual classifiers ------------------------

  # Classifier 1: best subset
  best_subset_out <- best_subset(train.data = cv_folds,
                                 verbose = TRUE)

  # Classifier 2: lasso
  lasso_out <- lasso(train.data = cv_folds,
                     verbose = TRUE)

  # Classifier 3: PCA

  # Classifier 4: GB

  # Classifier 5: SVM

}
