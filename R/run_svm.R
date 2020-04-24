#' Apply support vector machine classifier to MrP.
#'
#' \code{svm} is a wrapper function that applies the support vector machine
#' classifier to data provided by the user, evaluates prediction performance,
#' and chooses the best-performing model.
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
#' @param L2.reg Geographic region. A character scalar indicating the column
#'   name of the geographic region by which geographic units are grouped
#'   (L2.unit must be nested within L2.re).
#' @param loss.unit Loss function unit. A character-valued scalar indicating
#'   whether the loss should be evaluated at the level of individual respondents
#'   or the level of geographic units. Default is at the individual level.
#' @param loss.measure Loss function measure. A character-valued scalar
#'   indicating whether the loss should be measured by the mean squared error
#'   or the mean absolute error. Default is the MSE.
#' @param data Data for cross-validation. A list of k data.frames, one for
#'   each fold used in k-fold cross-validation.
#' @param kernel Kernel for SVM. A character string specifying the kernel to
#'   be used for SVM. The possible types are linear, polynomial, radial, and
#'   sigmoid. Default is radial.
#' @param error.fun
#' @param gamma.set Gamma parameter for SVM. This parameter is needed for all
#'   kernels except linear.
#' @param cost.set Cost parameter for SVM. This parameter specifies the cost of
#'   constraints violation.
#' @param k.folds Number of folds. An integer-valued scalar indicating the
#'   number of folds to be used for cross-validation. Defaults to the value of 5.
#' @param verbose Verbose output. A logical vector indicating whether or not
#'   verbose output should be printed.
#' @return
#' @examples

svm <- function(y, L1.x, L2.x, L2.unit, L2.reg,
                kernel = "radial", error.fun,
                gamma.set, cost.set, k.folds,
                data, verbose) {
  # Create model formula
  x <- paste(c(L1.x, L2.x, L2.unit, L2.reg), collapse = " + ")
  form <- as.formula(paste(y, " ~ ", x, sep = ""))

  # Prepare data
  data <- dplyr::bind_rows(data) %>%
    dplyr::mutate_at(.vars = y, as.factor)

  # Train and evaluate model using the supplied set of tuning parameters
  models <- svm_classifier(method = "svm",
                           form = form,
                           data = data,
                           kernel = kernel,
                           error.fun = error.fun,
                           probability = TRUE,
                           gamma.set = gamma.set,
                           cost.set = cost.set,
                           sampling = "cross",
                           cross = k.folds,
                           verbose = TRUE)

  # Extract the best model
  best_model <- models$best.model

  # Extract tuning parameters of best model
  out <- list(gamma = best_model$gamma,
              cost = best_model$cost)

  # Function output
  return(out)
}
