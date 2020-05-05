#' Apply support vector machine classifier to MrP.
#'
#' \code{run_svm} is a wrapper function that applies the support vector machine
#' classifier to data provided by the user, evaluates prediction performance,
#' and chooses the best-performing model.
#'
#' @param y Outcome variable. A character scalar containing the column name of
#'   the outcome variable in \code{survey}.
#' @param L1.x Individual-level covariates. A character vector containing the
#'   column names of the individual-level variables in \code{survey} and
#'   \code{census} used to predict outcome \code{y}. Note that geographic unit
#'   is specified in argument \code{L2.unit}.
#' @param L2.x Context-level covariates. A character vector containing the
#'   column names of the context-level variables in \code{survey} and
#'   \code{census} used to predict outcome \code{y}.
#' @param L2.unit Geographic unit. A character scalar containing the column
#'   name of the geographic unit in \code{survey} and \code{census} at which
#'   outcomes should be aggregated.
#' @param L2.reg Geographic region. A character scalar containing the column
#'   name of the geographic region in \code{survey} and \code{census} by which
#'   geographic units are grouped (\code{L2.unit} must be nested within
#'   \code{L2.reg}). Default is \code{NULL}.
#' @param loss.fun Loss function. A character-valued scalar indicating whether
#'   prediction loss should be measured by the mean squared error (\code{MSE})
#'   or the mean absolute error (\code{MAE}). Default is \code{MSE}.
#' @param kernel SVM kernel. A character-valued scalar specifying the kernel to
#'   be used by SVM. The possible values are \code{linear}, \code{polynomial},
#'   \code{radial}, and \code{sigmoid}. Default is \code{radial}.
#' @param loss.fun SVM loss function. If \code{NULL}, then SVM uses the
#'   misclassification error to measure the loss of categorical predictions and
#'   the mean squared error to measure the loss of numeric predictions. Default
#'   is \code{NULL}.
#' @param gamma SVM kernel parameter. A numeric vector whose values specify
#'   the gamma parameter in the SVM kernel. This parameter is needed for all
#'   kernel types except linear. Default is
#'   \eqn{c(0.3, 0.5, 0.55, 0.6, 0.65, 0.7, 0.8, 0.9, 1, 2, 3, 4)}.
#' @param cost SVM cost parameter. A numeric vector whose values specify the
#'   cost of constraints violation in SVM. Default is \eqn{c(1, 10)}.
#' @param data Data for cross-validation. A \code{list} of \eqn{k}
#'   \code{data.frames}, one for each fold to be used in \eqn{k}-fold
#'   cross-validation.
#' @param verbose Verbose output. A logical argument indicating whether or not
#'   verbose output should be printed. Default is \code{TRUE}.
#' @param cross The number cross-validation folds. A numeric scalar.
#' @return
#' @examples #not_yet

run_svm <- function(y, L1.x, L2.x, L2.unit, L2.reg,
                    kernel = "radial", loss.fun,
                    gamma, cost,
                    data, verbose) {

  # Create model formula
  x <- paste(c(L1.x, L2.x, L2.unit, L2.reg), collapse = " + ")
  form <- as.formula(paste(y, " ~ ", x, sep = ""))

  # Determine number of cross-validation folds provided
  k <- length(data)

  # Prepare data
  data <- dplyr::bind_rows(data) %>%
    dplyr::mutate_at(.vars = y, as.factor)

  # Train and evaluate model using the supplied set of tuning parameters
  models <- svm_classifier(method = "svm",
                           form = form,
                           data = data,
                           kernel = kernel,
                           error.fun = loss.fun,
                           probability = TRUE,
                           gamma = gamma,
                           cost = cost,
                           sampling = "cross",
                           cross = k,
                           verbose = TRUE)

  # Extract the best model
  best_model <- models$best.model

  # Extract tuning parameters of best model
  out <- list(gamma = best_model$gamma,
              cost = best_model$cost)

  # Function output
  return(out)
}
