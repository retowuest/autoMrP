#' SVM classifier
#'
#' \code{svm_classifier} applies support vector machine classification to a
#' data set.
#'
#' @param method Function. A character string specifying the name of the
#'   function to be tuned.
#' @param form Model formula. A two-sided linear formula describing
#'   the model to be fit, with the outcome on the LHS and the covariates
#'   separated by + operators on the RHS.
#' @param data Data. A data.frame containing the cross-validation data used to
#'   train and evaluate the model.
#' @param kernel Kernel for SVM. A character string specifying the kernel to
#'   be used for SVM. The possible types are linear, polynomial, radial, and
#'   sigmoid. Default is radial.
#' @param error.fun function returning the error measure to be minimized. It
#' takes two arguments: a vector of true values and a vector of predicted values.
#' If NULL, the misclassification error is used for categorical predictions and
#' the mean squared error for numeric predictions.
#' @param probability Probability predictions. A logical argument indicating
#'   whether the model should allow for probability predictions
#' @param svm.gamma Gamma parameter for SVM. This parameter is needed for all
#'   kernels except linear.
#' @param svm.cost Cost parameter for SVM. This parameter specifies the cost of
#'   constraints violation.
#' @param sampling Sampling scheme. A character string specifying the sampling
#'   scheme to be used. Possible values are cross, fix, and boot. Default is
#'   cross, which performs cross-validation.
#' @param cross The number cross-validation folds. A numeric scalar.
#' @param verbose Verbose output. A logical vector indicating whether or not
#'   verbose output should be printed.
#' @return The support vector machine model. An \code{\link[e1071]{tune}} object.
#' @examples \dontrun{
#' # Prepare data
#' survey_item <- dplyr::bind_rows(survey_item) %>%
#'   dplyr::mutate_at(.vars = y, as.factor)
#'
#' # run svm classifier
#' models <- svm_classifier(
#'   method = "svm",
#'   form = YES ~ L1x1 + L1x2 + state + region + L2.x1,
#'   data = survey_item,
#'   kernel = "radial",
#'   error.fun = "MSE",
#'   probability = TRUE,
#'   gamma = c(0.3, 0.5, 0.55, 0.6, 0.65, 0.7, 0.8, 0.9, 1, 2, 3, 4),
#'   cost = c(1, 10),
#'   sampling = "cross",
#'   cross = 5,
#'   verbose = TRUE)
#' }

svm_classifier <- function(method, form, data, kernel,
                           error.fun, probability,
                           svm.gamma, svm.cost,
                           sampling = "cross", cross,
                           verbose = c(TRUE, FALSE)) {
  # Train and evaluate model using the supplied set of tuning parameters
  if (isTRUE(verbose == TRUE)) {
    out <- e1071::tune(method = method,
                       train.x = form,
                       data = data,
                       kernel = kernel,
                       error.fun = error.fun,
                       probability = probability,
                       ranges = list(gamma = svm.gamma,
                                     cost = svm.cost),
                       tunecontrol = e1071::tune.control(sampling = sampling,
                                                         cross = cross))
  } else {
    out <- suppressMessages(suppressWarnings(
        e1071::tune(method = method,
                    train.x = form,
                    data = data,
                    kernel = kernel,
                    error.fun = error.fun,
                    probability = probability,
                    ranges = list(gamma = svm.gamma,
                                  cost = svm.cost),
                    tunecontrol = e1071::tune.control(sampling = sampling,
                                                      cross = cross))
      ))
  }

  # Function output
  return(out)
}
