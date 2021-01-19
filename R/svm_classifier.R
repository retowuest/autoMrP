#' SVM classifier
#'
#' \code{svm_classifier} applies support vector machine classification to a
#' data set.
#'
#' @param form Model formula. A two-sided linear formula describing
#'   the model to be fit, with the outcome on the LHS and the covariates
#'   separated by + operators on the RHS.
#' @param data Data. A data.frame containing the cross-validation data used to
#'   train and evaluate the model.
#' @param kernel Kernel for SVM. A character string specifying the kernel to
#'   be used for SVM. The possible types are linear, polynomial, radial, and
#'   sigmoid. Default is radial.
#' @param type svm can be used as a classification machine, as a regression machine, or for novelty detection. Depending of whether y is a factor or not, the default setting for type is C-classification or eps-regression, respectively, but may be overwritten by setting an explicit value. Valid options are: #' \enumerate{
#'   \item C-classification
#'   \item nu-classification
#'   \item one-classification (for novelty detection)
#'   \item eps-regression
#'   \item nu-regression
#' }
#' @param probability Probability predictions. A logical argument indicating
#'   whether the model should allow for probability predictions
#' @param svm.gamma Gamma parameter for SVM. This parameter is needed for all
#'   kernels except linear.
#' @param svm.cost Cost parameter for SVM. This parameter specifies the cost of
#'   constraints violation.
#' @param verbose Verbose output. A logical vector indicating whether or not
#'   verbose output should be printed.
#' @return The support vector machine model. An \code{\link[e1071]{svm}} object.

svm_classifier <- function(form, data, kernel, type, probability, svm.gamma,
                           svm.cost, verbose = c(TRUE, FALSE)) {

  # Train and evaluate model using the supplied set of tuning parameters
  if (isTRUE(verbose == TRUE)) {
    out <- e1071::svm(
      formula = form,
      data = data,
      type = type,
      kernel = kernel,
      gamma = svm.gamma,
      cost = svm.cost,
      probability = probability
    )
  } else {
    out <- suppressMessages(suppressWarnings(
      e1071::svm(
        formula = form,
        data = data,
        type = type,
        kernel = kernel,
        gamma = svm.gamma,
        cost = svm.cost,
        probability = probability
      )
    ))
  }

  # Function output
  return(out)
}
