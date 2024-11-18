#' KNN classifier
#'
#' \code{knn_classifier} applies k-nearest neighbors classification to a data
#' set.
#'
#' @inheritParams auto_MrP
#' @param form Model formula. A two-sided linear formula describing
#'   the model to be fit, with the outcome on the LHS and the covariates
#'   separated by + operators on the RHS.
#' @param data.train Training data. A \code{data.frame} containing the
#'   cross-validation data used to train the model.
#' @param data.valid Validation data. A \code{data.frame} containing the
#'   cross-validation data used to evaluate the model.
#' @param knn.k.value KNN number of neighbors. A positive integer-valued scalar
#'   specifying the number of neighbors to be considered in the KNN model.
#' @param knn.kernel KNN kernel. A character-valued scalar specifying the kernel
#'   to be used in the KNN model. The possible values are \code{rectangular}
#'   (which is standard unweighted KNN), \code{triangular}, \code{epanechnikov}
#'   (or beta(2,2)), \code{biweight} (or beta(3,3)), \code{triweight} (or
#'   beta(4,4)), \code{cos}, \code{inv}, \code{gaussian}, and \code{optimal}.
#' @param verbose Verbose output. A logical vector indicating whether or not
#'   verbose output should be printed.
#' @return The k-nearest neighbors model. A \code{\link[kknn]{kknn}} object.

knn_classifier <- function(
  y, form, data.train, data.valid, knn.k.value,
  knn.kernel = c(
    "rectangular", "triangular", "epanechnikov", "biweight", "triweight", "cos",
    "inv", "gaussian", "optimal"
  ),
  verbose = c(TRUE, FALSE)
) {

  # Train and evaluate model using the supplied value of k
  if (isTRUE(verbose)) {
    out <- kknn::kknn(
      formula = form,
      train = data.train,
      test = data.valid,
      k = knn.k.value,
      kernel = knn.kernel,
      scale = TRUE
    )
  } else {
    out <- suppressMessages(suppressWarnings(
      kknn::kknn(
        formula = form,
        train = data.train,
        test = data.valid,
        k = knn.k.value,
        kernel = knn.kernel,
        scale = TRUE
      )
    ))
  }

  # Function output
  return(out)
}
