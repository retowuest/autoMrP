#' Apply KNN classifier to MrP.
#'
#' \code{run_knn} is a wrapper function that applies the KNN classifier to data
#' provided by the user, evaluates prediction performance, and chooses the
#' best-performing model.
#'
#' @inheritParams auto_MrP
#' @param knn.k KNN number of neighbors. An integer-valued positive scalar
#'   specifying the number of neighbors to be considered in the KNN model.
#'   Default is \eqn{7}.
#' @param data Data for cross-validation. A \code{list} of \eqn{k}
#'   \code{data.frames}, one for each fold to be used in \eqn{k}-fold
#'   cross-validation.
#'
#' @return The tuned \code{knn.k} parameter. An integer scalar.

run_knn <- function(
    y, L1.x, L2.x, L2.unit, L2.reg, loss.unit, loss.fun,
    knn.k, data, verbose, cores
) {

}
