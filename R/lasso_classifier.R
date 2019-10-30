#' Lasso classifier
#'
#' \text{lasso_classifier} applies lasso classification to a data set.
#'
#' @param L2.fix Fixed effects. A two-sided linear formula describing
#'   the fixed effects part of the model, with the outcome on the LHS and
#'   the fixed effects separated by + operators on the RHS.
#' @param L1.re Random effects. A named list object, with the random effects
#'   providing the names of the list elements and ~ 1 being the list elements.
#' @param data.train Training data. A data.frame containing the training data
#'   used to train the model.
#' @param lambda Tuning parameter. Lambda is the penalty parameter that controls
#'   the shrinkage of fixed effects. Either a numeric vector of lambda values
#'   or a data.frame with two columns, the first containing the size by which
#'   lambda should increase and the second the upper threshold of lambda until
#'   which the step size applies.
#' @param model.family Model family. A variable indicating the model family
#'   to be used by glmmLasso. Defaults to binomial(link = "probit").
#' @param verbose Verbose output. A logical vector indicating whether or not
#'   verbose output should be printed.
#' @return
#' @examples

lasso_classifier <- function(L2.fix, L1.re, data.train, lambda,
                             model.family = binomial(link = "probit"),
                             verbose = c(TRUE, FALSE)) {
  # Train model on training data with lambda as tuning parameter
  if (verbose == TRUE) {
    out <- glmmLasso::glmmLasso(fix = L2.fix, rnd = L1.re,
                                data = data.train, lambda = lambda,
                                family = model.family,
                                switch.NR = FALSE, final.re = TRUE,
                                control = list(standarize = FALSE))
  } else {

  }

  # Function output
  return(out)
}
