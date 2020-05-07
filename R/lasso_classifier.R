#' Lasso classifier
#'
#' \code{lasso_classifier} applies lasso classification to a data set.
#'
#' @param L2.fix Fixed effects. A two-sided linear formula describing
#'   the fixed effects part of the model, with the outcome on the LHS and
#'   the fixed effects separated by + operators on the RHS.
#' @param L1.re Random effects. A named list object, with the random effects
#'   providing the names of the list elements and ~ 1 being the list elements.
#' @param data.train Training data. A data.frame containing the training data
#'   used to train the model.
#' @param lambda Tuning parameter. Lambda is the penalty parameter that controls
#'   the shrinkage of fixed effects.
#' @param model.family Model family. A variable indicating the model family
#'   to be used by glmmLasso. Defaults to binomial(link = "probit").
#' @param verbose Verbose output. A logical vector indicating whether or not
#'   verbose output should be printed.
#' @return A multilevel lasso model. An \code{\link[glmmLasso]{glmmLasso}} object.
#' @examples \dontrun{
#' m <- lasso_classifier(
#'   L2.fix = YES ~ L2.x1 + L2.x2,
#'   L1.re = list(L1x1 = ~1, L1x2 = ~1, state = ~1, region = ~1),
#'   data.train = survey_item,
#'   lambda = 5,
#'   model.family = binomial(link = "probit"),
#'   verbose = TRUE)
#' }

lasso_classifier <- function(L2.fix, L1.re, data.train,
                             lambda, model.family,
                             verbose = c(TRUE, FALSE)) {
  # Train model on training data with lambda as tuning parameter
  if (isTRUE(verbose == TRUE)) {
    out <- glmmLasso::glmmLasso(fix = L2.fix, rnd = L1.re,
                                data = data.train, lambda = lambda,
                                family = model.family,
                                switch.NR = FALSE, final.re = TRUE,
                                control = list(center = TRUE,
                                               standardize = TRUE))
  } else {
    out <- suppressMessages(suppressWarnings(
      glmmLasso::glmmLasso(fix = L2.fix, rnd = L1.re,
                           data = data.train, lambda = lambda,
                           family = model.family,
                           switch.NR = FALSE, final.re = TRUE,
                           control = list(center = TRUE,
                                          standardize = TRUE))
    ))
  }

  # Function output
  return(out)
}
