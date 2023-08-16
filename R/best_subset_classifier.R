#' Best subset classifier
#'
#' \code{best_subset_classifier} applies best subset classification to a data
#' set.
#'
#' @param model Multilevel model. A model formula describing the multilevel
#'   model to be estimated on the basis of the provided training data.
#' @param data.train Training data. A data.frame containing the training data
#'   used to train the model.
#' @param model.family Model family. A variable indicating the model family
#'   to be used by glmer. Defaults to binomial(link = "probit").
#' @param model.optimizer Optimization method. A character-valued scalar
#'   describing the optimization method to be used by glmer. Defaults to
#'   "bobyqa".
#' @param n.iter Iterations. A integer-valued scalar specifying the maximum
#'   number of function evaluations tried by the optimization method.
#' @param verbose Verbose output. A logical vector indicating whether or not
#'   verbose output should be printed.
#' @return The multilevel model. An \code{\link[lme4]{glmer}} object.

best_subset_classifier <- function(model, data.train,
                                   model.family, model.optimizer,
                                   n.iter, verbose = c(TRUE, FALSE)) {
  # Train model on training data
  if (isTRUE(verbose == TRUE)) {
    # optimizer
    if (model.optimizer == 'bobyqa'){
      out <- lme4::glmer(formula =  model,
                         data = data.train,
                         family = model.family,
                         lme4::glmerControl(
                           optimizer = model.optimizer,
                           optCtrl = list(maxfun = n.iter)))
    } else if (model.optimizer == 'nloptwrap') {
      out <- lme4::glmer(formula =  model,
                         data = data.train,
                         family = model.family,
                         lme4::glmerControl(
                           calc.derivs = FALSE,
                           optimizer = model.optimizer,
                           optCtrl = list(
                             method = "NLOPT_LN_NELDERMEAD",
                             starttests = TRUE, kkt = TRUE)))
    }
  } else {
    # optimizer
    if (model.optimizer == 'bobyqa') {
      out <- suppressMessages(suppressWarnings(
        lme4::glmer(formula =  model,
                    data = data.train,
                    family = model.family,
                    lme4::glmerControl(optimizer = model.optimizer,
                                       optCtrl = list(maxfun = n.iter)))
      ))
    } else if (model.optimizer == 'nloptwrap') {
      out <- suppressMessages(suppressWarnings(
       lme4::glmer(formula =  model,
                   data = data.train,
                   family = model.family,
                   lme4::glmerControl(
                     calc.derivs = FALSE,
                     optimizer = model.optimizer,
                     optCtrl = list(
                       method = "NLOPT_LN_NELDERMEAD",
                       starttests = TRUE, kkt = TRUE)))
      ))
    }
  }

  # Function output
  return(out)
}
