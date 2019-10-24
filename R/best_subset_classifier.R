#' Best subset classifier
#'
#' \text{best_subset_classifier} applies best subset classification to a data
#' set.
#'
#' @param model Multilevel model. A model formula describing the multilevel
#'   model to be estimated on the basis of the provided training data.
#' @param data.train Training data. A data.frame containing the training data
#'   used to train the model.
#' @param model.family Model family. A character-valued scalar indicating the
#'   model family to be used by glmer. Defaults to "binomial(link = 'probit')".
#' @param model.optimizer Optimization method. A character-valued scalar
#'   describing the optimization method to be used by glmer. Defaults to
#'   "bobyqa".
#' @param n.iter Iterations. A integer-valued scalar specifying the maximum
#'   number of function evaluations tried by the optimization method.
#' @return
#' @examples

best_subset_classifier <- function(model, data.train,
                                   model.family = "binomial(link = 'probit')",
                                   model.optimizer = "bobyqa",
                                   n.iter = 1000000) {
  # Train model on training data
  out <- lme4::glmer(model, data = data.train,
                     family = get(model.family),
                     glmerControl(optimizer = model.optimizer,
                                  optCtrl = list(maxfun = n.iter)))

  # Function output
  return(out)
}
