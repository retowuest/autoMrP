#' GB classifier
#'
#' \code{gb_classifier} applies gradient boosting classification to a data set.
#'
#' @param form Model formula. A two-sided linear formula describing
#'   the model to be fit, with the outcome on the LHS and the covariates
#'   separated by + operators on the RHS.
#' @param distribution Model distribution. A character string specifying the
#'   name of the distribution to be used.
#' @param data.train Training data. A data.frame containing the training data
#'   used to train the model.
#' @param n.trees Total number of trees. An integer-valued scalar specifying
#'   the total number of trees to be fit.
#' @param interaction.depth Interaction depth. An integer-valued scalar
#'   specifying the maximum depth of each tree.
#' @param n.minobsinnode Minimum number of observations in terminal nodes. An
#'   integer-valued scalar specifying the minimum number of observations in the
#'   terminal nodes of the trees.
#' @param shrinkage Learning rate. A numeric scalar specifying the shrinkage or
#'   learning rate applied to each tree in the expansion.
#' @param gbm_weights_xyz Weights. A numeric vector of weights to be applied to
#'   the training data. Note: Due to a bug in the gbm package, this argument
#'   is assigned to the global environemnt and then removed after the model is
#'   fit. Please make sure that you do not have an object named
#'   "gbm_weights_xyz" in your global environment.
#' @param verbose Verbose output. A logical vector indicating whether or not
#'   verbose output should be printed.
#' @return A gradient tree boosting model. A \code{\link[gbm]{gbm}} object.

gb_classifier <- function(
  form, distribution, data.train, n.trees, interaction.depth, n.minobsinnode,
  shrinkage, gbm_weights_xyz, verbose = c(TRUE, FALSE)) {

  # assign gbm_weights_xyz to parent environment
  # this is a temporary fix of a bug within the gbm package
  assign("gbm_weights_xyz", gbm_weights_xyz, envir = .GlobalEnv)

  # Train model on training data with number of total trees, interaction depth,
  # and learning rate as tuning parameters
  if (isTRUE(verbose == TRUE)) {
    out <- gbm::gbm(formula = form, distribution = distribution,
                    data = data.train, n.trees = n.trees,
                    interaction.depth = interaction.depth,
                    n.minobsinnode = n.minobsinnode,
                    shrinkage = shrinkage,
                    weights = gbm_weights_xyz,
                    train.fraction = 1, n.cores = 1,
                    keep.data = TRUE)
  } else {
    out <- suppressMessages(suppressWarnings(
      gbm::gbm(formula = form, distribution = distribution,
               data = data.train, n.trees = n.trees,
               interaction.depth = interaction.depth,
               n.minobsinnode = n.minobsinnode,
               shrinkage = shrinkage,
               weights = gbm_weights_xyz,
               train.fraction = 1, n.cores = 1,
               keep.data = TRUE)
    ))
  }

  # remove gbm_weights_xyz from parent environment
  # this is a temporary fix of a bug within the gbm package
  rm("gbm_weights_xyz", envir = .GlobalEnv)

  # Function output
  return(out)
}

#' GB classifier update
#'
#' \code{gb_classifier_update()} grows additional trees in gradient tree
#' boosting ensemble.
#'
#' @param object Gradient tree boosting output. A gbm object.
#' @param n.new.trees Number of additional trees to grow. A numeric scalar.
#' @param verbose Verbose output. A logical vector indicating whether or not
#'   verbose output should be printed.
#' @return An updated gradient tree boosting model.
#'   A \code{\link[gbm]{gbm.more}} object.

gb_classifier_update <- function(object, n.new.trees,
                                 verbose = c(TRUE, FALSE)) {
  # Train model on training data with number of total trees, interaction depth,
  # and learning rate as tuning parameters
  if (isTRUE(verbose == TRUE)) {
    out <- gbm::gbm.more(object = object,
                         n.new.trees = n.new.trees)
  } else {
    out <- suppressMessages(suppressWarnings(
      gbm::gbm.more(object = object,
                    n.new.trees = n.new.trees)
    ))
  }

  # Function output
  return(out)
}
