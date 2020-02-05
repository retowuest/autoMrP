#' Apply gradient boosting classifier to MrP.
#'
#' \code{gb} is a wrapper function that applies the gradient boosting classifier
#' to data provided by the user, evaluates prediction performance, and chooses
#' the best-performing model.
#'
#' @param y Outcome variable. A character scalar containing the column name of
#'   the outcome variable.
#' @param L1.x Individual-level covariates. A character vector of column names
#'   corresponding to the individual-level variables used to predict the outcome
#'   variable. Must include geographic unit.
#' @param L2.x Context-level covariates. A character vector of column names
#'   corresponding to the context-level variables used to predict the outcome
#'   variable.
#' @param L2.unit Geographic unit. A character scalar indicating the column
#'   name of the geographic unit at which outcomes should be aggregated.
#' @param L2.reg Geographic region. A character scalar indicating the column
#'   name of the geographic region by which geographic units are grouped
#'   (L2.unit must be nested within L2.re).
#' @param loss.unit Loss function unit. A character-valued scalar indicating
#'   whether the loss should be evaluated at the level of individual respondents
#'   or the level of geographic units. Default is at the individual level.
#' @param loss.measure Loss function measure. A character-valued scalar
#'   indicating whether the loss should be measured by the mean squared error
#'   or the mean absolute error. Default is the MSE.
#' @param data Data for cross-validation. A list of k data.frames, one for
#'   each fold used in k-fold cross-validation.
#' @param L2.unit.include Include L2.unit in GB. A logical argument indicating
#'   whether L2.unit is included in the GB models. Default is FALSE.
#' @param L2.reg.include Include L2.reg in GB. A logical argument indicating
#'   whether L2.reg is included in the GB models. Default is FALSE.
#' @param interaction.set Set of interaction depth values. An integer-valued
#'   vector whose values define the maximum depth of each tree. Interaction
#'   depth is used to tune the model.
#' @param shrinkage.set Learning rate. A numeric vector whose values define
#'   the learning rate or step-size reduction. Learning rate is used to tune
#'   the model. Values between 0.001 and 0.1 usually work, but a smaller
#'   learning rate typically requires more trees.
#' @param tree.start Initial total number of trees. An integer-valued scalar
#'   specifying the initial number of total trees. Default is 2.
#' @param tree.increase.set Increase in total number of trees. Either an
#'   integer-valued scalar specifying by how many trees the total number of
#'   trees is increased (until the maximum number of trees is reached) or an
#'   integer-valued vector of `length(shrinkage.set)` with each value being
#'   associated with a learning rate. Total number of trees is used to tune the
#'   model.
#' @param trees.max.set Maximum number of trees. Either an integer-valued
#'   scalar specifying the maximum number of trees or an integer-valued vector
#'   of `length(shrinkage.set)` with each value being associated with a
#'   learning rate and a number of tree increase.
#' @param n.minobsinnode Minimum number of observations in the terminal nodes.
#'   An integer-valued scalar specifying the minimum number of observations
#'   that each terminal node of the trees must contain. Default is 5.
#' @param verbose Verbose output. A logical vector indicating whether or not
#'   verbose output should be printed.
#' @return
#' @examples

gb <- function(y, L1.x, L2.x, L2.unit, L2.reg,
               L2.unit.include = c(TRUE, FALSE),
               L2.reg.include = c(TRUE, FALSE),
               loss.unit, loss.measure,
               interaction.set, shrinkage.set,
               tree.start, tree.increase.set, trees.max.set,
               n.minobsinnode, data, verbose) {
  # Evaluate inclusion of L2.unit
  if (isTRUE(L2.unit.include == FALSE)) {
    L2.unit <- NULL
  }

  # Evaluate inclusion of L2.reg
  if (isTRUE(L2.reg.include == FALSE)) {
    L2.reg <- NULL
  }

  # Create model formula
  x <- paste(c(L1.x, L2.x, L2.unit, L2.reg), collapse = " + ")
  form <- as.formula(paste(y, " ~ ", x, sep = ""))

  # Prepare data
  data <- lapply(data, function(k) {
    dplyr::select_at(k, c(y, L1.x, L2.x, L2.unit, L2.reg))
  })

  # Initialize number of total trees
  n_trees <- tree.start

  # Initialize counter for iterations
  iteration_no <- 0

  # Loop over interaction depth
  for (d in seq_along(interaction.set)) {
    # Set interaction depth
    depth <- interaction.set[d]

    if (length(tree.increase.set) == 1 & length(trees.max.set) == 1) {
      # Loop over learning rate
      for (s in seq_along(shrinkage.set)) {
        # Set learning rate
        shrinkage <- shrinkage.set[s]

        # Update counter for iterations
        iteration_no <- iteration_no + 1

        # Print tuning parameters
        if (isTRUE(verbose == TRUE)) {
          cat(paste("GB: Running interaction depth ", depth,
                    " and learning rate ", shrinkage, "\n",
                    "    (model no. ", iteration_no,
                    " -- no improvement evaluation)\n", sep = ""))
        }

        # Loop over each fold
        k_errors <- lapply(seq_along(data), function(k) {
          # Split data in training and validation sets
          data_train <- dplyr::bind_rows(data[-k])
          data_valid <- dplyr::bind_rows(data[k])

          # Convert individual-level, geographic unit, and geographic region
          # covariates to factor variables in training and validation sets
          data_train <- data_train %>%
            dplyr::mutate_at(.vars = c(L1.x, L2.unit, L2.reg), as.factor)

          data_valid <- data_valid %>%
            dplyr::mutate_at(.vars = c(L1.x, L2.unit, L2.reg), as.factor)

          # Train model using lambda on kth training set
          model_l <- gb_classifier(form = form,
                                   distribution = "bernoulli",
                                   data.train = data_train,
                                   n.trees = n_trees,
                                   interaction.depth = depth,
                                   n.minobsinnode = n.minobsinnode,
                                   shrinkage = shrinkage,
                                   verbose = verbose)

          # Use trained model to make predictions for kth validation set
          pred_l <- stats::predict(model_l, newdata = data_valid)

          # Evaluate predictions based on loss function
          perform_l <- loss_function(pred = pred_l, data.valid = data_valid,
                                     loss.unit = loss.unit,
                                     loss.measure = loss.measure,
                                     y = y, L2.unit = L2.unit)
        })


      }
    } else {
      # Combine learning rate, tree increase, and trees max in tuning set
      tuning_set <- dplyr::bind_cols(shrinkage.set = shrinkage.set,
                                     tree.increase.set = tree.increase.set,
                                     trees.max.set = trees.max.set)

      # Loop over tuning set
      for (s in seq_along(tuning_set)) {

      }
    }
  }
}
