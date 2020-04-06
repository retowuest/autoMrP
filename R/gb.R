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
#'   depth is used to tune the model. Default is c(1, 2, 3).
#' @param shrinkage.set Learning rate. A numeric vector whose values define
#'   the learning rate or step-size reduction. Learning rate is used to tune
#'   the model. Values between 0.001 and 0.1 usually work, but a smaller
#'   learning rate typically requires more trees. Default is
#'   c(0.04, 0.01, 0.008, 0.005, 0.001).
#' @param tree.start Initial total number of trees. An integer-valued scalar
#'   specifying the initial number of total trees. Default is 50.
#' @param tree.increase.set Increase in total number of trees. Either an
#'   integer-valued scalar specifying by how many trees the total number of
#'   trees is increased (until the maximum number of trees is reached) or an
#'   integer-valued vector of `length(shrinkage.set)` with each value being
#'   associated with a learning rate. Total number of trees is used to tune the
#'   model. Default is 50.
#' @param trees.max.set Maximum number of trees. Either an integer-valued
#'   scalar specifying the maximum number of trees or an integer-valued vector
#'   of `length(shrinkage.set)` with each value being associated with a
#'   learning rate and a number of tree increase. Default is 1000.
#' @param iterations.max Stopping rule. A numeric scalar specifying the
#'   maximum number of iterations without performance improvement the GB
#'   classifier runs before stopping. Default is 70.
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
               tree.start, tree.increase.set,
               trees.max.set, iterations.max,
               n.minobsinnode, data, verbose) {

  # for debugging
  #browser()

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

  # Initialize counter for iterations
  iteration_no <- 0

  # Loop over interaction depth
  out_d <- lapply(seq_along(interaction.set), function(d) {
    # Set interaction depth
    depth <- interaction.set[d]

    if (length(tree.increase.set) == 1 & length(trees.max.set) == 1) {
      # Loop over learning rate
      out_s <- lapply(seq_along(shrinkage.set), function(s, n_trees = tree.start) {
        # Set learning rate
        shrinkage <- shrinkage.set[s]

        # Update counter for iterations
        iteration_no <- iteration_no + 1

        # Print tuning parameters
        if (isTRUE(verbose == TRUE)) {
          cat(paste("GB: Running interaction depth ", depth,
                    ", learning rate ", shrinkage,
                    ", and number of total trees ", n_trees, "\n",
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

          # Train model using tuning parameters on kth training set
          model_l <- gb_classifier(form = form,
                                   distribution = "bernoulli",
                                   data.train = data_train,
                                   n.trees = n_trees,
                                   interaction.depth = depth,
                                   n.minobsinnode = n.minobsinnode,
                                   shrinkage = shrinkage,
                                   verbose = verbose)

          # Use trained model to make predictions for kth validation set
          pred_l <- gbm::predict.gbm(
            model_l, newdata = data_valid,
            n.trees = model_l$n.trees,
            type = "response")

          # Evaluate predictions based on loss function
          perform_l <- loss_function(pred = pred_l, data.valid = data_valid,
                                     loss.unit = loss.unit,
                                     loss.measure = loss.measure,
                                     y = y, L2.unit = L2.unit)

          # Function output
          return(list(perform_l = perform_l,
                      model_l = model_l))
        })

        # Mean over all k folds
        best_error <- mean(unlist(lapply(seq_along(k_errors),
                                         function(x) {k_errors[[x]]["perform_l"]})))

        # Initialize list of tuning parameters associated with currently best
        # error
        out <- list(n_trees = n_trees,
                    depth = depth,
                    shrinkage = shrinkage,
                    error = best_error,
                    models = lapply(seq_along(k_errors),
                                    function(x) {k_errors[[x]]["model_l"]}))

        # Initialize counter for iterations since last performance improvement
        iter_since_improv <- 0

        # Loop over number of total trees
        while (n_trees < trees.max.set) {
          # Set number of total trees
          n_trees <- n_trees + tree.increase.set

          # Update counter for iterations
          iteration_no <- iteration_no + 1

          # Print tuning parameters
          if (isTRUE(verbose == TRUE)) {
            cat(paste("GB: Running interaction depth ", depth,
                      ", learning rate ", shrinkage,
                      ", and number of total trees ", n_trees, "\n",
                      "    (model no. ", iteration_no,
                      " -- iterations w/o improvement: ",
                      iter_since_improv, ")\n", sep = ""))
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
            model_l <- gb_classifier_update(object = out$models[[k]]$model_l,
                                            n.new.trees = tree.increase.set,
                                            verbose = verbose)

            # Use trained model to make predictions for kth validation set
            pred_l <- gbm::predict.gbm(
              model_l, newdata = data_valid,
              n.trees = model_l$n.trees,
              type = "response")

            # Evaluate predictions based on loss function
            perform_l <- loss_function(pred = pred_l, data.valid = data_valid,
                                       loss.unit = loss.unit,
                                       loss.measure = loss.measure,
                                       y = y, L2.unit = L2.unit)

            # Function output
            return(list(perform_l = perform_l,
                        model_l = model_l))
          })

          # Mean over all k folds
          current_error <- mean(unlist(lapply(seq_along(k_errors),
                                              function(x) {k_errors[[x]]["perform_l"]})))

          # Check if current tuning parameters outperform tuning parameters
          # that were best so far
          if (current_error < best_error) {
            #if(verbose) cat(paste("Improvement on previous model \n"), sep = "")
            best_error <- current_error
            out <- list(n_trees = n_trees,
                        depth = depth,
                        shrinkage = shrinkage,
                        error = best_error,
                        models = lapply(seq_along(k_errors),
                                        function(x) {k_errors[[x]]["model_l"]}))
            iter_since_improv <- 0
          } else {
            iter_since_improv <- iter_since_improv + 1
          }

          # Break loop if maximum number of iterations without performance
          # improvement is reached
          if (!is.null(iterations.max)) {
            if (iter_since_improv > iterations.max) {
              break
            }
          }
        }
        out
      })
      # Function output
      return(out_s)
    } else {
      # Combine learning rate, tree increase, and trees max in tuning set
      tuning_set <- dplyr::bind_cols(shrinkage.set = shrinkage.set,
                                     tree.increase.set = tree.increase.set,
                                     trees.max.set = trees.max.set)

      # Loop over tuning set
      out_s <- lapply(seq_along(tuning_set$shrinkage.set), function(s, n_trees = tree.start) {
        # Set learning rate
        shrinkage <- tuning_set$shrinkage.set[s]

        # Update counter for iterations
        iteration_no <- iteration_no + 1

        # Print tuning parameters
        if (isTRUE(verbose == TRUE)) {
          cat(paste("GB: Running interaction depth ", depth,
                    ", learning rate ", shrinkage,
                    ", and number of total trees ", n_trees, "\n",
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
          pred_l <- gbm::predict.gbm(
            model_l, newdata = data_valid,
            n.trees = model_l$n.trees,
            type = "response")

          # Evaluate predictions based on loss function
          perform_l <- loss_function(pred = pred_l, data.valid = data_valid,
                                     loss.unit = loss.unit,
                                     loss.measure = loss.measure,
                                     y = y, L2.unit = L2.unit)

          # Function output
          return(list(perform_l = perform_l,
                      model_l = model_l))
        })

        # Mean over all k folds
        best_error <- mean(unlist(lapply(seq_along(k_errors),
                                         function(x) {k_errors[[x]]["perform_l"]})))

        # Initialize list of tuning parameters associated with currently best
        # error
        out <- list(n_trees = n_trees,
                    depth = depth,
                    shrinkage = shrinkage,
                    error = best_error,
                    models = lapply(seq_along(k_errors),
                                    function(x) {k_errors[[x]]["model_l"]}))

        # Initialize counter for iterations since last performance improvement
        iter_since_improv <- 0

        # Loop over number of total trees
        while (n_trees < tuning_set$trees.max.set[s]) {
          # Set number of total trees
          n_trees <- n_trees + tuning_set$tree.increase.set[s]

          # Update counter for iterations
          iteration_no <- iteration_no + 1

          # Print tuning parameters
          if (isTRUE(verbose == TRUE)) {
            cat(paste("GB: Running interaction depth ", depth,
                      ", learning rate ", shrinkage,
                      ", and number of total trees ", n_trees, "\n",
                      "    (model no. ", iteration_no,
                      " -- iterations w/o improvement: ",
                      iter_since_improv, ")\n", sep = ""))
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
            pred_l <- gbm::predict.gbm(
              model_l, newdata = data_valid,
              n.trees = model_l$n.trees,
              type = "response")

            # Evaluate predictions based on loss function
            perform_l <- loss_function(pred = pred_l, data.valid = data_valid,
                                       loss.unit = loss.unit,
                                       loss.measure = loss.measure,
                                       y = y, L2.unit = L2.unit)

            # Function output
            return(list(perform_l = perform_l,
                        model_l = model_l))
          })

          # Mean over all k folds
          current_error <- mean(unlist(lapply(seq_along(k_errors),
                                              function(x) {k_errors[[x]]["perform_l"]})))

          # Check if current tuning parameters outperform tuning parameters
          # that were best so far
          if (current_error < best_error) {
            best_error <- current_error
            out <- list(n_trees = n_trees,
                        depth = depth,
                        shrinkage = shrinkage,
                        error = best_error,
                        models = lapply(seq_along(k_errors),
                                        function(x) {k_errors[[x]]["model_l"]}))
            iter_since_improv <- 0
          } else {
            iter_since_improv <- iter_since_improv + 1
          }

          # Break loop if maximum number of iterations without performance
          # improvement is reached
          if (!is.null(iterations.max)) {
            if (iter_since_improv > iterations.max) {
              break
            }
          }
        }
        out
      })
      # Function output
      return(out_s)
    }
  })

  # Choose best-performing model
  tuning_grid <- expand.grid(d = seq_along(interaction.set),
                             s = shrinkage.set,
                             error = NA)

  for (i in 1:nrow(tuning_grid)) {
    tuning_grid$error[i] <- out_d[[tuning_grid$d[i]]][[tuning_grid$s[i]]]$error
  }

  min_e <- which.min(tuning_grid$error)

  out <- list(interaction_depth = tuning_grid$d[min_e],
              shrinkage = tuning_grid$s[min_e],
              n_trees = out_d[[tuning_grid$d[min_e]]][[tuning_grid$s[min_e]]]$n_trees)

  # Function output
  return(out)
}
