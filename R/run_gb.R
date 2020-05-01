#' Apply gradient boosting classifier to MrP.
#'
#' \code{run_gb} is a wrapper function that applies the gradient boosting
#' classifier to data provided by the user, evaluates prediction performance,
#' and chooses the best-performing model.
#'
#' @param y Outcome variable. A character scalar containing the column name of
#'   the outcome variable in \code{survey}.
#' @param L1.x Individual-level covariates. A character vector containing the
#'   column names of the individual-level variables in \code{survey} and
#'   \code{census} used to predict outcome \code{y}. Note that geographic unit
#'   is specified in argument \code{L2.unit}.
#' @param L2.x Context-level covariates. A character vector containing the
#'   column names of the context-level variables in \code{survey} and
#'   \code{census} used to predict outcome \code{y}.
#' @param L2.unit Geographic unit. A character scalar containing the column
#'   name of the geographic unit in \code{survey} and \code{census} at which
#'   outcomes should be aggregated.
#' @param L2.reg Geographic region. A character scalar containing the column
#'   name of the geographic region in \code{survey} and \code{census} by which
#'   geographic units are grouped (\code{L2.unit} must be nested within
#'   \code{L2.reg}). Default is \code{NULL}.
#' @param loss.unit Loss function unit. A character-valued scalar indicating
#'   whether performance loss should be evaluated at the level of individual
#'   respondents (\code{individuals}) or geographic units (\code{L2 units}).
#'   Default is \code{individuals}.
#' @param loss.fun Loss function. A character-valued scalar indicating whether
#'   prediction loss should be measured by the mean squared error (\code{MSE})
#'   or the mean absolute error (\code{MAE}). Default is \code{MSE}.
#' @param interaction.depth GB interaction depth. An integer-valued vector
#'   whose values specify the interaction depth of GB. The interaction depth
#'   defines the maximum depth of each tree grown (i.e., the maximum level of
#'   variable interactions). Default is \code{c(1, 2, 3)}.
#' @param shrinkage GB learning rate. A numeric vector whose values specify the
#'   learning rate or step-size reduction of GB. Values between \eqn{0.001}
#'   and \eqn{0.1} usually work, but a smaller learning rate typically requires
#'   more trees. Default is \code{c(0.04, 0.01, 0.008, 0.005, 0.001)}.
#' @param n.trees.init GB initial total number of trees. An integer-valued
#'   scalar specifying the initial number of total trees to fit by GB. Default
#'   is \eqn{50}.
#' @param n.trees.increase GB increase in total number of trees. An
#'   integer-valued scalar specifying by how many trees the total number of
#'   trees to fit should be increased (until \code{n.trees.max} is reached)
#'   or an integer-valued vector of length \code{length(shrinkage)} with each
#'   of its values being associated with a learning rate in \code{shrinkage}.
#'   Default is \eqn{50}.
#' @param n.trees.max GB maximum number of trees. An integer-valued scalar
#'   specifying the maximum number of trees to fit by GB or an integer-valued
#'   vector of length \code{length(shrinkage)} with each of its values being
#'   associated with a learning rate and an increase in the total number of
#'   trees. Default is \eqn{1000}.
#' @param n.iter GB number of iterations without improvement. A numeric scalar
#'   specifying the maximum number of iterations without performance
#'   improvement the algorithm runs before stopping. Default is \eqn{70}.
#' @param n.minobsinnode GB minimum number of observations in the terminal
#'   nodes. An integer-valued scalar specifying the minimum number of
#'   observations that each terminal node of the trees must contain. Default is
#'   \eqn{5}.
#' @param data Data for cross-validation. A \code{list} of \eqn{k}
#'   \code{data.frames}, one for each fold to be used in \eqn{k}-fold
#'   cross-validation.
#' @param verbose Verbose output. A logical argument indicating whether or not
#'   verbose output should be printed. Default is \code{TRUE}.
#' @return
#' @examples

run_gb <- function(y, L1.x, L2.x, L2.unit, L2.reg,
                   loss.unit, loss.fun,
                   interaction.depth, shrinkage,
                   n.trees.init, n.trees.increase,
                   n.trees.max, n.iter,
                   n.minobsinnode, data, verbose) {

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
  out_d <- lapply(seq_along(interaction.depth), function(d) {
    # Set interaction depth
    depth <- interaction.depth[d]

    if (length(n.trees.increase) == 1 & length(n.trees.max) == 1) {
      # Loop over learning rate
      out_s <- lapply(seq_along(shrinkage),
                      function(s, n_trees = n.trees.init) {
        # Set learning rate
        shrinkage_value <- shrinkage[s]

        # Update counter for iterations
        iteration_no <- iteration_no + 1

        # Print tuning parameters
        if (isTRUE(verbose)) {
          cat(paste("GB: Running interaction depth ", depth,
                    ", learning rate ", shrinkage_value,
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
                                   shrinkage = shrinkage_value,
                                   verbose = verbose)

          # Use trained model to make predictions for kth validation set
          pred_l <- gbm::predict.gbm(model_l, newdata = data_valid,
                                     n.trees = model_l$n.trees,
                                     type = "response")

          # Evaluate predictions based on loss function
          perform_l <- loss_function(pred = pred_l,
                                     data.valid = data_valid,
                                     loss.unit = loss.unit,
                                     loss.fun = loss.fun,
                                     y = y,
                                     L2.unit = L2.unit)

          # Function output
          return(list(perform_l = perform_l,
                      model_l = model_l))
        })

        # Mean over all k folds
        best_error <- mean(unlist(lapply(seq_along(k_errors),
                                         function(x) {k_errors[[x]]["perform_l"]})))

        # Initialize list of tuning parameters associated with the currently
        # best error
        out <- list(n_trees = n_trees,
                    depth = depth,
                    shrinkage = shrinkage_value,
                    error = best_error,
                    models = lapply(seq_along(k_errors),
                                    function(x) {k_errors[[x]]["model_l"]}))

        # Initialize counter for iterations since last performance improvement
        iter_since_improv <- 0

        # Loop over number of total trees
        while (n_trees < n.trees.max) {
          # Set number of total trees
          n_trees <- n_trees + n.trees.increase

          # Update counter for iterations
          iteration_no <- iteration_no + 1

          # Print tuning parameters
          if (isTRUE(verbose)) {
            cat(paste("GB: Running interaction depth ", depth,
                      ", learning rate ", shrinkage_value,
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
                                            n.new.trees = n.trees.increase,
                                            verbose = verbose)

            # Use trained model to make predictions for kth validation set
            pred_l <- gbm::predict.gbm(model_l, newdata = data_valid,
                                       n.trees = model_l$n.trees,
                                       type = "response")

            # Evaluate predictions based on loss function
            perform_l <- loss_function(pred = pred_l,
                                       data.valid = data_valid,
                                       loss.unit = loss.unit,
                                       loss.fun = loss.fun,
                                       y = y,
                                       L2.unit = L2.unit)

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
                        shrinkage = shrinkage_value,
                        error = best_error,
                        models = lapply(seq_along(k_errors),
                                        function(x) {k_errors[[x]]["model_l"]}))
            iter_since_improv <- 0
          } else {
            iter_since_improv <- iter_since_improv + 1
          }

          # Break loop if maximum number of iterations without performance
          # improvement is reached
          if (!is.null(n.iter)) {
            if (iter_since_improv > n.iter) {
              break
            }
          }
        }
        out
      })
      # Function output
      return(out_s)
    } else {
      # Combine shrinkage, n.trees.increase, and n.trees.max in tuning set
      tuning_set <- dplyr::bind_cols(shrinkage = shrinkage,
                                     n.trees.increase = n.trees.increase,
                                     n.trees.max = n.trees.max)

      # Loop over tuning set
      out_s <- lapply(seq_along(tuning_set$shrinkage),
                      function(s, n_trees = n.trees.init) {
        # Set learning rate
        shrinkage_value <- tuning_set$shrinkage[s]

        # Update counter for iterations
        iteration_no <- iteration_no + 1

        # Print tuning parameters
        if (isTRUE(verbose)) {
          cat(paste("GB: Running interaction depth ", depth,
                    ", learning rate ", shrinkage_value,
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
                                   shrinkage = shrinkage_value,
                                   verbose = verbose)

          # Use trained model to make predictions for kth validation set
          pred_l <- gbm::predict.gbm(model_l, newdata = data_valid,
                                     n.trees = model_l$n.trees,
                                     type = "response")

          # Evaluate predictions based on loss function
          perform_l <- loss_function(pred = pred_l,
                                     data.valid = data_valid,
                                     loss.unit = loss.unit,
                                     loss.fun = loss.fun,
                                     y = y,
                                     L2.unit = L2.unit)

          # Function output
          return(list(perform_l = perform_l,
                      model_l = model_l))
        })

        # Mean over all k folds
        best_error <- mean(unlist(lapply(seq_along(k_errors),
                                         function(x) {k_errors[[x]]["perform_l"]})))

        # Initialize list of tuning parameters associated with the currently
        # best error
        out <- list(n_trees = n_trees,
                    depth = depth,
                    shrinkage = shrinkage_value,
                    error = best_error,
                    models = lapply(seq_along(k_errors),
                                    function(x) {k_errors[[x]]["model_l"]}))

        # Initialize counter for iterations since last performance improvement
        iter_since_improv <- 0

        # Loop over number of total trees
        while (n_trees < tuning_set$n.trees.max[s]) {
          # Set number of total trees
          n_trees <- n_trees + tuning_set$n.trees.increase[s]

          # Update counter for iterations
          iteration_no <- iteration_no + 1

          # Print tuning parameters
          if (isTRUE(verbose)) {
            cat(paste("GB: Running interaction depth ", depth,
                      ", learning rate ", shrinkage_value,
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
                                     shrinkage = shrinkage_value,
                                     verbose = verbose)

            # Use trained model to make predictions for kth validation set
            pred_l <- gbm::predict.gbm(model_l, newdata = data_valid,
                                       n.trees = model_l$n.trees,
                                       type = "response")

            # Evaluate predictions based on loss function
            perform_l <- loss_function(pred = pred_l,
                                       data.valid = data_valid,
                                       loss.unit = loss.unit,
                                       loss.fun = loss.fun,
                                       y = y,
                                       L2.unit = L2.unit)

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
                        shrinkage = shrinkage_value,
                        error = best_error,
                        models = lapply(seq_along(k_errors),
                                        function(x) {k_errors[[x]]["model_l"]}))
            iter_since_improv <- 0
          } else {
            iter_since_improv <- iter_since_improv + 1
          }

          # Break loop if maximum number of iterations without performance
          # improvement is reached
          if (!is.null(n.iter)) {
            if (iter_since_improv > n.iter) {
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
  tuning_grid <- expand.grid(d = seq_along(interaction.depth),
                             s = shrinkage,
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
