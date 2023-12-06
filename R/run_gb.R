#' Apply gradient boosting classifier to MrP.
#'
#' \code{run_gb} is a wrapper function that applies the gradient boosting
#' classifier to data provided by the user, evaluates prediction performance,
#' and chooses the best-performing model.
#'
#' @inheritParams auto_MrP
#' @param L2.eval.unit Geographic unit for the loss function. A character scalar
#'   containing the column name of the geographic unit in \code{survey} and
#'   \code{census}.
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
#' @param n.minobsinnode GB minimum number of observations in the terminal
#'   nodes. An integer-valued scalar specifying the minimum number of
#'   observations that each terminal node of the trees must contain. Default is
#'   \eqn{5}.
#' @param data Data for cross-validation. A \code{list} of \eqn{k}
#'   \code{data.frames}, one for each fold to be used in \eqn{k}-fold
#'   cross-validation.
#' @param verbose Verbose output. A logical argument indicating whether or not
#'   verbose output should be printed. Default is \code{TRUE}.
#' @param cores The number of cores to be used. An integer indicating the number
#'   of processor cores used for parallel computing. Default is 1.
#' @return The tuned gradient boosting parameters. A list with three elements:
#'   \code{interaction_depth} contains the interaction depth parameter,
#'   \code{shrinkage} contains the learning rate, \code{n_trees} the number of
#'   trees to be grown.

run_gb <- function(y, L1.x, L2.x, L2.eval.unit, L2.unit, L2.reg,
                   loss.unit, loss.fun, interaction.depth, shrinkage,
                   n.trees.init, n.trees.increase, n.trees.max,
                   cores = cores, n.minobsinnode, data, verbose) {

  # Create model formula
  x <- paste(c(L1.x, L2.x, L2.unit, L2.reg), collapse = " + ")
  form <- as.formula(paste(y, " ~ ", x, sep = ""))

  # Prepare data
  data <- lapply(data, function(k) {
    dplyr::select_at(k, c(y, L1.x, L2.x, L2.eval.unit, L2.reg))
  })

  # Number of trees
  n_trees <- seq(from = n.trees.init, to = n.trees.max, by = n.trees.increase)

  # Search grid
  gb_grid <- expand.grid(interaction.depth, shrinkage, n_trees)
  names(gb_grid) <- c("depth", "shrinkage", "ntrees")

  ## tuning with 1) multiple cores; 2) a single core
  if (cores > 1){

    # 1) multiple cores
    grid_cells <- run_gb_mc(
      y = y, L1.x = L1.x, L2.eval.unit = L2.eval.unit, L2.unit = L2.unit,
      L2.reg = L2.reg, form = form, gb.grid = gb_grid,
      n.minobsinnode = n.minobsinnode,  loss.unit = loss.unit,
      loss.fun = loss.fun, data = data, cores = cores)
  } else{

    # 2) single core
    # loop over tuning grid
    grid_cells <- apply( gb_grid, 1, function(g) {

      # Set tuning parameters
      depth <- as.numeric(g["depth"])
      shrinkage_value <- as.numeric(g["shrinkage"])
      ntrees <- as.numeric(g["ntrees"])

      # Print tuning parameters
      if (isTRUE(verbose)) {
        cat(paste("GB: Running interaction depth ", depth,
                  ", learning rate ", shrinkage_value,
                  ", and number of total trees ", ntrees, "\n", sep = ""))
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
                                 n.trees = ntrees,
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
                                   L2.unit = L2.eval.unit)

      })

      # Mean over loss functions
      k_errors <- dplyr::bind_rows(k_errors) %>%
        dplyr::group_by(measure) %>%
        dplyr::summarise(value = mean(value), .groups = "drop") %>%
        dplyr::mutate(ntrees = ntrees, depth = depth, shrinkage = shrinkage_value)

    })
  }

  # Extract best tuning parameters
  grid_cells <- dplyr::bind_rows(grid_cells)
  best_params <- dplyr::slice(loss_score_ranking(score = grid_cells, loss.fun = loss.fun), 1)

  out <- list(interaction_depth = dplyr::pull(.data = best_params, var = depth),
              shrinkage = dplyr::pull(.data = best_params, var = shrinkage),
              n_trees = dplyr::pull(.data = best_params, var = ntrees))

  # Function output
  return(out)
}


###########################################################################
# Multicore tuning for GB -------------------------------------------------
###########################################################################
#' GB multicore tuning.
#'
#' \code{run_gb_mc} is called from within \code{run_gb}. It tunes using
#' multiple cores.
#'
#' @inheritParams run_gb
#' @param form The model formula. A formula object.
#' @param gb.grid The hyper-parameter search grid. A matrix of all
#'   hyper-parameter combinations.
#' @return The tuning parameter combinations and there associated loss function
#'   scores. A list.

run_gb_mc <- function(y, L1.x, L2.eval.unit, L2.unit, L2.reg, form, gb.grid,
                      n.minobsinnode, loss.unit, loss.fun, data, cores){

  # Binding for global variables
  g <- NULL
  `%>%` <- dplyr::`%>%`

  # Register cores
  cl <- multicore(cores = cores, type = "open", cl = NULL)

  # Train and evaluate each model
  grid_cells <- foreach::foreach(g = 1:nrow(gb.grid), .packages = 'autoMrP',
                                 .errorhandling = "pass") %dorng% {

    # Set tuning parameters
    depth <- as.numeric( gb.grid[g, "depth"] )
    shrinkage_value <- as.numeric( gb.grid[g, "shrinkage"] )
    ntrees <- as.numeric( gb.grid[g, "ntrees"] )

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
                               n.trees = ntrees,
                               interaction.depth = depth,
                               n.minobsinnode = n.minobsinnode,
                               shrinkage = shrinkage_value,
                               verbose = FALSE)

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
                                 L2.unit = L2.eval.unit)

    })

    # Mean over loss functions
    k_errors <- dplyr::bind_rows(k_errors) %>%
      dplyr::group_by(measure) %>%
      dplyr::summarise(value = mean(value), .groups = "drop") %>%
      dplyr::mutate(ntrees = ntrees, depth = depth, shrinkage = shrinkage_value)
  }

  # De-register cluster
  multicore(cores = cores, type = "close", cl = cl)

  # Function output
  return(grid_cells)
}
