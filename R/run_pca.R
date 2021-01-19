#' Apply PCA classifier to MrP.
#'
#' \code{run_pca} is a wrapper function that applies the PCA classifier to data
#' provided by the user, evaluates prediction performance, and chooses the
#' best-performing model.
#'
#' @inheritParams auto_MrP
#' @param data Data for cross-validation. A \code{list} of \eqn{k}
#'   \code{data.frames}, one for each fold to be used in \eqn{k}-fold
#'   cross-validation.
#'
#' @return A model formula of the winning best subset classifier model.

run_pca <- function(y, L1.x, L2.x, L2.unit, L2.reg,
                    loss.unit, loss.fun, data, cores,
                    verbose) {

  # List of all models to be evaluated
  models <- model_list_pca(
      y = y,
      L1.x = L1.x,
      L2.x = L2.x,
      L2.unit = L2.unit,
      L2.reg = L2.reg)

  # prallel tuning if cores > 1
  if( cores > 1 ){

    # Train all models in parallel
    m_errors <- run_best_subset_mc(
      verbose = verbose,
      models = models,
      data = data,
      loss.unit = loss.unit,
      loss.fun = loss.fun,
      y = y,
      L1.x = L1.x,
      L2.x = L2.x,
      L2.unit = L2.unit,
      L2.reg = L2.reg,
      cores = cores)
  } else{
    # Train and evaluate each model
    m_errors <- lapply(seq_along(models), function(m) {
      # Print model m
      if (isTRUE(verbose)) {
        M <- length(models)
        cat(paste("Best subset: Running model ", m,
                  " out of ", M, " models\n", sep = ""))
      }

      # Loop over each fold
      k_errors <- lapply(seq_along(data), function(k) {
        # Split data in training and validation sets
        data_train <- dplyr::bind_rows(data[-k])
        data_valid <- dplyr::bind_rows(data[k])

        # Train mth model on kth training set
        model_m <- best_subset_classifier(model = models[[m]],
                                          data.train = data_train,
                                          model.family = binomial(link = "probit"),
                                          model.optimizer = "bobyqa",
                                          n.iter = 1000000,
                                          verbose = verbose)

        # Use trained model to make predictions for kth validation set
        pred_m <- stats::predict(model_m, newdata = data_valid,
                                 type = "response", allow.new.levels = TRUE)

        # Evaluate predictions based on loss function
        perform_m <- loss_function(pred = pred_m,
                                   data.valid = data_valid,
                                   loss.unit = loss.unit,
                                   loss.fun = loss.fun,
                                   y = y,
                                   L2.unit = L2.unit)
      })

      # Mean over loss functions
      k_errors <- dplyr::bind_rows(k_errors) %>%
        dplyr::group_by(measure) %>%
        dplyr::summarise(value = mean(value), .groups = "drop") %>%
        dplyr::mutate(model = m)
    })
  }

  # Extract best tuning parameters
  grid_cells <- dplyr::bind_rows(m_errors)
  best_params <- dplyr::slice(loss_score_ranking(score = grid_cells, loss.fun = loss.fun), 1)

  # Choose best-performing model
  out <- models[[ dplyr::pull(.data = best_params, var = model) ]]

  # Function output
  return(out)

}
