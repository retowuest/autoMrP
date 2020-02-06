#' Apply best subset classifier to MrP.
#'
#' \code{best_subset} is a wrapper function that applies the best subset
#' classifier to a list of models provided by the user, evaluates the models'
#' prediction performance, and chooses the best-performing model.
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
#' @param verbose Verbose output. A logical vector indicating whether or not
#'   verbose output should be printed.
#' @return
#' @examples

best_subset <- function(y, L1.x, L2.x, L2.unit, L2.reg,
                        loss.unit, loss.measure,
                        data, verbose) {
  # List of all models to be evaluated
  models <- model_list(y = y,
                       L1.x = L1.x,
                       L2.x = L2.x,
                       L2.unit = L2.unit,
                       L2.reg = L2.reg)

  # Train and evaluate each model
  m_errors <- lapply(seq_along(models), function(m) {
    # Print model m
    if (isTRUE(verbose == TRUE)) {
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
      perform_m <- loss_function(pred = pred_m, data.valid = data_valid,
                                 loss.unit = loss.unit,
                                 loss.measure = loss.measure,
                                 y = y, L2.unit = L2.unit)
    })

    # Mean over all k folds
    mean(unlist(k_errors))
  })

  # Choose best-performing model
  min_m <- which.min(m_errors)
  out <- models[[min_m]]

  # Function output
  return(out)
}
