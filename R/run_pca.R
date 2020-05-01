#' Apply PCA classifier to MrP.
#'
#' \code{run_pca} is a wrapper function that applies the PCA classifier to data
#' provided by the user, evaluates prediction performance, and chooses the
#' best-performing model.
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
#' @param data Data for cross-validation. A \code{list} of \eqn{k}
#'   \code{data.frames}, one for each fold to be used in \eqn{k}-fold
#'   cross-validation.
#' @param verbose Verbose output. A logical argument indicating whether or not
#'   verbose output should be printed. Default is \code{TRUE}.
#' @return
#' @examples

run_pca <- function(y, L1.x, L2.x, L2.unit, L2.reg,
                    loss.unit, loss.fun,
                    data, verbose) {

  # List of all models to be evaluated
  models <- model_list_pca(y = y,
                           L1.x = L1.x,
                           L2.x = L2.x,
                           L2.unit = L2.unit,
                           L2.reg = L2.reg)

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

    # Mean over all k folds
    mean(unlist(k_errors))
  })

  # Choose best-performing model
  min.m <- which.min(m_errors)
  out <- models[[min.m]]

  # Function output
  return(out)
}
