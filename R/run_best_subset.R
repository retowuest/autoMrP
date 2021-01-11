#' Apply best subset classifier to MrP.
#'
#' \code{run_best_subset} is a wrapper function that applies the best subset
#' classifier to a list of models provided by the user, evaluates the models'
#' prediction performance, and chooses the best-performing model.
#'
#' @inheritParams auto_MrP
#' @param data Data for cross-validation. A \code{list} of \eqn{k}
#'   \code{data.frames}, one for each fold to be used in \eqn{k}-fold
#'   cross-validation.
#' @return A model formula of the winning best subset classifier model.

run_best_subset <- function(y, L1.x, L2.x, L2.unit, L2.reg,
                            loss.unit, loss.fun,
                            data, verbose, cores) {

  # List of all models to be evaluated
  models <- model_list(y = y,
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

################################################################################
#                Multicore tuning for best subset                              #
################################################################################
#' Best subset multicore tuning.
#'
#' \code{run_best_subset_mc} is called from within \code{run_best_subset}. It
#' tunes using multiple cores.
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
#' @param cores The number of cores to be used. An integer indicating the number
#'   of processor cores used for parallel computing. Default is 1.
#' @param models The models to perform best subset selection on. A list of model
#'   formulas.
#' @param verbose Verbose output. A logical argument indicating whether or not
#'   verbose output should be printed. Default is \code{TRUE}.
#' @return The cross-validation errors for all models. A list.
#' @examples \dontrun{
#' # not yet
#' }

run_best_subset_mc <- function(y, L1.x, L2.x, L2.unit, L2.reg,
                               loss.unit, loss.fun, data,
                               cores, models, verbose){

  # Binding for global variables
  m <- NULL

  # Register cores
  cl <- multicore(cores = cores, type = "open", cl = NULL)

  # Train and evaluate each model
  m_errors <- foreach::foreach(m = 1:length(models), .packages = 'autoMrP' ) %dorng% {

    # Loop over each fold
    k_errors <- lapply(seq_along(data), function(k) {
      # Split data in training and validation sets
      data_train <- dplyr::bind_rows(data[-k])
      data_valid <- dplyr::bind_rows(data[k])

      # Train mth model on kth training set
      model_m <- best_subset_classifier(
        model = models[[m]],
        data.train = data_train,
        model.family = binomial(link = "probit"),
        model.optimizer = "bobyqa",
        n.iter = 1000000,
        verbose = verbose)

      # Use trained model to make predictions for kth validation set
      pred_m <- stats::predict(
        model_m, newdata = data_valid,
        type = "response", allow.new.levels = TRUE)

      # Evaluate predictions based on loss function
      perform_m <- loss_function(
        pred = pred_m,
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
  }

  # De-register cluster
  multicore(cores = cores, type = "close", cl = cl)

  # Function output
  return(m_errors)
}
