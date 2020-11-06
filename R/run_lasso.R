#' Apply lasso classifier to MrP.
#'
#' \code{run_lasso} is a wrapper function that applies the lasso classifier to
#' data provided by the user, evaluates prediction performance, and chooses the
#' best-performing model.
#'
#' @inheritParams auto_MrP
#' @param lambda Lasso penalty parameter. A numeric \code{vector} of
#'   non-negative values or a \code{list} of two numeric vectors of equal size,
#'   with the first vector containing the step sizes by which the penalty
#'   parameter should increase and the second vector containing the upper
#'   thresholds of the intervals to which the step sizes apply. The penalty
#'   parameter controls the shrinkage of the context-level variables in the
#'   lasso model. Default is \code{1 / exp(- seq(from = -1, to = 4.5, length =
#'   100))}.
#' @param n.iter Lasso number of iterations without improvement. Either
#'   \code{NULL} or an integer-valued scalar specifying the maximum number of
#'   iterations without performance improvement the algorithm runs before
#'   stopping. Default is \eqn{70}.
#' @param data Data for cross-validation. A \code{list} of \eqn{k}
#'   \code{data.frames}, one for each fold to be used in \eqn{k}-fold
#'   cross-validation.
#'
#' @return The tuned lambda value. A numeric scalar.
#' @examples \dontrun{
#' # create list of cross-validation folds
#' cv_folds <- list(
#'   `1` = survey_item[1:200, ],
#'   `2` = survey_item[201:400, ],
#'   `3` = survey_item[401:1500, ])
#'
#' # run lasso classifier
#' lasso_out <- run_lasso(
#'   y = "YES",
#'   L1.x = c("L1x1", "L1x2"),
#'   L2.x = c("L2.x1", "L2.x2"),
#'   L2.unit = "state",
#'   L2.reg = "region",
#'   loss.unit = "individuals",
#'   loss.fun = "MSE",
#'   lambda = list(c(0.1, 0.3, 1), c(1, 10, 10000)),
#'   n.iter = 70,
#'   data = cv_folds,
#'   verbose = TRUE)
#' }

run_lasso <- function(y, L1.x, L2.x, L2.unit, L2.reg,
                      n.iter, loss.unit, loss.fun,
                      lambda, data, verbose, cores) {

  # Lasso search grid
  if ( is.null(lambda) ){
    lambda <- log_spaced(min = 0.1, max = 250, n = n.iter)
  }

  # Context-level fixed effects
  L2_fe <- paste(L2.x, collapse = " + ")
  L2_fe_form <- as.formula(paste(y, " ~ ", L2_fe, sep = ""))

  # Individual-level random effects as named list
  L1_re <- setNames(
    as.list(rep(c(~ 1), times = length(c(L1.x, L2.unit, L2.reg)))),
    c(L1.x, L2.unit, L2.reg))

  # Parallel processing
  if (cores > 1){
    lambda_errors <- run_lasso_mc_lambda(
      y = y, L1.x = L1.x, L2.x = L2.x, L2.unit = L2.unit, L2.reg = L2.reg,
      loss.unit = loss.unit, loss.fun = loss.fun, data = data,
      cores = cores, L2.fe.form = L2_fe_form, L1.re = L1_re,
      lambda = lambda)
  } else{

    # Train and evaluate each model
    lambda_errors <- lapply(seq_along(lambda), function(l) {

      # Print lambda value
      if (isTRUE(verbose)) {
        L <- length(lambda)
        cat(paste("Lasso: Running lambda w/ value ", lambda[l],
                  " (lambda ", l, " out of max. ",
                  L, " lambdas)\n", sep = ""))
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

        # Train model using lambda value on kth training set
        model_l <- lasso_classifier(L2.fix = L2_fe_form,
                                    L1.re = L1_re,
                                    data.train = data_train,
                                    lambda = lambda[l],
                                    model.family = binomial(link = "probit"),
                                    verbose = verbose)

        # Use trained model to make predictions for kth validation set
        pred_l <- stats::predict(model_l, newdata = data.frame(data_valid))

        # Evaluate predictions based on loss function
        perform_l <- loss_function(pred = pred_l, data.valid = data_valid,
                                   loss.unit = loss.unit,
                                   loss.fun = loss.fun,
                                   y = y, L2.unit = L2.unit)
      })

      # Mean over loss functions
      k_errors <- dplyr::bind_rows(k_errors) %>%
        dplyr::group_by(measure) %>%
        dplyr::summarise(value = mean(value), .groups = "drop") %>%
        dplyr::mutate(lambda = lambda[l] )
      })
  }
  # Extract best tuning parameters
  grid_cells <- dplyr::bind_rows(lambda_errors)
  best_params <- dplyr::slice(loss_score_ranking(score = grid_cells, loss.fun = loss.fun), 1)
}


################################################################################
#                Multicore tuning for lasso parallel across lambda values      #
################################################################################
#' Lasso multicore tuning.
#'
#' \code{run_lasso_mc_lambda} is called from within \code{run_lasso}. It
#' tunes using multiple cores.
#'
#' @inheritParams auto_MrP
#' @inheritParams run_lasso
#' @param L2.fe.form The fixed effects part of the Lasso classifier formula. The
#'   formula is inherited from \code{run_lasso}.
#' @param L1.re A list of random effects for the Lasso classifier formula. The
#'   formula is inherited from \code{run_lasso}.
#' @return The cross-validation errors for all models. A list.
#' @examples \dontrun{
#' # not yet
#' }

run_lasso_mc_lambda <- function(
  y, L1.x, L2.x, L2.unit, L2.reg,
  loss.unit, loss.fun, data,
  cores, L2.fe.form, L1.re, lambda){

  # Binding for global variables
  `%>%` <- dplyr::`%>%`
  l <- NULL

  # Register cores
  cl <- multicore(cores = cores, type = "open", cl = NULL)

  # Loop over each lambda value
  lambda_errors <- foreach::foreach(l = 1:length(lambda)) %dorng% {

    # Set lambda value to 0
    lambda_value <- lambda[l]

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

      # Train model using lambda value on kth training set
      model_l <- lasso_classifier(L2.fix = L2.fe.form,
                                  L1.re = L1.re,
                                  data.train = data_train,
                                  lambda = lambda_value,
                                  model.family = binomial(link = "probit"),
                                  verbose = FALSE)

      # Use trained model to make predictions for kth validation set
      pred_l <- stats::predict(model_l, newdata = data.frame(data_valid))

      # Evaluate predictions based on loss function
      perform_l <- loss_function(pred = pred_l, data.valid = data_valid,
                                 loss.unit = loss.unit,
                                 loss.fun = loss.fun,
                                 y = y, L2.unit = L2.unit)
    })

    # Mean over loss functions
    k_errors <- dplyr::bind_rows(k_errors) %>%
      dplyr::group_by(measure) %>%
      dplyr::summarise(value = mean(value), .groups = "drop") %>%
      dplyr::mutate(lambda = lambda[l] )

  }

  # De-register cluster
  multicore(cores = cores, type = "close", cl = cl)

  # Function output
  return(lambda_errors)

}
