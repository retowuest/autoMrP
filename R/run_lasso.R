#' Apply lasso classifier to MrP.
#'
#' \code{run_lasso} is a wrapper function that applies the lasso classifier to
#' data provided by the user, evaluates prediction performance, and chooses the
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
#' @param lambda Lasso penalty parameter. A numeric \code{vector} of
#'   non-negative values or a \code{list} of two numeric vectors of equal size,
#'   with the first vector containing the step sizes by which the penalty
#'   parameter should increase and the second vector containing the upper
#'   thresholds of the intervals to which the step sizes apply. The penalty
#'   parameter controls the shrinkage of the context-level variables in the
#'   lasso model. Default is \code{list(c(0.1, 0.3, 1), c(1, 10, 10000))}.
#' @param n.iter Lasso number of iterations without improvement. Either
#'   \code{NULL} or an integer-valued scalar specifying the maximum number of
#'   iterations without performance improvement the algorithm runs before
#'   stopping. Default is \eqn{70}.
#' @param data Data for cross-validation. A \code{list} of \eqn{k}
#'   \code{data.frames}, one for each fold to be used in \eqn{k}-fold
#'   cross-validation.
#' @param verbose Verbose output. A logical argument indicating whether or not
#'   verbose output should be printed. Default is \code{TRUE}.
#' @return
#' @examples

lasso <- function(y, L1.x, L2.x, L2.unit, L2.reg,
                  n.iter = NULL,
                  loss.unit, loss.fun,
                  lambda, iterations.max,
                  data, verbose) {

  # Context-level fixed effects
  L2_fe <- paste(L2.x, collapse = " + ")
  L2_fe_form <- as.formula(paste(y, " ~ ", L2_fe, sep = ""))

  # Individual-level random effects as named list
  L1_re <- setNames(as.list(rep(c(~ 1),
                                times = length(c(L1.x, L2.unit, L2.reg)))),
                    c(L1.x, L2.unit, L2.reg))

  # Train and evaluate each model
  if (is.vector(lambda)) {
    # Add value of 0 to lambda if not already included
    if (!0 %in% lambda) {
      lambda <- c(0, lambda)
    }

    # Set lambda value to 0
    lambda_value <- lambda[1]

    # Print lambda value
    if (isTRUE(verbose)) {
      L <- length(lambda)
      cat(paste("Lasso: Running lambda w/ value ", lambda_value,
                " (lambda ", 1, " out of max. ",
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
                                  lambda = as.numeric(lambda_value),
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

    # Mean over all k folds
    best_error <- mean(unlist(k_errors))

    # Initialize lambda value associated with currently best error
    lambda_out <- lambda_value

    # Initialize counter for iterations since last performance improvement
    iter_since_improv <- 0

    # Loop over lambda values in lambda
    for (l in 2:length(lambda)) {
      # Set lambda value
      lambda_value <- lambda[l]

      # Print lambda value
      if (isTRUE(verbose)) {
        L <- length(lambda)
        cat(paste("Lasso: Running lambda w/ value ", lambda_value,
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
                                    lambda = lambda_value,
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

      # Mean over all k folds
      current_error <- mean(unlist(k_errors))

      # Check if current lambda value outperforms the lambda value that was
      # best so far
      if (current_error < best_error) {
        best_error <- current_error
        lambda_out <- lambda_value
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

    # Function output
    return(lambda_out)
  } else {
    # Set lambda value to 0
    lambda_value <- 0

    # Initialize counter for lambda
    lambda_no <- 1

    # Print lambda value
    if (isTRUE(verbose)) {
      cat(paste("Lasso: Running lambda w/ value ", lambda_value,
                " (lambda no. ", lambda_no, " -- no improvement evaluation)\n",
                sep = ""))
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
                                  lambda = lambda_value,
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

    # Mean over all k folds
    best_error <- mean(unlist(k_errors))

    # Initialize lambda value associated with currently best error
    lambda_out <- lambda_value

    # Initialize counter for iterations since last performance improvement
    iter_since_improv <- 0

    # Loop over lambda values in lambda
    while(lambda_value < dplyr::last(lambda[, 2])) {
      # Set lambda value
      lambda_value <- round(lambda_value, digits = 10)
      lambda_value <- lambda_value +
        lambda[, 1][which(lambda_value < lambda[, 2])[1]]

      # Update counter for lambda
      lambda_no <- lambda_no + 1

      # Print lambda value
      if (isTRUE(verbose)) {
        cat(paste("Lasso: Running lambda w/ value ", lambda_value,
                  " (lambda no. ", lambda_no, " -- iterations w/o improvement: ",
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

        # Train model using lambda value on kth training set
        model_l <- lasso_classifier(L2.fix = L2_fe_form,
                                    L1.re = L1_re,
                                    data.train = data_train,
                                    lambda = lambda_value,
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

      # Mean over all k folds
      current_error <- mean(unlist(k_errors))

      # Check if current lambda value outperforms the lambda value that was
      # best so far
      if (current_error < best_error) {
        best_error <- current_error
        lambda_out <- lambda_value
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

    # Function output
    return(lambda_out)
  }
}
