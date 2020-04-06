#' Apply lasso classifier to MrP.
#'
#' \code{lasso} is a wrapper function that applies the lasso classifier to data
#' provided by the user, evaluates prediction performance, and chooses the
#' best-performing model.
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
#' @param lambda.set Set of tuning parameters. Lambda is the penalty parameter
#'   that controls the shrinkage of fixed effects. Either a numeric vector of
#'   lambda values or a data.frame with two columns, the first containing the
#'   size by which lambda should increase and the second the upper threshold of
#'   the interval of lambdas to which the step size applies. Default is
#'   data.frame(step_size = c(0.1, 0.3, 1), threshold = c(1, 10, 10000)).
#' @param iterations.max Stopping rule. A numeric scalar specifying the
#'   maximum number of iterations without performance improvement the
#'   algorithm runs before stopping. Default is 60
#' @param data Data for cross-validation. A list of k data.frames, one for
#'   each fold used in k-fold cross-validation.
#' @param verbose Verbose output. A logical vector indicating whether or not
#'   verbose output should be printed.
#' @return
#' @examples

lasso <- function(y, L1.x, L2.x, L2.unit, L2.reg,
                  loss.unit, loss.measure,
                  lambda.set, iterations.max,
                  data, verbose) {

  # Context-level fixed effects
  L2_fe <- paste(L2.x, collapse = " + ")
  L2_fe_form <- as.formula(paste(y, " ~ ", L2_fe, sep = ""))

  # Individual-level random effects as named list
  L1_re <- setNames(as.list(rep(c(~ 1), times = length(c(L1.x, L2.unit, L2.reg)))),
                    c(L1.x, L2.unit, L2.reg))

  # Train and evaluate each model
  if (is.vector(lambda.set)) {
    # Add value of 0 to lambda.set if not already included
    if (!0 %in% lambda.set) {
      lambda.set <- c(0, lambda.set)
    }

    # Set lambda to 0
    lambda <- lambda.set[1]

    # Print lambda
    if (isTRUE(verbose == TRUE)) {
      L <- length(lambda.set)
      cat(paste("Lasso: Running lambda w/ value ", lambda,
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

      # Train model using lambda on kth training set
      model_l <- lasso_classifier(L2.fix = L2_fe_form,
                                  L1.re = L1_re,
                                  data.train = data_train,
                                  lambda = lambda,
                                  model.family = binomial(link = "probit"),
                                  verbose = verbose)

      # Use trained model to make predictions for kth validation set
      pred_l <- stats::predict(model_l, newdata = data.frame(data_valid))

      # Evaluate predictions based on loss function
      perform_l <- loss_function(pred = pred_l, data.valid = data_valid,
                                 loss.unit = loss.unit,
                                 loss.measure = loss.measure,
                                 y = y, L2.unit = L2.unit)
    })

    # Mean over all k folds
    best_error <- mean(unlist(k_errors))

    # Initialize lambda associated with currently best error
    out <- lambda

    # Initialize counter for iterations since last performance improvement
    iter_since_improv <- 0

    # Loop over lambda values in lambda.set
    for (l in 2:length(lambda.set)) {
      # Set lambda value
      lambda <- lambda.set[l]

      # Print lambda
      if (isTRUE(verbose == TRUE)) {
        L <- length(lambda.set)
        cat(paste("Lasso: Running lambda w/ value ", lambda,
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

        # Train model using lambda on kth training set
        model_l <- lasso_classifier(L2.fix = L2_fe_form,
                                    L1.re = L1_re,
                                    data.train = data_train,
                                    lambda = lambda,
                                    model.family = binomial(link = "probit"),
                                    verbose = verbose)

        # Use trained model to make predictions for kth validation set
        pred_l <- stats::predict(model_l, newdata = data.frame(data_valid))

        # Evaluate predictions based on loss function
        perform_l <- loss_function(pred = pred_l, data.valid = data_valid,
                                   loss.unit = loss.unit,
                                   loss.measure = loss.measure,
                                   y = y, L2.unit = L2.unit)
      })

      # Mean over all k folds
      current_error <- mean(unlist(k_errors))

      # Check if current lambda outperforms lambda that was best so far
      if (current_error < best_error) {
        best_error <- current_error
        out <- lambda
        iter_since_improv <- 0
      } else {
        iter_since_improv <- iter_since_improv + 1
      }

      # Break loop if maximum number of iterations without performance
      # improvement is reached
      if (!is.null(iterations.max) & iter_since_improv > iterations.max) {
        break
      }
    }

    # Function output
    return(out)
  } else {
    # Set lambda to 0
    lambda <- 0

    # Initialize counter for lambda
    lambda_no <- 1

    # Print lambda
    if (isTRUE(verbose == TRUE)) {
      cat(paste("Lasso: Running lambda w/ value ", lambda,
                " (lambda no. ", lambda_no, " -- no improvement evaluation)\n", sep = ""))
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
      model_l <- lasso_classifier(L2.fix = L2_fe_form,
                                  L1.re = L1_re,
                                  data.train = data_train,
                                  lambda = lambda,
                                  model.family = binomial(link = "probit"),
                                  verbose = verbose)

      # Use trained model to make predictions for kth validation set
      pred_l <- stats::predict(model_l, newdata = data.frame(data_valid))

      # Evaluate predictions based on loss function
      perform_l <- loss_function(pred = pred_l, data.valid = data_valid,
                                 loss.unit = loss.unit,
                                 loss.measure = loss.measure,
                                 y = y, L2.unit = L2.unit)
    })

    # Mean over all k folds
    best_error <- mean(unlist(k_errors))

    # Initialize lambda associated with currently best error
    out <- lambda

    # Initialize counter for iterations since last performance improvement
    iter_since_improv <- 0

    # Loop over lambda values in lambda.set
    while(lambda < dplyr::last(lambda.set[, 2])) {
      # Set lambda value
      lambda <- round(lambda, digits = 10)
      lambda <- lambda + lambda.set[, 1][which(lambda < lambda.set[, 2])[1]]

      # Update counter for lambda
      lambda_no <- lambda_no + 1

      # Print lambda
      if (isTRUE(verbose == TRUE)) {
        cat(paste("Lasso: Running lambda w/ value ", lambda,
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

        # Train model using lambda on kth training set
        model_l <- lasso_classifier(L2.fix = L2_fe_form,
                                    L1.re = L1_re,
                                    data.train = data_train,
                                    lambda = lambda,
                                    model.family = binomial(link = "probit"),
                                    verbose = verbose)

        # Use trained model to make predictions for kth validation set
        pred_l <- stats::predict(model_l, newdata = data.frame(data_valid))

        # Evaluate predictions based on loss function
        perform_l <- loss_function(pred = pred_l, data.valid = data_valid,
                                   loss.unit = loss.unit,
                                   loss.measure = loss.measure,
                                   y = y, L2.unit = L2.unit)
      })

      # Mean over all k folds
      current_error <- mean(unlist(k_errors))

      # Check if current lambda outperforms lambda that was best so far
      if (current_error < best_error) {
        best_error <- current_error
        out <- lambda
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
    return(out)
  }
}
