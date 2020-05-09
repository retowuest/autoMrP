#' Apply support vector machine classifier to MrP.
#'
#' \code{run_svm} is a wrapper function that applies the support vector machine
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
#' @param loss.fun Loss function. A character-valued scalar indicating whether
#'   prediction loss should be measured by the mean squared error (\code{MSE})
#'   or the mean absolute error (\code{MAE}). Default is \code{MSE}.
#' @param kernel SVM kernel. A character-valued scalar specifying the kernel to
#'   be used by SVM. The possible values are \code{linear}, \code{polynomial},
#'   \code{radial}, and \code{sigmoid}. Default is \code{radial}.
#' @param loss.fun SVM loss function. If \code{NULL}, then SVM uses the
#'   misclassification error to measure the loss of categorical predictions and
#'   the mean squared error to measure the loss of numeric predictions. Default
#'   is \code{NULL}.
#' @param loss.unit Loss function unit. A character-valued scalar indicating
#'   whether performance loss should be evaluated at the level of individual
#'   respondents (\code{individuals}) or geographic units (\code{L2 units}).
#'   Default is \code{individuals}.
#' @param gamma SVM kernel parameter. A numeric vector whose values specify
#'   the gamma parameter in the SVM kernel. This parameter is needed for all
#'   kernel types except linear. Default is
#'   \eqn{c(0.3, 0.5, 0.55, 0.6, 0.65, 0.7, 0.8, 0.9, 1, 2, 3, 4)}.
#' @param cost SVM cost parameter. A numeric vector whose values specify the
#'   cost of constraints violation in SVM. Default is \eqn{c(1, 10)}.
#' @param data Data for cross-validation. A \code{list} of \eqn{k}
#'   \code{data.frames}, one for each fold to be used in \eqn{k}-fold
#'   cross-validation.
#' @param verbose Verbose output. A logical argument indicating whether or not
#'   verbose output should be printed. Default is \code{TRUE}.
#'
#' @return The support vector machine tuned parameters. A list.
#' @examples \dontrun{
#' # Create list of cross-validation folds
#' cv_folds <- list(
#'   `1` = survey_item[1:200, ],
#'   `2` = survey_item[201:400, ],
#'   `3` = survey_item[401:1500, ])
#'
#' # Run svm classifier
#' m <- run_svm(
#'   y = "YES",
#'   L1.x = c("L1x1", "L1x2"),
#'   L2.x = c("L2.x1", "L2.x2"),
#'   L2.unit = "state",
#'   L2.reg = "region",
#'   kernel = "radial",
#'   loss.fun = "MSE",
#'   loss.unit = "individuals",
#'   gamma = c(0.3, 0.1),
#'   cost = c(1, 50),
#'   data = cv_folds,
#'   verbose = TRUE)
#' }

run_svm <- function(y, L1.x, L2.x, L2.unit, L2.reg,
                    kernel = "radial", loss.fun,
                    loss.unit, gamma, cost, data,
                    verbose) {

  # Create model formula
  x <- paste(c(L1.x, L2.x, L2.unit, L2.reg), collapse = " + ")
  form <- as.formula(paste(y, " ~ ", x, sep = ""))

  # tuning parameter grid
  svm_grid <- expand.grid(gamma, cost)
  names(svm_grid) <- c("gamma", "cost")

  # loop over tuning grid
  grid_cells <- apply(svm_grid, 1, function(g) {

    # Set tuning parameters
    gamma_value <- as.numeric(g["gamma"])
    cost_value <- as.numeric(g["cost"])

    # Loop over each fold
    k_errors <- lapply(seq_along(data), function(k) {

      # Split data in training and validation sets and factorize DV
      data_train <- dplyr::bind_rows(data[-k]) %>%
        dplyr::mutate_at(.vars = y, as.factor)
      data_valid <- dplyr::bind_rows(data[k]) %>%
        dplyr::mutate_at(.vars = y, as.factor)

      # Svm classifier
      model_l <- svm_classifier(
        form = form,
        data = data_train,
        kernel = kernel,
        type = "C-classification",
        probability = TRUE,
        svm.gamma = gamma_value,
        svm.cost = cost_value,
        verbose = verbose
      )

      # Use trained model to make predictions for kth validation set
      pred_l <- predict(model_l, newdata = data.frame(data_valid),
                        probability = TRUE)
      pred_l <- as.numeric(attr(pred_l, "probabilities")[, "1"])

      # Transform factor DV to numeric for loss function
      data_valid <- data_valid %>%
        dplyr::mutate_at(.vars = y, function(x) as.numeric(levels(x))[x])

      # Evaluate predictions based on loss function
      perform_l <- loss_function(pred = pred_l, data.valid = data_valid,
                                 loss.unit = loss.unit,
                                 loss.fun = loss.fun,
                                 y = y, L2.unit = L2.unit)
    })

    # Mean over all k folds
    best_error <- mean(unlist(k_errors))

  })

  # Extract best tuning parameters
  out <- list(gamma = svm_grid[which.min(grid_cells), "gamma"],
              cost = svm_grid[which.min(grid_cells), "cost"])

  # Function output
  return(out)

}
