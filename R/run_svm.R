#' Apply support vector machine classifier to MrP.
#'
#' \code{run_svm} is a wrapper function that applies the support vector machine
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
#' @param loss.fun Loss function. A character-valued scalar indicating whether
#'   prediction loss should be measured by the mean squared error (\code{MSE})
#'   or the mean absolute error (\code{MAE}). Default is \code{MSE}.
#' @param kernel SVM kernel. A character-valued scalar specifying the kernel to
#'   be used by SVM. The possible values are \code{linear}, \code{polynomial},
#'   \code{radial}, and \code{sigmoid}. Default is \code{radial}.
#' @param gamma SVM kernel parameter. A numeric vector whose values specify the
#'   gamma parameter in the SVM kernel. This parameter is needed for all kernel
#'   types except linear. Default is a sequence with minimum = 1e-5, maximum =
#'   1e-1, and length = 20 that is equally spaced on the log-scale.
#' @param cost SVM cost parameter. A numeric vector whose values specify the
#'   cost of constraints violation in SVM. Default is a sequence with minimum =
#'   0.5, maximum = 10, and length = 5 that is equally spaced on the log-scale.
#' @param data Data for cross-validation. A \code{list} of \eqn{k}
#'   \code{data.frames}, one for each fold to be used in \eqn{k}-fold
#'   cross-validation.
#'
#' @return The support vector machine tuned parameters. A list.

run_svm <- function(y, L1.x, L2.x, L2.eval.unit, L2.unit, L2.reg,
                    kernel = "radial", loss.fun, loss.unit, gamma,
                    cost, data, verbose, cores) {

  # Create model formula
  x <- paste(c(L1.x, L2.x, L2.unit, L2.reg), collapse = " + ")
  form <- as.formula(paste(y, " ~ ", x, sep = ""))

  # Default Gamma values
  if( is.null(gamma) ){
    # SVM Gamma values
    gamma <- log_spaced(min = 1e-5, 1e-1, n = 20)
  }

  # Default Cost values
  if ( is.null(cost) ){
    cost <- log_spaced(min = 0.5, max = 10, n = 5)
  }

  # tuning parameter grid
  svm_grid <- expand.grid(gamma, cost, kernel)
  names(svm_grid) <- c("gamma", "cost", "kernel")

  # prallel tuning if cores > 1
  if( cores > 1 ){

    # Train all models in parallel
    grid_cells <- run_svm_mc(
      verbose = verbose,
      svm.grid = svm_grid,
      data = data,
      L2.eval.unit = L2.eval.unit,
      loss.unit = loss.unit,
      loss.fun = loss.fun,
      y = y,
      L1.x = L1.x,
      L2.x = L2.x,
      L2.unit = L2.unit,
      L2.reg = L2.reg,
      form = form,
      cores = cores)

  # Train all models sequentially
  } else {
    # loop over tuning grid
    grid_cells <- apply(svm_grid, 1, function(g) {

      # Set tuning parameters
      gamma_value <- as.numeric(g["gamma"])
      cost_value <- as.numeric(g["cost"])
      kernel_value <- as.character(g[["kernel"]])

      # Loop over each fold
      k_errors <- lapply(seq_along(data), function(k) {

        # Split data in training and validation sets and factorize DV
        data_train <- dplyr::bind_rows(data[-k]) %>%
          dplyr::mutate_at(.vars = y, as.factor) %>%
          dplyr::select( dplyr::all_of(c(y, L1.x, L2.x, L2.eval.unit, L2.reg)) ) %>%
          tidyr::drop_na()

        data_valid <- dplyr::bind_rows(data[k]) %>%
          dplyr::mutate_at(.vars = y, as.factor) %>%
          dplyr::select( dplyr::all_of(c(y, L1.x, L2.x, L2.eval.unit, L2.reg)) ) %>%
          tidyr::drop_na()

        # Svm classifier
        model_l <- svm_classifier(
          form = form,
          data = data_train,
          kernel = kernel_value,
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
                                   y = y, L2.unit = L2.eval.unit)
      })

      # Mean over loss functions
      k_errors <- dplyr::bind_rows(k_errors) %>%
        dplyr::group_by(measure) %>%
        dplyr::summarise(value = mean(value), .groups = "drop") %>%
        dplyr::mutate(gamma = gamma_value,
                      cost = cost_value,
                      kernel = kernel_value)

    })
  }

  # Extract best tuning parameters
  grid_cells <- dplyr::bind_rows(grid_cells)
  best_params <- dplyr::slice(loss_score_ranking(score = grid_cells, loss.fun = loss.fun), 1)

  out <- list(gamma =  dplyr::pull(.data = best_params, var = gamma),
              cost = dplyr::pull(.data = best_params, var = cost),
              kernel = dplyr::pull(.data = best_params, var = kernel))

  # Function output
  return(out)

}

################################################################################
#                     Multicore tuning for svm                                 #
################################################################################
#' SVM multicore tuning.
#'
#' \code{run_svm_mc} is called from within \code{run_svm}. It tunes using
#' multiple cores.
#'
#' @inheritParams run_svm
#' @param form The model formula. A formula object.
#' @param svm.grid The hyper-parameter search grid. A matrix of all
#'   hyper-parameter combinations.
#' @return The cross-validation errors for all models. A list.

run_svm_mc <- function(y, L1.x, L2.x, L2.eval.unit, L2.unit, L2.reg, form,
                       loss.unit, loss.fun, data, cores, svm.grid, verbose){

  # Binding for global variables
  g <- NULL
  `%>%` <- dplyr::`%>%`

  # Register cores
  cl <- multicore(cores = cores, type = "open", cl = NULL)

  # Train and evaluate each model
  grid_cells <- foreach::foreach(g = 1:nrow(svm.grid), .packages = 'autoMrP') %dorng% {

    # Set tuning parameters
    gamma_value <- as.numeric(svm.grid[g, "gamma"])
    cost_value <- as.numeric(svm.grid[g, "cost"])
    kernel_value <- svm.grid[g, "kernel"]

    # Loop over each fold
    k_errors <- lapply(seq_along(data), function(k) {

      # Split data in training and validation sets and factorize DV
      data_train <- dplyr::bind_rows(data[-k]) %>%
        dplyr::mutate_at(.vars = y, as.factor) %>%
        dplyr::select( dplyr::all_of(c(y, L1.x, L2.x, L2.eval.unit, L2.reg)) ) %>%
        tidyr::drop_na()

      data_valid <- dplyr::bind_rows(data[k]) %>%
        dplyr::mutate_at(.vars = y, as.factor) %>%
        dplyr::select( dplyr::all_of(c(y, L1.x, L2.x, L2.eval.unit, L2.reg)) ) %>%
        tidyr::drop_na()

      # Svm classifier
      model_l <- svm_classifier(
        form = form,
        data = data_train,
        kernel = kernel_value,
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
                                 y = y, L2.unit = L2.eval.unit)
    })

    # Mean over loss functions
    k_errors <- dplyr::bind_rows(k_errors) %>%
      dplyr::group_by(measure) %>%
      dplyr::summarise(value = mean(value), .groups = "drop") %>%
      dplyr::mutate(gamma = gamma_value,
                    cost = cost_value,
                    kernel = kernel_value)
  }

  # De-register cluster
  multicore(cores = cores, type = "close", cl = cl)

  # Function output
  return(grid_cells)
}
