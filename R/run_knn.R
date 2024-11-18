#' Apply KNN classifier to MrP.
#'
#' \code{run_knn} is a wrapper function that applies the KNN classifier to data
#' provided by the user, evaluates prediction performance, and chooses the
#' best-performing model.
#'
#' @inheritParams auto_MrP
#' @param knn.k.max KNN maximum number of neighbors. An positive integer-valued
#'   scalar specifying the maximum number of neighbors to be considered in the
#'   KNN model if \code{knn.k} is \code{NULL}. If \code{knn.k.max} is specified
#'   and \code{knn.k} is \code{NULL}, then the number of neighbors considered is
#'   the sequence \code{1:knn.k.max}. Default is \code{11}.
#' @param knn.k KNN number of neighbors. A \code{vector} of positive integers
#'   specifying the number of neighbors to be considered in the KNN model. If
#'   not \code{NULL}, \code{knn.k} takes precedence over \code{knn.k.max}.
#'   Default is \code{NULL}.
#' @param knn.kernel KNN kernel. A character-valued scalar specifying the kernel
#'   to be used in the KNN model. The possible values are \code{rectangular}
#'   (which is standard unweighted KNN), \code{triangular}, \code{epanechnikov}
#'   (or beta(2,2)), \code{biweight} (or beta(3,3)), \code{triweight} (or
#'   beta(4,4)), \code{cos}, \code{inv}, \code{gaussian}, and \code{optimal}.
#'   Default is \code{optimal}.
#' @param data Data for cross-validation. A \code{list} of \eqn{k}
#'   \code{data.frames}, one for each fold to be used in \eqn{k}-fold
#'   cross-validation.
#'
#' @return The tuned \code{knn.k} parameter. An integer-valued scalar.

run_knn <- function(
  y, L1.x, L2.x, L2.unit, L2.reg, loss.unit, loss.fun, knn.k.max, knn.k,
  knn.kernel, data, verbose, cores
) {

  # Create model formula
  x <- paste(c(L1.x, L2.x, L2.unit, L2.reg), collapse = " + ")
  form <- as.formula(paste(y, " ~ ", x, sep = ""))

  # ks search grid
  if (is.null(knn.k)) {
    ks <- 1:knn.k.max
  } else {
    ks <- sort(knn.k)
  }

  # Register cores
  cl <- multicore(cores = cores, type = "open", cl = NULL)

  # Loop over each value in ks
  knn_k_errors <- foreach::foreach(
    ki = seq_along(ks), .packages = "autoMrP",
    .export = c("deep_mrp_classifier", "loss_function", "contr.dummy"),
  ) %dorng% {

    # Set value for k
    knn_k_value <- ks[ki]

    # Loop over each fold
    k_errors <- lapply(seq_along(data), function(k) {

      # Split data in training and validation sets
      data_train <- dplyr::bind_rows(data[-k])
      data_valid <- dplyr::bind_rows(data[k])

      # Determine type of dependent variable
      if (
        data_train %>%
          dplyr::pull(!!y) %>%
          unique() %>%
          length() == 2
      ) {
        dv_type <- "binary"
      }

      if (dv_type == "binary") {
        # Convert the column named "YES" to a factor
        data_train <- data_train %>%
          dplyr::mutate(!!rlang::sym(y) := as.factor(!!rlang::sym(y)))
        data_valid <- data_train %>%
          dplyr::mutate(!!rlang::sym(y) := as.factor(!!rlang::sym(y)))

      }

      # Train model with kith value on kth training set and use trained model
      # to make predictions for kth validation set
      model_ki <- knn_classifier(
        y = y,
        form = form,
        data.train = data_train,
        data.valid = data_valid,
        knn.k.value = knn_k_value,
        knn.kernel = knn.kernel,
        verbose = verbose
      )

      # Get predictions for kth validation set
      pred_ki <- model_ki$fitted.values

      if (dv_type == "binary") {
        # class predictions to numeric
        pred_ki <- as.numeric(as.character(pred_ki))
        # y column in dv_valid back to numeric
        data_valid <- data_valid %>%
          dplyr::mutate(
            !!rlang::sym(y) := as.numeric(as.character(!!rlang::sym(y)))
          )
      }

      # Evaluate predictions based on loss function
      perform_ki <- loss_function(
        pred = pred_ki,
        data.valid = data_valid,
        loss.unit = loss.unit,
        loss.fun = loss.fun,
        y = y,
        L2.unit = L2.unit
      )
    })

    # Mean over loss functions
    k_errors <- dplyr::bind_rows(k_errors) %>%
      dplyr::group_by(measure) %>%
      dplyr::summarise(value = mean(value), .groups = "drop") %>%
      dplyr::mutate(knn.k.value = ks[ki])

    return(k_errors)
  }

  # De-register cluster
  multicore(cores = cores, type = "close", cl = cl)

  # Extract best tuning parameters
  knn_k_errors <- dplyr::bind_rows(knn_k_errors)
  best_params <- dplyr::slice(
    loss_score_ranking(score = knn_k_errors, loss.fun = loss.fun), 1
  )

  # Choose best-performing model
  out <- dplyr::pull(.data = best_params, var = knn.k.value)

  return(out)
}

################################################################################
#                     Multicore tuning for KNN                                 #
################################################################################
#' KNN multicore tuning.
#'
#' \code{run_knn_mc} is called from within \code{run_knn}. It tunes using
#' multiple cores.
#'
#' @inheritParams run_knn
#' @param form The model formula. A formula object.
#' @param ks The hyper-parameter search grid of all values for parameter k. A
#'   vector.
#' @return The cross-validation errors for all models. A list.

run_knn_mc <- function(
  y, L1.x, L2.x, L2.unit, L2.reg, form, loss.unit, loss.fun, data, cores,
  ks, knn.kernel, verbose
) {

  # Binding for global variables
  g <- NULL
  `%>%` <- dplyr::`%>%`



  # Function output
  return(knn_k_errors)
}
