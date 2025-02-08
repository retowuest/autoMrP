#' Apply deep mrp to the best subset classifier to MrP.
#'
#' \code{run_deep_bs} is a wrapper function that applies the bestsubset
#' classifier to a list of models provided by the user, evaluates the models'
#' prediction performance, and chooses the best-performing model. It differs
#' from \code{run_best_subset} in that it includes L1.x interactions.
#'
#' @inheritParams auto_MrP
#' @param data Data for cross-validation. A \code{list} of \eqn{k}
#'   \code{data.frames}, one for each fold to be used in \eqn{k}-fold
#'   cross-validation.
#' @return A model formula of the winning best subset classifier model.

run_deep_bs <- function(
  y, L1.x, L2.x, L2.unit, L2.reg, loss.unit, loss.fun, deep.splines, data,
  k.folds, verbose, cores
) {

  # Determine type of dependent variable
  if (
    data[[1]] %>%
      dplyr::pull(!!y) %>%
      unique() %>%
      length() == 2
  ) {
    dv_type <- "binary"
  } else {
    dv_type <- "linear"
  }

  # List of all models to be evaluated
  models <- model_list(
    y = y,
    L1.x = L1.x,
    L2.x = L2.x,
    L2.unit = L2.unit,
    L2.reg = L2.reg
  )

  # no nesting with deep interactions
  if (!is.null(L2.reg) && !is.null(L2.unit)) {
    models <- lapply(models, function(x) {
      # model formula to character
      m_form <- as.character(x)
      # replace (1 | region/state) with (1 | region) + (1 | state)
      m_form <- stringr::str_replace_all(
        string = m_form,
        pattern = sprintf("\\(1 \\| %s/%s\\)", L2.reg, L2.unit),
        replacement = sprintf("\\(1 | %s\\) + \\(1 | %s\\)", L2.unit, L2.reg)
      )
      # character to formula
      m_form <- as.formula(sprintf("%s%s%s", m_form[2], m_form[1], m_form[3]))
    })
  }

  # add interactions to the models
  models <- lapply(models, function(x) {

    # get all level 1 variables in the current model
    c_l1_x <- x %>%
      as.character() %>%
      .[3] %>%
      stringr::str_detect(pattern = L1.x)
    c_l1_x <- L1.x[c_l1_x]

    # generate all interactions of L1.x
    l1_comb <- unlist(lapply(2:length(c_l1_x), function(x) {
      apply(combn(L1.x, x), 2, paste, collapse = ".")
    }))

    # generate all interactions of L1.x with L2.unit
    l1_state <- paste(L1.x, L2.unit, sep = ".")

    # generate all interactions of L1.x with L2.reg
    if (!is.null(L2.reg)) {
      l1_region <- paste(L1.x, L2.reg, sep = ".")
    } else {
      l1_region <- NULL
    }

    # interactions
    add_interactions <- paste0(
      # interactions of L1x
      paste("(1 | ", l1_comb, ")", collapse = " + "), " + ",
      # interactions of L1x with L2.unit
      paste("(1 | ", l1_state, ")", collapse = " + "), " + ",
      # interactions of L1x with L2.reg
      if (any(!is.null(l1_region))) {
        paste("(1 | ", l1_region, ")", collapse = " + ")
      }
    )

    # remove trailing " + " from interactions
    add_interactions <- stringr::str_extract(
      string = add_interactions,
      pattern = "^.*\\)"
    )

    # character to formula
    add_interactions <- as.formula(paste("~ . +", add_interactions))

    # update formula with interactions
    x <- update(x, add_interactions)

    # add splines to context level variables
    if (deep.splines) {

      # formula to character
      char_form <- as.character(x)
      char_form <- sprintf("%s %s %s", char_form[2], char_form[1], char_form[3])

      # get all context level variables in the current model
      c_l2_x <- char_form %>%
        stringr::str_extract_all(pattern = "L2\\.x\\d+") %>%
        unlist()

      # replace in string
      for (i in seq_along(c_l2_x)) {
        char_form <- stringr::str_replace(
          string = char_form,
          pattern = c_l2_x[i],
          replacement = sprintf("v_s(%s)", c_l2_x[i])
        )
      }

      # character to formula
      x <- as.formula(char_form)

    }

    return(x)
  })

  # Register cores
  cl <- multicore(cores = cores, type = "open", cl = NULL)

  # loop over models
  m_errors <- foreach::foreach(
    m = seq_along(models), .packages = "autoMrP"
  ) %dorng% {

    `%>%` <- magrittr::`%>%`

    k_errors <- lapply(seq_len(k.folds), function(k) {

      # Split data in training and validation sets
      data_train <- dplyr::bind_rows(data[-k])
      data_valid <- dplyr::bind_rows(data[k])

      # Train mth model on kth training set
      model_m <- deep_mrp_classifier(
        form = models[[m]],
        y = y,
        data = data_train,
        verbose = TRUE
      )

      # predictions based on DV type (binary or continuous)
      if (dv_type == "binary") {
        # use trained model to make predictions for kth validation set
        pred_m <- vglmer::predict_MAVB(
          samples = 1000,
          model_m,
          newdata = data_valid,
          allow_missing_levels = TRUE
        )[["mean"]]

        # convert to response probabilities
        pred_m <- stats::plogis(pred_m)

      } else if (dv_type == "linear") {
        # Use trained model to make predictions for kth validation set
        pred_m <- predict(
          samples = 1000,
          object = model_m,
          newdata = data_valid,
          allow_missing_levels = TRUE
        )[["mean"]]
      }

      # evaluate predictions based on loss function
      perform_m <- loss_function(
        pred = pred_m,
        data.valid = data_valid,
        loss.unit = loss.unit,
        loss.fun = loss.fun,
        y = y,
        L2.unit = L2.unit
      )

      return(perform_m)
    })

    # Mean over loss functions
    k_errors <- dplyr::bind_rows(k_errors) %>%
      dplyr::group_by(measure) %>%
      dplyr::summarise(value = mean(value), .groups = "drop") %>%
      dplyr::mutate(model = m)

    return(k_errors)
  }

  # De-register cluster
  multicore(cores = cores, type = "close", cl = cl)

  # Extract best tuning parameters
  grid_cells <- dplyr::bind_rows(m_errors)
  best_params <- dplyr::slice(
    loss_score_ranking(
      score = grid_cells,
      loss.fun = loss.fun
    ), 1
  )

  # Choose best-performing model
  out <- models[[dplyr::pull(.data = best_params, var = model)]]


  # Function output
  return(out)

}