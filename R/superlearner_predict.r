#' Predict individual-level survey respones for new data.
#'
#' Predicts individual-level responses for new data based on a trained auto_MrP
#' object. This is uesful for combining autoMrP predictions with predictions
#' from other models into a superlearner.
#' @inheritParams auto_MrP
#' @param obj An autoMrP object.
#' @param data.test A new data frame with the same columns as the original data
#'  used to train autoMrP. Note: This data should not have been used to train
#'  the autoMrP model.
#' @return The individual-level predictions. A tibble. It contains the
#' individual-level predictions from all estimated classifiers, ebma, and
#'  stacking.
#' @export

superlearner_predict <- function(
  obj, data.test, L2.x.scale = TRUE
) {

  #---------------------------------------------------------------------------
  # pre-process data
  #---------------------------------------------------------------------------

  # Scale context-level variables in data.train and data.test
  if (isTRUE(L2.x.scale) && all(attr(obj, "evaluated_args")$L2.x != "")) {

    # scale context-level variables in data.train
    data.train <- dplyr::mutate_at(
      .tbl = attr(obj, "evaluated_args")$survey,
      .vars = attr(obj, "evaluated_args")$L2.x,
      .funs = function(x) {
        base::as.numeric(base::scale(x = x, center = TRUE, scale = TRUE))
      }
    )

    # scale context-level variables in data.test
    data.test <- dplyr::mutate_at(
      .tbl = data.test,
      .vars = attr(obj, "evaluated_args")$L2.x,
      .funs = function(x) {
        base::as.numeric(base::scale(x = x, center = TRUE, scale = TRUE))
      }
    )
  }

  # If not provided, compute the principal components of context-level variables
  if (
    is.null(attr(obj, "evaluated_args")$pcs) &&
      !is.null(attr(obj, "evaluated_args")$L2.x)
  ) {

    # Determine context-level covariates whose principal components are to be
    # computed
    if (is.null(attr(obj, "evaluated_args")$pca.L2.x)) {
      pca.L2.x <- attr(obj, "evaluated_args")$L2.x
    }

    # Compute principal components for data.train
    pca_out <- stats::prcomp(
      data.train[, pca.L2.x],
      retx = TRUE,
      center = TRUE,
      scale. = TRUE,
      tol = NULL
    )

    # Add PCs to survey data
    data.train <- data.train %>%
      dplyr::bind_cols(as.data.frame(pca_out$x))

    # Add PCs to census data
    data.test <- data.test %>%
      dplyr::bind_cols(
        as.data.frame(predict(object = pca_out, newdata = data.test))
      )
  }

  # Convert survey and census data to tibble
  data.train <- tibble::as_tibble(x = data.train)
  data.test <- tibble::as_tibble(x = data.test)

  # error if the dependent variable is not numeric
  if (
    !data.test %>%
      dplyr::pull(!!attr(obj, "evaluated_args")$y) %>%
      is.numeric(.)
  ) {
    stop("The dependent variable must be numeric.")
  }

  # type of dependent variable (binary or continuous)
  if (
    data.test %>%
      dplyr::pull(!!attr(obj, "evaluated_args")$y) %>%
      unique() %>%
      length() == 2
  ) {
    dv_type <- "binary"
  } else {
    dv_type <- "continuous"
  }

  # add deep interactions to newdata
  if (
    any(
      attr(obj, "evaluated_args")$best.subset.deep,
      attr(obj, "evaluated_args")$pca.deep,
      attr(obj, "evaluated_args")$deep.mrp
    )
  ) {

    # generate all interactions of L1.x
    l1_comb <- unlist(
      lapply(2:length(attr(obj, "evaluated_args")$L1.x), function(x) {
        apply(
          combn(attr(obj, "evaluated_args")$L1.x, x), 2, paste, collapse = "."
        )
      })
    )

    # generate all interactions of L1.x with L2.unit
    l1_state <- paste(
      attr(obj, "evaluated_args")$L1.x,
      attr(obj, "evaluated_args")$L2.unit, sep = "."
    )

    # generate all interactions of L1.x with L2.reg
    if (!is.null(attr(obj, "evaluated_args")$L2.reg)) {
      l1_region <- paste(
        attr(obj, "evaluated_args")$L1.x,
        attr(obj, "evaluated_args")$L2.reg, sep = "."
      )
    } else {
      l1_region <- NULL
    }

    # add the interactions to the data
    all_interactions <- c(l1_comb, l1_state, l1_region)

    # loop over all interactions for the data.train object
    x_data <- lapply(all_interactions, function(x) {

      # break down interaction components
      y <- stringr::str_extract(
        string = x,
        pattern = stringr::fixed(pattern = names(data.train))
      ) %>%
        .[!is.na(.)]

      # take each column of data and combine its values into a single string
      df_x <- data.train %>%
        dplyr::select({{y}}) %>%
        dplyr::rowwise() %>%
        dplyr::mutate({{x}} := paste(dplyr::c_across(
          dplyr::everything()
        ), collapse = "-")) %>%
        dplyr::ungroup() %>%
        dplyr::select(ncol(.))

      return(df_x)
    }) %>%
      dplyr::bind_cols()

    # combine survey and interactions
    data.train <- dplyr::bind_cols(data.train, x_data)
    rm(x_data)

    # loop over all interactions for the data.test object
    x_data <- lapply(all_interactions, function(x) {

      # break down interaction components
      y <- stringr::str_extract(
        string = x,
        pattern = stringr::fixed(pattern = names(data.test))
      ) %>%
        .[!is.na(.)]

      # take each column of data and combine its values into a single string
      df_x <- data.test %>%
        dplyr::select({{y}}) %>%
        dplyr::rowwise() %>%
        dplyr::mutate({{x}} := paste(dplyr::c_across(
          dplyr::everything()
        ), collapse = "-")) %>%
        dplyr::ungroup() %>%
        dplyr::select(ncol(.))

      return(df_x)
    }) %>%
      dplyr::bind_cols()

    # combine survey and interactions
    data.test <- dplyr::bind_cols(data.test, x_data)
    rm(x_data)
  } # end of add_deep_data

  #---------------------------------------------------------------------------
  # fit models to data.train and make predictions on data.test
  #---------------------------------------------------------------------------

  # best subset predictions
  if (attr(obj, "evaluated_args")$best.subset) {

    # Deep MrP
    if (attr(obj, "evaluated_args")$best.subset.deep) {

      # Train mth model on kth training set
      model_bs <- deep_mrp_classifier(
        form = attr(obj, "tuned_hyperparameters")$best_subset$form,
        y = attr(obj, "evaluated_args")$y,
        data = data.train,
        verbose = FALSE
      )

      # best subset predictions
      preds_bs <- vglmer::predict_MAVB(
        samples = 1000,
        model_bs,
        newdata = data.test,
        allow_missing_levels = TRUE
      )[["mean"]]

      # convert predictions to probabilities if binary DV
      if (dv_type == "binary") {
        preds_bs <- stats::plogis(preds_bs)
      }

    } else {

      # Train mth model on kth training set
      model_bs <- best_subset_classifier(
        model = attr(obj, "tuned_hyperparameters")$best_subset$form,
        y = attr(obj, "evaluated_args")$y,
        data.train = data.train,
        model.family = binomial(link = "probit"),
        model.optimizer = "bobyqa",
        n.iter = 1000000,
        verbose = verbose
      )

      # best subset predictions
      preds_bs <- stats::predict(
        model_bs, newdata = data.test, type = "response",
        allow.new.levels = TRUE
      )

    }
  } else {
    preds_bs <- NULL
  }

  # lasso predictions
  if (attr(obj, "evaluated_args")$lasso) {

    # lasso model on data.train
    model_l <- lasso_classifier(
      y = attr(obj, "evaluated_args")$y,
      L2.fix = attr(obj, "tuned_hyperparameters")$lasso$L2.fe.form,
      L1.re = attr(obj, "tuned_hyperparameters")$lasso$L1.re,
      data.train = data.train,
      lambda = attr(obj, "tuned_hyperparameters")$lasso$lasso.lambda,
      model.family = binomial(link = "probit"),
      verbose = FALSE
    )

    # Use trained model to make predictions for data.test
    preds_l <- stats::predict(model_l, newdata = data.frame(data.test))

  } else {
    preds_l <- NULL
  }

  # pca predictions
  if (attr(obj, "evaluated_args")$pca) {

    # deep MrP
    if (attr(obj, "evaluated_args")$pca.deep) {

      # Train mth model on kth training set
      model_pca <- deep_mrp_classifier(
        form = attr(obj, "tuned_hyperparameters")$pca$form,
        y = attr(obj, "evaluated_args")$y,
        data = data.train,
        verbose = FALSE
      )

      # pca predictions
      preds_pca <- vglmer::predict_MAVB(
        samples = 1000,
        model_pca,
        newdata = data.test,
        allow_missing_levels = TRUE
      )[["mean"]]

      # convert predictions to probabilities if binary DV
      if (dv_type == "binary") {
        preds_pca <- stats::plogis(preds_pca)
      }

    } else {

      # Train mth model on kth training set
      model_pca <- best_subset_classifier(
        model = attr(obj, "tuned_hyperparameters")$pca$form,
        y = attr(obj, "evaluated_args")$y,
        data.train = data.train,
        model.family = binomial(link = "probit"),
        model.optimizer = "bobyqa",
        n.iter = 1000000,
        verbose = verbose
      )

      # pca predictions
      preds_pca <- stats::predict(
        model_pca, newdata = data.test, type = "response",
        allow.new.levels = TRUE
      )
    }
  } else {
    preds_pca <- NULL
  }

  # GB predictions
  if (attr(obj, "evaluated_args")$gb) {

    # Train mth model on kth training set
    model_gb <- gb_classifier(
      y = attr(obj, "evaluated_args")$y,
      form = attr(obj, "tuned_hyperparameters")$gb$form,
      distribution = "bernoulli",
      data.train = data.train,
      n.trees = attr(obj, "tuned_hyperparameters")$gb$n.trees,
      interaction.depth =
        attr(obj, "tuned_hyperparameters")$gb$interaction.depth,
      n.minobsinnode = if (
        !is.null(attr(obj, "evaluated_args")$gb.n.minobsinnode)
      ) {
        attr(obj, "evaluated_args")$gb.n.minobsinnode
      } else {
        20 # default value
      },
      shrinkage = attr(obj, "tuned_hyperparameters")$gb$shrinkage,
      verbose = FALSE
    )

    # Use trained model to make predictions for data.test
    preds_gb <- if (dv_type == "binary") {
      stats::predict(
        model_gb,
        newdata = data.frame(data.test),
        n.trees = model_gb$n.trees,
        type = "response"
      )
    } else {
      stats::predict(
        model_gb,
        newdata = data.frame(data.test),
        n.trees = model_gb$n.trees
      )
    }
  } else {
    preds_gb <- NULL
  }

  # svm predictions
  if (attr(obj, "evaluated_args")$svm) {

    # Svm classifier
    model_svm <- svm_classifier(
      y = all.vars(attr(obj, "tuned_hyperparameters")$svm$form)[1],
      form = attr(obj, "tuned_hyperparameters")$svm$form,
      data = data.train %>%
        dplyr::mutate(
          !!all.vars(
            attr(obj, "tuned_hyperparameters")$svm$form
          )[1] := as.factor(
            !!rlang::sym(attr(obj, "evaluated_args")$y)
          )
        ),
      kernel = attr(obj, "tuned_hyperparameters")$svm$kernel,
      type = "C-classification",
      probability = TRUE,
      svm.gamma = attr(obj, "tuned_hyperparameters")$svm$gamma,
      svm.cost = attr(obj, "tuned_hyperparameters")$svm$cost,
      verbose = FALSE
    )

    preds_svm <- predict(
      model_svm,
      newdata = if (dv_type == "binary") {
        data.test %>%
          dplyr::mutate(
            !!all.vars(attr(obj, "tuned_hyperparameters")$svm$form)[1] :=
              as.factor(
                !!rlang::sym(attr(obj, "evaluated_args")$y)
              )
          )
      } else {
        data.frame(data.test)
      },
      probability = TRUE
    )
    if (!is.null(attr(preds_svm, "probabilities")[, "1"])) {
      preds_svm <- as.numeric(attr(preds_svm, "probabilities")[, "1"])
    }
  } else {
    preds_svm <- NULL
  }


  # knn predictions
  if (attr(obj, "evaluated_args")$knn) {

    # Train mth model on kth training set
    model_knn <- knn_classifier(
      y = attr(obj, "evaluated_args")$y,
      form = attr(obj, "tuned_hyperparameters")$knn$form,
      data.train = if (dv_type == "binary") {
        data.train %>%
          dplyr::mutate(!!attr(obj, "evaluated_args")$y := as.factor(
            !!rlang::sym(attr(obj, "evaluated_args")$y)
          ))
      } else {
        data.train
      },
      data.valid = if (dv_type == "binary") {
        data.test %>%
          dplyr::mutate(!!attr(obj, "evaluated_args")$y := as.factor(
            !!rlang::sym(attr(obj, "evaluated_args")$y)
          ))
      } else {
        data.test
      },
      knn.k.value = attr(obj, "tuned_hyperparameters")$knn$knn.k,
      knn.kernel = attr(obj, "tuned_hyperparameters")$knn$knn.kernel,
      verbose = FALSE
    )

    # Get predictions for data.test
    preds_knn <- if (dv_type == "binary") {
      kknn:::predict.kknn(
        model_knn,
        data.test %>%
          dplyr::mutate(!!attr(obj, "evaluated_args")$y := as.factor(
            !!rlang::sym(attr(obj, "evaluated_args")$y)
          )),
        type = "prob"
      )[, "1"]
    } else {
      model_knn$fit
    }
  } else {
    preds_knn <- NULL
  }

  # MrP predictions
  if (attr(obj, "evaluated_args")$mrp) {

    # Train mth model on kth training set
    model_mrp <- best_subset_classifier(
      model = attr(obj, "tuned_hyperparameters")$mrp,
      y = attr(obj, "evaluated_args")$y,
      data.train = data.train,
      model.family = binomial(link = "probit"),
      model.optimizer = "bobyqa",
      n.iter = 1000000,
      verbose = FALSE
    )

    # MrP predictions
    preds_mrp <- stats::predict(
      model_mrp, newdata = data.test,
      type = "response", allow.new.levels = TRUE
    )
  } else {
    preds_mrp <- NULL
  }

  # deep MrP predictions
  if (attr(obj, "evaluated_args")$deep.mrp) {

    # Train mth model on kth training set
    model_deep <- deep_mrp_classifier(
      form = attr(obj, "tuned_hyperparameters")$deep,
      y = attr(obj, "evaluated_args")$y,
      data = data.train,
      verbose = FALSE
    )

    # deep MrP predictions
    preds_deep <- vglmer::predict_MAVB(
      samples = 1000,
      model_deep,
      newdata = data.test,
      allow_missing_levels = TRUE
    )[["mean"]]

    # convert predictions to probabilities if binary DV
    if (dv_type == "binary") {
      preds_deep <- stats::plogis(preds_deep)
    }
  } else {
    preds_deep <- NULL
  }

  # combine classifier predictions to a tibble
  preds <- tibble::tibble(
    best_subset = preds_bs,
    lasso = if (!is.null(preds_l)) {
      as.numeric(preds_l)
    } else {
      NULL
    },
    pca = preds_pca,
    gb = preds_gb,
    svm = preds_svm,
    knn = preds_knn,
    mrp = preds_mrp,
    deep_mrp = preds_deep
  )

  # add EBMA predictions
  if (!is.null(obj$weights)) {
    preds$ebma <- as.numeric(as.matrix(preds) %*% as.matrix(obj$weights))
  }

  # add stacking predictions
  if (!is.null(obj$stacking)) {

    # non-negative least squares
    preds$stack_nnls <- tryCatch(
      {
        as.numeric(
          as.matrix(preds[, names(obj$classifiers)[-1]]) %*%
            as.matrix(obj$stacking_weights$stack_nnls)
        )
      },
      error = function(e) {
        NULL
      }
    )

    # optim with non-negative constraints
    preds$stack_optim <- tryCatch(
      {
        as.numeric(
          as.matrix(preds[, names(obj$classifiers)[-1]]) %*%
            as.matrix(obj$stacking_weights$stack_optim)
        )
      },
      error = function(e) {
        NULL
      }
    )

    # stack quadratic programming
    preds$stack_qp <- tryCatch(
      {
        as.numeric(
          as.matrix(preds[, names(obj$classifiers)[-1]]) %*%
            as.matrix(obj$stacking_weights$stack_qp)
        )
      },
      error = function(e) {
        NULL
      }
    )

    # stack ornstein
    preds$stack_ornstein <- tryCatch(
      {
        as.numeric(
          as.matrix(preds[, names(obj$classifiers)[-1]]) %*%
            as.matrix(obj$stacking_weights$stack_ornstein)
        )
      },
      error = function(e) {
        NULL
      }
    )

    # stack of stacks
    stacking_procs <- names(obj$stacking)[
      !names(obj$stacking) %in% c("stack_of_stacks", "stack_of_stacks_ebma")
    ][-1]

    preds$stack_of_stacks <- tryCatch(
      {
        as.numeric(
          as.matrix(preds[, stacking_procs]) %*%
            as.matrix(obj$stacking_weights$stack_of_stacks)
        )
      },
      error = function(e) {
        NULL
      }
    )

    # stack of stacks EBMA
    preds$stack_of_stacks_ebma <- tryCatch(
      {
        as.numeric(
          as.matrix(preds[, c("stack_of_stacks", "ebma")]) %*%
            as.matrix(obj$stacking_weights$stack_of_stacks_ebma)
        )
      },
      error = function(e) {
        NULL
      }
    )
  }

  return(preds)

}