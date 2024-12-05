# generates out of sample predictions from the tuned classifiers
get_predictions <- function(
  y, L1.x, L2.x, L2.unit, L2.reg, best.subset.opt, best.subset.deep,
  lasso.opt, lasso.L2.x, pca.opt, pca.deep, gb.opt, svm.opt, svm.L2.reg,
  svm.L2.unit, svm.L2.x, knn.opt, knn.L2.reg, knn.L2.unit,
  knn.L2.x, mrp.include, n.minobsinnode, L2.unit.include, L2.reg.include,
  kernel, knn.kernel, mrp.L2.x, deep.mrp, deep.splines, deep.L2.unit, deep.L2.x,
  deep.L2.reg, data, ebma.fold, verbose, cv.sampling, k.folds = k.folds,
  all_data = TRUE
) {

  message("this only works for binary dependent variables for now")

  #------------------------------------------------------------
  # dependent variable
  #------------------------------------------------------------
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

  #------------------------------------------------------------
  # data
  #------------------------------------------------------------
  # split data into k without holdout fold
  if (all_data) {
    cv_data <- dplyr::bind_rows(data) %>%
      dplyr::bind_rows(ebma.fold)
  } else {
    cv_data <- dplyr::bind_rows(data)
  }

  # add factor version of y for svm
  if (!is.null(svm.opt)) {
    # Prepare data
    cv_data <- cv_data %>%
      dplyr::mutate(y_svm := as.factor(!!rlang::sym(y))) %>%
      dplyr::relocate(y_svm, .after = !!rlang::sym(y))
  }

  #------------------------------------------------------------
  # model formulas
  #------------------------------------------------------------
  # set lasso formula
  if (!is.null(lasso.L2.x)) {

    # Determine context-level covariates
    if (is.null(lasso.L2.x)) {
      lasso.L2.x <- L2.x
    }

    # Context-level fixed effects
    L2_fe <- paste(lasso.L2.x, collapse = " + ")
    L2_fe_form <- as.formula(paste(y, " ~ ", L2_fe, sep = ""))

    # Individual-level random effects as named list
    L1_re <- setNames(
      as.list(
        rep(
          c(~ 1),
          times = length(
            c(L1.x, L2.unit, L2.reg)
          )
        )
      ), c(L1.x, L2.unit, L2.reg)
    )
  }

  # set the gb formula
  if (!is.null(gb.opt)) {

    # Evaluate inclusion of L2.unit
    if (isTRUE(L2.unit.include == TRUE)) {
      L2.unit.gb <- L2.unit
    } else {
      L2.unit.gb <- NULL
    }

    # Evaluate inclusion of L2.reg
    if (isTRUE(L2.reg.include == TRUE)) {
      L2.reg.gb <- L2.reg
    } else {
      L2.reg.gb <- NULL
    }

    # Create model formula
    x <- paste(c(L1.x, L2.x, L2.unit.gb, L2.reg.gb), collapse = " + ")
    form_gb <- as.formula(paste(y, " ~ ", x, sep = ""))
  }

  # svm
  if (!is.null(svm.opt)) {

    # Determine context-level covariates
    if (is.null(svm.L2.x)) {
      svm.L2.x <- L2.x
    }

    # Evaluate inclusion of L2.unit in GB
    if (isTRUE(svm.L2.unit)) {
      svm.L2.unit <- L2.unit
    } else {
      svm.L2.unit <- NULL
    }

    # Evaluate inclusion of L2.reg in GB
    if (isTRUE(svm.L2.reg)) {
      svm.L2.reg <- L2.reg
    } else {
      svm.L2.reg <- NULL
    }

    # Create model formula
    x <- paste(c(L1.x, svm.L2.x, svm.L2.unit, svm.L2.reg), collapse = " + ")
    form_svm <- as.formula(paste("y_svm ~ ", x, sep = ""))
  }

  # KNN
  if (!is.null(knn.opt)) {

    # Determine context-level covariates
    if (is.null(knn.L2.x)) {
      knn.L2.x <- L2.x
    }

    # Evaluate inclusion of L2.unit in KNN
    if (isTRUE(knn.L2.unit)) {
      knn.L2.unit <- L2.unit
    } else {
      knn.L2.unit <- NULL
    }

    # Evaluate inclusion of L2.reg in KNN
    if (isTRUE(knn.L2.reg)) {
      knn.L2.reg <- L2.reg
    } else {
      knn.L2.reg <- NULL
    }

    # Create model formula
    x <- paste(c(L1.x, knn.L2.x, knn.L2.unit, knn.L2.reg), collapse = " + ")
    form_knn <- as.formula(paste(y, " ~ ", x, sep = ""))
  }

  # mrp
  if (isTRUE(mrp.include == TRUE)) {
    # Create model formula
    # Individual-level random effects
    L1_re <- paste(paste("(1 | ", L1.x, ")", sep = ""), collapse = " + ")

    # Geographic unit or geographic unit-geographic region random effects
    if (is.null(L2.reg)) {
      L2_re <- paste("(1 | ", L2.unit, ")", sep = "")
    } else {
      L2_re <- paste(
        paste("(1 | ", L2.reg, "/", L2.unit, ")", sep = ""),
        collapse = " + "
      )
    }

    # Combine all random effects
    all_re <- paste(c(L1_re, L2_re), collapse = " + ")

    # Context-level fixed effects
    if (is.null(mrp.L2.x)) {
      L2_fe <- paste(L2.x, collapse = " + ")
    } else {
      if (length(mrp.L2.x) == 1) {
        if (mrp.L2.x == "empty") {
          L2_fe <- ""
        } else {
          L2_fe <- paste(mrp.L2.x, collapse = " + ")
        }
      } else {
        L2_fe <- paste(mrp.L2.x, collapse = " + ")
      }
    }

    # std MrP formula
    form_mrp <- as.formula(paste(y, " ~ ", L2_fe, " + ", all_re, sep = ""))

    # Individual-level random effects as named list
    L1_re <- setNames(
      as.list(
        rep(
          c(~ 1),
          times = length(
            c(L1.x, L2.unit, L2.reg)
          )
        )
      ), c(L1.x, L2.unit, L2.reg)
    )
  }

  # deep mrp
  if (deep.mrp) {

    # custom L2x
    if (!is.null(deep.L2.x)) {
      L2.x <- deep.L2.x
    }

    # custom L2.unit
    if (isTRUE(deep.L2.unit)) {
      deep.L2.unit <- L2.unit
    } else {
      deep.L2.unit <- NULL
    }

    # custom L2.reg
    if (isTRUE(deep.L2.reg)) {
      deep.L2.reg <- L2.reg
    } else {
      deep.L2.reg <- NULL
    }

    # generate all interactions of L1.x
    l1_comb <- unlist(lapply(2:length(L1.x), function(x) {
      apply(combn(L1.x, x), 2, paste, collapse = ".")
    }))

    # generate all interactions of L1.x with L2.unit
    if (!is.null(deep.L2.unit)) {
      l1_state <- paste(L1.x, L2.unit, sep = ".")
    } else {
      l1_state <- NULL
    }

    # generate all interactions of L1.x with L2.reg
    if (!is.null(deep.L2.reg)) {
      l1_region <- paste(L1.x, L2.reg, sep = ".")
    } else {
      l1_region <- NULL
    }

    # model formula
    if (!deep.splines) {
      form <- paste0(
        # DV
        y, " ~ ",
        # L2x
        paste(L2.x, collapse = " + "), " + ",
        # individual level variables
        paste("(1 | ", L1.x, ")", collapse = " + "), " + ",
        # interactions of L1x
        paste("(1 | ", l1_comb, ")", collapse = " + "), " + ",
        # interactions of L1x with L2.unit
        if (any(!is.null(l1_state))) {
          paste("(1 | ", l1_state, ")", collapse = " + ")
        }, " + ",
        # interactions of L1x with L2.reg
        if (any(!is.null(l1_region))) {
          paste("(1 | ", l1_region, ")", collapse = " + ")
        }
      )
    } else {
      form <- paste0(
        # DV
        y, " ~ ",
        # L2x
        paste0("v_s(", L2.x, ")", collapse = " + "), " + ",
        # individual level variables
        paste("(1 | ", L1.x, ")", collapse = " + "), " + ",
        # interactions of L1x
        paste("(1 | ", l1_comb, ")", collapse = " + "), " + ",
        # interactions of L1x with L2.unit
        if (any(!is.null(l1_state))) {
          paste("(1 | ", l1_state, ")", collapse = " + ")
        }, " + ",
        # interactions of L1x with L2.reg
        if (any(!is.null(l1_region))) {
          paste("(1 | ", l1_region, ")", collapse = " + ")
        }
      )
    }

    # match everything from the beginning to the last ")"
    # for example if the string ends in +
    form <- stringr::str_extract(string = form, pattern = "^.*\\)")
  }

  #------------------------------------------------------------
  # cross-validation folds
  #------------------------------------------------------------

  # get folds
  cv_folds <- cv_folding(
    data = cv_data,
    L2.unit = L2.unit,
    k.folds = k.folds,
    cv.sampling = cv.sampling
  )

  #------------------------------------------------------------
  # cross-validation loop
  #------------------------------------------------------------
  # generate out of sample predictions from the best classifiers
  k_preds <- lapply(seq_len(k.folds), function(k) {

    # Split data in training and validation sets
    data_train <- dplyr::bind_rows(cv_folds[-k])
    data_valid <- dplyr::bind_rows(cv_folds[k])

    # 1) best subset
    if (!is.null(best.subset.opt)) {

      # fit best subset model
      if (best.subset.deep) {
        best_subset <-  deep_mrp_classifier(
          form = best.subset.opt,
          y = y,
          data = data_train,
          verbose = TRUE
        )
        # predictions based on DV type (binary or continuous)
        if (dv_type == "binary") {
          # use trained model to make predictions for kth validation set
          bs_preds <- vglmer::predict_MAVB(
            samples = 1000,
            best_subset,
            newdata = data_valid,
            allow_missing_levels = TRUE
          )[["mean"]]

          # convert to response probabilities
          bs_preds <- stats::plogis(bs_preds)

        } else if (dv_type == "linear") {
          # Use trained model to make predictions for kth validation set
          bs_preds <- predict(
            samples = 1000,
            object = best_subset,
            newdata = data_valid,
            allow_missing_levels = TRUE
          )[["mean"]]
        }
      } else {
        best_subset <- best_subset_classifier(
          y = y,
          model = best.subset.opt,
          data.train = data_train,
          model.family = binomial(link = "probit"),
          model.optimizer = "bobyqa",
          n.iter = 1000000,
          verbose = verbose
        )
        if (dv_type == "binary") {
          # predict on validation set
          bs_preds <- stats::predict(
            object = best_subset,
            newdata = data_valid,
            allow.new.levels = TRUE,
            type = "response"
          )
        } else {
          # predict on validation set
          bs_preds <- stats::predict(
            object = best_subset,
            newdata = data_valid,
            allow.new.levels = TRUE
          )
        }
      }

    }

    # 2) lasso
    if (!is.null(lasso.opt)) {

      # fit lasso model
      lasso <- lasso_classifier(
        y = y,
        L2.fix = L2_fe_form,
        L1.re = L1_re,
        data.train = data_train,
        lambda = lasso.opt,
        model.family = binomial(link = "probit"),
        verbose = verbose
      )

      # predict on validation set
      lasso_preds <- predict_glmmLasso(
        m = lasso,
        lasso.L2.x = lasso.L2.x,
        L2.unit = L2.unit,
        L2.reg = L2.reg,
        L1.x = L1.x,
        census = data_valid,
        type = "response" # detects if binary or continuous
      )
    }

    # 3) pca
    if (!is.null(pca.opt)) {

      # fit pca model
      if (pca.deep) {
        # call deep_mrp_classifier
        pca <-  deep_mrp_classifier(
          form = pca.opt,
          y = y,
          data = data_train,
          verbose = TRUE
        )
        # predictions based on DV type (binary or continuous)
        if (dv_type == "binary") {
          # use trained model to make predictions for kth validation set
          pca_preds <- vglmer::predict_MAVB(
            samples = 1000,
            pca,
            newdata = data_valid,
            allow_missing_levels = TRUE
          )[["mean"]]

          # convert to response probabilities
          pca_preds <- stats::plogis(pca_preds)

        } else if (dv_type == "linear") {
          # Use trained model to make predictions for kth validation set
          pca_preds <- predict(
            samples = 1000,
            object = pca,
            newdata = data_valid,
            allow_missing_levels = TRUE
          )[["mean"]]
        }
      } else {
        # call best_subset_classifier
        pca <- best_subset_classifier(
          y = y,
          model = pca.opt,
          data.train = data_train,
          model.family = binomial(link = "probit"),
          model.optimizer = "bobyqa",
          n.iter = 1000000,
          verbose = verbose
        )

        # predict on validation set
        pca_preds <- stats::predict(
          object = pca,
          newdata = data_valid,
          allow.new.levels = TRUE,
          type = "response"
        )
      }
    }

    # 4) gbm
    if (!is.null(gb.opt)) {

      # fit gbm model
      gb <- gb_classifier(
        y = y,
        form = form_gb,
        distribution = "bernoulli",
        data.train = data_train,
        n.trees = gb.opt$n_trees,
        interaction.depth = gb.opt$interaction_depth,
        n.minobsinnode = n.minobsinnode,
        shrinkage = gb.opt$shrinkage,
        verbose = verbose
      )

      # predict on validation set
      gb_preds <- gbm::predict.gbm(
        object = gb,
        newdata = data_valid,
        n.trees = gb.opt$n_trees,
        type = "response"
      )

    }

    # 5) svm
    if (!is.null(svm.opt)) {

      # fit svm model
      svm <- svm_classifier(
        y = "y_svm",
        form = form_svm,
        data = data_train,
        kernel = svm.opt$kernel,
        type = "C-classification",
        probability = TRUE,
        svm.gamma = svm.opt$gamma,
        svm.cost = svm.opt$cost,
        verbose = verbose
      )

      # predict on validation set
      if (
        data_train %>%
          dplyr::pull(var = "y_svm") %>%
          unique() %>%
          length() > 2
      ) {
        svm_preds <- stats::predict(
          object = svm,
          newdata = data_valid,
          type = "response"
        )
      } else {
        svm_preds <- attr(stats::predict(
          object = svm,
          newdata = data_valid,
          probability = TRUE
        ), "probabilities")[, "1"]
      }

    }

    # 6) knn
    if (!is.null(knn.opt)) {

      # fit knn model
      knn <- knn_classifier(
        y = y,
        form = form_knn,
        data.train = if (dv_type == "binary") {
          data_train %>%
            dplyr::mutate(!!rlang::sym(y) := as.factor(!!rlang::sym(y)))
        } else {
          data_train
        },
        data.valid = data_valid,
        knn.k.value = knn.opt,
        knn.kernel = knn.kernel,
        verbose = verbose,
        contr.dummy = kknn:::contr.dummy
      )

      # predictions on validation set
      knn_preds <- if (dv_type == "binary") {
        kknn:::predict.kknn(knn, type = "prob")[, "1"]
      } else {
        knn$fit
      }

    }

    # 7) mrp
    if (isTRUE(mrp.include)) {

      # fit mrp model
      mrp <- best_subset_classifier(
        y = y,
        model = form_mrp,
        data.train = data_train,
        model.family = binomial(link = "probit"),
        model.optimizer = "bobyqa",
        n.iter = 1000000,
        verbose = verbose
      )

      # predict on validation set
      mrp_preds <- stats::predict(
        object = mrp,
        newdata = data_valid,
        allow.new.levels = TRUE,
        type = "response"
      )

    }

    # 8) deep mrp
    if (deep.mrp) {

      # fit deep mrp model
      deep_mrp <- deep_mrp_classifier(
        y = y,
        form = form,
        data = data_train,
        verbose = verbose
      )

      # Determine type of dependent variable
      if (
        data_train %>%
          dplyr::pull(!!y) %>%
          unique() %>%
          length() == 2
      ) {
        dv_type <- "binary"
      } else {
        dv_type <- "continuous"
      }

      # binary or continuous DV
      if (dv_type == "binary") {

        # predict on validation set for binary DV
        deep_preds <- vglmer::predict_MAVB(
          samples = 1000,
          deep_mrp,
          newdata = data_valid,
          allow_missing_levels = TRUE
        )[["mean"]]

        # convert to response probabilities
        deep_preds <- stats::plogis(deep_preds)

      } else if (dv_type == "continuous") {

        # predict on validation set for continuous DV
        deep_preds <- predict(
          samples = 1000,
          object = deep_mrp,
          newdata = data_valid,
          allow_missing_levels = TRUE
        )[["mean"]]

      }
    }

    # combine predictions into one table
    preds <- tibble::tibble(
      y = data_valid[[y]],
      !!rlang::sym(L2.unit) := data_valid[[L2.unit]]
    )
    if (exists("bs_preds")) {
      preds <- dplyr::mutate(
        .data = preds, best_subset = bs_preds
      )
    }
    if (exists("lasso_preds")) {
      preds <- dplyr::mutate(
        .data = preds, lasso = as.numeric(lasso_preds)
      )
    }
    if (exists("pca_preds")) {
      preds <- dplyr::mutate(.data = preds, pca = pca_preds)
    }
    if (exists("gb_preds")) {
      preds <- dplyr::mutate(.data = preds, gb = gb_preds)
    }
    if (exists("svm_preds")) {
      preds <- dplyr::mutate(.data = preds, svm = svm_preds)
    }
    if (exists("knn_preds")) {
      preds <- dplyr::mutate(.data = preds, knn = knn_preds)
    }
    if (exists("mrp_preds")) {
      preds <- dplyr::mutate(.data = preds, mrp = mrp_preds)
    }
    if (exists("deep_preds")) {
      preds <- dplyr::mutate(.data = preds, deep_mrp = deep_preds)
    }

    # return predictions table
    return(preds)

  }) # end of loop over k folds

  # combine predictions into one table
  preds <- dplyr::bind_rows(k_preds)

  # return model formulas
  best_classifiers <- list(
    best_subset = if (!is.null(pca.opt)) {
      list(
        form = best.subset.opt
      )
    } else {
      NULL
    },
    lasso = if (!is.null(lasso.opt)) {
      list(
        L2.fe.form = L2_fe_form,
        L1.re = L1_re,
        lasso.lambda = lasso.opt
      )
    } else {
      NULL
    },
    pca = if (!is.null(pca.opt)) {
      list(
        form = pca.opt
      )
    } else {
      NULL
    },
    gb = if (!is.null(gb.opt)) {
      list(
        form = form_gb,
        interaction.depth = gb.opt$interaction_depth,
        shrinkage = gb.opt$shrinkage,
        n.trees = gb.opt$n_trees
      )
    } else {
      NULL
    },
    svm = if (!is.null(svm.opt)) {
      list(
        form = form_svm,
        gamma = svm.opt$gamma,
        cost = svm.opt$cost,
        kernel = svm.opt$kernel
      )
    } else {
      NULL
    },
    knn = if (!is.null(knn.opt)) {
      list(
        form = form_knn,
        knn.k = knn.opt,
        knn.kernel = knn.kernel
      )
    } else {
      NULL
    },
    mrp = if (mrp.include) {
      form_mrp
    } else {
      NULL
    },
    deep_mrp = if (deep.mrp) {
      form
    } else {
      NULL
    }
  )

  return(list(preds = preds, best_classifiers = best_classifiers))
}
