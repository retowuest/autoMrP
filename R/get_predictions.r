# generates out of sample predictions from the tuned classifiers
get_predictions <- function(
  y, L1.x, L2.x, L2.unit, L2.reg, best.subset.opt, lasso.opt, lasso.L2.x,
  pca.opt, gb.opt, svm.opt, svm.L2.reg, svm.L2.unit, svm.L2.x, mrp.include,
  n.minobsinnode, L2.unit.include, L2.reg.include, kernel, mrp.L2.x, deep.mrp,
  deep.L2.x, deep.L2.reg, deep.splines, data, ebma.fold, verbose,
  cv.sampling, k.folds = k.folds, all_data = TRUE
) {

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
      dplyr::mutate(
        dplyr::across(dplyr::all_of(y), ~as.factor(.), .names = "y_svm")
      ) %>%
      dplyr::relocate(y_svm, .after = {{y}})
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
    l1_state <- paste(L1.x, L2.unit, sep = ".")

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
        paste("(1 | ", l1_state, ")", collapse = " + "), " + ",
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
        paste("(1 | ", l1_state, ")", collapse = " + "), " + ",
        # interactions of L1x with L2.reg
        if (any(!is.null(l1_region))) {
          paste("(1 | ", l1_region, ")", collapse = " + ")
        }
      )
    }

    # match everything from the beginning to the last ")"
    # for example if the string ends in +
    form <- stringr::str_extract(string = form, pattern = "^.*\\)")

    # add the interactions to the data
    all_interactions <- c(l1_comb, l1_state, l1_region)

    # loop over all interactions for data object
    x_data <- lapply(all_interactions, function(x) {

      # break down interaction components
      y <- stringr::str_extract(
        string = x,
        pattern = stringr::fixed(pattern = names(cv_data))
      ) %>%
        .[!is.na(.)]

      # take each column of data and combine its values into a single string
      df_x <- cv_data %>%
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

    # combine data and interactions
    cv_data <- dplyr::bind_cols(cv_data, x_data)
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
      best_subset <- best_subset_classifier(
        y = y,
        model = best.subset.opt,
        data.train = data_train,
        model.family = binomial(link = "probit"),
        model.optimizer = "bobyqa",
        n.iter = 1000000,
        verbose = verbose
      )

      # predict on validation set
      bs_preds <- stats::predict(
        object = best_subset,
        newdata = data_valid,
        allow.new.levels = TRUE,
        type = "response"
      )
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

      # fit gbm model
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

    # 6) mrp
    if (isTRUE(mrp.include == TRUE)) {

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

    # 7) deep mrp
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

  return(preds)
}