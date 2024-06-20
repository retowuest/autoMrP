#' Apply post-stratification to classifiers.
#'
#' @inheritParams auto_MrP
#' @param best.subset.opt Optimal tuning parameters from best subset selection
#'   classifier. A list returned by \code{run_best_subset()}.
#' @param lasso.opt Optimal tuning parameters from lasso classifier A list
#'   returned by \code{run_lasso()}.
#' @param pca.opt Optimal tuning parameters from best subset selection with
#'   principal components classifier A list returned by \code{run_pca()}.
#' @param gb.opt Optimal tuning parameters from gradient tree boosting
#'   classifier A list returned by \code{run_gb()}.
#' @param svm.opt Optimal tuning parameters from support vector machine
#'   classifier A list returned by \code{run_svm()}.
#' @param mrp.include Whether to run MRP classifier. A logical argument
#'   indicating whether the standard MRP classifier should be used for
#'   predicting outcome \code{y}. Passed from \code{autoMrP()} argument
#'   \code{mrp}.
#' @param n.minobsinnode GB minimum number of observations in the terminal
#'   nodes. An integer-valued scalar specifying the minimum number of
#'   observations that each terminal node of the trees must contain. Passed from
#'   \code{autoMrP()} argument \code{gb.n.minobsinnode}.
#' @param L2.unit.include GB L2.unit. A logical argument indicating whether
#'   \code{L2.unit} should be included in the GB classifier. Passed from
#'   \code{autoMrP()} argument \code{gb.L2.unit}.
#' @param L2.reg.include A logical argument indicating whether \code{L2.reg}
#'   should be included in the GB classifier. Passed from \code{autoMrP()}
#'   argument \code{GB L2.reg}.
#' @param kernel SVM kernel. A character-valued scalar specifying the kernel to
#'   be used by SVM. The possible values are \code{linear}, \code{polynomial},
#'   \code{radial}, and \code{sigmoid}. Passed from \code{autoMrP()} argument
#'   \code{svm.kernel}.
#' @param data A data.frame containing the survey data used in classifier
#'   training.
#' @param ebma.fold A data.frame containing the data not used in classifier
#'   training.

post_stratification <- function(
  y, L1.x, L2.x, L2.unit, L2.reg, best.subset.opt, lasso.opt, lasso.L2.x,
  pca.opt, gb.opt, svm.opt, svm.L2.reg, svm.L2.unit, svm.L2.x, mrp.include,
  n.minobsinnode, L2.unit.include, L2.reg.include, kernel, mrp.L2.x,
  data, ebma.fold, census, verbose, deep.mrp, stacking, deep.L2.x, deep.L2.reg,
  deep.splines
) {

  # globals
  lasso <- NULL
  pca <- NULL
  gb <- NULL
  svm <- NULL
  mrp <- NULL
  best_subset <- NULL

  # Copy L2.unit b/c it is needed twice but might be reset depending on call
  L2_unit <- L2.unit

  # post-stratification without level 2 variables
  if (all(L2.x == "")) {
    L2.x <- NULL
  }

  # model container for EBMA
  models <- list()

  # data that is used for models that will not enter EBMA
  no_ebma_data <- dplyr::bind_rows(data, ebma.fold)

  # remove missing values
  data <- tidyr::drop_na(
    data = data, dplyr::all_of(c(y, L1.x, L2.x, L2.unit, L2.reg))
  )
  no_ebma_data <- tidyr::drop_na(
    data = no_ebma_data, dplyr::all_of(c(y, L1.x, L2.x, L2.unit, L2.reg))
  )

  # ----- Fit optimal model and make prediction for individual classifiers -----

  # Classifier 1: Best Subset
  if (!is.null(best.subset.opt)) {

    # Fit optimal model for EBMA
    best_subset_opt_ebma <- best_subset_classifier(
      y = y,
      model = best.subset.opt,
      data.train = data,
      model.family = binomial(link = "probit"),
      model.optimizer = "bobyqa",
      n.iter = 1000000,
      verbose = verbose
    )

    # Fit optimal model for post-stratification w/o EBMA
    best_subset_opt_poststrat_only <- best_subset_classifier(
      y = y,
      model = best.subset.opt,
      data.train = no_ebma_data,
      model.family = binomial(link = "probit"),
      model.optimizer = "bobyqa",
      n.iter = 1000000,
      verbose = verbose
    )

    # post-stratification
    bs_preds <- census %>%
      dplyr::mutate(
        best_subset = stats::predict(
          object = best_subset_opt_poststrat_only,
          newdata = .,
          allow.new.levels = TRUE,
          type = "response"
        )
      ) %>%
      dplyr::group_by(!! rlang::sym(L2.unit)) %>%
      dplyr::summarize(best_subset = stats::weighted.mean(
        x = best_subset, w = prop
      ), .groups = "keep") %>%
      dplyr::ungroup()

    # individual level predictions for EBMA
    bs_ind <- stats::predict(object = best_subset_opt_ebma, type = "response")

    # individual level predictions for stacking
    if (stacking) {
      # ddd best subset predictions to the census
      census <- census %>%
        dplyr::mutate(
          best_subset = stats::predict(
            object = best_subset_opt_poststrat_only,
            newdata = .,
            allow.new.levels = TRUE
          )
        )
      # individual level predictions for stacking on the base data
      bs_stacking <- stats::predict(
        object = best_subset_opt_poststrat_only, type = "response"
      )
    }

    # model for EBMA
    models$best_subset <- best_subset_opt_ebma

  }

  # Classifier 2: Lasso
  if (!is.null(lasso.opt)) {

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

    # Fit optimal model for EBMA
    lasso_opt_ebma <- lasso_classifier(
      y = y,
      L2.fix = L2_fe_form,
      L1.re = L1_re,
      data.train = data,
      lambda = lasso.opt,
      model.family = binomial(link = "probit"),
      verbose = verbose
    )

    # Fit optimal model for post-stratification w/o EBMA
    lasso_opt_poststrat_only <- lasso_classifier(
      y = y,
      L2.fix = L2_fe_form,
      L1.re = L1_re,
      data.train = no_ebma_data,
      lambda = lasso.opt,
      model.family = binomial(link = "probit"),
      verbose = verbose
    )

    # predictions for post-stratification only (no EBMA)
    lasso_ests <- predict_glmmLasso(
      m = lasso_opt_poststrat_only,
      lasso.L2.x = lasso.L2.x,
      L2.unit = L2.unit,
      L2.reg = L2.reg,
      L1.x = L1.x,
      census = census
    )

    # post-stratification
    lasso_preds <- census %>%
      dplyr::mutate(lasso = lasso_ests) %>%
      dplyr::group_by(!! rlang::sym(L2.unit)) %>%
      dplyr::summarize(
        lasso = stats::weighted.mean(x = lasso, w = prop), .groups = "keep"
      ) %>%
      dplyr::ungroup()

    # individual level predictions for EBMA
    lasso_ind <- stats::predict(object = lasso_opt_ebma, type = "response")

    # individual level predictions for stacking
    if (stacking) {

      # add lasso predictions to census
      census <- census %>%
        dplyr::mutate(
          lasso = predict_glmmLasso(
            m = lasso_opt_poststrat_only,
            lasso.L2.x = lasso.L2.x,
            L2.unit = L2.unit,
            L2.reg = L2.reg,
            L1.x = L1.x,
            census = census,
            type = "predictors"
          )
        )

      # individual level predictions for stacking
      lasso_stacking <- stats::predict(
        object = lasso_opt_poststrat_only, type = "response"
      )
    }

    # model for EBMA
    models$lasso <- lasso_opt_ebma

  }

  # Classifier 3: PCA
  if (!is.null(pca.opt)) {
    # Fit optimal model for EBMA
    pca_opt_ebma <- best_subset_classifier(
      y = y,
      model = pca.opt,
      data.train = data,
      model.family = binomial(link = "probit"),
      model.optimizer = "bobyqa",
      n.iter = 1000000,
      verbose = verbose
    )

    # Fit optimal model for for post-stratification w/o EBMA
    pca_opt_poststrat_only <- best_subset_classifier(
      y = y,
      model = pca.opt,
      data.train = no_ebma_data,
      model.family = binomial(link = "probit"),
      model.optimizer = "bobyqa",
      n.iter = 1000000,
      verbose = verbose
    )

    # post-stratification
    pca_preds <- census %>%
      dplyr::mutate(pca = stats::predict(
        object = pca_opt_poststrat_only,
        newdata = .,
        allow.new.levels = TRUE,
        type = "response"
      )) %>%
      dplyr::group_by(!! rlang::sym(L2.unit)) %>%
      dplyr::summarize(
        pca = stats::weighted.mean(x = pca, w = prop), .groups = "keep"
      ) %>%
      dplyr::ungroup()

    # individual level predictions for EBMA
    pca_ind <- stats::predict(object = pca_opt_ebma, type = "response")

    # individual level predictions for stacking
    if (stacking) {
      # add pca predictions to census
      census <- census %>%
        dplyr::mutate(pca = stats::predict(
          object = pca_opt_poststrat_only,
          newdata = .,
          allow.new.levels = TRUE
        ))
      # individual level predictions for stacking
      pca_stacking <- stats::predict(
        object = pca_opt_poststrat_only, type = "response"
      )
    }

    # model for EBMA
    models$pca <- pca_opt_ebma

  }

  # Classifier 4: GB
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

    # Fit optimal model for EBMA
    gb_opt_ebma <- gb_classifier(
      y = y,
      form = form_gb,
      distribution = "bernoulli",
      data.train = data,
      n.trees = gb.opt$n_trees,
      interaction.depth = gb.opt$interaction_depth,
      n.minobsinnode = n.minobsinnode,
      shrinkage = gb.opt$shrinkage,
      verbose = verbose
    )

    # Fit optimal model for for post-stratification w/o EBMA
    gb_opt_poststrat_only <- gb_classifier(
      y = y,
      form = form_gb,
      distribution = "bernoulli",
      data.train = no_ebma_data,
      n.trees = gb.opt$n_trees,
      interaction.depth = gb.opt$interaction_depth,
      n.minobsinnode = n.minobsinnode,
      shrinkage = gb.opt$shrinkage,
      verbose = verbose
    )

    # post-stratification
    gb_preds <- census %>%
      dplyr::mutate(gb = gbm::predict.gbm(
        object = gb_opt_poststrat_only,
        newdata = .,
        n.trees = gb.opt$n_trees,
        type = "response"
      )) %>%
      dplyr::group_by(!! rlang::sym(L2.unit)) %>%
      dplyr::summarize(
        gb = stats::weighted.mean(x = gb, w = prop), .groups = "keep"
      ) %>%
      dplyr::ungroup()

    # individual level predictions for EBMA
    gb_ind <- gbm::predict.gbm(
      object = gb_opt_ebma, n.trees = gb.opt$n_trees, type = "response"
    )

    # individual level predictions for stacking
    if (stacking) {

      # add gb predictions to census
      census <- census %>%
        dplyr::mutate(gb = gbm::predict.gbm(
          object = gb_opt_poststrat_only,
          newdata = .,
          n.trees = gb.opt$n_trees,
        ))

      # individual level predictions for stacking
      gb_stacking <- gbm::predict.gbm(
        object = gb_opt_poststrat_only,
        n.trees = gb.opt$n_trees,
        type = "response"
      )
    }

    # model for EBMA
    models$gb <- gb_opt_ebma

  }

  # Classifier 5: SVM
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

    # Prepare data
    svm_data <- data %>%
      dplyr::mutate_at(.vars = y, .funs = list(y_svm = ~as.factor(.))) %>%
      dplyr::select(
        y_svm, dplyr::all_of(c(y, L1.x, svm.L2.x, L2.unit, L2.reg))
      )

    svm_data_no_ebma <- data %>%
      dplyr::bind_rows(ebma.fold) %>%
      dplyr::mutate_at(.vars = y, .funs = list(y_svm = ~as.factor(.))) %>%
      dplyr::select(
        y_svm, dplyr::all_of(c(y, L1.x, svm.L2.x, L2.unit, L2.reg))
      )

    # Create model formula
    x <- paste(c(L1.x, svm.L2.x, svm.L2.unit, svm.L2.reg), collapse = " + ")
    form_svm <- as.formula(paste("y_svm ~ ", x, sep = ""))

    # Fit optimal model for EBMA
    svm_opt_ebma <- svm_classifier(
      y = "y_svm",
      form = form_svm,
      data = svm_data,
      kernel = svm.opt$kernel,
      type = "C-classification",
      probability = TRUE,
      svm.gamma = svm.opt$gamma,
      svm.cost = svm.opt$cost,
      verbose = verbose
    )

    # Fit optimal model for post-stratification w/o EBMA
    svm_opt_poststrat_only <- svm_classifier(
      y = "y_svm",
      form = form_svm,
      data = svm_data_no_ebma,
      kernel = kernel,
      type = "C-classification",
      probability = TRUE,
      svm.gamma = svm.opt$gamma,
      svm.cost = svm.opt$cost,
      verbose = verbose
    )

    # post-stratification step 1: only L2 units in training data
    svm_preds <- census %>%
      dplyr::filter(
        dplyr::pull(.data = census, !!L2.unit) %in%
          dplyr::pull(.data = svm_data_no_ebma, !!L2.unit)
      )

    # post-stratification step 2: predictions depending the DV
    if (
      svm_data_no_ebma %>%
        dplyr::pull(var = "y_svm") %>%
        unique() %>%
        length() > 2
    ) {
      # continuous DV
      svm_preds <- svm_preds %>%
        dplyr::mutate(
          svm = stats::predict(
            object = svm_opt_poststrat_only,
            newdata = .,
            allow.new.levels = TRUE,
            type = "response"
          )
        )
    } else {
      # binary DV
      svm_preds <- svm_preds %>%
        dplyr::mutate(
          svm = attr(stats::predict(
            object = svm_opt_poststrat_only,
            newdata = .,
            allow.new.levels = TRUE,
            probability = TRUE
          ), "probabilities")[, "1"]
        )
    }

    # add svm predictions to census if stacking
    if (stacking) {

      # post-stratification step 1: only L2 units in training data
      svm_stack_census <- census %>%
        dplyr::filter(
          dplyr::pull(.data = census, !!L2.unit) %in%
            dplyr::pull(.data = svm_data_no_ebma, !!L2.unit)
        )

      # stacking svm predictions on the linear predictor scale
      if (
        svm_stack_census %>%
          dplyr::pull(var = "y_svm") %>%
          unique() %>%
          length() > 2
      ) {
        # continuous DV
        census_stack <- svm_stack_census %>%
          dplyr::mutate(
            svm = stats::predict(
              object = svm_opt_poststrat_only,
              newdata = .,
              allow.new.levels = TRUE
            )
          )
      } else {
        # binary DV
        svm_preds <- svm_preds %>%
          dplyr::mutate(
            svm = attr(stats::predict(
              object = svm_opt_poststrat_only,
              newdata = .,
              allow.new.levels = TRUE,
              probability = FALSE
            ), "probabilities")[, "1"]
          )
      }
      census <- svm_preds
    }

    # post-stratification step 3: weighted mean
    svm_preds <- svm_preds %>%
      dplyr::group_by(!!rlang::sym(L2.unit)) %>%
      dplyr::summarize(
        svm = stats::weighted.mean(x = svm, w = prop), .groups = "keep"
      ) %>%
      dplyr::ungroup()

    # individual level predictions for EBMA
    if (
      svm_data_no_ebma %>%
        dplyr::pull(var = "y_svm") %>%
        unique() %>%
        length() > 2
    ) {
      svm_ind <- stats::predict(
        object = svm_opt_ebma,
        newdata = svm_data,
        type = "response"
      )
    } else {
      svm_ind <- attr(stats::predict(
        object = svm_opt_ebma,
        newdata = svm_data,
        probability = TRUE
      ), "probabilities")[, "1"]
    }

    # individual level predictions for stacking
    if (stacking) {
      # individual level predictions for EBMA
      if (
        svm_data_no_ebma %>%
          dplyr::pull(var = "y_svm") %>%
          unique() %>%
          length() > 2
      ) {
        svm_stacking <- stats::predict(
          object = svm_opt_poststrat_only,
          newdata = svm_data_no_ebma,
          type = "response"
        )
      } else {
        svm_stacking <- attr(stats::predict(
          object = svm_opt_poststrat_only,
          newdata = svm_data_no_ebma,
          allow.new.levels = TRUE,
          probability = TRUE
        ), "probabilities")[, "1"]
      }

    }

    # model for EBMA
    models$svm <- svm_opt_ebma
  }

  # Classifier 6: MRP
  # Fit model
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

    # fit model for EBMA
    mrp_model_ebma <- best_subset_classifier(
      y = y,
      model = form_mrp,
      data.train = data,
      model.family = binomial(link = "probit"),
      model.optimizer = "bobyqa",
      n.iter = 1000000,
      verbose = verbose
    )

    # fit MrP model for post-stratification only
    mrp_model_poststrat_only <- best_subset_classifier(
      y = y,
      model = form_mrp,
      data.train = no_ebma_data,
      model.family = binomial(link = "probit"),
      model.optimizer = "bobyqa",
      n.iter = 1000000,
      verbose = verbose
    )

    # post-stratification
    mrp_preds <- census %>%
      dplyr::mutate(
        mrp = stats::predict(
          object = mrp_model_poststrat_only,
          newdata = .,
          allow.new.levels = TRUE,
          type = "response"
        )
      ) %>%
      dplyr::group_by(!! rlang::sym(L2.unit)) %>%
      dplyr::summarize(
        mrp = stats::weighted.mean(x = mrp, w = prop), .groups = "keep"
      ) %>%
      dplyr::ungroup()

    # individual level predictions for EBMA
    mrp_ind <- stats::predict(object = mrp_model_ebma, type = "response")

    # individual level predictions for stacking
    if (stacking) {

      # add mrp predictions to census
      census <- census %>%
        dplyr::mutate(
          mrp = stats::predict(
            object = mrp_model_poststrat_only,
            newdata = .,
            allow.new.levels = TRUE,
            type = "response"
          )
        )

      # individual level predictions for stacking
      mrp_stacking <- stats::predict(
        object = mrp_model_poststrat_only,
        newdata = no_ebma_data,
        type = "response"
      )
    }


    # model for EBMA
    models$mrp <- mrp_model_ebma
  }

  # Classifier 7: deep MRP
  # Fit model
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
        pattern = stringr::fixed(pattern = names(data))
      ) %>%
        .[!is.na(.)]

      # take each column of data and combine its values into a single string
      df_x <- data %>%
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
    data <- dplyr::bind_cols(data, x_data)

    # loop over all interactions for data no_ebma_data object
    x_no_ebma <- lapply(all_interactions, function(x) {

      # break down interaction components
      y <- stringr::str_extract(
        string = x,
        pattern = stringr::fixed(pattern = names(no_ebma_data))
      ) %>%
        .[!is.na(.)]

      # take each column of data and combine its values into a single string
      df_x <- no_ebma_data %>%
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
    no_ebma_data <- dplyr::bind_cols(no_ebma_data, x_no_ebma)

    # run deep mrp model for ebma
    deep_mrp_model <- deep_mrp_classifier(
      y = y,
      form = form,
      data = data,
      verbose = verbose
    )

    # run deep mrp model for postratification only
    deep_mrp_model_poststrat_only <- deep_mrp_classifier(
      y = y,
      form = form,
      data = no_ebma_data,
      verbose = verbose
    )

    # loop over all interactions for census data
    x_census <- lapply(all_interactions, function(x) {

      # break down interaction components
      y <- stringr::str_extract(
        string = x,
        pattern = stringr::fixed(pattern = names(census))
      ) %>%
        .[!is.na(.)]

      # take each column of data and combine its values into a single string
      df_x <- census %>%
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
    deep_census <- dplyr::bind_cols(census, x_census)

    # Determine type of dependent variable
    if (
      data %>%
        dplyr::pull(!!y) %>%
        unique() %>%
        length() == 2
    ) {
      dv_type <- "binary"
    } else {
      dv_type <- "linear"
    }

    # binary or continuous DV
    if (dv_type == "binary") {
      # predictions for post-stratification only (no EBMA)
      pred_d <- vglmer::predict_MAVB(
        samples = 1000,
        deep_mrp_model_poststrat_only,
        newdata = deep_census,
        allow_missing_levels = TRUE
      )[["mean"]]

      # convert to response probabilities
      pred_d <- stats::plogis(pred_d)

    } else if (dv_type == "linear") {
      # predictions for post-stratification only (no EBMA)
      pred_d <- predict(
        samples = 1000,
        object = deep_mrp_model_poststrat_only,
        newdata = deep_census,
        allow_missing_levels = TRUE
      )[["mean"]]
    }

    # individual level predictions for stacking
    if (stacking) {

      # add deep mrp predictions to census
      census <- census %>%
        dplyr::mutate(deep_mrp = pred_d)

      # data for stacking
      deep_data_no_ebma <- data %>%
        dplyr::bind_rows(ebma.fold) %>%
        dplyr::select(
          dplyr::all_of(c(y, L1.x, svm.L2.x, L2.unit, L2.reg))
        )

      # loop over all interactions for data_stack
      x2_data <- lapply(all_interactions, function(x) {

        # break down interaction components
        y <- stringr::str_extract(
          string = x,
          pattern = stringr::fixed(pattern = names(deep_data_no_ebma))
        ) %>%
          .[!is.na(.)]

        # take each column of data and combine its values into a single string
        df_x <- deep_data_no_ebma %>%
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
      x2_data <- dplyr::bind_cols(deep_data_no_ebma, x2_data)
      rm(deep_data_no_ebma)

      if (dv_type == "binary") {
        deep_stacking <- vglmer::predict_MAVB(
          samples = 1000,
          deep_mrp_model_poststrat_only,
          newdata = x2_data,
          allow_missing_levels = TRUE
        )[["mean"]]

        # convert to response probabilities
        deep_stacking <- stats::plogis(deep_stacking)

      } else if (dv_type == "linear") {

        # predictions for stacking
        deep_stacking <- predict(
          samples = 1000,
          object = deep_mrp_model_poststrat_only,
          newdata = x2_data,
          allow_missing_levels = TRUE
        )[["mean"]]
      }

    }

    # post-stratification
    deep_preds <- deep_census %>%
      dplyr::mutate(deep_mrp = pred_d) %>%
      dplyr::group_by(!!rlang::sym(L2.unit)) %>%
      dplyr::summarize(
        deep_mrp = stats::weighted.mean(
          x = deep_mrp, w = prop
        ), .groups = "keep"
      )

    # binary or continuous DV
    if (dv_type == "binary") {
      # individual level predictions for EBMA
      deep_mrp_ind <- vglmer::predict_MAVB(
        samples = 1000,
        deep_mrp_model,
        newdata = data,
        allow_missing_levels = TRUE
      )[["mean"]]

      # convert response to probabilities
      deep_mrp_ind <- stats::plogis(deep_mrp_ind)
    } else if (dv_type == "linear") {
      # individual level predictions for EBMA
      deep_mrp_ind <- predict(
        samples = 1000,
        object = deep_mrp_model,
        newdata = data,
        allow_missing_levels = TRUE
      )[["mean"]]
    }

    # model for EBMA
    models$deep_mrp <- deep_mrp_model

  } # end of deep.mrp

  # stacking:
  if (stacking) {

    # stacking data
    data_stack <- tibble::tibble(
      y = no_ebma_data[[y]]
    )
    if (exists("bs_stacking")) {
      data_stack <- dplyr::mutate(.data = data_stack, best_subset = bs_stacking)
    }
    if (exists("lasso_stacking")) {
      data_stack <- dplyr::mutate(
        .data = data_stack, lasso = as.numeric(lasso_stacking)
      )
    }
    if (exists("pca_stacking")) {
      data_stack <- dplyr::mutate(.data = data_stack, pca = pca_stacking)
    }
    if (exists("gb_stacking")) {
      data_stack <- dplyr::mutate(.data = data_stack, gb = gb_stacking)
    }
    if (exists("svm_stacking")) {
      data_stack <- dplyr::mutate(.data = data_stack, svm = svm_stacking)
    }
    if (exists("mrp_stacking")) {
      data_stack <- dplyr::mutate(.data = data_stack, mrp = mrp_stacking)
    }
    if (exists("deep_stacking")) {
      data_stack <- dplyr::mutate(.data = data_stack, deep_mrp = deep_stacking)
    }

    # get stacking weights
    stack_weights(s_data = data_stack)


    # stack model
    stack_out <- glm(
      formula = form_stack,
      data = data_stack,
      family = binomial(link = "probit")
    )

    # stacking weights
    stack_weights <- stack_out$coefficients

  }

  # --------------------------- combine l2 level predictions ------------------

  # tibble of L2 units
  L2_preds <- dplyr::select(census, one_of(L2.unit)) %>%
    dplyr::distinct()
  # add existing classifiers
  if (exists("bs_preds")) {
    L2_preds <- dplyr::left_join(
      x = L2_preds,
      y = bs_preds,
      by = L2.unit
    )
  }
  if (exists("lasso_preds")) {
    L2_preds <- dplyr::left_join(
      x = L2_preds,
      y = lasso_preds,
      by = L2.unit
    )
  }
  if (exists("pca_preds")) {
    L2_preds <- dplyr::left_join(
      x = L2_preds,
      y = pca_preds,
      by = L2.unit
    )
  }
  if (exists("gb_preds")) {
    L2_preds <- dplyr::left_join(
      x = L2_preds,
      y = gb_preds,
      by = L2.unit
    )
  }
  if (exists("svm_preds")) {
    L2_preds <- dplyr::left_join(
      x = L2_preds,
      y = svm_preds,
      by = L2.unit
    )
  }
  if (exists("mrp_preds")) {
    L2_preds <- dplyr::left_join(
      x = L2_preds,
      y = mrp_preds,
      by = L2.unit
    )
  }
  if (exists("deep_preds")) {
    L2_preds <- dplyr::left_join(
      x = L2_preds,
      y = deep_preds,
      by = L2.unit
    )
  }

  ###########################
  # current
  # stacking combination
  if (stacking) {
    L2_stacking <- cbind(
      !!colnames(L2_preds)[1] == L2_preds[, 1], # L2 ids
      rep(1, nrow(L2_preds))
    )
  }

  # individual predictions for EBMA
  L1_preds <- data %>%
    dplyr::select(one_of(y)) %>%
    tidyr::drop_na() %>%
    # add best subset
    dplyr::mutate(
      best_subset = if (exists("bs_ind")) {
        as.numeric(bs_ind)
      } else {
        NA
      },
      # add lasso
      lasso = if (exists("lasso_ind")) {
        as.numeric(lasso_ind)
      } else {
        NA
      },
      # add pca
      pca = if (exists("pca_ind")) {
        as.numeric(pca_ind)
      } else {
        NA
      },
      # add gb
      gb = if (exists("gb_ind")) {
        as.numeric(gb_ind)
      } else {
        NA
      },
      # add svm
      svm = if (exists("svm_ind")) {
        as.numeric(svm_ind)
      } else {
        NA
      },
      # add MrP
      mrp = if (exists("mrp_ind")) {
        as.numeric(mrp_ind)
      } else {
        NA
      },
      # add deep MrP
      deep_mrp = if (exists("deep_mrp_ind")) {
        as.numeric(deep_mrp_ind)
      } else {
        NA
      }
    )

  # remove NA's
  L1_preds <- L1_preds[, apply(
    X = L1_preds, MARGIN = 2, FUN = function(x) {
      all(!is.na(x))
    }
  )]

  # Function output
  return(
    ps = list(
      predictions = list(
        Level1 = L1_preds,
        Level2 = L2_preds
      ),
      models = models
    )
  )
}
