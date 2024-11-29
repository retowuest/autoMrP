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
#' @param knn.opt Optimal tuning parameters from k-nearest neighbors
#'   classifier A list returned by \code{knn_svm()}.
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
  y, L1.x, L2.x, L2.unit, L2.reg, best.subset.opt, best.subset.deep,
  lasso.opt, lasso.L2.x, pca.opt, pca.deep, gb.opt, svm.opt, knn.opt,
  svm.L2.reg, svm.L2.unit, svm.L2.x, knn.L2.reg, knn.L2.unit, knn.L2.x,
  mrp.include, n.minobsinnode, L2.unit.include, L2.reg.include, kernel,
  knn.kernel, mrp.L2.x, data, ebma.fold, census, verbose, deep.mrp,
  deep.splines, deep.L2.x, deep.L2.unit, deep.L2.reg
) {

  # globals
  lasso <- NULL
  pca <- NULL
  gb <- NULL
  svm <- NULL
  knn <- NULL
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

  # Determine type of dependent variable
  if (
    no_ebma_data %>%
      dplyr::pull(!!y) %>%
      unique() %>%
      length() == 2
  ) {
    dv_type <- "binary"
  } else {
    dv_type <- "linear"
  }

  # ----- Fit optimal model and make prediction for individual classifiers -----

  # Classifier 1: Best Subset
  if (!is.null(best.subset.opt)) {

    # deep MrP or not
    if (best.subset.deep) {

      # Fit optimal model for EBMA
      best_subset_opt_ebma <- deep_mrp_classifier(
        form = best.subset.opt,
        y = y,
        data = data,
        verbose = TRUE
      )

      # Fit optimal model for post-stratification w/o EBMA
      best_subset_opt_poststrat_only <- deep_mrp_classifier(
        form = best.subset.opt,
        y = y,
        data = no_ebma_data,
        verbose = TRUE
      )

      if (dv_type == "binary") {

        # post-stratification
        bs_preds <- census %>%
          dplyr::mutate(
            best_subset = vglmer::predict_MAVB(
              samples = 1000,
              best_subset_opt_poststrat_only,
              newdata = census,
              allow_missing_levels = TRUE
            )[["mean"]] %>%
              stats::plogis()
          ) %>%
          dplyr::group_by(!! rlang::sym(L2.unit)) %>%
          dplyr::summarize(best_subset = stats::weighted.mean(
            x = best_subset, w = prop
          ), .groups = "keep") %>%
          dplyr::ungroup()

        # individual level predictions for EBMA
        bs_ind <- vglmer::predict_MAVB(
          samples = 1000,
          best_subset_opt_ebma,
          newdata = data,
          allow_missing_levels = TRUE
        )[["mean"]] %>%
          stats::plogis()

        # model for EBMA
        models$best_subset <- best_subset_opt_ebma

      } else {

        # post-stratification
        bs_preds <- census %>%
          dplyr::mutate(
            best_subset = vglmer::predict_MAVB(
              samples = 1000,
              best_subset_opt_poststrat_only,
              newdata = census,
              allow_missing_levels = TRUE
            )[["mean"]]
          ) %>%
          dplyr::group_by(!! rlang::sym(L2.unit)) %>%
          dplyr::summarize(best_subset = stats::weighted.mean(
            x = best_subset, w = prop
          ), .groups = "keep") %>%
          dplyr::ungroup()

        # individual level predictions for EBMA
        bs_ind <- vglmer::predict_MAVB(
          samples = 1000,
          best_subset_opt_ebma,
          newdata = data,
          allow_missing_levels = TRUE
        )[["mean"]]

        # model for EBMA
        models$best_subset <- best_subset_opt_ebma
      }

    } else {

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

      # model for EBMA
      models$best_subset <- best_subset_opt_ebma
    }
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

    # model for EBMA
    models$lasso <- lasso_opt_ebma

  }

  # Classifier 3: PCA
  if (!is.null(pca.opt)) {

    # deep MrP or not
    if (pca.deep) {

      # Fit optimal model for EBMA
      pca_opt_ebma <- deep_mrp_classifier(
        form = pca.opt,
        y = y,
        data = data,
        verbose = TRUE
      )

      # Fit optimal model for for post-stratification w/o EBMA
      pca_opt_poststrat_only <- deep_mrp_classifier(
        form = pca.opt,
        y = y,
        data = no_ebma_data,
        verbose = TRUE
      )

      # dependent variable type
      if (dv_type == "binary") {

        # post-stratification
        pca_preds <- census %>%
          dplyr::mutate(
            pca = vglmer::predict_MAVB(
              samples = 1000,
              pca_opt_poststrat_only,
              newdata = census,
              allow_missing_levels = TRUE
            )[["mean"]] %>%
              stats::plogis()
          ) %>%
          dplyr::group_by(!! rlang::sym(L2.unit)) %>%
          dplyr::summarize(
            pca = stats::weighted.mean(x = pca, w = prop), .groups = "keep"
          ) %>%
          dplyr::ungroup()

        # individual level predictions for EBMA
        pca_ind <- vglmer::predict_MAVB(
          samples = 1000,
          pca_opt_ebma,
          newdata = data,
          allow_missing_levels = TRUE
        )[["mean"]] %>%
          stats::plogis()

      } else {

        # post-stratification
        pca_preds <- census %>%
          dplyr::mutate(
            pca = vglmer::predict_MAVB(
              samples = 1000,
              pca_opt_poststrat_only,
              newdata = census,
              allow_missing_levels = TRUE
            )[["mean"]]
          ) %>%
          dplyr::group_by(!! rlang::sym(L2.unit)) %>%
          dplyr::summarize(
            pca = stats::weighted.mean(x = pca, w = prop), .groups = "keep"
          ) %>%
          dplyr::ungroup()

        # individual level predictions for EBMA
        pca_ind <- vglmer::predict_MAVB(
          samples = 1000,
          pca_opt_ebma,
          newdata = data,
          allow_missing_levels = TRUE
        )[["mean"]]

      }

      # model for EBMA
      models$pca <- pca_opt_ebma

    } else {

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

      # model for EBMA
      models$pca <- pca_opt_ebma
    }
  }

  # Classifier 4: GB
  if (!is.null(gb.opt)) {

    # Evaluate inclusion of L2.unit
    if (isTRUE(L2.unit.include)) {
      L2.unit.gb <- L2.unit
    } else {
      L2.unit.gb <- NULL
    }

    # Evaluate inclusion of L2.reg
    if (isTRUE(L2.reg.include)) {
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

    # model for EBMA
    models$gb <- gb_opt_ebma

  }

  # Classifier 5: SVM
  if (!is.null(svm.opt)) {

    # Determine context-level covariates
    if (is.null(svm.L2.x)) {
      svm.L2.x <- L2.x
    }

    # Evaluate inclusion of L2.unit in SVM
    if (isTRUE(svm.L2.unit)) {
      svm.L2.unit <- L2.unit
    } else {
      svm.L2.unit <- NULL
    }

    # Evaluate inclusion of L2.reg in SVM
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

    # model for EBMA
    models$svm <- svm_opt_ebma
  }

  # Classifier 6: KNN
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

    # Fit optimal model for EBMA
    knn_opt_ebma <- knn_classifier(
      y = y,
      form = form_knn,
      data.train = if (dv_type == "binary") {
        data %>%
          dplyr::mutate(!!rlang::sym(y) := as.factor(!!rlang::sym(y)))
      } else {
        data
      },
      data.valid = data,
      knn.k.value = knn.opt,
      knn.kernel = knn.kernel,
      verbose = verbose
    )

    # Fit optimal model for post-stratification w/o EBMA
    knn_opt_poststrat_only <- knn_classifier(
      y = y,
      form = form_knn,
      data.train = if (dv_type == "binary") {
        no_ebma_data %>%
          dplyr::mutate(!!rlang::sym(y) := as.factor(!!rlang::sym(y)))
      } else {
        no_ebma_data
      },
      data.valid = census,
      knn.k.value = knn.opt,
      knn.kernel = knn.kernel,
      verbose = verbose
    )

    # post-stratification
    knn_preds <- census %>%
      dplyr::mutate(
        knn = if (dv_type == "binary") {
          kknn:::predict.kknn(
            object = knn_opt_poststrat_only,
            type = "prob"
          )[, "1"]
        } else {
          knn_opt_poststrat_only$fit
        }
      ) %>%
      dplyr::group_by(!! rlang::sym(L2.unit)) %>%
      dplyr::summarize(
        knn = stats::weighted.mean(x = knn, w = prop), .groups = "keep"
      ) %>%
      dplyr::ungroup()

    # individual level predictions for EBMA
    knn_ind <- if (dv_type == "binary") {
      kknn:::predict.kknn(
        object = knn_opt_ebma,
        type = "prob"
      )[, "1"]
    } else {
      knn_opt_ebma$fit
    }

    # model for EBMA
    models$knn <- knn_opt_ebma

  }

  # Classifier 7: MRP
  # Fit model
  if (isTRUE(mrp.include == TRUE)) {

    mrp_start_time <- Sys.time()

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

    # model for EBMA
    models$mrp <- mrp_model_ebma
    mrp_end_time <- Sys.time()
    mrp_runtime <- difftime(
      time1 = mrp_end_time, time2 = mrp_start_time, units = "mins"
    )
  } else {
    mrp_runtime <- NULL
  }

  # Classifier 8: deep MRP
  # Fit model
  if (deep.mrp) {

    deep_start_time <- Sys.time()

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
        newdata = census,
        allow_missing_levels = TRUE
      )[["mean"]]

      # convert to response probabilities
      pred_d <- stats::plogis(pred_d)

    } else if (dv_type == "linear") {
      # predictions for post-stratification only (no EBMA)
      pred_d <- predict(
        samples = 1000,
        object = deep_mrp_model_poststrat_only,
        newdata = census,
        allow_missing_levels = TRUE
      )[["mean"]]
    }

    # post-stratification
    deep_preds <- census %>%
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
    deep_end_time <- Sys.time()
    deep_runtime <- difftime(
      time1 = deep_end_time, time2 = deep_start_time, units = "mins"
    )

  } else {
    deep_runtime <- NULL
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
  if (exists("knn_preds")) {
    L2_preds <- dplyr::left_join(
      x = L2_preds,
      y = knn_preds,
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
      # add knn
      knn = if (exists("knn_ind")) {
        as.numeric(knn_ind)
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
      models = models,
      mrp_runtime = mrp_runtime,
      deep_mrp_runtime = deep_runtime
    )
  )
}
