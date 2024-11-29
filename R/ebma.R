#' Bayesian Ensemble Model Averaging EBMA
#'
#' \code{ebma} tunes EBMA and generates weights for classifier averaging.
#'
#' @inheritParams auto_MrP
#' @param ebma.fold New data for EBMA tuning. A list containing the the data
#'   that must not have been used in classifier training.
#' @param pc.names Principal Component Variable names. A character vector
#'   containing the names of the context-level principal components variables.
#' @param post.strat Post-stratification results. A list containing the best
#'   models for each of the tuned classifiers, the individual level predictions
#'   on the data classifier trainig data and the post-stratified context-level
#'   predictions.
#' @param n.draws EBMA number of samples. An integer-valued scalar specifying
#'   the number of bootstrapped samples to be drawn from the EBMA fold and used
#'   for tuning EBMA. Default is \eqn{100}. Passed on from \code{ebma.n.draws}.
#' @param tol EBMA tolerance. A numeric vector containing the tolerance values
#'   for improvements in the log-likelihood before the EM algorithm stops
#'   optimization. Values should range at least from \eqn{0.01} to \eqn{0.001}.
#'   Default is \code{c(0.01, 0.005, 0.001, 0.0005, 0.0001, 0.00005, 0.00001)}.
#'   Passed on from \code{ebma.tol}.
#' @param best.subset.opt Tuned best subset parameters. A list returned from
#'   \code{run_best_subset()}.
#' @param pca.opt Tuned best subset with principal components parameters. A list
#'   returned from \code{run_pca()}.
#' @param lasso.opt Tuned lasso parameters. A list returned from
#'   \code{run_lasso()}.
#' @param gb.opt Tuned gradient tree boosting parameters. A list returned from
#'   \code{run_gb()}.
#' @param svm.opt Tuned support vector machine parameters. A list returned from
#'   \code{run_svm()}.
#' @param knn.opt Tuned k-nearest neighbors parameters. A list returned from
#'   \code{run_knn()}.

ebma <- function(
  ebma.fold, cv.data, y, L1.x, L2.x, L2.unit, L2.reg, pc.names, post.strat,
  n.draws, tol, best.subset.opt, best.subset.deep, pca.opt, pca.deep,
  lasso.opt, gb.opt, svm.opt, knn.opt, verbose, cores, preds_all, knn.kernel
) {

  # Run EBMA if at least two classifiers selected
  if (
    sum(
      unlist(lapply(X = post.strat$models, FUN = function(x) !is.null(x)))
    ) > 1
  ) {

    if (verbose) {
      cli::cli_progress_step(
        "Performing ensemble Bayesian model averaging"
      )
    }

    # dependent variable type
    if (
      ebma.fold %>%
        dplyr::pull(!!y) %>%
        unique() %>%
        length() == 2
    ) {
      dv_type <- "binary"
    } else {
      dv_type <- "linear"
    }

    # EBMA w/o L2.x variables
    if (all(L2.x == "")) L2.x <- NULL

    # Models
    model_bs <- post.strat$models$best_subset
    model_pca <- post.strat$models$pca
    model_lasso <- post.strat$models$lasso
    model_gb <- post.strat$models$gb
    model_svm <- post.strat$models$svm
    model_knn <- post.strat$models$knn
    model_mrp <- post.strat$models$mrp
    model_deep <- post.strat$models$deep

    # CV data to NULL if not needed
    if (is.null(model_knn)) {
      cv.data <- NULL
    } else {
      if (dv_type == "binary") {
        cv.data <- cv.data %>%
          dplyr::mutate(!!rlang::sym(y) := as.factor(!!rlang::sym(y)))
      }
    }

    # Training predictions
    train_preds <- post.strat$predictions$Level1 %>%
      dplyr::select(-dplyr::one_of(y))

    # Training set outcomes
    train_y <- dplyr::pull(.data = post.strat$predictions$Level1, var = y)

    # Distribute clusters over tolerance values or n.draws
    if (length(tol) <= n.draws * 3) {

      final_model_weights <- ebma_mc_draws(
        train.preds = train_preds,
        train.y = train_y,
        ebma.fold = ebma.fold,
        cv.data = cv.data,
        y = y,
        L1.x = L1.x,
        L2.x = L2.x,
        L2.unit = L2.unit,
        L2.reg = L2.reg,
        pc.names = pc.names,
        model.bs = model_bs,
        model.pca = model_pca,
        model.lasso = model_lasso,
        model.gb = model_gb,
        model.svm = model_svm,
        model.knn = model_knn,
        model.mrp = model_mrp,
        model.deep = model_deep,
        tol = tol,
        n.draws = n.draws,
        cores = cores,
        preds_all = preds_all,
        post.strat = post.strat,
        dv_type = dv_type,
        best.subset.deep = best.subset.deep,
        pca.deep = pca.deep,
        knn.opt = knn.opt,
        knn.kernel = knn.kernel
      )
      # unlist
      ebma_preds <- final_model_weights$ebma_preds
      final_model_weights <- final_model_weights$final_model_weights

    } else {

      final_model_weights <- ebma_mc_tol(
        train.preds = train_preds,
        train.y = train_y,
        ebma.fold = ebma.fold,
        cv.data = cv.data,
        y = y,
        L1.x = L1.x,
        L2.x = L2.x,
        L2.unit = L2.unit,
        L2.reg = L2.reg,
        pc.names = pc.names,
        model.bs = model_bs,
        model.pca = model_pca,
        model.lasso = model_lasso,
        model.gb = model_gb,
        model.svm = model_svm,
        model.mrp = model_mrp,
        tol = tol,
        n.draws = n.draws,
        cores = cores,
        preds_all = preds_all,
        post.strat = post.strat,
        dv_type = dv_type,
        best.subset.deep = best.subset.deep,
        pca.deep = pca.deep,
        knn.opt = knn.opt,
        knn.kernel = knn.kernel
      )
      # unlist
      ebma_preds <- final_model_weights$ebma_preds
      final_model_weights <- final_model_weights$final_model_weights
    }

    # weighted average
    model_preds <- as.matrix(
      post.strat$predictions$Level2[, names(final_model_weights)]
    )
    w_avg <- as.numeric(model_preds %*% final_model_weights)

    # deal with missings
    if (any(is.na(w_avg))) {

      # missing row indexes
      m_idx <- which(is.na(w_avg))

      # which classifiers contain missings
      NA_classifiers <- which(apply(model_preds, 2, function(x) any(is.na(x))))

      # readjust weights without that classifier
      NA_adj_weights <- final_model_weights[-NA_classifiers] / sum(
        final_model_weights[-NA_classifiers]
      )

      # remove columns with NAs
      no_Na_preds <- model_preds[
        m_idx, # rows
        apply(X = model_preds, MARGIN = 2, FUN = function(x) {
          all(!is.na(x))
        }) # columns
      ]

      # predictions for rows with NAs on at least 1 classifier
      no_NA_avg <- as.numeric(no_Na_preds %*% NA_adj_weights)

      # replace predictions
      w_avg[m_idx] <- no_NA_avg
    }

    # L2 preds object
    L2_preds <- dplyr::tibble(
      !! rlang::sym(L2.unit) := dplyr::pull(
        .data = post.strat$predictions$Level2, var = L2.unit
      ),
      ebma = w_avg
    )

    # function output
    return(
      list(
        ebma = L2_preds,
        classifiers = post.strat$predictions$Level2,
        weights = final_model_weights,
        individual_level_predictions = ebma_preds
      )
    )

  } else {
    if (verbose) {
      message("\n Skipping EBMA (only 1 classifier selected) \n")
    }

    # function output
    return(
      list(
        ebma = "EBMA step skipped (only 1 classifier run)",
        classifiers = post.strat$predictions$Level2
      )
    )
  }
}


################################################################################
#     Multicore tuning for EBMA - parallel over tolerance                      #
################################################################################
#' EBMA multicore tuning - parallelises over tolerance values.
#'
#' \code{ebma_mc_tol} is called from within \code{ebma}. It tunes using
#' multiple cores.
#'
#' @inheritParams auto_MrP
#' @inheritParams ebma
#' @param train.preds Predictions of classifiers on the classifier training
#'   data. A tibble.
#' @param train.y Outcome variable of the classifier training data. A numeric
#'   vector.
#' @param ebma.fold The data used for EBMA tuning. A tibble.
#' @param model.bs The tuned model from the multilevel regression with best
#'   subset selection classifier. An \code{\link[lme4]{glmer}} object.
#' @param model.pca The tuned model from the multilevel regression with
#'   principal components as context-level predictors classifier. An
#'   \code{\link[lme4]{glmer}} object.
#' @param model.lasso The tuned model from the multilevel regression with L1
#'   regularization classifier. A \code{\link[glmmLasso]{glmmLasso}} object.
#' @param model.gb The tuned model from the gradient boosting classifier. A
#'   \code{\link[gbm]{gbm}} object.
#' @param model.svm The tuned model from the support vector machine classifier.
#'   An \code{\link[e1071]{svm}} object.
#' @param model.mrp The standard MrP model. An \code{\link[lme4]{glmer}} object
#' @param tol The tolerance values used for EBMA. A numeric vector.
#' @param dv_type The type of the depenedent variable. A character string.
#'   Either "binary" or "linear".
#' @return The classifier weights. A numeric vector.
#' @examples \dontrun{
#' # not yet
#' }

ebma_mc_tol <- function(
  train.preds, train.y, ebma.fold, cv.data, y, L1.x, L2.x, L2.unit, L2.reg,
  pc.names, model.bs, best.subset.deep, model.pca, pca.deep, model.lasso,
  model.gb, model.svm, model.knn, model.mrp, model.deep, tol, n.draws, cores, preds_all, post.strat, dv_type, knn.opt, knn.kernel
) {

  # Binding for global variables
  `%>%` <- dplyr::`%>%`
  idx.tol <- NULL

  # list of draws
  draws <- vector("list", n.draws)

  # loop over Ndraws with equal obs/state
  for (idx.Ndraws in 1:n.draws) {

    # Determine number per group to sample
    n_per_group <- as.integer(nrow(ebma.fold) / length(
      levels(ebma.fold[[L2.unit]])
    ))

    # Test set with n_per_group persons per state (with resampling)
    test <- ebma.fold %>%
      dplyr::group_by_at(.vars = L2.unit) %>%
      dplyr::sample_n(n_per_group, replace = TRUE) %>%
      dplyr::ungroup() %>%
      dplyr::mutate_at(.vars = c(L1.x, L2.unit, L2.reg), .funs = as.factor) %>%
      tidyr::drop_na()

    # predict outcomes in test set
    test_preds <- dplyr::tibble(
      best_subset = if (!is.null(model.bs)) {
        # regular best subset without level 1 interactions
        if (!best.subset.deep) {
          predict(
            object = model.bs,
            newdata = test,
            type = "response",
            allow.new.levels = TRUE
          )
        } else {
          # best subset with level 1 interactions
          if (dv_type == "binary") {
            vglmer::predict_MAVB(
              samples = 1000,
              model.bs,
              newdata = test,
              allow_missing_levels = TRUE
            )[["mean"]] %>%
              stats::plogis()
          } else if (dv_type == "linear") {
            vglmer::predict_MAVB(
              samples = 1000,
              model.bs,
              newdata = test,
              allow_missing_levels = TRUE
            )[["mean"]]
          }
        }
      } else {
        NA
      },
      pca = if (!is.null(model.pca)) {
        # regular pca without level 1 interactions
        if (!pca.deep) {
          predict(
            object = model.pca,
            newdata = test,
            type = "response",
            allow.new.levels = TRUE
          )
        } else {
          # pca with level 1 interactions
          if (dv_type == "binary") {
            vglmer::predict_MAVB(
              samples = 1000,
              model.pca,
              newdata = test,
              allow_missing_levels = TRUE
            )[["mean"]] %>%
              stats::plogis()
          } else {
            vglmer::predict_MAVB(
              samples = 1000,
              model.pca,
              newdata = test,
              allow_missing_levels = TRUE
            )[["mean"]]
          }
        }
      } else {
        NA
      },
      lasso = if (!is.null(model.lasso)) {
        predict_glmmLasso(
          census = test,
          m = model.lasso,
          L1.x = L1.x,
          lasso.L2.x = L2.x,
          L2.unit = L2.unit,
          L2.reg = L2.reg
        )
      } else {
        NA
      },
      gb = if (!is.null(model.gb)) {
        gbm::predict.gbm(
          object = model.gb,
          newdata = test,
          n.trees = model.gb$n.trees,
          type = "response"
        )
      } else {
        NA
      },
      svm = if (!is.null(model.svm)) {
        # binary DV
        if (length(unique(test[[y]])) == 2) {
          as.numeric(
            attr(
              predict(
                object = model.svm,
                newdata = test,
                probability = TRUE
              ), "probabilities"
            )[, "1"]
          )
        } else {
          # continuous DV
          predict(
            object = model.svm,
            newdata = test
          )
        }
      } else {
        NA
      },
      knn = if (!is.null(model.knn)) {
        if (dv_type == "binary") {
          c_pred <- knn_classifier(
            y = y,
            form = formula(model.knn$terms),
            data.train = cv.data,
            data.valid = test,
            knn.k.value = knn.opt,
            knn.kernel = knn.kernel
          )
          c_pred$prob[, "1"]
        } else {
          c_pred <- knn_classifier(
            y = y,
            form = formula(model.knn$terms),
            data.train = cv.data,
            data.valid = test,
            knn.k.value = knn.opt,
            knn.kernel = knn.kernel
          )
          c_pred$fit
        }
      } else {
        NA
      },
      mrp = if (!is.null(model.mrp)) {
        predict(
          object = model.mrp,
          newdata = test,
          type = "response",
          allow.new.levels = TRUE
        )
      } else {
        NA
      },
      deep = if (!is.null(model.deep)) {
        # binary DV
        if (length(unique(test[[y]])) == 2) {
          # predictions for post-stratification only (no EBMA)
          pred_d <- vglmer::predict_MAVB(
            samples = 1000,
            model.deep,
            newdata = test,
            allow_missing_levels = TRUE
          )[["mean"]]
          # convert to response probabilities
          pred_d <- stats::plogis(pred_d)
        } else {
          # continuous DV
          pred_d <- predict(
            samples = 1000,
            object = model.deep,
            newdata = test,
            allow_missing_levels = TRUE
          )[["mean"]]
        }
      } else {
        NA
      }
    )

    # remove NA's
    test_preds <- test_preds[, apply(
      X = test_preds, MARGIN = 2, FUN = function(x) {
        all(!is.na(x))
      }
    )]

    # Outcome on the test
    test_y <- dplyr::select(.data = test, dplyr::one_of(y))

    # Register cores
    cl <- multicore(cores = cores, type = "open", cl = NULL)

    # Loop over tolerance values
    ebma_tune <- foreach::foreach(
      idx.tol = seq_along(tol),
      .packages = c("glmmLasso", "e1071", "gbm", "vglmer", "kknn")
    ) %dorng% {

      # EBMA
      forecast.data <- suppressWarnings(
        EBMAforecast::makeForecastData(
          .predCalibration = data.frame(train.preds),
          .outcomeCalibration = as.numeric(unlist(train.y)),
          .predTest = data.frame(test_preds),
          .outcomeTest = as.numeric(unlist(test_y))
        )
      )

      # vector of initial model weights
      # Note: On some Mac versions the model weights do not sum to exactly
      # 1 for repeating decimal weights
      W <- rep(
        x = 1 / ncol(forecast.data@predCalibration),
        times = ncol(forecast.data@predCalibration)
      )
      W[length(W)] <- 1 - sum(W[-length(W)])

      # calibrate EBMA ensemble
      forecast.out <- EBMAforecast::calibrateEnsemble(
        forecast.data,
        model = "normal",
        useModelParams = FALSE,
        tol = tol[idx.tol],
        W = W
      )

      # mse
      mse_collector <- mean((
        as.numeric(unlist(test_y)) - as.numeric(
          attributes(forecast.out)$predTest[, 1, 1]
        )
      )^2)

      # model weights
      weights_box <- matrix(
        data = attributes(forecast.out)$modelWeights,
        nrow = 1,
        ncol = ncol(train.preds)
      )

      return(list(
        MSEs = mse_collector,
        weights = weights_box,
        EBMA_model = forecast.out
      ))
    }

    # De-register cluster
    multicore(cores = cores, type = "close", cl = cl)

    # MSEs for each tolerance value
    MSEs <- unlist(lapply(ebma_tune, `[[`, 1))

    # weights
    weights <- lapply(ebma_tune, `[[`, 2)
    weights <- matrix(
      data = unlist(weights),
      ncol = ncol(train.preds),
      byrow = TRUE,
      dimnames = list(
        c(paste0("Tol_", tol)),
        c(names(train.preds))
      )
    )

    # models for each tolerance value
    EBMA_models <- lapply(ebma_tune, `[[`, 3)


    # Store MSEs and weights in current draw
    draws[[idx.Ndraws]] <- list(
      MSEs = MSEs, weights = weights, EBMA_models = EBMA_models
    )
  }

  # Container of best weights per draw
  wgt <- matrix(
    data = NA,
    nrow = n.draws,
    ncol = ncol(draws[[1]][["weights"]])
  )

  # container of best models per draw
  EBMA_models <- vector("list", n.draws)

  # generate EBMA predictions stacking
  train_preds <- as.matrix(preds_all[, -c(1, 2)])

  # Select best weights per draw
  for (idx_w in seq_along(draws)) {

    # Tolerance with lowest MSE in draw idx_w
    best_tolerance <- which(draws[[idx_w]]$MSEs == min(draws[[idx_w]]$MSEs))

    if (length(best_tolerance) > 1) {
      wgt[idx_w, ] <- apply(draws[[idx_w]]$weights[best_tolerance, ], 2, mean)
    } else {
      wgt[idx_w, ] <- draws[[idx_w]]$weights[best_tolerance, ]
    }

    # loop over best models
    for (idx_model in seq_along(draws[[idx_w]]$EBMA_models[best_tolerance])) {

      # predictions from the EBMA model on the individual level data
      ebma_preds <- EBMAforecast::EBMApredict(
        EBMAmodel = draws[[idx_w]]$EBMA_models[best_tolerance][[idx_model]],
        Predictions = train_preds,
        Outcome = as.numeric(preds_all$y)
      )

      # extract predictions
      ebma_preds <- ebma_preds@predTest
      if (idx_model == 1) {
        individual_preds <- as.matrix(ebma_preds[, "EBMA", 1])
      } else {
        individual_preds <- cbind(
          individual_preds,
          as.matrix(ebma_preds[, "EBMA", 1])
        )
      }
    }
  }

  # Average model weights
  final_model_weights <- apply(wgt, 2, mean)
  names(final_model_weights) <- names(train.preds)

  # average EBMA individual level predictions
  ebma_preds <- apply(individual_preds, 1, mean)

  # L2 preds object
  L2_preds <- dplyr::tibble(
    !! rlang::sym(L2.unit) := dplyr::pull(
      .data = post.strat$predictions$Level2, var = L2.unit
    ),
    ebma = w_avg
  )

  # Function output
  return(
    list(
      final_model_weights = final_model_weights,
      ebma_preds = tibble::tibble(
        y = as.numeric(preds_all$y),
        !!rlang::sym(L2.unit) := preds_all %>% dplyr::pull(var = L2.unit),
        ebma_preds = ebma_preds
      )
    )
  )

}


################################################################################
#             Multicore tuning for EBMA - parallel over the draws              #
################################################################################
#' EBMA multicore tuning - parallelises over draws.
#'
#' \code{ebma_mc_draws} is called from within \code{ebma}. It tunes using
#' multiple cores.
#'
#' @inheritParams auto_MrP
#' @inheritParams ebma
#' @inheritParams ebma_mc_tol
#' @return The classifier weights. A numeric vector.

ebma_mc_draws <- function(
  train.preds, train.y, ebma.fold, cv.data, y, L1.x, L2.x, L2.unit, L2.reg,
  pc.names, model.bs, model.pca, model.lasso, model.gb, model.svm, model.knn,
  model.mrp, model.deep, tol, n.draws, cores, preds_all, post.strat, dv_type,
  best.subset.deep, pca.deep, knn.kernel, knn.opt
) {

  # Binding for global variables
  `%>%` <- dplyr::`%>%`

  # Register cores
  cl <- multicore(cores = cores, type = "open", cl = NULL)

  ebma_tune <- foreach::foreach(
    idx.Ndraws = 1:n.draws,
    .packages = c("glmmLasso", "e1071", "gbm", "vglmer", "kknn"),
    .export = "predict_glmmLasso"
  ) %dorng% {

    # Determine number per group to sample
    n_per_group <- as.integer(
      nrow(ebma.fold) / length(levels(ebma.fold[[L2.unit]]))
    )

    # Test set with n_per_group persons per state (with resampling)
    test <- ebma.fold %>%
      dplyr::group_by_at(.vars = L2.unit) %>%
      dplyr::sample_n(n_per_group, replace = TRUE) %>%
      dplyr::ungroup() %>%
      dplyr::mutate_at(.vars = c(L1.x, L2.unit, L2.reg), .funs = as.factor) %>%
      tidyr::drop_na()

    # predict outcomes in test set
    test_preds <- dplyr::tibble(
      best_subset = if (!is.null(model.bs)) {
        # regular best subset without level 1 interactions
        if (!best.subset.deep) {
          predict(
            object = model.bs,
            newdata = test,
            type = "response",
            allow.new.levels = TRUE
          )
        } else {
          # best subset with level 1 interactions
          if (dv_type == "binary") {
            vglmer::predict_MAVB(
              samples = 1000,
              model.bs,
              newdata = test,
              allow_missing_levels = TRUE
            )[["mean"]] %>%
              stats::plogis()
          } else if (dv_type == "linear") {
            vglmer::predict_MAVB(
              samples = 1000,
              model.bs,
              newdata = test,
              allow_missing_levels = TRUE
            )[["mean"]]
          }
        }
      } else {
        NA
      },
      pca = if (!is.null(model.pca)) {
        # regular pca without level 1 interactions
        if (!pca.deep) {
          predict(
            object = model.pca,
            newdata = test,
            type = "response",
            allow.new.levels = TRUE
          )
        } else {
          # pca with level 1 interactions
          if (dv_type == "binary") {
            vglmer::predict_MAVB(
              samples = 1000,
              model.pca,
              newdata = test,
              allow_missing_levels = TRUE
            )[["mean"]] %>%
              stats::plogis()
          } else {
            vglmer::predict_MAVB(
              samples = 1000,
              model.pca,
              newdata = test,
              allow_missing_levels = TRUE
            )[["mean"]]
          }
        }
      } else {
        NA
      },
      lasso = if (!is.null(model.lasso)) {
        predict_glmmLasso(
          census = test,
          m = model.lasso,
          L1.x = L1.x,
          lasso.L2.x = L2.x,
          L2.unit = L2.unit,
          L2.reg = L2.reg
        )
      } else {
        NA
      },
      gb = if (!is.null(model.gb)) {
        gbm::predict.gbm(
          object = model.gb,
          newdata = test,
          n.trees = model.gb$n.trees,
          type = "response"
        )
      } else {
        NA
      },
      svm = if (!is.null(model.svm)) {
        # binary DV
        if (length(unique(test[[y]])) == 2) {
          as.numeric(
            attr(
              predict(
                object = model.svm,
                newdata = test,
                probability = TRUE
              ), "probabilities"
            )[, "1"]
          )
        } else {
          # continuous DV
          predict(
            object = model.svm,
            newdata = test
          )
        }
      } else {
        NA
      },
      knn = if (!is.null(model.knn)) {
        if (dv_type == "binary") {
          c_pred <- knn_classifier(
            y = y,
            form = formula(model.knn$terms),
            data.train = cv.data,
            data.valid = test,
            knn.k.value = knn.opt,
            knn.kernel = knn.kernel
          )
          c_pred$prob[, "1"]
        } else {
          c_pred <- knn_classifier(
            y = y,
            form = formula(model.knn$terms),
            data.train = cv.data,
            data.valid = test,
            knn.k.value = knn.opt,
            knn.kernel = knn.kernel
          )
          c_pred$fit
        }
      } else {
        NA
      },
      mrp = if (!is.null(model.mrp)) {
        predict(
          object = model.mrp,
          newdata = test,
          type = "response",
          allow.new.levels = TRUE
        )
      } else {
        NA
      },
      deep = if (!is.null(model.deep)) {
        # binary DV
        if (length(unique(test[[y]])) == 2) {
          # predictions for post-stratification only (no EBMA)
          pred_d <- vglmer::predict_MAVB(
            samples = 1000,
            model.deep,
            newdata = test,
            allow_missing_levels = TRUE
          )[["mean"]]
          # convert to response probabilities
          pred_d <- stats::plogis(pred_d)
        } else {
          # continuous DV
          pred_d <- predict(
            samples = 1000,
            object = model.deep,
            newdata = test,
            allow_missing_levels = TRUE
          )[["mean"]]
        }
      } else {
        NA
      }
    )

    # remove NA's
    test_preds <- test_preds[, apply(
      X = test_preds, MARGIN = 2, FUN = function(x) {
        all(!is.na(x))
      }
    )]

    # Outcome on the test
    test_y <- dplyr::select(.data = test, dplyr::one_of(y))

    # Container of MSEs per tolerance value and draw combination
    mse_collector <- NA

    # Container of model weights over tolerances in current draw
    weights_box <- matrix(
      data = NA,
      nrow = length(tol),
      ncol = ncol(train.preds),
      dimnames = list(
        c(paste0("Tol_", tol)),
        c(names(train.preds))
      )
    )

    # container of models
    model_box <- vector("list", length(tol))

    # Loop over tolerance values
    for (idx.tol in seq_along(tol)) {

      # EBMA
      forecast.data <- suppressWarnings(
        EBMAforecast::makeForecastData(
          .predCalibration = data.frame(train.preds),
          .outcomeCalibration = as.numeric(unlist(train.y)),
          .predTest = data.frame(test_preds),
          .outcomeTest = as.numeric(unlist(test_y))
        )
      )

      # vector of initial model weights
      # Note: On some Mac versions the model weights do not sum to exactly
      # 1 for repeating decimal weights
      W <- rep(
        x = 1 / ncol(forecast.data@predCalibration),
        times = ncol(forecast.data@predCalibration)
      )
      W[length(W)] <- 1 - sum(W[-length(W)])

      # calibrate EBMA ensemble
      forecast.out <- EBMAforecast::calibrateEnsemble(
        forecast.data,
        model = "normal",
        useModelParams = FALSE,
        tol = tol[idx.tol],
        W = W
      )

      # mse
      mse_collector[idx.tol] <- mean((
        as.numeric(unlist(test_y)) - as.numeric(
          attributes(forecast.out)$predTest[, 1, 1]
        )
      )^2)

      # model weights
      weights_box[idx.tol, ] <- attributes(forecast.out)$modelWeights

      # model box
      model_box[[idx.tol]] <- forecast.out
    }

    return(list(
      MSEs = mse_collector,
      weights = weights_box,
      EBMA_model = model_box
    ))

  }

  # De-register cluster
  multicore(cores = cores, type = "close", cl = cl)

  # Container of best weights per draw
  wgt <- matrix(
    data = NA,
    nrow = n.draws,
    ncol = ncol(ebma_tune[[1]][["weights"]])
  )

  # generate EBMA predictions stacking
  train_preds <- as.matrix(preds_all[, -c(1, 2)])

  # Select best weights per draw
  for (idx_w in seq_along(ebma_tune)) {

    # Tolerance with lowest MSE in draw idx_w
    best_tolerance <- which(
      ebma_tune[[idx_w]]$MSEs == min(ebma_tune[[idx_w]]$MSEs)
    )

    if (length(best_tolerance) > 1) {
      wgt[idx_w, ] <- apply(
        ebma_tune[[idx_w]]$weights[best_tolerance, ], 2, mean
      )
    } else {
      wgt[idx_w, ] <- ebma_tune[[idx_w]]$weights[best_tolerance, ]
    }

    # loop over best models
    for (
      idx_model in seq_along(ebma_tune[[idx_w]]$EBMA_model[best_tolerance])
    ) {

      # predictions from the EBMA model on the individual level data
      ebma_preds <- EBMAforecast::EBMApredict(
        EBMAmodel = ebma_tune[[idx_w]]$EBMA_model[best_tolerance][[idx_model]],
        Predictions = train_preds,
        Outcome = as.numeric(preds_all$y)
      )

      # extract predictions
      ebma_preds <- ebma_preds@predTest
      if (idx_model == 1) {
        individual_preds <- as.matrix(ebma_preds[, "EBMA", 1])
      } else {
        individual_preds <- cbind(
          individual_preds,
          as.matrix(ebma_preds[, "EBMA", 1])
        )
      }
    }
  }

  # Average model weights
  final_model_weights <- apply(wgt, 2, mean)
  names(final_model_weights) <- names(train.preds)

  # average EBMA individual level predictions
  ebma_preds <- apply(individual_preds, 1, mean)

  # Function output
  return(
    list(
      final_model_weights = final_model_weights,
      ebma_preds = tibble::tibble(
        y = as.numeric(preds_all$y),
        !!rlang::sym(L2.unit) := preds_all %>% dplyr::pull(var = L2.unit),
        ebma_preds = ebma_preds
      )
    )
  )
}
