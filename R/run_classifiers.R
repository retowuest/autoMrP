#' Optimal individual classifiers
#'
#' \code{run_classifiers} tunes classifiers, post-stratifies, and carries out
#' EBMA.
#'
#' @inheritParams auto_MrP
#' @param cv.folds Data for cross-validation. A \code{list} of \eqn{k}
#'   \code{data.frames}, one for each fold to be used in \eqn{k}-fold
#'   cross-validation.
#' @param cv.data A data.frame containing the survey data used in classifier
#'   training.
#' @param ebma.fold A data.frame containing the data not used in classifier
#'   training.
#' @param pc.names A character vector of the principal component variable names
#'   in the data.

run_classifiers <- function(
  y, L1.x, L2.x, mrp.L2.x, L2.unit, L2.reg, L2.x.scale, pcs, pc.names, folds,
  bin.proportion, bin.size, cv.folds, cv.data, ebma.fold, census, ebma.size,
  ebma.n.draws, k.folds, cv.sampling, loss.unit, loss.fun, best.subset,
  lasso, pca, gb, svm, knn, mrp, deep.mrp, best.subset.L2.x, lasso.L2.x, pca.L2.x,
  gb.L2.x, svm.L2.x, knn.L2.x, gb.L2.unit, gb.L2.reg, svm.L2.unit, svm.L2.reg,
  knn.L2.unit, knn.L2.reg, deep.splines, lasso.lambda, lasso.n.iter,
  gb.interaction.depth, gb.shrinkage, gb.n.trees.init, gb.n.trees.increase,
  gb.n.trees.max, gb.n.minobsinnode, svm.kernel, svm.gamma, svm.cost,
  knn.k.max, knn.k, knn.kernel, ebma.tol, cores, verbose
) {

  # Classifier 1: Best Subset
  if (isTRUE(best.subset)) {

    # get start time
    best_subset_start_time <- Sys.time()

    if (verbose) {
      cli::cli_progress_step(
        "Tuning multilevel regression with best subset selection classifier"
      )
    }

    # Determine context-level covariates
    if (is.null(best.subset.L2.x)) {
      best.subset.L2.x <- L2.x
    }

    # interactions of L1.x yes/no
    if (isTRUE(deep.mrp)) {
      # Run classifier with L1.x interactions
      best_subset_out <- run_deep_bs(
        y = y,
        L1.x = L1.x,
        L2.x = best.subset.L2.x,
        L2.unit = L2.unit,
        L2.reg = L2.reg,
        deep.splines = deep.splines,
        loss.unit = loss.unit,
        loss.fun = loss.fun,
        k.folds = k.folds,
        data = cv.folds,
        verbose = verbose,
        cores = cores
      )
    } else {
      # Run classifier without L1.x interactions
      best_subset_out <- run_best_subset(
        y = y,
        L1.x = L1.x,
        L2.x = best.subset.L2.x,
        L2.unit = L2.unit,
        L2.reg = L2.reg,
        loss.unit = loss.unit,
        loss.fun = loss.fun,
        data = cv.folds,
        verbose = verbose,
        cores = cores
      )
    }

    # get end time
    best_subset_end_time <- Sys.time()

    # best subset runtime
    best_subset_runtime <- difftime(
      time1 = best_subset_end_time, time2 = best_subset_start_time,
      units = "mins"
    )

  } else {
    best_subset_out <- NULL
    best_subset_runtime <- NULL
  }

  # Classifier 2: Lasso

  # message if Lasso is set to TRUE but no context level variables provided
  if (isTRUE(lasso) && is.null(L2.x)) {
    if (verbose) {
      message("Lasso requires L2.x variables to be specified. Skipping Lasso.")
    }
  }

  if (isTRUE(lasso) && !is.null(L2.x)) {

    # get start time
    lasso_start_time <- Sys.time()

    if (verbose) {
      cli::cli_progress_step(
        "Tuning multilevel regression with L1 regularization"
      )
    }

    # Determine context-level covariates
    if (is.null(lasso.L2.x)) {
      lasso.L2.x <- L2.x
    }

    # Run classifier
    lasso_out <- run_lasso(
      y = y,
      L1.x = L1.x,
      L2.x = lasso.L2.x,
      L2.unit = L2.unit,
      L2.reg = L2.reg,
      loss.unit = loss.unit,
      loss.fun = loss.fun,
      lambda = lasso.lambda,
      n.iter = lasso.n.iter,
      data = cv.folds,
      verbose = verbose,
      cores = cores
    )

    # get end time
    lasso_end_time <- Sys.time()

    # lasso runtime
    lasso_runtime <- difftime(
      time1 = lasso_end_time, time2 = lasso_start_time, units = "mins"
    )

  } else {
    lasso_out <- NULL
    lasso_runtime <- NULL
  }

  # Classifier 3: PCA

  # message if pca is TRUE but no level 2 variables or pc_names provided
  if (isTRUE(pca) && is.null(pca.L2.x)) {
    message(
      "PCA requires that L2.x variables are specified or alternatively",
      " that the pcs argument is specified."
    )
  }
  if (isTRUE(pca) && !is.null(pca.L2.x)) {

    # get start time
    pca_start_time <- Sys.time()

    if (verbose) {
      cli::cli_progress_step(
        paste0(
          "Tuning multilevel regression with principal components as context",
          " level variables"
        )
      )
    }

    # interactions of L1.x yes/no
    if (isTRUE(deep.mrp)) {
      # Run classifier with L1.x interactions
      pca_out <- run_deep_pca(
        y = y,
        L1.x = L1.x,
        L2.x = pc.names,
        L2.unit = L2.unit,
        L2.reg = L2.reg,
        deep.splines = deep.splines,
        loss.unit = loss.unit,
        loss.fun = loss.fun,
        data = cv.folds,
        verbose = verbose,
        cores = cores
      )
    } else {
      # run classifier without L1.x interactions
      pca_out <- run_pca(
        y = y,
        L1.x = L1.x,
        L2.x = pc.names,
        L2.unit = L2.unit,
        L2.reg = L2.reg,
        loss.unit = loss.unit,
        loss.fun = loss.fun,
        data = cv.folds,
        verbose = verbose,
        cores = cores
      )
    }

    # get end time
    pca_end_time <- Sys.time()

    # pca runtime
    pca_runtime <- difftime(
      time1 = pca_end_time, time2 = pca_start_time, units = "mins"
    )

  } else {
    pca_out <- NULL
    pca_runtime <- NULL
  }

  # Classifier 4: GB
  if (gb) {

    # get start time
    gb_start_time <- Sys.time()

    if (verbose) {
      cli::cli_progress_step("Tuning gradient tree boosting")
    }

    # Determine context-level covariates
    if (is.null(gb.L2.x)) {
      gb.L2.x <- L2.x
    }

    # Evaluate inclusion of L2.unit in GB
    if (gb.L2.unit) {
      gb.L2.unit <- L2.unit
    } else {
      gb.L2.unit <- NULL
    }

    # Evaluate inclusion of L2.reg in GB
    if (gb.L2.reg) {
      gb.L2.reg <- L2.reg
    } else {
      gb.L2.reg <- NULL
    }

    # Run classifier
    gb_out <- run_gb(
      y = y,
      L1.x = L1.x,
      L2.x = gb.L2.x,
      L2.eval.unit = L2.unit,
      L2.unit = gb.L2.unit,
      L2.reg = gb.L2.reg,
      loss.unit = loss.unit,
      loss.fun = loss.fun,
      interaction.depth = gb.interaction.depth,
      shrinkage = gb.shrinkage,
      n.trees.init = gb.n.trees.init,
      n.trees.increase = gb.n.trees.increase,
      n.trees.max = gb.n.trees.max,
      n.minobsinnode = gb.n.minobsinnode,
      data = cv.folds,
      cores = cores,
      verbose = verbose
    )

    # get end time
    gb_end_time <- Sys.time()

    # gb runtime
    gb_runtime <- difftime(
      time1 = gb_end_time, time2 = gb_start_time, units = "mins"
    )

  } else {
    gb_out <- NULL
    gb_runtime <- NULL
  }

  # Classifier 5: SVM
  if (isTRUE(svm)) {

    # get start time
    svm_start_time <- Sys.time()

    if (verbose) {
      cli::cli_progress_step("Tuning support vector machine")
    }

    # Determine context-level covariates
    if (is.null(svm.L2.x)) {
      svm.L2.x <- L2.x
    }

    # Evaluate inclusion of L2.unit in SVM
    if (svm.L2.unit) {
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

    # Run classifier
    svm_out <- run_svm(
      y = y,
      L1.x = L1.x,
      L2.x = svm.L2.x,
      L2.eval.unit = L2.unit,
      L2.unit = svm.L2.unit,
      L2.reg = svm.L2.reg,
      kernel = svm.kernel,
      loss.fun = loss.fun,
      loss.unit = loss.unit,
      gamma = svm.gamma,
      cost = svm.cost,
      data = cv.folds,
      verbose = verbose,
      cores = cores
    )

    # get end time
    svm_end_time <- Sys.time()

    # svm runtime
    svm_runtime <- difftime(
      time1 = svm_end_time, time2 = svm_start_time, units = "mins"
    )

  } else {
    svm_out <- NULL
    svm_runtime <- NULL
  }

  # Classifier 6: KNN
  if (isTRUE(knn)) {

    # get start time
    knn_start_time <- Sys.time()

    if (verbose) {
      cli::cli_progress_step("Tuning k-nearest neighbors")
    }

    # Determine context-level covariates
    if (is.null(knn.L2.x)) {
      knn.L2.x <- L2.x
    }

    # Evaluate inclusion of L2.unit in KNN
    if (knn.L2.unit) {
      knn.L2.unit <- L2.unit
    } else {
      knn.L2.unit <- NULL
    }

    # Evaluate inclusion of L2.reg in KNN
    if (knn.L2.reg) {
      knn.L2.reg <- L2.reg
    } else {
      knn.L2.reg <- NULL
    }

    # Run classifier
    knn_out <- run_knn(
      y = y,
      L1.x = L1.x,
      L2.x = knn.L2.x,
      L2.unit = knn.L2.unit,
      L2.reg = knn.L2.reg,
      loss.unit = loss.unit,
      loss.fun = loss.fun,
      knn.k.max = knn.k.max,
      knn.k = knn.k,
      knn.kernel = knn.kernel,
      data = cv.folds,
      verbose = verbose,
      cores = cores
    )

    # get end time
    knn_end_time <- Sys.time()

    # knn runtime
    knn_runtime <- difftime(
      time1 = knn_end_time, time2 = knn_start_time, units = "mins"
    )

  } else {
    knn_out <- NULL
    knn_runtime <- NULL
  }

  # Individual level predictions for all data -------------------------------

  if (verbose) {
    cli::cli_progress_step(
      "Generating out of sample predictions from tuned classifiers"
    )
  }

  # get start time
  preds_all_start_time <- Sys.time()

  preds_all <- suppressWarnings(
    suppressMessages(
      get_predictions(
        y = y,
        L1.x = L1.x,
        L2.x = L2.x,
        L2.unit = L2.unit,
        L2.reg = L2.reg,
        best.subset.opt = best_subset_out,
        lasso.opt = lasso_out,
        lasso.L2.x = lasso.L2.x,
        pca.opt = pca_out,
        gb.opt = gb_out,
        svm.opt = svm_out,
        svm.L2.reg = svm.L2.reg,
        svm.L2.unit = svm.L2.unit,
        svm.L2.x = svm.L2.x,
        knn.opt = knn_out,
        knn.L2.reg = knn.L2.reg,
        knn.L2.unit = knn.L2.unit,
        knn.L2.x = knn.L2.x,
        mrp.include = mrp,
        n.minobsinnode = gb.n.minobsinnode,
        L2.unit.include = gb.L2.unit,
        L2.reg.include = gb.L2.reg,
        kernel = svm.kernel,
        knn.kernel = knn.kernel,
        mrp.L2.x = mrp.L2.x,
        deep.mrp = deep.mrp,
        deep.splines = deep.splines,
        data = cv.folds,
        ebma.fold = ebma.fold,
        verbose = verbose,
        cv.sampling = cv.sampling,
        k.folds = k.folds,
        all_data = TRUE
      )
    )
  )

  # get end time
  preds_all_end_time <- Sys.time()

  # preds_all runtime
  preds_all_runtime <- difftime(
    time1 = preds_all_end_time, time2 = preds_all_start_time, units = "mins"
  )

  # Post-stratification -----------------------------------------------------

  if (verbose) {
    cli::cli_progress_step("Post-stratification")
  }

  # get start time
  ps_start_time <- Sys.time()

  ps_out <- post_stratification(
    y = y,
    L1.x = L1.x,
    L2.x = L2.x,
    L2.unit = L2.unit,
    L2.reg = L2.reg,
    best.subset.opt = best_subset_out,
    lasso.opt = lasso_out,
    lasso.L2.x = lasso.L2.x,
    pca.opt = pca_out,
    gb.opt = gb_out,
    svm.opt = svm_out,
    svm.L2.reg = svm.L2.reg,
    svm.L2.unit = svm.L2.unit,
    svm.L2.x = svm.L2.x,
    mrp.include = mrp,
    n.minobsinnode = gb.n.minobsinnode,
    L2.unit.include = gb.L2.unit,
    L2.reg.include = gb.L2.reg,
    kernel = svm.kernel,
    mrp.L2.x = mrp.L2.x,
    deep.mrp = deep.mrp,
    deep.splines = deep.splines,
    data = cv.data,
    ebma.fold = ebma.fold,
    census = census,
    verbose = verbose
  )

  # get start time
  ps_end_time <- Sys.time()

  # ps runtime
  ps_runtime <- difftime(
    time1 = ps_end_time, time2 = ps_start_time, units = "mins"
  )

  # EBMA --------------------------------------------------------------------

  # get start time
  ebma_start_time <- Sys.time()

  ebma_out <- ebma(
    ebma.fold = ebma.fold,
    y = y,
    L1.x = L1.x,
    L2.x = L2.x,
    L2.unit = L2.unit,
    L2.reg = L2.reg,
    post.strat = ps_out,
    n.draws = ebma.n.draws,
    tol = ebma.tol,
    best.subset.opt = best_subset_out,
    pca.opt = pca_out,
    lasso.opt = lasso_out,
    gb.opt = gb_out,
    svm.opt = svm_out,
    deep.mrp = deep.mrp,
    pc.names = pc.names,
    verbose = verbose,
    cores = cores,
    preds_all = preds_all
  )

  # get end time
  ebma_end_time <- Sys.time()

  # ebma runtime
  ebma_runtime <- difftime(
    time1 = ebma_end_time, time2 = ebma_start_time, units = "mins"
  )

  # Stacking  ----------------------------------------------------------------

  # get start time
  stack_start_time <- Sys.time()

  if (verbose) {
    cli::cli_progress_step("Stacking")
  }

  # get stacking weights
  stack_out <- autoMrP:::stacking_weights(
    preds = preds_all, ebma_out = ebma_out, L2.unit = L2.unit,
    k.folds = k.folds, cores = cores
  )

  # apply stacking weights
  ebma_out <- apply_stack_weights(
    ebma_out = ebma_out,
    stack_out = stack_out,
    L2.unit = L2.unit,
    y = y,
    preds_all = preds_all
  )

  # get end time
  stack_end_time <- Sys.time()

  # stack runtime
  stack_runtime <- difftime(
    time1 = stack_end_time, time2 = stack_start_time, units = "mins"
  )

  # Detailed runtime ---------------------------------------------------------
  runtime_detailed <- tibble::tibble(
    best_subset = best_subset_runtime,
    lasso = lasso_runtime,
    pca = pca_runtime,
    gb = gb_runtime,
    svm = svm_runtime,
    individual_level_predictions = preds_all_runtime,
    post_stratification = ps_runtime,
    ebma = ebma_runtime,
    stacking = stack_runtime
  )
  ebma_out$runtime <- runtime_detailed

  return(ebma_out)
}
