#' Optimal individual classifiers
#'
#' \code{run_classifiers} tunes classifiers, post-stratifies and carries out EMBA.
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

run_classifiers <- function(y, L1.x, L2.x, mrp.L2.x, L2.unit, L2.reg,
                            L2.x.scale, pcs, pc.names, folds, bin.proportion,
                            bin.size, cv.folds, cv.data, ebma.fold, census, ebma.size,
                            ebma.n.draws, k.folds, cv.sampling, loss.unit,
                            loss.fun, best.subset, lasso, pca, gb, svm, mrp,
                            forward.select, best.subset.L2.x,
                            lasso.L2.x, pca.L2.x, gb.L2.x, svm.L2.x,
                            gb.L2.unit, gb.L2.reg, svm.L2.unit, svm.L2.reg,
                            lasso.lambda, lasso.n.iter, gb.interaction.depth,
                            gb.shrinkage, gb.n.trees.init,
                            gb.n.trees.increase, gb.n.trees.max,
                            gb.n.minobsinnode, svm.kernel,
                            svm.gamma, svm.cost, ebma.tol, cores, verbose) {

  # Classifier 1: Best Subset
  if (isTRUE(best.subset)) {
    if (verbose) {
      message("Starting multilevel regression with best subset selection classifier tuning")
    }

    # Determine context-level covariates
    if (is.null(best.subset.L2.x)) {
      best.subset.L2.x <- L2.x
    }

    # Run classifier
    best_subset_out <- run_best_subset(y = y,
                                       L1.x = L1.x,
                                       L2.x = best.subset.L2.x,
                                       L2.unit = L2.unit,
                                       L2.reg = L2.reg,
                                       loss.unit = loss.unit,
                                       loss.fun = loss.fun,
                                       data = cv.folds,
                                       verbose = verbose,
                                       cores = cores)
  } else {
    best_subset_out <- NULL
  }

  # Classifier 2: Lasso

  # message if Lasso is set to TRUE but no context level variables provided
  if (isTRUE(lasso) & is.null(L2.x)) {
    if (verbose) {
      message('Lasso requires L2.x variables to be specified. Skipping Lasso.')
    }
  }

  if (isTRUE(lasso) & !is.null(L2.x)) {

    if (verbose) {
      message("Starting multilevel regression with L1 regularization tuning")
    }

    # Determine context-level covariates
    if (is.null(lasso.L2.x)) {
      lasso.L2.x <- L2.x
    }

    # Run classifier
    lasso_out <- run_lasso(y = y,
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
                           cores = cores)
  } else {
    lasso_out <- NULL
  }

  # Classifier 3: PCA

  # message if pca is TRUE but no level 2 variables or pc_names provided
  if (isTRUE(pca) & is.null(pca.L2.x)) {
    message(paste0('PCA requires that L2.x variables are specified or alternatively',
                   ' that the pcs argument is specified.'))
  }
  if (isTRUE(pca) & !is.null(pca.L2.x)) {

    if (verbose) {
      message("Starting multilevel regression with principal components as context level variables tuning")
    }

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
      cores = cores)

  } else {
    pca_out <- NULL
  }

  # Classifier 4: GB
  if (isTRUE(gb)) {

    if (verbose) {
      message("Starting gradient tree boosting tuning")
    }

    # Determine context-level covariates
    if (is.null(gb.L2.x)) {
      gb.L2.x <- L2.x
    }

    # Evaluate inclusion of L2.unit in GB
    if (isTRUE(gb.L2.unit)) {
      gb.L2.unit <- L2.unit
    } else {
      gb.L2.unit <- NULL
    }

    # Evaluate inclusion of L2.reg in GB
    if (isTRUE(gb.L2.reg)) {
      gb.L2.reg <- L2.reg
    } else {
      gb.L2.reg <- NULL
    }

    # Run classifier
    gb_out <- run_gb(y = y,
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
                     verbose = verbose)
  } else {
    gb_out <- NULL
  }

  # Classifier 5: SVM
  if ( isTRUE(svm) ) {

    if (verbose) {
      message("Starting support vector machine tuning")
    }

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
      cores = cores)
  } else {
    svm_out <- NULL
  }


  # Post-stratification -----------------------------------------------------

  if (verbose) {
    message("Starting post-stratification")
  }

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
    data = cv.data,
    ebma.fold = ebma.fold,
    census = census,
    verbose = verbose
  )


  # EBMA --------------------------------------------------------------------


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
    pc.names = pc.names,
    verbose = verbose,
    cores = cores
  )
}
