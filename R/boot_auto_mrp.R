boot_auto_mrp <- function(y, L1.x, L2.x, mrp.L2.x, L2.unit, L2.reg,
                          L2.x.scale, pcs, folds, bin.proportion,
                          bin.size, survey, census, ebma.size,
                          k.folds, cv.sampling, loss.unit, loss.fun,
                          best.subset, lasso, pca, gb, svm, mrp,
                          forward.select, best.subset.L2.x,
                          lasso.L2.x, pca.L2.x, gb.L2.x, svm.L2.x,
                          gb.L2.unit, gb.L2.reg, lasso.lambda,
                          lasso.n.iter, gb.interaction.depth,
                          gb.shrinkage, gb.n.trees.init,
                          gb.n.trees.increase, gb.n.trees.max,
                          gb.n.iter, gb.n.minobsinnode, svm.kernel,
                          svm.gamma, svm.cost, ebma.tol, seed){

  # Binding for global variables
  `%>%` <- dplyr::`%>%`

  # Register cores
  cl <- multicore(cores = cores, type = "open", cl = NULL)

  # Bootstrap iterations
  boot_out <- foreach::foreach(idx_boot = 1:boot.iter, .packages = "autoMrP") %dopar% {

    # Bootstrapped survey sample
    boot_sample <- dplyr::sample_n(tbl = survey, size = nrow(survey), replace = TRUE)

    # # ------------------------------- Create folds -------------------------------
    #
    # if (is.null(folds)) {
    #   # EBMA hold-out fold
    #   ebma.size <- round(nrow(survey) * ebma.size, digits = 0)
    #
    #   if(ebma.size>0){
    #     ebma_folding_out <- ebma_folding(data = survey,
    #                                      L2.unit = L2.unit,
    #                                      ebma.size = ebma.size)
    #     ebma_fold <- ebma_folding_out$ebma_fold
    #     cv_data <- ebma_folding_out$cv_data
    #   } else{
    #     ebma_fold <- NULL
    #     cv_data <- survey
    #   }
    #
    #   # K folds for cross-validation
    #   cv_folds <- cv_folding(data = cv_data,
    #                          L2.unit = L2.unit,
    #                          k.folds = k.folds,
    #                          cv.sampling = cv.sampling)
    # } else {
    #
    #   if (ebma.size > 0){
    #     # EBMA hold-out fold
    #     ebma_fold <- survey %>%
    #       dplyr::filter_at(dplyr::vars(dplyr::one_of(folds)),
    #                        dplyr::any_vars(. == k.folds + 1))
    #   }
    #
    #   # K folds for cross-validation
    #   cv_data <- survey %>%
    #     dplyr::filter_at(dplyr::vars(dplyr::one_of(folds)),
    #                      dplyr::any_vars(. != k.folds + 1))
    #
    #   cv_folds <- cv_data %>%
    #     dplyr::group_split(.data[[folds]])
    # }
    #
    # # ---------------------- Optimal individual classifiers ----------------------
    #
    # # Classifier 1: Best Subset
    # if (isTRUE(best.subset)) {
    #
    #   message("Starting multilevel regression with best subset selection classifier tuning")
    #
    #   # Determine context-level covariates
    #   if (is.null(best.subset.L2.x)) {
    #     best.subset.L2.x <- L2.x
    #   }
    #
    #   # Run classifier
    #   best_subset_out <- run_best_subset(y = y,
    #                                      L1.x = L1.x,
    #                                      L2.x = best.subset.L2.x,
    #                                      L2.unit = L2.unit,
    #                                      L2.reg = L2.reg,
    #                                      loss.unit = loss.unit,
    #                                      loss.fun = loss.fun,
    #                                      data = cv_folds,
    #                                      verbose = verbose,
    #                                      cores = cores)
    # } else {
    #   best_subset_out <- NULL
    # }

    # Estimate on 1 sample in autoMrP
    boot_mrp <- auto_MrP(
      survey = boot_sample,
      ebma.n.draws = 1,
      uncertainty = FALSE,
      verbose = FALSE,
      cores = 1,
      y = y,
      L1.x = L1.x,
      L2.x = L2.x,
      mrp.L2.x = mrp.L2.x,
      L2.unit = L2.unit,
      L2.reg = L2.reg,
      L2.x.scale = L2.x.scale,
      pcs = pcs,
      folds = folds,
      bin.proportion = bin.proportion,
      bin.size = bin.size,
      census = census,
      ebma.size = ebma.size,
      k.folds = k.folds,
      cv.sampling = cv.sampling,
      loss.unit = loss.unit,
      loss.fun = loss.fun,
      best.subset = best.subset,
      lasso = lasso,
      pca = pca,
      gb = gb,
      svm = svm,
      mrp = mrp,
      forward.select = forward.select,
      best.subset.L2.x = best.subset.L2.x,
      lasso.L2.x = lasso.L2.x,
      pca.L2.x = pca.L2.x,
      gb.L2.x = gb.L2.x,
      svm.L2.x = svm.L2.x,
      gb.L2.unit = gb.L2.unit,
      gb.L2.reg = gb.L2.reg,
      lasso.lambda = lasso.lambda,
      lasso.n.iter = lasso.n.iter,
      gb.interaction.depth = gb.interaction.depth,
      gb.shrinkage = gb.shrinkage,
      gb.n.trees.init = gb.n.trees.init,
      gb.n.trees.increase = gb.n.trees.increase,
      gb.n.trees.max = gb.n.trees.max,
      gb.n.iter = gb.n.iter,
      gb.n.minobsinnode = gb.n.minobsinnode,
      svm.kernel = svm.kernel,
      svm.gamma = svm.gamma,
      svm.cost = svm.cost,
      ebma.tol = ebma.tol,
      seed = seed
    )

  }

  # De-register cluster
  multicore(cores = cores, type = "close", cl = cl)

  return(boot_out)

}
