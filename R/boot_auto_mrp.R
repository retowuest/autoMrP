#' Bootstrappinng wrapper for auto_mrp
#'
#' \code{boot_auto_mrp} estimates uncertainty for auto_mrp via botstrapping.
#'
#' @inheritParams auto_MrP
#' @param pc.names A character vector of the principal component variable names
#'   in the data.

boot_auto_mrp <- function(y, L1.x, L2.x, mrp.L2.x, L2.unit, L2.reg,
                          L2.x.scale, pcs, folds, bin.proportion,
                          bin.size, survey, census, ebma.size,
                          k.folds, cv.sampling, loss.unit, loss.fun,
                          best.subset, lasso, pca, gb, svm, mrp,
                          forward.select, best.subset.L2.x,
                          lasso.L2.x, pca.L2.x, pc.names, gb.L2.x, svm.L2.x,
                          svm.L2.unit, svm.L2.reg,gb.L2.unit, gb.L2.reg,
                          lasso.lambda, lasso.n.iter, gb.interaction.depth,
                          gb.shrinkage, gb.n.trees.init,
                          gb.n.trees.increase, gb.n.trees.max,
                          gb.n.minobsinnode, svm.kernel,
                          svm.gamma, svm.cost, ebma.tol,
                          boot.iter, cores){

  # Binding for global variables
  `%>%` <- dplyr::`%>%`
  idx_boot <- NULL

  # Register cores
  cl <- multicore(cores = cores, type = "open", cl = NULL)

  # Bootstrap iterations
  boot_out <- foreach::foreach(idx_boot = 1:boot.iter, .packages = "autoMrP") %dorng% {

    boot_fun(
      y = y, L1.x = L1.x, L2.x = L2.x, mrp.L2.x = mrp.L2.x,
      L2.unit = L2.unit, L2.reg = L2.reg, pcs = pcs,
      folds = folds, survey = survey, census = census, k.folds = k.folds,
      cv.sampling = cv.sampling, ebma.size = ebma.size,
      loss.unit = loss.unit, loss.fun = loss.fun,
      best.subset = best.subset, lasso = lasso, pca = pca,
      gb = gb, svm = svm, mrp = mrp, forward.select = forward.select,
      best.subset.L2.x = best.subset.L2.x,
      lasso.L2.x = lasso.L2.x, pca.L2.x = pca.L2.x, pc.names = pc.names,
      gb.L2.x = gb.L2.x, svm.L2.x = svm.L2.x, svm.L2.unit = svm.L2.unit,
      svm.L2.reg = svm.L2.reg, gb.L2.unit = gb.L2.unit, gb.L2.reg = gb.L2.reg,
      lasso.lambda = lasso.lambda, lasso.n.iter = lasso.n.iter,
      gb.interaction.depth = gb.interaction.depth,
      gb.shrinkage = gb.shrinkage,
      gb.n.trees.init = gb.n.trees.init,
      gb.n.trees.increase = gb.n.trees.increase,
      gb.n.trees.max = gb.n.trees.max,
      gb.n.minobsinnode = gb.n.minobsinnode,
      svm.kernel = svm.kernel, svm.gamma = svm.gamma,
      svm.cost = svm.cost, ebma.tol = ebma.tol, cores = cores,
      verbose = verbose)
  } # end of foreach loop

  # Median and standard deviation of EBMA estimates
  if ( !any(boot_out[[1]]$ebma == "EBMA step skipped (only 1 classifier run)") ){
    ebma <- base::do.call(base::rbind, base::do.call(base::rbind, boot_out)[,"ebma"] ) #%>%
    #dplyr::group_by(.dots = list(L2.unit)) %>%
    #dplyr::summarise(median = median(ebma, na.rm = TRUE),
    #                 sd = sd(ebma, na.rm = TRUE), .groups = "drop")

    # weights
    weights <- base::do.call(base::rbind, base::do.call(base::rbind, boot_out)[,"weights"] ) %>%
      dplyr::as_tibble() %>%
      dplyr::select(
        contains("best_subset"),
        contains("pca"),
        contains("lasso"),
        contains("gb"),
        contains("svm"),
        contains("mrp")
      )

  } else {
    ebma <- "EBMA step skipped (only 1 classifier run)"
    weights <- NULL
  }

  # Median and standard deviations for classifier estimates
  classifiers <- base::do.call(base::rbind, base::do.call(base::rbind, boot_out)[,"classifiers"] ) %>%
    dplyr::select(
      one_of(L2.unit),
      contains("best_subset"),
      contains("pca"),
      contains("lasso"),
      contains("gb"),
      contains("svm"),
      contains("mrp")
    )

  #dplyr::group_by(.dots = list(L2.unit)) %>%
  #dplyr::summarise_all(.funs = c(median = median, sd = sd), na.rm = TRUE) %>%
  #dplyr::select(
  #  state,
  #  contains("best_subset"),
  #  contains("pca"),
  #  contains("lasso"),
  #  contains("gb"),
  #  contains("svm"),
  #  contains("mrp")
  #)

  #dplyr::summarise_all(.funs = c(median = median, sd = sd), na.rm = TRUE) %>%
  #dplyr::select(
  #  contains("best_subset"),
  #  contains("pca"),
  #  contains("lasso"),
  #  contains("gb"),
  #  contains("svm"),
  #  contains("mrp")
  #)

  if ( !is.null(weights) ) {
    boot_out <- list(ebma = ebma, classifiers = classifiers, weights = weights)
  } else {
    boot_out <- list(ebma = ebma, classifiers = classifiers)
  }

  # De-register cluster
  multicore(cores = cores, type = "close", cl = cl)

  return(boot_out)

}
