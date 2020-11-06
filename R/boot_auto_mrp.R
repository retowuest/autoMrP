#' Bootstrappinng wrapper for auto_mrp
#'
#' \code{boot_auto_mrp} estimates uncertainty for auto_mrp via botstrapping.
#'
#' @inheritParams auto_MrP

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
                          gb.n.minobsinnode, svm.kernel,
                          svm.gamma, svm.cost, ebma.tol, seed,
                          boot.iter, cores){

  # Binding for global variables
  `%>%` <- dplyr::`%>%`
  idx_boot <- NULL

  # Register cores
  cl <- multicore(cores = cores, type = "open", cl = NULL)

  # Bootstrap iterations
  boot_out <- foreach::foreach(idx_boot = 1:boot.iter,
                               .errorhandling = "remove",
                               .packages = "autoMrP") %dorng% {

    # Bootstrapped survey sample
    #boot_sample <- dplyr::slice_sample(.data = survey, n = base::nrow(survey), replace = TRUE)

    boot_sample <- survey %>%
      dplyr::group_by( !! rlang::sym(L2.unit) ) %>%
      tidyr::nest() %>%
      dplyr::mutate(data = purrr::map(data, function(x){
        data = dplyr::slice_sample(.data = x, n = nrow(x), replace = TRUE)
      })) %>%
      tidyr::unnest(data) %>%
      dplyr::ungroup()

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
      gb.n.minobsinnode = gb.n.minobsinnode,
      svm.kernel = svm.kernel,
      svm.gamma = svm.gamma,
      svm.cost = svm.cost,
      ebma.tol = ebma.tol,
      seed = seed,
      boot.iter = NULL
    )
                               }

  # Median and standard deviation of EBMA estimates
  if ( !any(boot_out[[1]]$ebma == "EBMA step skipped (only 1 classifier run)") ){
    ebma <- base::do.call(base::rbind, base::do.call(base::rbind, boot_out)[,"ebma"] ) #%>%
    #dplyr::group_by(.dots = list(L2.unit)) %>%
    #dplyr::summarise(median = median(ebma, na.rm = TRUE),
    #                 sd = sd(ebma, na.rm = TRUE), .groups = "drop")
  } else {
    ebma <- "EBMA step skipped (only 1 classifier run)"
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

  #dplyr::summarise_all(.funs = c(median = median, sd = sd), na.rm = TRUE) %>%
  #dplyr::select(
  #  contains("best_subset"),
  #  contains("pca"),
  #  contains("lasso"),
  #  contains("gb"),
  #  contains("svm"),
  #  contains("mrp")
  #)

  boot_out <- list(ebma = ebma, classifiers = classifiers, weights = weights)

  # De-register cluster
  multicore(cores = cores, type = "close", cl = cl)

  return(boot_out)

}
