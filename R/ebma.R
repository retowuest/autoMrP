ebma <- function(ebma.fold, y, L1.x, L2.x, L2.unit, L2.reg, post.strat,
                 n.draws, tol, best.subset.opt, pca.opt, lasso.opt,
                 gb.opt, svm.opt, verbose){

  # run EBMA if at least two classifiers selected
  if (sum(unlist(lapply(X = post.strat$models, FUN = function(x) !is.null(x)))) > 1){

    # models
    model_bs <- post.strat$models$best_subset
    model_pca <- post.strat$models$pca
    model_lasso <- post.strat$models$lasso
    model_gb <- post.strat$models$gb
    model_svm <- post.strat$models$svm
    model_mrp <- post.strat$models$mrp

    # training predictions
    train_preds <- post.strat$predictions$Level1 %>%
      dplyr::select(-one_of(y))

    # training set outcomes
    train_y <- dplyr::select(.data = post.strat$predictions$Level1, one_of(y))

    # container to store the MSE on the test folds
    mse_collector <- matrix(NA, n.draws, length(tol))

    # container for model weights for each draw and tolerance value
    weights_box <- array(NA, dim = c(n.draws, ncol(train_preds), length(tol)))

    # counter for verbose screen output
    counter <- 0

    # loop over tolerance values
    for (idx.tol in 1:length(tol)){

      # loop over Ndraws wit equal obs/state
      for (idx.Ndraws in 1:n.draws){

        # increase counter
        counter <- counter +1

        # sample with replacement equal obs/state
        test <- ebma.fold %>%
          dplyr::group_by(state) %>%
          dplyr::mutate(n_L2 = dplyr::n_groups(ebma.fold)) %>%
          dplyr::sample_n(as.integer(nrow(ebma.fold) /  n_L2), replace = TRUE) %>%
          dplyr::ungroup() %>%
          dplyr::mutate_at(.vars = c( L1.x, L2.unit, L2.reg), .funs = as.factor)

        # predict outcomes in test set
        test_preds <- dplyr::tibble(
          best_subset = if(!is.null(model_bs)){
            predict(object = model_bs, newdata = test, type = "response", allow.new.levels = TRUE)
          } else{NA},
          pca = if(!is.null(model_pca)){
            predict(object = model_pca, newdata = test, type = "response", allow.new.levels = TRUE)
          } else{NA},
          lasso = if(!is.null(model_lasso)){
            as.numeric(predict(object = model_lasso, newdata = data.frame(test), type = "response"))
          } else{NA},
          gb = if(!is.null(model_gb)){
            gbm::predict.gbm(object = model_gb, newdata = test, n.trees = model_gb$n.trees, type = "response")
          } else{NA},
          svm = if(!is.null(model_svm)){
            as.numeric(attr(predict(object = model_svm, newdata = test, probability = TRUE),"probabilities")[,"1"])
          } else{NA},
          mrp = if(!is.null(model_mrp)){
            predict(object = model_mrp, newdata = test, type = "response", allow.new.levels = TRUE)
          } else{NA}
        )
        # remove NA's
        test_preds <- test_preds[,apply(X = test_preds, MARGIN = 2, FUN = function(x){
          all(!is.na(x))})]

        # outcome on the test
        test_y <- dplyr::select(.data = test, one_of(y))

        # EBMA
        forecast.data <- suppressWarnings(
          EBMAforecast::makeForecastData(
            .predCalibration = data.frame(train_preds),
            .outcomeCalibration = as.numeric(unlist(train_y)),
            .predTest = data.frame(test_preds),
            .outcomeTest = as.numeric(unlist(test_y)))
        )

        forecast.out <- EBMAforecast::calibrateEnsemble(
          forecast.data,
          model = "normal",
          useModelParams = FALSE,
          tol = tol[idx.tol])

        # mse
        mse_collector[idx.Ndraws, idx.tol] <- mean(( as.numeric(unlist(test_y)) - as.numeric( attributes(forecast.out)$predTest[,1,1]))^2)

        # model weights
        weights_box[idx.Ndraws, , idx.tol] <- attributes(forecast.out)$modelWeights

        # progress
        if (verbose) cat(paste("\n","EBMA: ", round(counter / (length(tol) * n.draws),2)*100, "% done",sep=""))

      }
    }

    # which tolerance value minimizes the mse on the test set
    best_tolerance <- apply(mse_collector, 1, function(x) which.min(x))

    # container of best model weights
    weights_mat <- matrix(data = NA, nrow = n.draws, ncol = ncol(train_preds))

    # model weights; rows = observations, columns = model weights, layers = tolerance values
    if (length(tol)>1){
      for (idx.tol in 1:length(best_tolerance)){
        weights_mat[idx.tol, ] <- weights_box[idx.tol, ,][ ,best_tolerance[idx.tol]]
      }
    } else{
      weights_mat <- weights_box[, , 1]
    }

    # average model weights
    final_model_weights <- apply(weights_mat, 2, mean)
    names(final_model_weights) <- names(train_preds)

    # weighted average
    w_avg <- as.numeric(as.matrix(post.strat$predictions$Level2[,names(final_model_weights)]) %*% final_model_weights)

    # L2 preds object
    L2_preds <- dplyr::tibble(
      state = post.strat$predictions$Level2$state,
      ebma = w_avg
    )

    # function output
    return(list(ebma = L2_preds, classifiers = post.strat$predictions$Level2))

  } else{
   cat("\n Skipping EBMA (only 1 classifier selected) \n")

    # function output
    return(list(ebma = "EBMA step skipped (only 1 classifier run)", classifiers = post.strat$predictions$Level2))
  }
}
