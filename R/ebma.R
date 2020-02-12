ebma <- function(ebma.fold, L1.x, L2.x, L2.unit, L2.reg, post.strat,
                 Ndraws, tol.values, best.subset, pca, lasso, gb,
                 svm.out, verbose){

  # models
  model_bs <- post.strat$models$best_subset
  model_pca <- post.strat$models$pca
  model_lasso <- post.strat$models$lasso
  model_gb <- post.strat$models$gb
  model_svm <- post.strat$models$svm

  # training predictions
  train_preds <- post.strat$predictions$Level1 %>%
    dplyr::select(-y)

  # training set outcomes
  train_y <- dplyr::select(.data = post.strat$predictions$Level1, y)

  # container to store the MSE on the test folds
  mse_collector <- matrix(NA, Ndraws, length(tol.values))

  # container for model weights for each draw and tolerance value
  weights_box <- array(NA, dim = c(Ndraws, ncol(train_preds), length(tol.values)))

  # counter for verbose screen output
  counter <- 0

  # loop over tolerance values
  for (idx.tol in 1:length(tol.values)){

    # loop over Ndraws wit equal obs/state
    for (idx.Ndraws in 1:Ndraws){

      # increase counter
      counter <- counter +1

      # number of states in the sample
      dplyr::n_groups(x = test)

      # sample with replacement equal obs/state
      test <- ebma.fold %>%
        dplyr::group_by(state) %>%
        dplyr::mutate(n_L2 = dplyr::n_groups(test)) %>%
        dplyr::sample_n(as.integer(nrow(ebma.fold) /  n_L2), replace = TRUE) %>%
        dplyr::ungroup() %>%
        dplyr::mutate_at(.vars = c( L1.x, L2.unit, L2.reg), .funs = as.factor)

      # predict outcomes in test set
      test_preds <- dplyr::tibble(
        best_subset = predict(object = model_bs, newdata = test, type = "response", allow.new.levels = TRUE),
        pca = predict(object = model_pca, newdata = test, type = "response", allow.new.levels = TRUE),
        lasso = as.numeric(predict(object = model_lasso, newdata = data.frame(test), type = "response")),
        gb = gbm::predict.gbm(object = model_gb, newdata = test, n.trees = model_gb$n.trees, type = "response"),
        svm = attr(predict(object = model_svm, newdata = test, probability = TRUE),"probabilities")[,"1"]
      )

      # outcome on the test
      test_y <- dplyr::select(.data = test, one_of(y))

      # EBMA
      forecast.data <- EBMAforecast::makeForecastData(
        .predCalibration = data.frame(train_preds),
        .outcomeCalibration = as.numeric(unlist(train_y)),
        .predTest = data.frame(test_preds),
        .outcomeTest = as.numeric(unlist(test_y)))

      forecast.out <- EBMAforecast::calibrateEnsemble(
        forecast.data,
        model = "normal",
        useModelParams = FALSE,
        tol = tol.values[idx.tol])

      # mse
      mse_collector[idx.Ndraws, idx.tol] <- mean(( as.numeric(unlist(test_y)) - as.numeric( attributes(forecast.out)$predTest[,1,1]))^2)

      # model weights
      weights_box[idx.Ndraws, , idx.tol] <- attributes(forecast.out)$modelWeights

      # progress
      if (verbose) cat(paste("\n","EBMA: ", round(counter / (length(tol.values) * Ndraws),2)*100, "% done",sep=""))

    }
  }

  # which tolerance value minimizes the mse on the test set
  best_tolerance <- apply(mse_collector, 1, function(x) which.min(x))

  # container of best model weights
  weights_mat <- matrix(data = NA, nrow = Ndraws, ncol = ncol(train_preds))

  # model weights; rows = observations, columns = model weights, layers = tolerance values
  for (idx.tol in 1:length(best_tolerance)){
    weights_mat[idx.tol, ] <- weights_box[idx.tol, ,][ ,best_tolerance[idx.tol]]
  }

  # average model weights
  final_model_weights <- apply(weights_mat, 2, mean)
  names(final_model_weights) <- c("BestSubset", "PCA", "Lasso", "GB", "SVM")

  # weighted average
  L2_preds <- post.strat$predictions$Level2 %>%
    dplyr::mutate(ebma = as.numeric(cbind(best_subset, pca, lasso, gb, svm) %*% final_model_weights))

  # function output
  return(L2_preds)
}
