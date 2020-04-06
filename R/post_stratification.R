#' Apply post-stratification to classifiers.

post_stratification <- function(data, census, y, L1.x, L2.x, L2.unit, L2.reg,
                                best.subset, pca, lasso, gb, n.minobsinnode,
                                L2.unit.include, L2.reg.include, svm.out,
                                kernel = "radial", verbose){

  # copy argument b/c it is needed twice but might be reset depending on call
  # see lines: 59-63
  L2_unit <- L2.unit

  # bind together survey sample data
  data <- dplyr::bind_rows(data) %>%
    dplyr::mutate_at(.vars = c(L1.x, L2.unit, L2.reg), .funs = as.factor)

  # factors in census
  census <- census %>%
    dplyr::ungroup() %>%
    dplyr::mutate_at(.vars = c(L1.x, L2.unit, L2.reg), .funs = as.factor)

  ## fit best model from each classifier
  # 1) multilevel model with best subset
  model_bs <- best_subset_classifier(
    model = best.subset,
    data.train = data,
    model.family = binomial(link = "probit"),
    model.optimizer = "bobyqa",
    n.iter = 1000000,
    verbose = verbose)

  # 2) multilevel model with best subset of principal components
  model_pca <- best_subset_classifier(
    model = pca,
    data.train = data,
    model.family = binomial(link = "probit"),
    model.optimizer = "bobyqa",
    n.iter = 1000000,
    verbose = verbose)

  # 3) multilevel model with L1 regularization (lasso)
  # Context-level fixed effects
  L2_fe <- paste(L2.x, collapse = " + ")
  L2_fe_form <- as.formula(paste(y, " ~ ", L2_fe, sep = ""))
  # Individual-level random effects as named list
  L1_re <- setNames(as.list(rep(c(~ 1), times = length(c(L1.x, L2.unit, L2.reg)))),
                    c(L1.x, L2.unit, L2.reg))
  model_l <- lasso_classifier(L2.fix = L2_fe_form,
                              L1.re = L1_re,
                              data.train = data,
                              lambda = lasso,
                              model.family = binomial(link = "probit"),
                              verbose = verbose)
  # temporary: predict census data outcome line by line for Lasso
  #lasso_p <- NA
  #for (idx.census in 1:nrow(census)){
  #  lasso_p[idx.census] <- predict(model_l, newdata = data.frame(census[idx.census,]), type = "response")
  #}
  lasso.census <- census
  lasso.census[,y] <- 1
  lasso_p <- predict(model_l, newdata = data.frame(lasso.census), type = "response")

  # 4) boosting
  # Evaluate inclusion of L2.unit
  if (isTRUE(L2.unit.include == FALSE)) {
    L2.unit <- NULL
  }

  # Evaluate inclusion of L2.reg
  if (isTRUE(L2.reg.include == FALSE)) {
    L2.reg <- NULL
  }
  # Create model formula
  x <- paste(c(L1.x, L2.x, L2.unit, L2.reg), collapse = " + ")
  form <- as.formula(paste(y, " ~ ", x, sep = ""))
  # call boosting classifier
  model_gb <- gb_classifier(
    form = form,
    distribution = "bernoulli",
    data.train = data,
    n.trees = gb$n_trees,
    interaction.depth = gb$interaction_depth,
    n.minobsinnode = n.minobsinnode,
    shrinkage = gb$shrinkage,
    verbose = verbose)

  # 5) SVM

  # Create model formula
  x <- paste(c(L1.x, L2.x, L2.unit, L2.reg), collapse = " + ")
  form <- as.formula(paste(y, " ~ ", x, sep = ""))

  # Prepare data
  data <- dplyr::mutate_at(.tbl = data, .vars = y, as.factor)

  # estimate model with tuned set of parameters
  model_svm <- e1071::svm(
    formula = form,
    data = data,
    type = "C-classification",
    cost = svm.out$cost,
    gamma = svm.out$gamma,
    kernel = "radial",
    scale = FALSE,
    probability = TRUE)
  
  # convert factor DV back to numeric
  data <- dplyr::mutate_at(.tbl = data, .vars = y, as.numeric)
  data[,y] <- data[,y] - 1

  # post-stratification: 1) predict, 2) weighted mean by state
  L2_preds <- census %>%
    dplyr::ungroup() %>%
    dplyr::mutate(best_subset = predict(object = model_bs, newdata = census, allow.new.levels = TRUE, type = "response"),
                  pca = predict(object = model_pca, newdata = census, allow.new.levels = TRUE, type = "response"),
                  lasso = lasso_p,
                  gb = gbm::predict.gbm(object = model_gb, newdata = census, n.trees = gb$n_trees, type = "response"),
                  svm = attr(predict(object = model_svm, newdata = census, probability = TRUE),"probabilities")[,"1"]) %>%
    dplyr::group_by(.dots = list(L2_unit)) %>%
    dplyr::summarize(best_subset = weighted.mean(x = best_subset, w = prop),
                     pca = weighted.mean(x = pca, w = prop),
                     lasso = weighted.mean(x = lasso, w = prop),
                     gb = weighted.mean(x = gb, w = prop),
                     svm = weighted.mean(x = svm, w = prop))
  
  
  # individual predictions for EBMA
  L1_preds <- dplyr::tibble(
    best_subset = predict(object = model_bs, type = "response"),
    pca = predict(object = model_pca, type = "response"),
    lasso = as.numeric(predict(object = model_l, type = "response")),
    gb = gbm::predict.gbm(object = model_gb, n.trees = gb$n_trees, type = "response"),
    svm = attr(predict(object = model_svm, newdata = data, probability = TRUE),"probabilities")[,"1"]
  )
  L1_preds <- L1_preds %>%
    dplyr::mutate(!!y := as.numeric(unlist(data[,y]))) %>%
    dplyr::select(one_of(y), "best_subset", "pca", "lasso", "gb", "svm")

  # Function output
  return(ps = list(
    predictions = list(Level1 = L1_preds, Level2 = L2_preds),
    models = list(best_subset = model_bs,
                  pca = model_pca,
                  lasso = model_l,
                  gb = model_gb,
                  svm = model_svm)))
}
