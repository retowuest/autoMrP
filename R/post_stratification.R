#' Apply post-stratification to classifiers.

post_stratification <- function(data, census, L1.x, L2.x, L2.unit, L2.reg, best.subset, pca, lasso, verbose){
  
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
  lasso_p <- NA 
  for (idx.census in 1:nrow(census)){
    lasso_p[idx.census] <- predict(model_l, newdata = data.frame(census[idx.census,]), type = "response")
  }
  
  # post-stratification: 1) predict, 2) weighted mean by state
  L2_preds <- census %>%
    dplyr::mutate(best_subset = predict(object = model_bs, newdata = census, allow.new.levels = TRUE, type = "response"),
                  pca = predict(object = model_pca, newdata = census, allow.new.levels = TRUE, type = "response"),
                  lasso = lasso_p) %>%
    dplyr::group_by(.dots = list(L2.unit)) %>%
    dplyr::summarize(best_subset = weighted.mean(x = best_subset, w = prop),
                     pca = weighted.mean(x = pca, w = prop),
                     lasso = weighted.mean(x = lasso, w = prop))
  
  # individual predictions for EBMA
  L1_preds <- dplyr::tibble(
    y = data$y,
    best_subest = predict(object = model_bs, type = "response"),
    pca = predict(object = model_pca, type = "response"),
    lasso = as.numeric(predict(object = model_l, type = "response"))
  ) 
  
  # Function output  
  return(ps = list(
    predictions = list(Level1 = L1_preds, Level2 = L2_preds),
    models = list(best_subset = model_bs,
                  pca = model_pca,
                  lasso = model_l)))
}
