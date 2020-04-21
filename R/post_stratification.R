#' Apply post-stratification to classifiers.

post_stratification <- function(y, L1.x, L2.x, L2.unit, L2.reg,
                                best.subset.opt, lasso.opt,
                                pca.opt, gb.opt, svm.opt,
                                mrp.include, n.minobsinnode,
                                L2.unit.include, L2.reg.include,
                                kernel, L2.x.mrp, data, census,
                                verbose){

  # Copy L2.unit b/c it is needed twice but might be reset depending on call
  # See lines: 59-63
  L2_unit <- L2.unit

  # ----- Fit optimal model and make prediction for individual classifiers -----

  # Classifier 1: Best Subset
  if (!is.null(best.subset.opt)) {
    # Fit optimal model
    best_subset_opt <- best_subset_classifier(model = best.subset.opt,
                                              data.train = data,
                                              model.family = binomial(link = "probit"),
                                              model.optimizer = "bobyqa",
                                              n.iter = 1000000,
                                              verbose = verbose)
  }

  # Classifier 2: Lasso
  if (!is.null(lasso.opt)) {
    # Context-level fixed effects
    L2_fe <- paste(L2.x, collapse = " + ")
    L2_fe_form <- as.formula(paste(y, " ~ ", L2_fe, sep = ""))

    # Individual-level random effects as named list
    L1_re <- setNames(as.list(rep(c(~ 1), times = length(c(L1.x, L2.unit, L2.reg)))),
                      c(L1.x, L2.unit, L2.reg))

    # Fit optimal model
    lasso_opt <- lasso_classifier(L2.fix = L2_fe_form,
                                  L1.re = L1_re,
                                  data.train = data,
                                  lambda = lasso.opt,
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
  }

  # Classifier 3: PCA
  if (!is.null(pca.opt)) {
    # Fit optimal model
    pca_opt <- best_subset_classifier(model = pca.opt,
                                      data.train = data,
                                      model.family = binomial(link = "probit"),
                                      model.optimizer = "bobyqa",
                                      n.iter = 1000000,
                                      verbose = verbose)
  }

  # Classifier 4: GB
  if (!is.null(gb.opt)) {
    # Evaluate inclusion of L2.unit
    if (isTRUE(L2.unit.include == TRUE)) {
      L2.unit.gb <- L2.unit
    } else {
      L2.unit.gb <- NULL
    }

    # Evaluate inclusion of L2.reg
    if (isTRUE(L2.reg.include == TRUE)) {
      L2.reg.gb <- L2.reg
    } else {
      L2.reg.gb <- NULL
    }

    # Create model formula
    x <- paste(c(L1.x, L2.x, L2.unit.gb, L2.reg.gb), collapse = " + ")
    form_gb <- as.formula(paste(y, " ~ ", x, sep = ""))

    # Fit optimal model
    gb_opt <- gb_classifier(form = form_gb,
                            distribution = "bernoulli",
                            data.train = data,
                            n.trees = gb.opt$n_trees,
                            interaction.depth = gb.opt$interaction_depth,
                            n.minobsinnode = n.minobsinnode,
                            shrinkage = gb.opt$shrinkage,
                            verbose = verbose)
  }

  # Classifier 5: SVM
  if (!is.null(svm.opt)) {
    # Prepare data
    data <- data %>%
      dplyr::mutate_at(.vars = y, .funs = list(y_svm = ~as.factor(.)))

    # Create model formula
    x <- paste(c(L1.x, L2.x, L2.unit, L2.reg), collapse = " + ")
    form_svm <- as.formula(paste("y_svm ~ ", x, sep = ""))

    # Fit optimal model
    if (isTRUE(verbose == TRUE)) {
      svm_opt <- e1071::svm(formula = form_svm,
                            data = data,
                            scale = FALSE,
                            type = "C-classification",
                            kernel = kernel,
                            gamma = svm.opt$gamma,
                            cost = svm.opt$cost,
                            probability = TRUE)
    } else {
      svm_opt <- suppressMessages(suppressWarnings(
        e1071::svm(formula = form_svm,
                   data = data,
                   scale = FALSE,
                   type = "C-classification",
                   kernel = kernel,
                   gamma = svm.opt$gamma,
                   cost = svm.opt$cost,
                   probability = TRUE)
      ))
    }
  }

  # Classifier 6: MRP
  # Create model formula
  # Individual-level random effects
  L1_re <- paste(paste("(1 | ", L1.x, ")", sep = ""), collapse = " + ")

  # Geographic unit or geographic unit-geographic region random effects
  if (is.null(L2.reg)) {
    L2_re <- paste("(1 | ", L2.unit, ")", sep = "")
  } else {
    L2_re <- paste(paste("(1 | ", L2.reg, "/", L2.unit, ")", sep = ""),
                   collapse = " + ")
  }

  # Combine all random effects
  all_re <- paste(c(L1_re, L2_re), collapse = " + ")

  # Context-level fixed effects
  L2_fe <- paste(L2.x.mrp, collapse = " + ")

  # Empty model
  form_mrp <- as.formula(paste(y, " ~ ", L2_fe, " + ", all_re, sep = ""))

  # Fit optimal model
  if (isTRUE(mrp.include == TRUE)) {
    if (isTRUE(verbose == TRUE)) {
      mrp_opt <- lme4::glmer(formula = form_mrp,
                             data = data,
                             family = binomial(link = "probit"),
                             lme4::glmerControl(optimizer = "bobyqa",
                                                optCtrl = list(maxfun = 1000000)))
    } else {
      mrp_opt <- suppressMessages(suppressWarnings(
        lme4::glmer(formula = form_mrp,
                    data = data.train,
                    family = binomial(link = "probit"),
                    lme4::glmerControl(optimizer = "bobyqa",
                                       optCtrl = list(maxfun = 1000000)))
      ))
    }
  }

  # --------------------------- Post-stratification ----------------------------

  # post-stratification: 1) predict, 2) weighted mean by state
  L2_preds <- census %>%
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
  return(ps = list(predictions = list(Level1 = L1_preds, Level2 = L2_preds),
                   models = list(best_subset = model_bs,
                                 pca = model_pca,
                                 lasso = model_l,
                                 gb = model_gb,
                                 svm = model_svm)))
}
