#' Apply post-stratification to classifiers.

post_stratification <- function(y, L1.x, L2.x, L2.unit, L2.reg,
                                best.subset.opt, lasso.opt,
                                pca.opt, gb.opt, svm.opt,
                                mrp.include, n.minobsinnode,
                                L2.unit.include, L2.reg.include,
                                kernel, mrp.L2.x, data, census,
                                verbose){

  #browser()

  # Copy L2.unit b/c it is needed twice but might be reset depending on call
  # See lines: 59-63
  L2_unit <- L2.unit

  # model container for EBMA
  models <- list()

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

    # post-stratification
    bs_preds <- census %>%
      dplyr::mutate(best_subset = stats::predict(object = best_subset_opt, newdata = census, allow.new.levels = TRUE, type = "response")) %>%
      dplyr::group_by(.dots = list(L2_unit)) %>%
      dplyr::summarize(best_subset = weighted.mean(x = best_subset, w = prop))

    # individual level predictions for EBMA
    bs_ind <- stats::predict(object = best_subset_opt, type = "response")

    # model for EBMA
    models$best_subset <- best_subset_opt

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

    lasso.census <- census
    lasso.census[,y] <- 1
    lasso_p <- stats::predict(lasso_opt, newdata = data.frame(lasso.census), type = "response")

    # # post-stratification
    lasso_preds <- census %>%
      dplyr::mutate(lasso = lasso_p) %>%
      dplyr::group_by(.dots = list(L2_unit)) %>%
      dplyr::summarize(lasso = weighted.mean(x = lasso, w = prop))

    # individual level predictions for EBMA
    lasso_ind <- stats::predict(object = lasso_opt, type = "response")

    # model for EBMA
    models$lasso <- lasso_opt

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

    # post-stratification
    pca_preds <- census %>%
      dplyr::mutate(pca = stats::predict(object = pca_opt, newdata = census, allow.new.levels = TRUE, type = "response")) %>%
      dplyr::group_by(.dots = list(L2_unit)) %>%
      dplyr::summarize(pca = weighted.mean(x = pca, w = prop))

    # individual level predictions for EBMA
    pca_ind <- stats::predict(object = pca_opt, type = "response")

    # model for EBMA
    models$pca <- pca_opt


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

    # post-stratification
    gb_preds <- census %>%
      dplyr::mutate(gb = gbm::predict.gbm(object = gb_opt, newdata = census, n.trees = gb.opt$n_trees, type = "response")) %>%
      dplyr::group_by(.dots = list(L2_unit)) %>%
      dplyr::summarize(gb = weighted.mean(x = gb, w = prop))

    # individual level predictions for EBMA
    gb_ind <- gbm::predict.gbm(object = gb_opt, n.trees = gb.opt$n_trees, type = "response")

    # model for EBMA
    models$gb <- gb_opt

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

    # post-stratification
    svm_preds <- census %>%
      dplyr::mutate(svm = attr(stats::predict(object = svm_opt, newdata = census, probability = TRUE),"probabilities")[,"1"]) %>%
      dplyr::group_by(.dots = list(L2_unit)) %>%
      dplyr::summarize(svm = weighted.mean(x = svm, w = prop))

    # individual level predictions for EBMA
    attr(stats::predict(object = svm_opt, newdata = data, probability = TRUE),"probabilities")[,"1"]

    # individual level predictions for EBMA
    svm_ind <- attr(stats::predict(object = svm_opt, newdata = data, probability = TRUE),"probabilities")[,"1"]

    # model for EBMA
    models$svm <- svm_opt

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
  L2_fe <- paste(mrp.L2.x, collapse = " + ")

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

  # --------------------------- combine l2 level predictions ----------------------------

  # tibble of L2 units
  L2_preds <- dplyr::select(census, one_of(L2.unit)) %>%
    dplyr::distinct()
  # add existing classifiers
  if(exists("bs_preds")) L2_preds <- dplyr::left_join(x = L2_preds, y = bs_preds, by = L2.unit)
  if(exists("lasso_preds")) L2_preds <- dplyr::left_join(x = L2_preds, y = lasso_preds, by = L2.unit)
  if(exists("pca_preds")) L2_preds <- dplyr::left_join(x = L2_preds, y = pca_preds, by = L2.unit)
  if(exists("gb_preds")) L2_preds <- dplyr::left_join(x = L2_preds, y = gb_preds, by = L2.unit)
  if(exists("svm_preds")) L2_preds <- dplyr::left_join(x = L2_preds, y = svm_preds, by = L2.unit)


  # individual predictions for EBMA
  L1_preds <- data %>%
    dplyr::select(one_of(y)) %>%
    # add best subset
    dplyr::mutate(
      best_subset = if(exists("bs_ind")){
        as.numeric(bs_ind)
      } else{
        NA
        },
      # add lasso
      lasso = if(exists("lasso_ind")){
        as.numeric(lasso_ind)
        } else{
          NA
          },
      # add pca
      pca = if(exists("pca_ind")){
        as.numeric(pca_ind)
      } else{
        NA
      },
      # add gb
      gb = if(exists("gb_ind")){
        as.numeric(gb_ind)
      } else{
        NA
      },
      # add svm
      svm = if(exists("svm_ind")){
        as.numeric(svm_ind)
      } else{
        NA
      }
    )
  # remove NA's
  L1_preds <- L1_preds[,apply(X = L1_preds, MARGIN = 2, FUN = function(x){
    all(!is.na(x))})]

  # Function output
  return(ps = list(predictions = list(Level1 = L1_preds, Level2 = L2_preds),
                   models = models))
}
