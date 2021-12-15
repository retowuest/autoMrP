#' Apply post-stratification to classifiers.
#'
#' @inheritParams auto_MrP
#' @param best.subset.opt Optimal tuning parameters from best subset selection
#'   classifier. A list returned by \code{run_best_subset()}.
#' @param lasso.opt Optimal tuning parameters from lasso classifier A list
#'   returned by \code{run_lasso()}.
#' @param pca.opt Optimal tuning parameters from best subset selection with
#'   principal components classifier A list returned by \code{run_pca()}.
#' @param gb.opt Optimal tuning parameters from gradient tree boosting
#'   classifier A list returned by \code{run_gb()}.
#' @param svm.opt Optimal tuning parameters from support vector machine
#'   classifier A list returned by \code{run_svm()}.
#' @param mrp.include Whether to run MRP classifier. A logical argument
#'   indicating whether the standard MRP classifier should be used for
#'   predicting outcome \code{y}. Passed from \code{autoMrP()} argument
#'   \code{mrp}.
#' @param n.minobsinnode GB minimum number of observations in the terminal
#'   nodes. An integer-valued scalar specifying the minimum number of
#'   observations that each terminal node of the trees must contain. Passed from
#'   \code{autoMrP()} argument \code{gb.n.minobsinnode}.
#' @param L2.unit.include GB L2.unit. A logical argument indicating whether
#'   \code{L2.unit} should be included in the GB classifier. Passed from
#'   \code{autoMrP()} argument \code{gb.L2.unit}.
#' @param L2.reg.include A logical argument indicating whether \code{L2.reg}
#'   should be included in the GB classifier. Passed from \code{autoMrP()}
#'   argument \code{GB L2.reg}.
#' @param kernel SVM kernel. A character-valued scalar specifying the kernel to
#'   be used by SVM. The possible values are \code{linear}, \code{polynomial},
#'   \code{radial}, and \code{sigmoid}. Passed from \code{autoMrP()} argument
#'   \code{svm.kernel}.
#' @param data A data.frame containing the survey data used in classifier
#'   training.
#' @param ebma.fold A data.frame containing the data not used in classifier
#'   training.

post_stratification <- function(y, L1.x, L2.x, L2.unit, L2.reg,
                                best.subset.opt, lasso.opt, lasso.L2.x,
                                pca.opt, gb.opt, svm.opt, svm.L2.reg,
                                svm.L2.unit, svm.L2.x, mrp.include,
                                n.minobsinnode, L2.unit.include,
                                L2.reg.include, kernel, mrp.L2.x,
                                data, ebma.fold, census, verbose){

  # globals
  lasso <- NULL
  pca <- NULL
  gb <- NULL
  svm <- NULL
  mrp <- NULL
  best_subset <- NULL

  # Copy L2.unit b/c it is needed twice but might be reset depending on call
  L2_unit <- L2.unit

  # post-stratification without level 2 variables
  if (all(L2.x == "")) L2.x <- NULL

  # model container for EBMA
  models <- list()

  # data that is used for models that will not enter EBMA
  no_ebma_data <- dplyr::bind_rows(data, ebma.fold)

  # remove missing values
  data <- tidyr::drop_na(data = data, dplyr::all_of(c(y, L1.x, L2.x, L2.unit, L2.reg)))
  no_ebma_data <- tidyr::drop_na(data = no_ebma_data, dplyr::all_of(c(y, L1.x, L2.x, L2.unit, L2.reg)))

  # ----- Fit optimal model and make prediction for individual classifiers -----

  # Classifier 1: Best Subset
  if (!is.null(best.subset.opt)) {

    # Fit optimal model for EBMA
    best_subset_opt_ebma <- best_subset_classifier(
      model = best.subset.opt,
      data.train = data,
      model.family = binomial(link = "probit"),
      model.optimizer = "bobyqa",
      n.iter = 1000000,
      verbose = verbose)

    # Fit optimal model for post-stratification w/o EBMA
    best_subset_opt_poststrat_only <- best_subset_classifier(
      model = best.subset.opt,
      data.train = no_ebma_data,
      model.family = binomial(link = "probit"),
      model.optimizer = "bobyqa",
      n.iter = 1000000,
      verbose = verbose)

    # post-stratification
    bs_preds <- census %>%
      dplyr::mutate(best_subset = stats::predict(object = best_subset_opt_poststrat_only,
                                                 newdata = ., allow.new.levels = TRUE,
                                                 type = "response")) %>%
      dplyr::group_by( !! rlang::sym(L2.unit)) %>%
      dplyr::summarize(best_subset = stats::weighted.mean(x = best_subset, w = prop), .groups = "keep")

    # individual level predictions for EBMA
    bs_ind <- stats::predict(object = best_subset_opt_ebma, type = "response")

    # model for EBMA
    models$best_subset <- best_subset_opt_ebma

  }

  # Classifier 2: Lasso
  if (!is.null(lasso.opt)) {

    # Determine context-level covariates
    if (is.null(lasso.L2.x)) {
      lasso.L2.x <- L2.x
      }

    # Context-level fixed effects
    L2_fe <- paste(lasso.L2.x, collapse = " + ")
    L2_fe_form <- as.formula(paste(y, " ~ ", L2_fe, sep = ""))

    # Individual-level random effects as named list
    L1_re <- setNames(as.list(rep(c(~ 1), times = length(c(L1.x, L2.unit, L2.reg)))),
                      c(L1.x, L2.unit, L2.reg))

    # Fit optimal model for EBMA
    lasso_opt_ebma <- lasso_classifier(
      L2.fix = L2_fe_form,
      L1.re = L1_re,
      data.train = data,
      lambda = lasso.opt,
      model.family = binomial(link = "probit"),
      verbose = verbose)

    # Fit optimal model for post-stratification w/o EBMA
    lasso_opt_poststrat_only <- lasso_classifier(
      L2.fix = L2_fe_form,
      L1.re = L1_re,
      data.train = no_ebma_data,
      lambda = lasso.opt,
      model.family = binomial(link = "probit"),
      verbose = verbose)

    # predictions for post-stratification only (no EBMA)
    lasso_ests <- predict_glmmLasso(
      m = lasso_opt_poststrat_only,
      lasso.L2.x = lasso.L2.x,
      L2.unit = L2.unit,
      L2.reg = L2.reg,
      L1.x = L1.x,
      census = census)

    # post-stratification
    lasso_preds <- census %>%
      dplyr::mutate(lasso = lasso_ests) %>%
      dplyr::group_by( !! rlang::sym(L2.unit)) %>%
      dplyr::summarize(lasso = stats::weighted.mean(x = lasso, w = prop), .groups = "keep")

    # individual level predictions for EBMA
    lasso_ind <- stats::predict(object = lasso_opt_ebma, type = "response")

    # model for EBMA
    models$lasso <- lasso_opt_ebma

  }

  # Classifier 3: PCA
  if (!is.null(pca.opt)) {
    # Fit optimal model for EBMA
    pca_opt_ebma <- best_subset_classifier(
      model = pca.opt,
      data.train = data,
      model.family = binomial(link = "probit"),
      model.optimizer = "bobyqa",
      n.iter = 1000000,
      verbose = verbose)

    # Fit optimal model for for post-stratification w/o EBMA
    pca_opt_poststrat_only <- best_subset_classifier(
      model = pca.opt,
      data.train = no_ebma_data,
      model.family = binomial(link = "probit"),
      model.optimizer = "bobyqa",
      n.iter = 1000000,
      verbose = verbose)

    # post-stratification
    pca_preds <- census %>%
      dplyr::mutate(pca = stats::predict(object = pca_opt_poststrat_only,
                                         newdata = .,
                                         allow.new.levels = TRUE,
                                         type = "response")) %>%
      dplyr::group_by( !! rlang::sym(L2.unit)) %>%
      dplyr::summarize(pca = stats::weighted.mean(x = pca, w = prop), .groups = "keep")

    # individual level predictions for EBMA
    pca_ind <- stats::predict(object = pca_opt_ebma, type = "response")

    # model for EBMA
    models$pca <- pca_opt_ebma

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

    # Fit optimal model for EBMA
    gb_opt_ebma <- gb_classifier(
      form = form_gb,
      distribution = "bernoulli",
      data.train = data,
      n.trees = gb.opt$n_trees,
      interaction.depth = gb.opt$interaction_depth,
      n.minobsinnode = n.minobsinnode,
      shrinkage = gb.opt$shrinkage,
      verbose = verbose)

    # Fit optimal model for for post-stratification w/o EBMA
    gb_opt_poststrat_only <- gb_classifier(
      form = form_gb,
      distribution = "bernoulli",
      data.train = no_ebma_data,
      n.trees = gb.opt$n_trees,
      interaction.depth = gb.opt$interaction_depth,
      n.minobsinnode = n.minobsinnode,
      shrinkage = gb.opt$shrinkage,
      verbose = verbose)

    # post-stratification
    gb_preds <- census %>%
      dplyr::mutate(gb = gbm::predict.gbm(object = gb_opt_poststrat_only,
                                          newdata = .,
                                          n.trees = gb.opt$n_trees,
                                          type = "response")) %>%
      dplyr::group_by( !! rlang::sym(L2.unit)) %>%
      dplyr::summarize(gb = stats::weighted.mean(x = gb, w = prop), .groups = "keep")

    # individual level predictions for EBMA
    gb_ind <- gbm::predict.gbm(object = gb_opt_ebma, n.trees = gb.opt$n_trees, type = "response")

    # model for EBMA
    models$gb <- gb_opt_ebma

  }

  # Classifier 5: SVM
  if (!is.null(svm.opt)) {

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

    # Prepare data
    svm_data <- data %>%
      dplyr::mutate_at(.vars = y, .funs = list(y_svm = ~as.factor(.))) %>%
      dplyr::select( y_svm, dplyr::all_of(c(y, L1.x, svm.L2.x, L2.unit, L2.reg)) )

    svm_data_no_ebma <- data %>%
      dplyr::bind_rows(ebma.fold) %>%
      dplyr::mutate_at(.vars = y, .funs = list(y_svm = ~as.factor(.))) %>%
      dplyr::select( y_svm, dplyr::all_of(c(y, L1.x, svm.L2.x, L2.unit, L2.reg)) )

    # Create model formula
    x <- paste(c(L1.x, svm.L2.x, svm.L2.unit, svm.L2.reg), collapse = " + ")
    form_svm <- as.formula(paste("y_svm ~ ", x, sep = ""))

    # Fit optimal model for EBMA
    svm_opt_ebma <- svm_classifier(
      form = form_svm,
      data = svm_data,
      kernel = svm.opt$kernel,
      type = "C-classification",
      probability = TRUE,
      svm.gamma = svm.opt$gamma,
      svm.cost = svm.opt$cost,
      verbose = verbose
    )

    # Fit optimal model for post-stratification w/o EBMA
    svm_opt_poststrat_only <- svm_classifier(
      form = form_svm,
      data = svm_data_no_ebma,
      kernel = kernel,
      type = "C-classification",
      probability = TRUE,
      svm.gamma = svm.opt$gamma,
      svm.cost = svm.opt$cost,
      verbose = verbose
    )

    # post-stratification
    svm_preds <- census %>%
      dplyr::filter(dplyr::pull(.data = census, !!L2.unit) %in% dplyr::pull(.data = svm_data_no_ebma, !!L2.unit)) %>%
      dplyr::mutate(svm = attr(stats::predict(object = svm_opt_poststrat_only,
                                              newdata = ., allow.new.levels = TRUE,
                                              probability = TRUE),"probabilities")[,"1"]) %>%
      dplyr::group_by( !! rlang::sym(L2.unit)) %>%
      dplyr::summarize(svm = stats::weighted.mean(x = svm, w = prop), .groups = "keep")

    # individual level predictions for EBMA
    svm_ind <- attr(stats::predict(object = svm_opt_ebma, newdata = svm_data, probability = TRUE),"probabilities")[,"1"]

    # model for EBMA
    models$svm <- svm_opt_ebma

  }

  # Classifier 6: MRP
  # Fit model
  if (isTRUE(mrp.include == TRUE)) {

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
    if(is.null(mrp.L2.x)){
      L2_fe <- paste(L2.x, collapse = " + ")
    } else{
      if(length(mrp.L2.x) == 1){
        if(mrp.L2.x == "empty"){
          L2_fe <- ""
        } else {
          L2_fe <- paste(mrp.L2.x, collapse = " + ")
        }
      } else{
        L2_fe <- paste(mrp.L2.x, collapse = " + ")
      }
    }

    # std MrP formula
    form_mrp <- as.formula(paste(y, " ~ ", L2_fe, " + ", all_re, sep = ""))

    # fit model for EBMA
    mrp_model_ebma <- best_subset_classifier(
      model = form_mrp,
      data.train = data,
      model.family = binomial(link = "probit"),
      model.optimizer = "bobyqa",
      n.iter = 1000000,
      verbose = verbose)
    # mrp_model_ebma <- best_subset_classifier(
    #   model = form_mrp,
    #   data.train = data,
    #   model.family = binomial(link = "probit"),
    #   model.optimizer = "nloptwrap",
    #   n.iter = NULL,
    #   verbose = verbose)

    # fit MrP model for post-stratification only
    mrp_model_poststrat_only <- best_subset_classifier(
      model = form_mrp,
      data.train = no_ebma_data,
      model.family = binomial(link = "probit"),
      model.optimizer = "bobyqa",
      n.iter = 1000000,
      verbose = verbose)
    # mrp_model_poststrat_only <- best_subset_classifier(
    #   model = form_mrp,
    #   data.train = no_ebma_data,
    #   model.family = binomial(link = "probit"),
    #   model.optimizer = "nloptwrap",
    #   n.iter = NULL,
    #   verbose = verbose)

    # post-stratification
    mrp_preds <- census %>%
      dplyr::mutate(mrp = stats::predict(object = mrp_model_poststrat_only,
                                         newdata = ., allow.new.levels = TRUE,
                                         type = "response")) %>%
      dplyr::group_by( !! rlang::sym(L2.unit)) %>%
      dplyr::summarize(mrp = stats::weighted.mean(x = mrp, w = prop), .groups = "keep")

    # individual level predictions for EBMA
    mrp_ind <- stats::predict(object = mrp_model_ebma, type = "response")

    # model for EBMA
    models$mrp <- mrp_model_ebma
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
  if(exists("mrp_preds")) L2_preds <- dplyr::left_join(x = L2_preds, y = mrp_preds, by = L2.unit)


  # individual predictions for EBMA
  L1_preds <- data %>%
    dplyr::select(one_of(y)) %>%
    tidyr::drop_na() %>%
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
      },
      # add MrP
      mrp = if(exists("mrp_ind")){
        as.numeric(mrp_ind)
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
