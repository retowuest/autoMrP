#' Improve MrP through ensemble learning.
#'
#' \code{auto_MrP} improves the prediction performance of multilevel regression
#' with post-stratification (MrP) by combining a number of machine learning
#' methods through ensemble bayesian model averaging (EBMA).
#'
#' @param y Outcome variable. A character scalar containing the column name of
#'   the outcome variable.
#' @param L1.x Individual-level covariates. A character vector of column names
#'   corresponding to the individual-level variables used to predict the outcome
#'   variable. Note that geographic unit is specified in argument `L2.unit`.
#' @param L2.x Context-level covariates. A character vector of column names
#'   corresponding to the context-level variables used to predict the outcome
#'   variable.
#' @param L2.unit Geographic unit. A character scalar indicating the column
#'   name of the geographic unit at which outcomes should be aggregated.
#' @param L2.reg Geographic region. A character scalar indicating the column
#'   name of the geographic region by which geographic units are grouped
#'   (`L2.unit` must be nested within `L2.reg`). Default is NULL.
#' @param survey Survey data. A data.frame whose column names include `y`,
#'    `L1.x`, `L2.x`, `L2.unit`, and, if specified, `L2.reg`.
#' @param census Census data. A data.frame whose column names include `L1.x`,
#'   `L2.x`, `L2.unit`, and, if specified, `L2.reg`.
#' @param L2.x.scale Scale context-level covariates. A logical argument
#'   specifying whether context-level covariates should be normalized. Default
#'   is TRUE. Note that if L2.x.scale is FALSE, then context-level covariates
#'   should be normalized prior to calling auto_MrP.
#' @param bin.proportion Proportion of ideal types. A character scalar
#'   indicating the column name of the variable in census containing the
#'   propotion of individuals of an ideal type in a geographic unit. Default is
#'   NULL. Note: Not needed if `bin.size` is provided.
#' @param bin.size Bin size for ideal types. A character scalar indicating the
#'   column name of the variable in census containing the bin size for ideal
#'   types in a geographic unit. Default is NULL. Note: Not needed if
#'   `bin.proportion` is provided.
#' @param uncertainty Uncertainty estimates. A logical argument indicating
#'   whether uncertainty estimates should be computed. Default is FALSE.
#' @param best.subset Best subset classifier. A logical argument indicating
#'   whether best subset classifier should be used for prediction of the
#'   outcome. Default is TRUE.
#' @param lasso Lasso classifier. A logical argument indicating whether lasso
#'   classifier should be used for prediction of the outcome. Default is TRUE.
#' @param pca pca classifier. A logical argument indicating whether pca
#'   classifier should be used for prediction of the outcome. Default is TRUE.
#' @param gb gb classifier. A logical argument indicating whether gb classifier
#'   should be used for prediction of the outcome. Default is TRUE.
#' @param svm svm classifier. A logical argument indicating whether svm
#'   classifier should be used for prediction of the outcome. Default is TRUE.
#' @param mrp mrp classifier. A logical argument indicating whether standard
#'   mrp classifier should be used for prediction of the outcome. Default is TRUE.
#' @param forward.selection Apply forward selection for mulilevel model with
#'   post-stratification classifier instead of best subset selection. A logical
#'   argument indicating whether to apply forward selection instead of best subset
#'   selection or not. Default is FALSE. Note: With more than 8 context level
#'   variables, forward selection is recommended.
#' @param ebma.size Size of EBMA hold-out fold. A rational number in the open
#'   unit interval indicating the share of respondents to be contained in the
#'   EBMA hold-out fold. If left unspecified (NULL), then ebma.size is set to
#'   1/3 of the survey sample size. Default is NULL.
#' @param k.folds Number of folds. An integer-valued scalar indicating the
#'   number of folds to be used for cross-validation. Defaults to the value of 5.
#' @param custom.folds Custom folds. A character scalar indicating the column
#'   name of the variable in survey data that specifies the fold to which an
#'   observation should be allocated. The variable should contain integers
#'   running from 1 to k + 1, where k is the number of folds used in
#'   cross-validation. Value k + 1 refers to the ebma fold. Default is NULL.
#' @param custom.pc Custom pcs. A character vector of column names corresponding
#'   to the principal components of the context-level variables in `survey` and
#'   `census`. Note that the columns containing the principal components in
#'   `survey` and `census` must be named identically. Default is NULL.
#' @param cv.sampling Sampling method. A character-valued scalar indicating
#'   whether sampling in the creation of cross-validation folds should be done
#'   by respondents or geographic units. Default is by geographic units.
#' @param loss.unit Loss function unit. A character-valued scalar indicating
#'   whether the loss should be evaluated at the level of individual respondents
#'   or the level of geographic units. Default is at the individual level.
#' @param loss.measure Loss function measure. A character-valued scalar
#'   indicating whether the loss should be measured by the mean squared error
#'   or the mean absolute error. Default is the MSE.
#' @param lasso.lambda.set Set of tuning parameters. Lambda is the penalty
#'   parameter that controls the shrinkage of fixed effects. Either a numeric
#'   vector of lambda values or a data.frame with two columns, the first
#'   containing the size by which lambda should increase and the second the
#'   upper threshold of the interval of lambdas to which the step size applies.
#' @param lasso.iterations.max Stopping rule. A numeric scalar specifying the
#'   maximum number of iterations without performance improvement the
#'   algorithm runs before stopping. Default is NULL.
#' @param gb.L2.unit.include Include L2.unit in GB. A logical argument indicating
#'   whether L2.unit is included in the GB models. Default is FALSE.
#' @param gb.L2.reg.include Include L2.reg in GB. A logical argument indicating
#'   whether L2.reg is included in the GB models. Default is FALSE.
#' @param gb.interaction.set Set of interaction depth values. An integer-valued
#'   vector whose values define the maximum depth of each tree. Interaction
#'   depth is used to tune the model. Default is Default is c(1, 2, 3).
#' @param gb.shrinkage.set Learning rate. A numeric vector whose values define
#'   the learning rate or step-size reduction. Learning rate is used to tune
#'   the model. Values between 0.001 and 0.1 usually work, but a smaller
#'   learning rate typically requires more trees. Default is
#'   c(0.04, 0.01, 0.008, 0.005, 0.001).
#' @param gb.tree.start Initial total number of trees. An integer-valued scalar
#'   specifying the initial number of total trees. Default is 50.
#' @param gb.tree.increase.set Increase in total number of trees. Either an
#'   integer-valued scalar specifying by how many trees the total number of
#'   trees is increased (until the maximum number of trees is reached) or an
#'   integer-valued vector of `length(gb.shrinkage.set)` with each value being
#'   associated with a learning rate. Total number of trees is used to tune the
#'   model. Default is 50.
#' @param gb.trees.max.set Maximum number of trees. Either an integer-valued
#'   scalar specifying the maximum number of trees or an integer-valued vector
#'   of `length(gb.shrinkage.set)` with each value being associated with a
#'   learning rate and a number of tree increase. Default is 1000.
#' @param gb.iterations.max Stopping rule. A numeric scalar specifying the
#'   maximum number of iterations without performance improvement the GB
#'   classifier runs before stopping. Default is 70.
#' @param gb.n.minobsinnode Minimum number of observations in the terminal nodes.
#'   An integer-valued scalar specifying the minimum number of observations
#'   that each terminal node of the trees must contain. Default is 5.
#' @param svm.kernel Kernel for SVM. A character scalar specifying the kernel to
#'   be used for SVM. The possible types are linear, polynomial, radial, and
#'   sigmoid. Default is radial.
#' @param svm.error.fun. Loss function. Defaults to misclassification error for
#'   factor dependent variables and MSE for numeric dependent variables.
#' @param svm.gamma.set Gamma parameter for SVM. This parameter is needed for
#'   all kernels except linear. Default is
#'   c(0.3, 0.5, 0.55, 0.6, 0.65, 0.7, 0.8, 0.9, 1, 2, 3, 4).
#' @param svm.cost.set Cost parameter for SVM. This parameter specifies the cost
#'   of constraints violation. Default is c(1, 10).
#' @param ebma.n.draws Number of ebma samples. An integer-valued scalar
#'   specifying the number of bootstrapped samples to be drawn from the ebma
#'   fold and used for tuning ebma. Default is 100.
#' @param ebma.tol.values Tolerance for ebma. A numeric vector containing the
#'   tolerance values for improvements in the log-likelihood before the em
#'   algorithm will stop optimization. Values should range at least from 0.01
#'   to 0.001. Default is c(0.01, 0.005, 0.001, 0.0005, 0.0001, 0.00005, 0.00001).
#' @param seed Seed. An integer-valued scalar to control random number
#'   generation. If left unspecified (NULL), then seed is set to 546213978.
#' @param verbose Verbose output. A logical argument indicating whether or not
#'   verbose output should be printed. Default is TRUE.
#' @return
#' @keywords MRP multilevel regression post-stratification machine learning
#'   EBMA ensemble bayesian model averaging
#' @examples
#' @export

auto_MrP <- function(y, L1.x, L2.x,
                     L2.unit, L2.reg = NULL,
                     survey, census,
                     L2.x.scale = TRUE,
                     bin.proportion = NULL,
                     bin.size = NULL,
                     uncertainty = FALSE,
                     best.subset = TRUE,
                     lasso = TRUE,
                     pca = TRUE,
                     gb = TRUE,
                     svm = TRUE,
                     mrp = FALSE,
                     forward.selection = FALSE,
                     ebma.size = NULL,
                     k.folds = 5,
                     custom.folds = NULL,
                     custom.pc = NULL,
                     cv.sampling = "L2 units",
                     loss.unit = "individual",
                     loss.measure = "mse",
                     lasso.lambda.set = data.frame(step_size = c(0.1, 0.3, 1),
                                                 threshold = c(1, 10, 10000)),
                     lasso.iterations.max = 70,
                     gb.L2.unit.include = FALSE,
                     gb.L2.reg.include = FALSE,
                     gb.interaction.set = c(1, 2, 3),
                     gb.shrinkage.set = c(0.04, 0.01, 0.008, 0.005, 0.001),
                     gb.tree.start = 50,
                     gb.tree.increase.set = 50,
                     gb.trees.max.set = 1000,
                     gb.iterations.max = 70,
                     gb.n.minobsinnode = 5,
                     svm.kernel = "radial",
                     svm.error.fun = "MSE",
                     svm.gamma.set = c(0.3, 0.5, 0.55, 0.6, 0.65, 0.7, 0.8, 0.9, 1, 2, 3, 4),
                     svm.cost.set = c(1, 10),
                     ebma.n.draws = 100,
                     ebma.tol.values = c(0.01, 0.005, 0.001,
                                         0.0005, 0.0001, 0.00005, 0.00001),
                     seed = NULL,
                     verbose = TRUE) {

  # Set seed
  if (is.null(seed)) {
    set.seed(546213978)
  } else {
    set.seed(seed)
  }

  # Error and warning checks
  if (!all(L1.x %in% colnames(survey))) {
    stop(paste("Individual-level variable(s) '",
               L1.x[which(!(L1.x %in% colnames(survey)))],
               "' is/are not in your survey data.", sep = ""))
  }

  if (!all(L1.x %in% colnames(census))) {
    stop(paste("Individual-level variable(s) '",
               L1.x[which(!(L1.x %in% colnames(census)))],
               "' is/are not in your census data.", sep = ""))
  }

  if (!all(L2.x %in% colnames(survey))) {
    stop(paste("Context-level variable(s) '",
               L2.x[which(!(L2.x %in% colnames(survey)))],
               "' is/are not in your survey data.", sep = ""))
  }

  if (!all(L2.x %in% colnames(census))) {
    stop(paste("Context-level variable(s) '",
               L2.x[which(!(L2.x %in% colnames(census)))],
               "' is/are not in your census data.", sep = ""))
  }

  if (!(y %in% colnames(survey))) {
    stop(paste("Outcome '", y,
               "' is not in your survey data.", sep = ""))
  }

  if (!(L2.unit %in% colnames(survey))) {
    stop(paste("The geographic unit '", L2.unit,
               "' is not in your survey data.", sep = ""))
  }

  if (!(L2.unit %in% colnames(census))) {
    stop(paste("The geographic unit '", L2.unit,
               "' is not in your census data.", sep = ""))
  }

  if (!is.null(L2.reg)) {
    if (!(L2.reg %in% colnames(survey))) {
      stop(paste("The geographic region '", L2.reg,
                 "' is not in your survey data.", sep = ""))
    }

    if (!(L2.reg %in% colnames(census))) {
      stop(paste("The geographic region '", L2.reg,
                 "' is not in your census data.", sep = ""))
    }

    if (any(unlist(lapply(dplyr::group_split(survey, .data[[L2.unit]]),
                          function(x) length(unique(x[[L2.reg]])))) > 1)) {
      stop(paste("The geographic unit(s) '",
                 which(unlist(lapply(dplyr::group_split(survey, .data[[L2.unit]]),
                                     function(x) length(unique(x[[L2.reg]])))) > 1),
                 "' is/are nested in multiple regions in your survey data."))
    }

    if (any(unlist(lapply(dplyr::group_split(census, .data[[L2.unit]]),
                          function(x) length(unique(x[[L2.reg]])))) > 1)) {
      stop(paste("The geographic unit(s) '",
                 which(unlist(lapply(dplyr::group_split(census, .data[[L2.unit]]),
                                     function(x) length(unique(x[[L2.reg]])))) > 1),
                 "' is/are nested in multiple regions in your census data."))
    }
  }

  if (is.null(ebma.size)) {
    ebma.size <- round(nrow(survey) / 3, digits = 0)
  } else if (is.numeric(ebma.size) & ebma.size > 0 & ebma.size < 1) {
    ebma.size <- round(nrow(survey) * ebma.size, digits = 0)
  } else {
    stop("ebma.size must be a rational number in the open unit interval.")
  }

  if (!((is.integer(k.folds) | all(as.integer(k.folds) == k.folds)) &
        length(k.folds) == 1)) {
    stop("k.folds must be an integer number.")
  }

  if (!cv.sampling %in% c("respondents", "L2 units")) {
    stop("cv.sampling must take either the value 'respondents' or 'L2 units'.")
  }

  if (!(is.vector(lasso.lambda.set) | is.data.frame(lasso.lambda.set))) {
    stop(paste("lasso.lambda.set must be either a numeric vector or a data.frame ",
               "with two columns, one for step size increase and the other ",
               "for the upper threshold of the interval of lambdas to which ",
               "the step size applies", sep = ""))
  }

  if (!(is.null(lasso.iterations.max) | (is.numeric(lasso.iterations.max) &
                                         length(lasso.iterations.max) == 1))) {
    stop("lasso.iterations.max must be either a numeric scalar or NULL.")
  }

  if (!(is.integer(gb.interaction.set) |
        all(as.integer(gb.interaction.set) == gb.interaction.set))) {
    stop("gb.interaction.set must be an integer-valued vector.")
  }

  if (!is.numeric(gb.shrinkage.set)) {
    stop("gb.shrinkage.set must be a numeric vector")
  } else if (min(gb.shrinkage.set) < 0.001 | max(gb.shrinkage.set) > 0.1) {
    warning("gb.shrinkage.set should have values lying between 0.001 and 0.1.")
  }

  if (!((is.integer(gb.tree.start) |
         all(as.integer(gb.tree.start) == gb.tree.start)) &
        length(gb.tree.start) == 1)) {
    stop("gb.tree.start must be an integer-valued scalar.")
  }

  if (!(is.integer(gb.tree.increase.set) |
        all(as.integer(gb.tree.increase.set) == gb.tree.increase.set))) {
    stop("gb.tree.increase.set must be an integer-valued scalar or vector.")
  } else if (length(gb.tree.increase.set) > 1 &
             length(gb.tree.increase.set) != length(gb.shrinkage.set)) {
    stop(paste("gb.tree.increase.set must be either a scalar or a vector of ",
               "size `length(gb.shrinkage.set)`.", sep = ""))
  }

  if (!(is.integer(gb.trees.max.set) |
        all(as.integer(gb.trees.max.set) == gb.trees.max.set))) {
    stop("gb.trees.max.set must be an integer-valued scalar or vector.")
  } else if (length(gb.trees.max.set) > 1 &
             length(gb.trees.max.set) != length(gb.shrinkage.set)) {
    stop(paste("gb.trees.max.set must be either a scalar or a vector of size ",
               "`length(gb.shrinkage.set)`.", sep = ""))
  }

  if (!is.null(custom.folds)) {
    if (!(is.double(survey %>% dplyr::select_at(.vars = custom.folds) %>% dplyr::pull()) |
          is.integer(survey %>% dplyr::select_at(.vars = custom.folds) %>% dplyr::pull()))) {
      stop("The variable specifying folds must of type integer or double.")
    }

    if (!all(survey %>% dplyr::select_at(.vars = custom.folds) %>%
             dplyr::pull() %>% unique() %>% sort() == 1:(k.folds + 1))) {
      stop("custom.folds must contain integers ranging from 1 to `k.folds` + 1.")
    }
  }

  # ------------------------------- Prepare data -------------------------------

  # If not provided in census data, calculate bin size and bin proportion for
  # each ideal type in a geographic unit
  if (is.null(bin.proportion)) {
    if (is.null(bin.size)) {
      census <- census %>%
        dplyr::group_by(.dots = c(L1.x, L2.unit)) %>%
        dplyr::summarise(n = dplyr::n())
    } else {
      census$n <- census[[bin.size]]
    }
    census <- census %>%
      dplyr::group_by(.dots = L2.unit) %>%
      dplyr::mutate(prop = n / sum(n))
  } else {
    census <- census %>%
      dplyr::rename(prop = one_of(bin.proportion))
  }

  if (is.null(custom.pc)) {
    # Compute principal components for survey data
    pca_out <- stats::prcomp(survey[, L2.x],
                             retx = TRUE,
                             center = TRUE,
                             scale. = TRUE,
                             tol = NULL)

    # Add PCs to survey data
    survey <- survey %>%
      dplyr::bind_cols(as.data.frame(pca_out$x))

    # Add PCs to census data
    pc_names <- colnames(pca_out$x)

    census <- census %>%
      dplyr::left_join(unique(survey %>% dplyr::select(all_of(L2.unit),
                                                       all_of(pc_names))),
                       by = L2.unit)
  }

  # Scale context-level variables in survey and census data
  if (isTRUE(L2.x.scale)) {
    survey[, L2.x] <- scale(survey[, L2.x], center = TRUE, scale = TRUE)
    census[, L2.x] <- scale(census[, L2.x], center = TRUE, scale = TRUE)
  }

  # Convert survey and census data to tibble
  survey <- tibble::as_tibble(x = survey)
  census <- tibble::as_tibble(x = census)

  # ------------------------------- Create folds -------------------------------

  if (is.null(custom.folds)) {
    # EBMA hold-out fold
    ebma_folding_out <- ebma_folding(data = survey,
                                     L2.unit = L2.unit,
                                     ebma.size = ebma.size)

    ebma_fold <- ebma_folding_out$ebma_fold
    cv_data <- ebma_folding_out$cv_data

    # K folds for cross-validation
    cv_folds <- cv_folding(data = cv_data,
                           L2.unit = L2.unit,
                           k.folds = k.folds,
                           cv.sampling = cv.sampling)
  } else {
    # EBMA hold-out fold
    ebma_fold <- survey %>%
      dplyr::filter_at(dplyr::vars(dplyr::one_of(custom.folds)),
                       dplyr::any_vars(. == k.folds + 1))

    # K folds for cross-validation
    cv_data <- survey %>%
      dplyr::filter_at(dplyr::vars(dplyr::one_of(custom.folds)),
                       dplyr::any_vars(. != k.folds + 1))

    cv_folds <- cv_data %>%
      dplyr::group_split(.data[[custom.folds]])
  }

  # ------------------------ Run individual classifiers ------------------------

  # Classifier 1: Best Subset
  best_subset_out <- best_subset(y = y,
                                 L1.x = L1.x,
                                 L2.x = L2.x,
                                 L2.unit = L2.unit,
                                 L2.reg = L2.reg,
                                 loss.unit = loss.unit,
                                 loss.measure = loss.measure,
                                 data = cv_folds,
                                 verbose = verbose)

  # Classifier 2: Lasso
  lasso_out <- lasso(y = y,
                     L1.x = L1.x,
                     L2.x = L2.x,
                     L2.unit = L2.unit,
                     L2.reg = L2.reg,
                     loss.unit = loss.unit,
                     loss.measure = loss.measure,
                     lambda.set = lasso.lambda.set,
                     iterations.max = lasso.iterations.max,
                     data = cv_folds,
                     verbose = verbose)

  # Classifier 3: PCA
  pca_out <- pca(y = y,
                 L1.x = L1.x,
                 L2.x = pc_names,
                 L2.unit = L2.unit,
                 L2.reg = L2.reg,
                 loss.unit = loss.unit,
                 loss.measure = loss.measure,
                 data = cv_folds,
                 verbose = verbose)

  # Classifier 4: GB
  gb_out <- gb(y = y,
               L1.x = L1.x,
               L2.x = L2.x,
               L2.unit = L2.unit,
               L2.reg = L2.reg,
               L2.unit.include = gb.L2.unit.include,
               L2.reg.include = gb.L2.reg.include,
               loss.unit = loss.unit,
               loss.measure = loss.measure,
               interaction.set = gb.interaction.set,
               shrinkage.set = gb.shrinkage.set,
               tree.start = gb.tree.start,
               tree.increase.set = gb.tree.increase.set,
               trees.max.set = gb.trees.max.set,
               iterations.max = gb.iterations.max,
               n.minobsinnode = gb.n.minobsinnode,
               data = cv_folds,
               verbose = verbose)

  # Classifier 5: SVM
  svm_out <- svm(y = y,
                 L1.x = L1.x,
                 L2.x = L2.x,
                 L2.unit = L2.unit,
                 L2.reg = L2.reg,
                 kernel = svm.kernel,
                 error.fun = svm.error.fun,
                 gamma.set = svm.gamma.set,
                 cost.set = svm.cost.set,
                 k.folds = k.folds,
                 data = cv_folds,
                 verbose = verbose)

  # --------------------------- Post-stratification ----------------------------

  ps_out <- post_stratification(
    data = cv_folds,
    census = census,
    y = y,
    L1.x = L1.x,
    L2.x = L2.x,
    L2.unit = L2.unit,
    L2.reg = L2.reg,
    best.subset = best_subset_out,
    pca = pca_out,
    lasso = lasso_out,
    gb = gb_out,
    n.minobsinnode = gb.n.minobsinnode,
    L2.unit.include = gb.L2.unit.include,
    L2.reg.include = gb.L2.reg.include,
    svm.out = svm_out,
    kernel = svm.kernel,
    verbose = verbose)

  # ----------------------------------- EBMA -----------------------------------

  ebma_out <- ebma(
    ebma.fold = ebma_fold,
    y = y,
    L1.x = L1.x,
    L2.x = L2.x,
    L2.unit = L2.unit,
    L2.reg = L2.reg,
    post.strat = ps_out,
    Ndraws = ebma.n.draws,
    tol.values = ebma.tol.values,
    best.subset = best_subset_out,
    pca = pca_out,
    lasso = lasso_out,
    gb = gb_out,
    svm.out = svm_out,
    verbose = verbose)

  # ---------------------------------- Output ----------------------------------

  return(ebma_out)
}
