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
#'   variable. Must include geographic unit.
#' @param L2.x Context-level covariates. A character vector of column names
#'   corresponding to the context-level variables used to predict the outcome
#'   variable.
#' @param L2.unit Geographic unit. A character scalar indicating the column
#'   name of the geographic unit at which outcomes should be aggregated.
#' @param L2.reg Geographic region. A character scalar indicating the column
#'   name of the geographic region by which geographic units are grouped
#'   (`L2.unit` must be nested within `L2.reg`). Default is NULL.
#' @param survey Survey data. A data.frame containing the y and x column names.
#' @param census Census data. A data.frame containing the x column names.
#' @param bin.size Bin size for ideal types. A character vector indicating the
#'   column name of the variable in census containing the bin size for ideal
#'   types in a geographic unit. Default is NULL.
#' @param ebma.size Size of EBMA hold-out fold. A rational number in the open
#'   unit interval indicating the share of respondents to be contained in the
#'   EBMA hold-out fold. If left unspecified (NULL), then ebma.size is set to
#'   1/3 of the survey sample size. Default is NULL.
#' @param k.folds Number of folds. An integer-valued scalar indicating the
#'   number of folds to be used for cross-validation. Defaults to the value of 5.
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
#'   depth is used to tune the model.
#' @param gb.shrinkage.set Learning rate. A numeric vector whose values define
#'   the learning rate or step-size reduction. Learning rate is used to tune
#'   the model. Values between 0.001 and 0.1 usually work, but a smaller
#'   learning rate typically requires more trees.
#' @param gb.tree.start Initial total number of trees. An integer-valued scalar
#'   specifying the initial number of total trees. Default is 50.
#' @param gb.tree.increase.set Increase in total number of trees. Either an
#'   integer-valued scalar specifying by how many trees the total number of
#'   trees is increased (until the maximum number of trees is reached) or an
#'   integer-valued vector of `length(gb.shrinkage.set)` with each value being
#'   associated with a learning rate. Total number of trees is used to tune the
#'   model.
#' @param gb.trees.max.set Maximum number of trees. Either an integer-valued
#'   scalar specifying the maximum number of trees or an integer-valued vector
#'   of `length(gb.shrinkage.set)` with each value being associated with a
#'   learning rate and a number of tree increase.
#' @param gb.iterations.max Stopping rule. A numeric scalar specifying the
#'   maximum number of iterations without performance improvement the GB
#'   classifier runs before stopping. Default is NULL.
#' @param gb.n.minobsinnode Minimum number of observations in the terminal nodes.
#'   An integer-valued scalar specifying the minimum number of observations
#'   that each terminal node of the trees must contain. Default is 5.
#' @param svm.kernel Kernel for SVM. A character string specifying the kernel to
#'   be used for SVM. The possible types are linear, polynomial, radial, and
#'   sigmoid. Default is radial.
#' @param svm.error.fun
#' @param svm.gamma.set Gamma parameter for SVM. This parameter is needed for
#'   all kernels except linear.
#' @param svm.cost.set Cost parameter for SVM. This parameter specifies the cost
#'   of constraints violation.
#' @param ebma_n_draws
#' @param ebma_tol_values
#' @param seed Seed. An integer-valued scalar to control random number
#'   generation. If left unspecified (NULL), then seed is set to 12345.
#' @param verbose Verbose output. A logical argument indicating whether or not
#'   verbose output should be printed. Default is TRUE.
#' @return
#' @keywords MRP multilevel regression post-stratification machine learning
#'   EBMA ensemble bayesian model averaging
#' @examples
#' @export

auto_MrP <- function(y, L1.x, L2.x, L2.unit, L2.reg = NULL, survey, census,
                     bin.size = NULL, ebma.size = NULL, k.folds = 5,
                     cv.sampling = "L2 units", loss.unit = "individual",
                     loss.measure = "mse", lasso.lambda.set,
                     lasso.iterations.max = NULL, gb.L2.unit.include = FALSE,
                     gb.L2.reg.include = FALSE, gb.interaction.set,
                     gb.shrinkage.set, gb.tree.start = 50,
                     gb.tree.increase.set, gb.trees.max.set,
                     gb.iterations.max = NULL, gb.n.minobsinnode = 5,
                     svm.kernel, svm.error.fun, svm.gamma.set, svm.cost.set,
                     ebma_n_draws, ebma_tol_values,
                     seed = NULL, verbose = TRUE) {
  # Set seed
  if (is.null(seed)) {
    set.seed(12345)
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

  # ------------------------------- Prepare data -------------------------------

  # If not provided in census data, calculate bin size for each ideal type in
  # a geographic unit
  if (is.null(bin.size)) {
    census <- census %>%
      dplyr::group_by(.dots = c(L1.x, L2.unit)) %>%
      dplyr::summarise(n = dplyr::n())
  } else {
    census$n <- census[[bin.size]]
  }

  # In census data, calculate bin proportion for each ideal type in a
  # geographic unit
  census <- census %>%
    dplyr::group_by(.dots = L2.unit) %>%
    dplyr::mutate(prop = n / sum(n))

  # Scale context-level variables in survey and census data
  survey[, L2.x] <- scale(survey[, L2.x], center = TRUE, scale = TRUE)
  census[, L2.x] <- scale(census[, L2.x], center = TRUE, scale = TRUE)

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
    dplyr::left_join(unique(dplyr::select(survey, L2.unit, pc_names)),
                     by = L2.unit)

  # ------------------------------- Create folds -------------------------------

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

  ps_out <- post_stratification(data = cv_folds,
                                census = census,
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

  ebma_out <- ebma(ebma.fold = ebma_fold,
                   L1.x = L1.x,
                   L2.x = L2.x,
                   L2.unit = L2.unit,
                   L2.reg = L2.reg,
                   post.strat = ps_out,
                   Ndraws = ebma_n_draws,
                   tol.values = ebma_tol_values,
                   best.subset = best_subset_out,
                   pca = pca_out,
                   lasso = lasso_out,
                   gb = gb_out,
                   svm.out = svm_out,
                   verbose = verbose)

  # ---------------------------------- Output ----------------------------------

  return(ebma_out)
}
