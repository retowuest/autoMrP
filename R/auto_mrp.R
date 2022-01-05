#' Improve MrP through ensemble learning.
#'
#' This package improves the prediction performance of multilevel
#' regression with post-stratification (MrP) by combining a number of machine
#' learning methods through ensemble Bayesian model averaging (EBMA).
#'
#' @param y Outcome variable. A character vector containing the column names of
#'   the outcome variable. A character scalar containing the column name of
#'   the outcome variable in \code{survey}.
#' @param L1.x Individual-level covariates. A character vector containing the
#'   column names of the individual-level variables in \code{survey} and
#'   \code{census} used to predict outcome \code{y}. Note that geographic unit
#'   is specified in argument \code{L2.unit}.
#' @param L2.x Context-level covariates. A character vector containing the
#'   column names of the context-level variables in \code{survey} and
#'   \code{census} used to predict outcome \code{y}. Note: For the empty MrP
#'   model, set \code{L2.x = NULL} and \code{mrp.L2.x = ""}.
#' @param L2.unit Geographic unit. A character scalar containing the column
#'   name of the geographic unit in \code{survey} and \code{census} at which
#'   outcomes should be aggregated.
#' @param L2.reg Geographic region. A character scalar containing the column
#'   name of the geographic region in \code{survey} and \code{census} by which
#'   geographic units are grouped (\code{L2.unit} must be nested within
#'   \code{L2.reg}). Default is \code{NULL}.
#' @param L2.x.scale Scale context-level covariates. A logical argument
#'   indicating whether the context-level covariates should be normalized.
#'   Default is \code{TRUE}. Note that if set to \code{FALSE}, then the
#'   context-level covariates should be normalized prior to calling
#'   \code{auto_MrP()}.
#' @param pcs Principal components. A character vector containing the column
#'   names of the principal components of the context-level variables in
#'   \code{survey} and \code{census}. Default is \code{NULL}.
#' @param folds EBMA and cross-validation folds. A character scalar containing
#'   the column name of the variable in \code{survey} that specifies the fold
#'   to which an observation is allocated. The variable should contain integers
#'   running from \eqn{1} to \eqn{k + 1}, where \eqn{k} is the number of
#'   cross-validation folds. Value \eqn{k + 1} refers to the EBMA fold. Default
#'   is \code{NULL}. \emph{Note:} if \code{folds} is \code{NULL}, then
#'   \code{ebma.size}, \code{k.folds}, and \code{cv.sampling} must be specified.
#' @param bin.proportion Proportion of ideal types. A character scalar
#'   containing the column name of the variable in \code{census} that indicates
#'   the proportion of individuals by ideal type and geographic unit. Default is
#'   \code{NULL}. \emph{Note:} if \code{bin.proportion} is \code{NULL}, then
#'   \code{bin.size} must be specified.
#' @param bin.size Bin size of ideal types. A character scalar containing the
#'   column name of the variable in \code{census} that indicates the bin size of
#'   ideal types by geographic unit. Default is \code{NULL}. \emph{Note:}
#'   ignored if \code{bin.proportion} is provided, but must be specified
#'   otherwise.
#' @param survey Survey data. A \code{data.frame} whose column names include
#'   \code{y}, \code{L1.x}, \code{L2.x}, \code{L2.unit}, and, if specified,
#'   \code{L2.reg}, \code{pcs}, and \code{folds}.
#' @param census Census data. A \code{data.frame} whose column names include
#'   \code{L1.x}, \code{L2.x}, \code{L2.unit}, if specified, \code{L2.reg} and
#'   \code{pcs}, and either \code{bin.proportion} or \code{bin.size}.
#' @param ebma.size EBMA fold size. A number in the open unit interval
#'   indicating the proportion of respondents to be allocated to the EBMA fold.
#'   Default is \eqn{1/3}. \emph{Note:} ignored if \code{folds} is provided, but
#'   must be specified otherwise.
#' @param cores The number of cores to be used. An integer indicating the number
#'   of processor cores used for parallel computing. Default is 1.
#' @param k.folds Number of cross-validation folds. An integer-valued scalar
#'   indicating the number of folds to be used in cross-validation. Default is
#'   \eqn{5}. \emph{Note:} ignored if \code{folds} is provided, but must be
#'   specified otherwise.
#' @param cv.sampling Cross-validation sampling method. A character-valued
#'   scalar indicating whether cross-validation folds should be created by
#'   sampling individual respondents (\code{individuals}) or geographic units
#'   (\code{L2 units}). Default is \code{L2 units}. \emph{Note:} ignored if
#'   \code{folds} is provided, but must be specified otherwise.
#' @param loss.unit Loss function unit. A character-valued scalar indicating
#'   whether performance loss should be evaluated at the level of individual
#'   respondents (\code{individuals}), geographic units (\code{L2 units}) or at
#'   both levels. Default is \code{c("individuals", "L2 units")}. With multiple
#'   loss units, parameters are ranked for each loss unit and the loss unit with
#'   the lowest rank sum is chosen. Ties are broken according to the order in
#'   the search grid.
#' @param loss.fun Loss function. A character-valued scalar indicating whether
#'   prediction loss should be measured by the mean squared error (\code{MSE}),
#'   the mean absolute error (\code{MAE}), binary cross-entropy
#'   (\code{cross-entropy}), mean squared false error (\code{msfe}), the f1
#'   score (\code{f1}), or a combination thereof. Default is \code{c("MSE",
#'   "cross-entropy","msfe", "f1")}. With multiple loss functions, parameters
#'   are ranked for each loss function and the parameter combination with the
#'   lowest rank sum is chosen. Ties are broken according to the order in the
#'   search grid.
#' @param best.subset Best subset classifier. A logical argument indicating
#'   whether the best subset classifier should be used for predicting outcome
#'   \code{y}. Default is \code{TRUE}.
#' @param lasso Lasso classifier. A logical argument indicating whether the
#'   lasso classifier should be used for predicting outcome \code{y}. Default is
#'   \code{TRUE}.
#' @param pca PCA classifier. A logical argument indicating whether the PCA
#'   classifier should be used for predicting outcome \code{y}. Default is
#'   \code{TRUE}.
#' @param gb GB classifier. A logical argument indicating whether the GB
#'   classifier should be used for predicting outcome \code{y}. Default is
#'   \code{TRUE}.
#' @param svm SVM classifier. A logical argument indicating whether the SVM
#'   classifier should be used for predicting outcome \code{y}. Default is
#'   \code{TRUE}.
#' @param mrp MRP classifier. A logical argument indicating whether the standard
#'   MRP classifier should be used for predicting outcome \code{y}. Default is
#'   \code{FALSE}.
#' @param oversampling Over sample to create balance on the dependent variable.
#'   A logical argument. Default is \code{FALSE}.
#' @param forward.select Forward selection classifier. A logical argument
#'   indicating whether to use forward selection rather than best subset
#'   selection. Default is \code{FALSE}. \emph{Note:} forward selection is
#'   recommended if there are more than \eqn{8} context-level variables.
#'   \emph{Note:} forward selection is not implemented yet.
#' @param best.subset.L2.x Best subset context-level covariates. A character
#'   vector containing the column names of the context-level variables in
#'   \code{survey} and \code{census} to be used by the best subset classifier.
#'   If \code{NULL} and \code{best.subset} is set to \code{TRUE}, then best
#'   subset uses the variables specified in \code{L2.x}. Default is \code{NULL}.
#' @param lasso.L2.x Lasso context-level covariates. A character vector
#'   containing the column names of the context-level variables in
#'   \code{survey} and \code{census} to be used by the lasso classifier. If
#'   \code{NULL} and \code{lasso} is set to \code{TRUE}, then lasso uses the
#'   variables specified in \code{L2.x}. Default is \code{NULL}.
#' @param pca.L2.x PCA context-level covariates. A character vector containing
#'   the column names of the context-level variables in \code{survey} and
#'   \code{census} whose principal components are to be used by the PCA
#'   classifier. If \code{NULL} and \code{pca} is set to \code{TRUE}, then PCA
#'   uses the principal components of the variables specified in \code{L2.x}.
#'   Default is \code{NULL}.
#' @param gb.L2.x GB context-level covariates. A character vector containing the
#'   column names of the context-level variables in \code{survey} and
#'   \code{census} to be used by the GB classifier. If \code{NULL} and \code{gb}
#'   is set to \code{TRUE}, then GB uses the variables specified in \code{L2.x}.
#'   Default is \code{NULL}.
#' @param svm.L2.x SVM context-level covariates. A character vector containing
#'   the column names of the context-level variables in \code{survey} and
#'   \code{census} to be used by the SVM classifier. If \code{NULL} and
#'   \code{svm} is set to \code{TRUE}, then SVM uses the variables specified in
#'   \code{L2.x}. Default is \code{NULL}.
#' @param mrp.L2.x MRP context-level covariates. A character vector containing
#'   the column names of the context-level variables in \code{survey} and
#'   \code{census} to be used by the MRP classifier. The character vector
#'   \emph{empty} if no context-level variables should be used by the MRP
#'   classifier. If \code{NULL} and \code{mrp} is set to \code{TRUE}, then MRP
#'   uses the variables specified in \code{L2.x}. Default is \code{NULL}. Note:
#'   For the empty MrP model, set \code{L2.x = NULL} and \code{mrp.L2.x = ""}.
#' @param gb.L2.unit GB L2.unit. A logical argument indicating whether
#'   \code{L2.unit} should be included in the GB classifier. Default is
#'   \code{FALSE}.
#' @param gb.L2.reg GB L2.reg. A logical argument indicating whether
#'   \code{L2.reg} should be included in the GB classifier. Default is
#'   \code{FALSE}.
#' @param svm.L2.unit SVM L2.unit. A logical argument indicating whether
#'   \code{L2.unit} should be included in the SVM classifier. Default is
#'   \code{FALSE}.
#' @param svm.L2.reg SVM L2.reg. A logical argument indicating whether
#'   \code{L2.reg} should be included in the SVM classifier. Default is
#'   \code{FALSE}.
#' @param lasso.lambda Lasso penalty parameter. A numeric \code{vector} of
#'   non-negative values. The penalty parameter controls the shrinkage of the
#'   context-level variables in the lasso model. Default is a sequence with
#'   minimum 0.1 and maximum 250 that is equally spaced on the log-scale. The
#'   number of values is controlled by the \code{lasso.n.iter} parameter.
#' @param lasso.n.iter Lasso number of lambda values. An integer-valued scalar
#'   specifying the number of lambda values to search over. Default is \eqn{100}.
#'   \emph{Note:} Is ignored if a vector of \code{lasso.lambda} values is
#'   provided.
#' @param gb.interaction.depth GB interaction depth. An integer-valued vector
#'   whose values specify the interaction depth of GB. The interaction depth
#'   defines the maximum depth of each tree grown (i.e., the maximum level of
#'   variable interactions). Default is \code{c(1, 2, 3)}.
#' @param gb.shrinkage GB learning rate. A numeric vector whose values specify
#'   the learning rate or step-size reduction of GB. Values between \eqn{0.001}
#'   and \eqn{0.1} usually work, but a smaller learning rate typically requires
#'   more trees. Default is \code{c(0.04, 0.01, 0.008, 0.005, 0.001)}.
#' @param gb.n.trees.init GB initial total number of trees. An integer-valued
#'   scalar specifying the initial number of total trees to fit by GB. Default
#'   is \eqn{50}.
#' @param gb.n.trees.increase GB increase in total number of trees. An
#'   integer-valued scalar specifying by how many trees the total number of
#'   trees to fit should be increased (until \code{gb.n.trees.max} is reached).
#'   Default is \eqn{50}.
#' @param gb.n.trees.max GB maximum number of trees. An integer-valued scalar
#'   specifying the maximum number of trees to fit by GB. Default is \eqn{1000}.
#' @param gb.n.minobsinnode GB minimum number of observations in the terminal
#'   nodes. An integer-valued scalar specifying the minimum number of
#'   observations that each terminal node of the trees must contain. Default is
#'   \eqn{20}.
#' @param svm.kernel SVM kernel. A character-valued scalar specifying the kernel
#'   to be used by SVM. The possible values are \code{linear}, \code{polynomial},
#'   \code{radial}, and \code{sigmoid}. Default is \code{radial}.
#' @param svm.gamma SVM kernel parameter. A numeric vector whose values specify
#'   the gamma parameter in the SVM kernel. This parameter is needed for all
#'   kernel types except linear. Default is a sequence with minimum = 1e-5,
#'   maximum = 1e-1, and length = 20 that is equally spaced on the log-scale.
#' @param svm.cost SVM cost parameter. A numeric vector whose values specify the
#'   cost of constraints violation in SVM. Default is a sequence with minimum =
#'   0.5, maximum = 10, and length = 5 that is equally spaced on the log-scale.
#' @param ebma.n.draws EBMA number of samples. An integer-valued scalar
#'   specifying the number of bootstrapped samples to be drawn from the EBMA
#'   fold and used for tuning EBMA. Default is \eqn{100}.
#' @param ebma.tol EBMA tolerance. A numeric vector containing the
#'   tolerance values for improvements in the log-likelihood before the EM
#'   algorithm stops optimization. Values should range at least from \eqn{0.01}
#'   to \eqn{0.001}. Default is
#'   \code{c(0.01, 0.005, 0.001, 0.0005, 0.0001, 0.00005, 0.00001)}.
#' @param uncertainty Uncertainty estimates. A logical argument indicating
#'   whether uncertainty estimates should be computed. Default is \code{FALSE}.
#' @param boot.iter Number of bootstrap iterations. An integer argument
#'   indicating the number of bootstrap iterations to be computed. Will be
#'   ignored unless \code{uncertainty = TRUE}. Default is \code{200} if
#'   \code{uncertainty = TRUE} and \code{NULL} if \code{uncertainty = FALSE}.
#' @param seed Seed. Either \code{NULL} or an integer-valued scalar controlling
#'   random number generation. If \code{NULL}, then the seed is set to
#'   \eqn{546213978}. Default is \code{NULL}.
#' @param verbose Verbose output. A logical argument indicating whether or not
#'   verbose output should be printed. Default is \code{FALSE}.
#' @return The context-level predictions. A list with two elements. The first
#'   element, \code{EBMA}, contains the post-stratified ensemble bayesian model
#'   avaeraging (EBMA) predictions. The second element, \code{classifiers},
#'   contains the post-stratified predictions from all estimated classifiers.
#' @keywords MRP multilevel regression post-stratification machine learning
#'   EBMA ensemble Bayesian model averaging
#' @examples
#' # An MrP model without machine learning
#' m <- auto_MrP(
#'   y = "YES",
#'   L1.x = c("L1x1"),
#'   L2.x = c("L2.x1", "L2.x2"),
#'   L2.unit = "state",
#'   bin.proportion = "proportion",
#'   survey = taxes_survey,
#'   census = taxes_census,
#'   ebma.size = 0,
#'   cores = max_cores,
#'   best.subset = FALSE,
#'   lasso = FALSE,
#'   pca = FALSE,
#'   gb = FALSE,
#'   svm = FALSE,
#'   mrp = TRUE
#' )
#'
#' # summarize and plot results
#' summary(m)
#' plot(m)
#'
#' # The empty MrP model
#' m <- auto_MrP(
#'   y = "YES",
#'   L1.x = "L1x1",
#'   L2.x = NULL,
#'   mrp.L2.x = "",
#'   L2.unit = "state",
#'   bin.proportion = "proportion",
#'   survey = taxes_survey,
#'   census = taxes_census,
#'   ebma.size = 0,
#'   cores = 1,
#'   best.subset = FALSE,
#'   lasso = FALSE,
#'   pca = FALSE,
#'   gb = FALSE,
#'   svm = FALSE,
#'   mrp = TRUE
#'   )
#'
#' \donttest{
#' # MrP model only:
#' mrp_out <- auto_MrP(
#'   y = "YES",
#'   L1.x = c("L1x1", "L1x2", "L1x3"),
#'   L2.x = c("L2.x1", "L2.x2", "L2.x3", "L2.x4", "L2.x5", "L2.x6"),
#'   L2.unit = "state",
#'   L2.reg = "region",
#'   bin.proportion = "proportion",
#'   survey = taxes_survey,
#'   census = taxes_census,
#'   ebma.size = 0,
#'   best.subset = FALSE,
#'   lasso = FALSE,
#'   pca = FALSE,
#'   gb = FALSE,
#'   svm = FALSE,
#'   mrp = TRUE
#' )
#'
#' # Predictions through machine learning
#'
#' # detect number of available cores
#' max_cores <- parallel::detectCores()
#'
#' # autoMrP with machine learning
#' ml_out <- auto_MrP(
#'   y = "YES",
#'   L1.x = c("L1x1", "L1x2", "L1x3"),
#'   L2.x = c("L2.x1", "L2.x2", "L2.x3", "L2.x4", "L2.x5", "L2.x6"),
#'   L2.unit = "state",
#'   L2.reg = "region",
#'   bin.proportion = "proportion",
#'   survey = taxes_survey,
#'   census = taxes_census,
#'   gb.L2.reg = TRUE,
#'   svm.L2.reg = TRUE,
#'   cores = max_cores
#'   )
#'   }
#' @export
#' @importFrom stats as.formula binomial predict setNames weighted.mean median sd
#' @importFrom utils combn
#' @importFrom dplyr %>%
#' @importFrom rlang .data
#' @importFrom foreach %dopar%
#' @importFrom doRNG %dorng%


auto_MrP <- function(y, L1.x, L2.x, L2.unit, L2.reg = NULL, L2.x.scale = TRUE, pcs = NULL,
                     folds = NULL, bin.proportion = NULL, bin.size = NULL, survey, census,
                     ebma.size = 1/3, cores = 1, k.folds = 5, cv.sampling = "L2 units",
                     loss.unit = c("individuals", "L2 units"),
                     loss.fun = c("msfe", "cross-entropy", "f1", "MSE"),
                     best.subset = TRUE, lasso = TRUE, pca = TRUE, gb = TRUE, svm = TRUE,
                     mrp = FALSE, oversampling = FALSE, forward.select = FALSE,
                     best.subset.L2.x = NULL, lasso.L2.x = NULL, pca.L2.x = NULL,
                     gb.L2.x = NULL, svm.L2.x = NULL, mrp.L2.x = NULL, gb.L2.unit = TRUE,
                     gb.L2.reg = FALSE, svm.L2.unit = TRUE, svm.L2.reg = FALSE,
                     lasso.lambda = NULL,
                     lasso.n.iter = 100,
                     gb.interaction.depth = c(1, 2, 3),
                     gb.shrinkage = c(0.04, 0.01, 0.008, 0.005, 0.001),
                     gb.n.trees.init = 50,
                     gb.n.trees.increase = 50,
                     gb.n.trees.max = 1000,
                     gb.n.minobsinnode = 20,
                     svm.kernel = c("radial"),
                     svm.gamma = NULL,
                     svm.cost = NULL,
                     ebma.n.draws = 100,
                     ebma.tol = c(0.01, 0.005, 0.001, 0.0005, 0.0001, 0.00005, 0.00001),
                     seed = NULL,
                     verbose = FALSE,
                     uncertainty = FALSE,
                     boot.iter = NULL) {


# Error checks ------------------------------------------------------------

  # Call to function doing the error checks
  error_checks(y = y,
               L1.x = L1.x,
               L2.x = L2.x,
               L2.unit = L2.unit,
               L2.reg = L2.reg,
               L2.x.scale = L2.x.scale,
               pcs = pcs,
               folds = folds,
               bin.proportion = bin.proportion,
               bin.size = bin.size,
               survey = survey,
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
               gb.L2.x = gb.L2.x,
               svm.L2.x = svm.L2.x,
               mrp.L2.x = mrp.L2.x,
               gb.L2.unit = gb.L2.unit,
               gb.L2.reg = gb.L2.reg,
               lasso.lambda = lasso.lambda,
               lasso.n.iter = lasso.n.iter,
               uncertainty = uncertainty,
               boot.iter = boot.iter)


# Prepare data ------------------------------------------------------------

    # Coerce individual-level variables and geographic variables to factors in
    # survey and census data
    survey <- survey %>%
      dplyr::mutate_at(.vars = c(L1.x, L2.unit, L2.reg), .funs = as.factor)

    census <- census %>%
      dplyr::mutate_at(.vars = c(L1.x, L2.unit, L2.reg), .funs = as.factor)

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

    # If not provided in survey and census data, compute the principal components
    # of context-level variables
    if (is.null(pcs) & !is.null(L2.x)) {

      # Determine context-level covariates whose principal components are to be
      # computed
      if (is.null(pca.L2.x)) {
        pca.L2.x <- L2.x
      }

      # Compute principal components for survey data
      pca_out <- stats::prcomp(survey[, pca.L2.x],
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
    } else {
      pc_names <- pcs
    }

    # Scale context-level variables in survey and census data
    if (isTRUE(L2.x.scale) & all(L2.x != "")) {

      # scale context-level variables in survey
      survey <- dplyr::mutate_at(
        .tbl = survey,
        .vars = L2.x,
        .funs = function(x) base::as.numeric(base::scale(x = x, center = TRUE, scale = TRUE)))

      # scale context-level variables in census
      census <- dplyr::mutate_at(
        .tbl = census,
        .vars = L2.x,
        .funs = function(x) base::as.numeric(base::scale(x = x, center = TRUE, scale = TRUE)))
    }

    # Convert survey and census data to tibble
    survey <- tibble::as_tibble(x = survey)
    census <- tibble::as_tibble(x = census)

    # Random over-sampling
    if ( isTRUE(oversampling) ){
      add_rows <- survey %>%
        dplyr::group_by( .dots = L2.unit ) %>%
        tidyr::nest() %>%
        dplyr::mutate(os = purrr::map(data, function( x ){
          n <- nrow(x)
          os <- dplyr::group_by(.data = x, !! rlang::sym(y) )
          y_1 <- sum(dplyr::pull(.data = os, var = !! rlang::sym(y)))
          y_0 <- n - y_1
          if (y_1 > 0 & y_0 > 0){
            y_needed <- ifelse(test = y_1 > y_0, yes = 0, no = 1)
            n_needed <- ifelse(test = y_needed == 0, yes = y_1 - y_0, no = y_0 - y_1)
            os <- dplyr::filter(.data = os, !! rlang::sym(y) == y_needed )
            os <- dplyr::slice_sample(.data = os, replace = TRUE, n = n_needed)
          }
          return(os)
        })) %>%
        tidyr::unnest(os) %>%
        dplyr::ungroup()
      survey <- dplyr::bind_rows(survey, add_rows)
    }

# No bootstrapping --------------------------------------------------------

    if (!uncertainty){

# Create folds ------------------------------------------------------------

    if (is.null(folds)) {

      # EBMA hold-out fold
      ebma.size <- round(nrow(survey) * ebma.size, digits = 0)

      if(ebma.size>0){
        ebma_folding_out <- ebma_folding(data = survey,
                                         L2.unit = L2.unit,
                                         ebma.size = ebma.size)
        ebma_fold <- ebma_folding_out$ebma_fold
        cv_data <- ebma_folding_out$cv_data
      } else{
        ebma_fold <- NULL
        cv_data <- survey
      }

      # K folds for cross-validation
      cv_folds <- cv_folding(data = cv_data,
                             L2.unit = L2.unit,
                             k.folds = k.folds,
                             cv.sampling = cv.sampling)
    } else {

      if (ebma.size > 0){
        # EBMA hold-out fold
        ebma_fold <- survey %>%
          dplyr::filter_at(dplyr::vars(dplyr::one_of(folds)),
                           dplyr::any_vars(. == k.folds + 1))
      }

      # K folds for cross-validation
      cv_data <- survey %>%
        dplyr::filter_at(dplyr::vars(dplyr::one_of(folds)),
                         dplyr::any_vars(. != k.folds + 1))

      cv_folds <- cv_data %>%
        dplyr::group_split(.data[[folds]])
    }

# Optimal individual classifiers ------------------------------------------

      ebma_out <- run_classifiers(
        y = y, L1.x = L1.x, L2.x = L2.x, mrp.L2.x = mrp.L2.x,
        L2.unit = L2.unit, L2.reg = L2.reg, pcs = pcs,
        folds = folds, cv.folds = cv_folds, cv.data = cv_data,
        ebma.fold = ebma_fold, census = census, k.folds = k.folds,
        cv.sampling = cv.sampling,
        loss.unit = loss.unit, loss.fun = loss.fun,
        best.subset = best.subset, lasso = lasso, pca = pca,
        gb = gb, svm = svm, mrp = mrp, forward.select = forward.select,
        best.subset.L2.x = best.subset.L2.x,
        lasso.L2.x = lasso.L2.x, pca.L2.x = pca.L2.x,
        gb.L2.x = gb.L2.x, svm.L2.x = svm.L2.x,
        gb.L2.unit = gb.L2.unit, gb.L2.reg = gb.L2.reg,
        lasso.lambda = lasso.lambda, lasso.n.iter = lasso.n.iter,
        gb.interaction.depth = gb.interaction.depth,
        gb.shrinkage = gb.shrinkage,
        gb.n.trees.init = gb.n.trees.init,
        gb.n.trees.increase = gb.n.trees.increase,
        gb.n.trees.max = gb.n.trees.max,
        gb.n.minobsinnode = gb.n.minobsinnode,
        svm.kernel = svm.kernel, svm.gamma = svm.gamma,
        svm.cost = svm.cost, ebma.tol = ebma.tol, cores = cores,
        verbose = verbose)

# Boostrapping wrapper ----------------------------------------------------

  } else{

    if (is.null(boot.iter)){
      boot.iter <- 200
    }

    ebma_out <- boot_auto_mrp(
      y = y, L1.x = L1.x, L2.x = L2.x, mrp.L2.x = mrp.L2.x,
      L2.unit = L2.unit, L2.reg = L2.reg, L2.x.scale = L2.x.scale,
      pcs = pcs, folds = folds, bin.proportion = bin.proportion,
      bin.size = bin.size, survey = survey, census = census,
      ebma.size = ebma.size, k.folds = k.folds,
      cv.sampling = cv.sampling, loss.unit = loss.unit,
      loss.fun = loss.fun, best.subset = best.subset,
      lasso = lasso, pca = pca, gb = gb, svm = svm, mrp = mrp,
      forward.select = forward.select,
      best.subset.L2.x = best.subset.L2.x,
      lasso.L2.x = lasso.L2.x, pca.L2.x = pca.L2.x,
      gb.L2.x = gb.L2.x, svm.L2.x = svm.L2.x,
      gb.L2.unit = gb.L2.unit, gb.L2.reg = gb.L2.reg,
      lasso.lambda = lasso.lambda, lasso.n.iter = lasso.n.iter,
      gb.interaction.depth = gb.interaction.depth,
      gb.shrinkage = gb.shrinkage,
      gb.n.trees.init = gb.n.trees.init,
      gb.n.trees.increase = gb.n.trees.increase,
      gb.n.trees.max = gb.n.trees.max,
      gb.n.minobsinnode = gb.n.minobsinnode,
      svm.kernel = svm.kernel, svm.gamma = svm.gamma,
      svm.cost = svm.cost, ebma.tol = ebma.tol,
      boot.iter = boot.iter, cores = cores)
  }


# autoMrP function output -------------------------------------------------

  class(ebma_out) <- c("autoMrP", "list")
  class(ebma_out$ebma) <- c("autoMrP", "ensemble", class(ebma_out$ebma))
  class(ebma_out$classifiers) <- c("autoMrP", "classifiers", class(ebma_out$classifiers))
  if ("weights" %in% names(ebma_out)){
    class(ebma_out$weights) <- c("autoMrP", "weights", class(ebma_out$weights))
  } else{
    ebma_out$weights <- "EBMA step skipped (only 1 classifier run)"
    class(ebma_out$weights) <- c("autoMrP", "weights", class(ebma_out$weights))
  }
  return(ebma_out)
}
