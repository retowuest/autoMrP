################################################################################
#            Function to check if arguments to auto_MrP() are valid            #
################################################################################

#' Catches user input errors
#'
#' \code{error_checks()} checks for incorrect data entry in \code{autoMrP()}
#' call.
#'
#' @param y Outcome variable. A character vector containing the column names of
#'   the outcome variable.
#' @param L1.x #' @param L1.x Individual-level covariates. A character vector
#'   containing the column names of the individual-level variables in
#'   \code{survey} and \code{census} used to predict outcome \code{y}. Note that
#'   geographic unit is specified in argument \code{L2.unit}.
#' @param L2.x Context-level covariates. A character vector containing the
#'   column names of the context-level variables in \code{survey} and
#'   \code{census} used to predict outcome \code{y}.
#' @param L2.unit Geographic unit. A character scalar containing the column name
#'   of the geographic unit in \code{survey} and \code{census} at which outcomes
#'   should be aggregated.
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
#'   the column name of the variable in \code{survey} that specifies the fold to
#'   which an observation is allocated. The variable should contain integers
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
#'   respondents (\code{individuals}) or geographic units (\code{L2 units}).
#'   Default is \code{individuals}.
#' @param loss.fun Loss function. A character-valued scalar indicating whether
#'   prediction loss should be measured by the mean squared error (\code{MSE})
#'   or the mean absolute error (\code{MAE}). Default is \code{MSE}.
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
#' @param forward.select Forward selection classifier. A logical argument
#'   indicating whether to use forward selection rather than best subset
#'   selection. Default is \code{FALSE}. \emph{Note:} forward selection is
#'   recommended if there are more than \eqn{8} context-level variables.
#' @param best.subset.L2.x Best subset context-level covariates. A character
#'   vector containing the column names of the context-level variables in
#'   \code{survey} and \code{census} to be used by the best subset classifier.
#'   If \code{NULL} and \code{best.subset} is set to \code{TRUE}, then best
#'   subset uses the variables specified in \code{L2.x}. Default is \code{NULL}.
#' @param lasso.L2.x Lasso context-level covariates. A character vector
#'   containing the column names of the context-level variables in \code{survey}
#'   and \code{census} to be used by the lasso classifier. If \code{NULL} and
#'   \code{lasso} is set to \code{TRUE}, then lasso uses the variables specified
#'   in \code{L2.x}. Default is \code{NULL}.
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
#'   uses the variables specified in \code{L2.x}. Default is \code{NULL}.
#' @param gb.L2.unit GB L2.unit. A logical argument indicating whether
#'   \code{L2.unit} should be included in the GB classifier. Default is
#'   \code{FALSE}.
#' @param gb.L2.reg GB L2.reg. A logical argument indicating whether
#'   \code{L2.reg} should be included in the GB classifier. Default is
#'   \code{FALSE}.
#' @param lasso.lambda Lasso penalty parameter. A numeric \code{vector} of
#'   non-negative values or a \code{list} of two numeric vectors of equal size,
#'   with the first vector containing the step sizes by which the penalty
#'   parameter should increase and the second vector containing the upper
#'   thresholds of the intervals to which the step sizes apply. The penalty
#'   parameter controls the shrinkage of the context-level variables in the
#'   lasso model. Default is \code{list(c(0.1, 0.3, 1), c(1, 10, 10000))}.
#' @param lasso.n.iter Lasso number of iterations without improvement. Either
#'   \code{NULL} or an integer-valued scalar specifying the maximum number of
#'   iterations without performance improvement the algorithm runs before
#'   stopping. Default is \eqn{70}.

error_checks <- function(y, L1.x, L2.x, L2.unit, L2.reg, L2.x.scale, pcs,
                         folds, bin.proportion, bin.size, survey, census,
                         ebma.size, k.folds, cv.sampling, loss.unit, loss.fun,
                         best.subset, lasso, pca, gb, svm, mrp, forward.select,
                         best.subset.L2.x, lasso.L2.x, gb.L2.x, svm.L2.x,
                         mrp.L2.x, gb.L2.unit, gb.L2.reg, lasso.lambda,
                         lasso.n.iter) {


  # Check if y is a character scalar
  if (!(is.character(y) & length(y) == 1)) {
    stop(paste("The argument 'y', specifying the outcome variable, must be a",
               " character scalar.", sep = ""))
  }

  # Check if y is in survey data
  if (!(y %in% colnames(survey))) {
    stop(paste("Outcome '", y,
               "' is not in your survey data.", sep = ""))
  }

  # Check if L1.x is a character vector
  if (!is.character(L1.x)) {
    stop(paste("The argument 'L1.x', specifying the individual-level variables",
               " to be used to predict y, must be a character vector.",
               sep = ""))
  }

  # Check if L1.x is in survey data
  if (!all(L1.x %in% colnames(survey))) {
    stop(cat(paste("Individual-level variable '",
                   L1.x[which(!(L1.x %in% colnames(survey)))],
                   "', specified in argument 'L1.x', is not in your survey",
                   " data.\n", sep = ""), sep = ""))
  }

  # Check if L1.x is in census data
  if (!all(L1.x %in% colnames(census))) {
    stop(cat(paste("Individual-level variable '",
                   L1.x[which(!(L1.x %in% colnames(census)))],
                   "', specified in argument 'L1.x', is not in your census",
                   " data.\n", sep = ""), sep = ""))
  }

  # Check if L2.x is a character vector
  if (!is.character(L2.x)) {
    stop(paste("The argument 'L2.x', specifying the context-level variables to",
               " be used to predict y, must be a character vector.", sep = ""))
  }

  # Check if L2.x is in survey data
  if (!all(L2.x %in% colnames(survey))) {
    stop(cat(paste("Context-level variable '",
                   L2.x[which(!(L2.x %in% colnames(survey)))],
                   "', specified in argument 'L2.x', is not in your survey",
                   " data.\n", sep = ""), sep = ""))
  }

  # Check if L2.x is in census data
  if (!all(L2.x %in% colnames(census))) {
    stop(cat(paste("Context-level variable '",
                   L2.x[which(!(L2.x %in% colnames(census)))],
                   "', specified in argument 'L2.x', is not in your census",
                   " data.\n", sep = ""), sep = ""))
  }

  # Check if L2.unit is a character scalar
  if (!(is.character(L2.unit) & length(L2.unit) == 1)) {
    stop(paste("The argument 'L2.unit', specifying the geographic unit,",
               " must be a character scalar.", sep = ""))
  }

  # Check if L2.unit is in survey data
  if (!(L2.unit %in% colnames(survey))) {
    stop(paste("The geographic unit '", L2.unit,
               "' is not in your survey data.", sep = ""))
  }

  # Check if L2.unit is in census data
  if (!(L2.unit %in% colnames(census))) {
    stop(paste("The geographic unit '", L2.unit,
               "' is not in your census data.", sep = ""))
  }

  # Check if L2.reg is NULL
  if (!is.null(L2.reg)) {
    # Check if L2.reg is a character scalar
    if (!(is.character(L2.reg) & length(L2.reg) == 1)) {
      stop(paste("The argument 'L2.reg', specifying the geographic region,",
                 " must be a character scalar.", sep = ""))
    }

    # Check if L2.reg is in survey data
    if (!(L2.reg %in% colnames(survey))) {
      stop(paste("The geographic region '", L2.reg,
                 "' is not in your survey data.", sep = ""))
    }

    # Check if L2.reg is in census data
    if (!(L2.reg %in% colnames(census))) {
      stop(paste("The geographic region '", L2.reg,
                 "' is not in your census data.", sep = ""))
    }

    # Check if each geographic unit is nested in only one geographic region
    # in survey data
    if (any(unlist(lapply(dplyr::group_split(survey, .data[[L2.unit]]),
                          function(x) length(unique(x[[L2.reg]])))) > 1)) {
      stop(cat(paste("The geographic unit '",
                     L2.unit[which(unlist(lapply(dplyr::group_split(survey, .data[[L2.unit]]),
                                                 function(x) length(unique(x[[L2.reg]])))) > 1)],
                     "' is nested in multiple regions in your survey data.\n"), sep = ""))
    }

    # Check if each geographic unit is nested in only one geographic region
    # in census data
    if (any(unlist(lapply(dplyr::group_split(census, .data[[L2.unit]]),
                          function(x) length(unique(x[[L2.reg]])))) > 1)) {
      stop(cat(paste("The geographic unit '",
                     L2.unit[which(unlist(lapply(dplyr::group_split(census, .data[[L2.unit]]),
                                                 function(x) length(unique(x[[L2.reg]])))) > 1)],
                     "' is nested in multiple regions in your census data.\n"), sep = ""))
    }
  }

  # Check if L2.x.scale is logical
  if (!is.logical(L2.x.scale)) {
    stop(paste("The logical argument 'L2.x.scale' must be either TRUE or",
               " FALSE.", sep = ""))
  }

  # Check if folds is NULL
  if (!is.null(folds)) {
    # Check if folds is a character scalar
    if (!(is.character(folds) & length(folds) == 1)) {
      stop(paste("The argument 'folds', specifying the fold to which each",
                 " observation is to be allocated, must be a character scalar.",
                 sep = ""))
    }

    # Check if folds is in survey data
    if (!(folds %in% colnames(survey))) {
      stop(paste("Fold variable '", folds,
                 "' is not in your survey data.", sep = ""))
    }

    # Check if folds contains a sequence of integer numbers
    folds_var <- survey %>%
      dplyr::select_at(.vars = folds) %>%
      dplyr::pull() %>%
      unique() %>%
      sort()

    if (isFALSE(all(dplyr::near(folds_var, as.integer(folds_var))))) {
      stop(paste("Fold variable '", folds,
                 "' must contain integer numbers only.", sep = ""))
    }

    if (!(folds_var == 1:max(folds_var))) {
      stop(paste("Fold variable '", folds,
                 "' must contain a sequence of integers running from 1 to ",
                 max(folds_var), ".", sep = ""))
    }
  } else {
    # Check if ebma.size is NULL
    if (is.null(ebma.size)) {
      stop(paste("If argument 'folds' is NULL, then argument 'ebma.size' must",
                 " be specified.", sep = ""))
    } else {
      # Check if ebma.size is a proportion in the open unit interval
      if (!(is.numeric(ebma.size) & ebma.size > 0 & ebma.size < 1)) {
        stop(paste("The argument 'ebma.size', specifying the share of",
                   " respondents to be allocated to the EBMA fold, must take a",
                   " number in the open unit interval.", sep = ""))
      }
    }

    # Check if k.folds is NULL
    if (is.null(k.folds)) {
      stop(paste("If argument 'folds' is NULL, then argument 'k.folds' must",
                 " be specified.", sep = ""))
    } else {
      # Check if k.folds is an integer-valued scalar
      if (!(dplyr::near(k.folds, as.integer(k.folds)) &
            length(k.folds) == 1)) {
        stop(paste("The argument 'k.folds', specifying the number of folds to",
                   " be used in cross-validation, must be an integer-valued",
                   " scalar.",
                   sep = ""))
      } else {
        # Check if k.folds is less than or equal to the number of survey
        # respondents
        if (k.folds > nrow(survey)) {
          stop(paste("The argument 'k.folds', specifying the number of folds",
                     " to be used in cross-validation, cannot be larger than",
                     " the number of survey respondents, ", nrow(survey), ".",
                     sep = ""))
        }
      }
    }

    # Check if cv.sampling is NULL
    if (is.null(cv.sampling)) {
      stop(paste("If argument 'folds' is NULL, then argument 'cv.sampling'",
                 " must be specified.", sep = ""))
    } else {
      # Check if cv.sampling is either "individuals" or "L2 units"
      if (!cv.sampling %in% c("individuals", "L2 units")) {
        stop(paste("The argument 'cv.sampling', specifying the sampling method",
                   " to be used for cross-validation, must be either",
                   " 'individuals' or 'L2 units'.", sep = ""))
      }
    }
  }

  # Check if bin.proportion is NULL
  if (!is.null(bin.proportion)) {
    # Check if bin.proportion is a character scalar
    if (!(is.character(bin.proportion) & length(bin.proportion) == 1)) {
      stop(paste("The argument 'bin.proportion', specifying the variable that",
                 " indicates the proportion of ideal types in the census data,",
                 " must be a character scalar.", sep = ""))
    }

    # Check if bin.proportion is in census data
    if (!(bin.proportion %in% colnames(census))) {
      stop(paste("Variable '", bin.proportion,
                 "', indicating the proportion of ideal types, is not in your",
                 " census data.", sep = ""))
    }

    # Check if bin.proportion is a proportion
    bin_proportion_var <- census %>%
      dplyr::select_at(.vars = bin.proportion) %>%
      dplyr::pull() %>%
      unique()

    if (!(is.numeric(bin_proportion_var) &
          min(bin_proportion_var) >= 0 &
          max(bin_proportion_var) <= 1)) {
      stop(paste("Variable '", bin.proportion,
                 "', indicating the proportion of ideal types, can only take",
                 " values lying in the unit interval.", sep = ""))
    }
  } else {
    # Check if bin.size is NULL
    if (!is.null(bin.size)) {
      # Check if bin.size is a character scalar
      if (!(is.character(bin.size) & length(bin.size) == 1)) {
        stop(paste("The argument 'bin.size', specifying the variable that",
                   " indicates the bin size of ideal types in the census data,",
                   " must be a character scalar.", sep = ""))
      }

      # Check if bin.size is in census data
      if (!(bin.size %in% colnames(census))) {
        stop(paste("Variable '", bin.size,
                   "', indicating the bin size of ideal types, is not in your",
                   " census data.", sep = ""))
      }

      # Check if bin.size contains only non-negative numbers
      bin_size_var <- census %>%
        dplyr::select_at(.vars = bin.size) %>%
        dplyr::pull() %>%
        unique()

      if (!is.numeric(bin_size_var)) {
        stop(paste("Variable '", bin.size,
                   "', indicating the bin size of ideal types, must be numeric.",
                   sep = ""))
      }

      if (min(bin_size_var) < 0) {
        stop(paste("Variable '", bin.size,
                   "', indicating the bin size of ideal types, can only take",
                   " non-negative values.", sep = ""))
      }
    } else {
      stop(paste("Either argument 'bin.proportion' or argment 'bin.size' must",
                 " be specified to perform post-stratification.", sep = ""))
    }
  }

  # Check if survey data is provided as a data.frame
  if (is.null(survey)) {
    stop(paste("Argument 'survey' cannot be NULL. Please provide survey data.",
               sep = ""))
  } else {
    if (!is.data.frame(survey)) {
      stop(paste("The argument 'survey', specifying the survey data,",
                 " must be a data.frame.", sep = ""))
    }
  }

  # Check if census data is provided as a data.frame
  if (is.null(census)) {
    stop(paste("Argument 'census' cannot be NULL. Please provide census data.",
               sep = ""))
  } else {
    if (!is.data.frame(census)) {
      stop(paste("The argument 'census', specifying the census data,",
                 " must be a data.frame.", sep = ""))
    }
  }

  # Check if loss.unit is either "individuals" or "L2 units"
  if (!loss.unit %in% c("individuals", "L2 units")) {
    stop(paste("The argument 'loss.unit', specifying the level at which to",
               " evaluate prediction performance, must be either",
               " 'individuals' or 'L2 units'.", sep = ""))
  }

  # Check if loss.fun is either "MSE" or "MAE"
  if (!loss.fun %in% c("MSE", "MAE")) {
    stop(paste("The argument 'loss.fun', specifying the loss function used",
               " to measure prediction performance, must be either",
               " 'MSE' or 'MAE'.", sep = ""))
  }

  # Check if best.subset is logical
  if (is.logical(best.subset)) {
    # Check if best.subset is TRUE
    if (isTRUE(best.subset)) {
      # Check if best.subset.L2.x is NULL
      if (!is.null(best.subset.L2.x)) {
        # Check if best.subset.L2.x is a character vector
        if (!is.character(best.subset.L2.x)) {
          stop(paste("The argument 'best.subset.L2.x', specifying the context-level",
                     " variables to be used by the best subset classifier, must be",
                     " a character vector.", sep = ""))
        }

        # Check if best.subset.L2.x is in survey data
        if (!all(best.subset.L2.x %in% colnames(survey))) {
          stop(cat(paste("Context-level variable '",
                         best.subset.L2.x[which(!(best.subset.L2.x %in% colnames(survey)))],
                         "', specified in argument 'best.subset.L2.x' to be used by the",
                         " best subset classifier, is not in your survey data.", sep = ""),
                   sep = ""))
        }

        # Check if best.subset.L2.x is in census data
        if (!all(best.subset.L2.x %in% colnames(census))) {
          stop(cat(paste("Context-level variable '",
                         best.subset.L2.x[which(!(best.subset.L2.x %in% colnames(census)))],
                         "', specified in argument 'best.subset.L2.x' to be used by the",
                         " best subset classifier, is not in your census data.", sep = ""),
                   sep = ""))
        }
      }
    } else {
      # Check if best.subset.L2.x is NULL
      if (!is.null(best.subset.L2.x)) {
        warning(paste("The argument 'best.subset.L2.x', specifying the context-level",
                      " variables to be used by the best subset classifier, will be",
                      " ignored because 'best.subset' is set to FALSE.", sep = ""))
      }
    }
  } else {
    stop(paste("The logical argument 'best.subset', indicating whether the",
               " best subset classifier is to be used for predicting y,",
               " must be either TRUE or FALSE.", sep = ""))
  }

  # Check if lasso is logical
  if (is.logical(lasso)) {
    # Check if lasso is TRUE
    if (isTRUE(lasso)) {
      # Check if lasso.L2.x is NULL
      if (!is.null(lasso.L2.x)) {
        # Check if lasso.L2.x is a character vector
        if (!is.character(lasso.L2.x)) {
          stop(paste("The argument 'lasso.L2.x', specifying the context-level",
                     " variables to be used by the lasso classifier, must be",
                     " a character vector.", sep = ""))
        }

        # Check if lasso.L2.x is in survey data
        if (!all(lasso.L2.x %in% colnames(survey))) {
          stop(cat(paste("Context-level variable '",
                         lasso.L2.x[which(!(lasso.L2.x %in% colnames(survey)))],
                         "', specified in argument 'lasso.L2.x' to be used by the",
                         " lasso classifier, is not in your survey data.", sep = ""),
                   sep = ""))
        }

        # Check if lasso.L2.x is in census data
        if (!all(lasso.L2.x %in% colnames(census))) {
          stop(cat(paste("Context-level variable '",
                         lasso.L2.x[which(!(lasso.L2.x %in% colnames(census)))],
                         "', specified in argument 'lasso.L2.x' to be used by the",
                         " lasso classifier, is not in your census data.", sep = ""),
                   sep = ""))
        }
      }

      # Check if lasso.lambda is a list
      if (is.list(lasso.lambda)) {
        # Check if lasso.lambda is a list of two numeric vectors of equal size
        if (!(length(lasso.lambda) == 2 &
              all(sapply(lasso.lambda, function(x) {is.numeric(x)})) &
              length(unique(sapply(lasso.lambda, function(x) {length(x)}))) == 1)) {
          stop(paste("If provided as a list, argument 'lasso.lambda' must be a",
                     " list of two numeric vectors of equal size, with the",
                     " first vector containing the step sizes by which the",
                     " penalty parameter should increase and the second vector",
                     " containing the upper thresholds of the intervals to",
                     " which the step sizes apply.",
                     sep = ""))
        } else {
          # Check if step sizes do not exceed upper thresholds of intervals to
          # which they should be applied
          if (any(lasso.lambda[[1]] > lasso.lambda[[2]])) {
            stop(paste("If argument 'lasso.lambda' is specified as a list of",
                       " two vectors, the value in the first vector indicating",
                       " the step size cannot exceed the corresponding value",
                       " in the second vector indicating the upper threshold",
                       " of the interval to which the step size should apply.",
                       sep = ""))
          } else {
            # Check if lasso.lambda contains only non-negative values
            if (min(unlist(lasso.lambda)) < 0) {
              stop(paste("The argument 'lasso.lambda' can only take",
                         " non-negative values.", sep = ""))
            }
          }
        }
      } else {
        # Check if lasso.lambda is a numeric vector
        if (is.numeric(lasso.lambda)) {
          # Check if lasso.lambda contains only non-negative values
          if (min(lasso.lambda) < 0) {
            stop(paste("The argument 'lasso.lambda' can only take",
                       " non-negative values.", sep = ""))
          }
        } else {
          stop(paste("The argument 'lasso.lambda' must be either a numeric",
                     " vector of non-negative values or a list of two numeric",
                     " vectors of equal size, with the first vector containing",
                     " the step sizes by which the penalty parameter should",
                     " increase and the second vector containing the upper",
                     " thresholds of the intervals to which the step sizes apply.",
                     sep = ""))
        }
      }

      # Check if lasso.n.iter is NULL
      if (is.null(lasso.n.iter)) {

      }

      # Check if lasso.n.iter is an integer-valued scalar
      if (!(dplyr::near(lasso.n.iter, as.integer(lasso.n.iter)) &
            length(lasso.n.iter) == 1)) {
        stop(paste("The argument 'lasso.n.iter', specifying the number of folds to",
                   " be used in cross-validation, must be an integer-valued",
                   " scalar.",
                   sep = ""))
      }
    } else {
      # Check if lasso.L2.x is NULL
      if (!is.null(lasso.L2.x)) {
        warning(paste("The argument 'lasso.L2.x', specifying the context-level",
                      " variables to be used by the lasso classifier, will be",
                      " ignored because 'lasso' is set to FALSE.", sep = ""))
      }
    }
  } else {
    stop(paste("The logical argument 'lasso', indicating whether the lasso",
               " classifier is to be used for predicting y,",
               " must be either TRUE or FALSE.", sep = ""))
  }

  # Check if pca is logical
  if (is.logical(pca)) {
    # Check if pca is TRUE
    if (isTRUE(pca)) {
      # Check if pcs is NULL
      if (!is.null(pcs)) {
        # Check if pcs is a character vector
        if (!is.character(pcs)) {
          stop(paste("The argument 'pcs', specifying the principal components of",
                     " the context-level variables, must be a character vector.",
                     sep = ""))
        }

        # Check if pcs is in survey data
        if (!all(pcs %in% colnames(survey))) {
          stop(cat(paste("Principal component '",
                         pcs[which(!(pcs %in% colnames(survey)))],
                         "', specified in argument 'pcs', is not in your survey",
                         " data.\n", sep = ""), sep = ""))
        }

        # Check if pcs is in census data
        if (!all(pcs %in% colnames(census))) {
          stop(cat(paste("Principal component '",
                         pcs[which(!(pcs %in% colnames(census)))],
                         "', specified in argument 'pcs', is not in your census",
                         " data.\n", sep = ""), sep = ""))
        }
      } else{
        # Check if pcs are not specified but column names contain "PC" followed by at least one number
        if (is.null(pcs)){
          if (any(grepl(pattern = "PC[0-9]?", x = names(survey)))){
            stop(paste("Survey contains the column names: ",
                       paste(names(survey)[grepl(pattern = "PC[0-9]?", x = names(survey))], collapse = ", "),
                       ". These must be specified in the argument 'pcs' or removed from survey or renamed in survey.", sep = ""))
          }
          if (any(grepl(pattern = "PC[0-9]?", x = names(census)))){
            stop(paste("Census contains the column names: ",
                       paste(names(census)[grepl(pattern = "PC[0-9]?", x = names(census))], collapse = ", "),
                       ". These must be specified in the argument 'pcs' or removed from census or renamed in census.", sep = ""))
          }
        }
      }
    } else {
      # Check if pcs is NULL
      if (!is.null(pcs)) {
        warning(paste("The argument 'pcs', specifying the principal components",
                      " of the context-level variables, will be ignored because",
                      " 'pca' is set to FALSE.", sep = ""))
      }
    }
  } else {
    stop(paste("The logical argument 'pca', indicating whether the PCA",
               " classifier is to be used for predicting y,",
               " must be either TRUE or FALSE.", sep = ""))
  }

  # Check if gb is logical
  if (is.logical(gb)) {
    # Check if gb is TRUE
    if (isTRUE(gb)) {
      # Check if gb.L2.x is NULL
      if (!is.null(gb.L2.x)) {
        # Check if gb.L2.x is a character vector
        if (!is.character(gb.L2.x)) {
          stop(paste("The argument 'gb.L2.x', specifying the context-level",
                     " variables to be used by the GB classifier, must be",
                     " a character vector.", sep = ""))
        }

        # Check if gb.L2.x is in survey data
        if (!all(gb.L2.x %in% colnames(survey))) {
          stop(cat(paste("Context-level variable '",
                         gb.L2.x[which(!(gb.L2.x %in% colnames(survey)))],
                         "', specified in argument 'gb.L2.x' to be used by the GB",
                         " classifier, is not in your survey data.", sep = ""),
                   sep = ""))
        }

        # Check if gb.L2.x is in census data
        if (!all(gb.L2.x %in% colnames(census))) {
          stop(cat(paste("Context-level variable '",
                         gb.L2.x[which(!(gb.L2.x %in% colnames(census)))],
                         "', specified in argument 'gb.L2.x' to be used by the GB",
                         " classifier, is not in your census data.", sep = ""),
                   sep = ""))
        }
      }

      # Check if gb.L2.unit is logical
      if (!is.logical(gb.L2.unit)) {
        stop(paste("The logical argument 'gb.L2.unit', indicating whether",
                   " 'L2.unit' should be included in the GB classifier must be",
                   " either TRUE or FALSE.", sep = ""))
      }

      # Check if gb.L2.reg is logical
      if (!is.logical(gb.L2.reg)) {
        stop(paste("The logical argument 'gb.L2.reg', indicating whether",
                   " 'L2.reg' should be included in the GB classifier must be",
                   " either TRUE or FALSE.", sep = ""))
      }
    } else {
      # Check if gb.L2.x is NULL
      if (!is.null(gb.L2.x)) {
        warning(paste("The argument 'gb.L2.x', specifying the context-level",
                      " variables to be used by the GB classifier, will be",
                      " ignored because 'gb' is set to FALSE.", sep = ""))
      }

      # Check if gb.L2.unit has a value other than the default
      if (!isFALSE(gb.L2.unit)) {
        stop(paste("The argument 'gb.L2.unit', indicating whether 'L2.unit'",
                   " should be included in the GB classifier, will be",
                   " ignored because 'gb' is set to FALSE.", sep = ""))
      }

      # Check if gb.L2.reg has a value other than the default
      if (!isFALSE(gb.L2.reg)) {
        stop(paste("The argument 'gb.L2.reg', indicating whether 'L2.reg'",
                   " should be included in the GB classifier, will be",
                   " ignored because 'gb' is set to FALSE.", sep = ""))
      }
    }
  } else {
    stop(paste("The logical argument 'gb', indicating whether the GB",
               " classifier is to be used for predicting y,",
               " must be either TRUE or FALSE.", sep = ""))
  }

  # Check if svm is logical
  if (is.logical(svm)) {
    # Check if svm is TRUE
    if (isTRUE(svm)) {
      # Check if svm.L2.x is NULL
      if (!is.null(svm.L2.x)) {
        # Check if svm.L2.x is a character vector
        if (!is.character(svm.L2.x)) {
          stop(paste("The argument 'svm.L2.x', specifying the context-level",
                     " variables to be used by the SVM classifier, must be",
                     " a character vector.", sep = ""))
        }

        # Check if svm.L2.x is in survey data
        if (!all(svm.L2.x %in% colnames(survey))) {
          stop(cat(paste("Context-level variable '",
                         svm.L2.x[which(!(svm.L2.x %in% colnames(survey)))],
                         "', specified in argument 'svm.L2.x' to be used by the",
                         " SVM classifier, is not in your survey data.", sep = ""),
                   sep = ""))
        }

        # Check if svm.L2.x is in census data
        if (!all(svm.L2.x %in% colnames(census))) {
          stop(cat(paste("Context-level variable '",
                         svm.L2.x[which(!(svm.L2.x %in% colnames(census)))],
                         "', specified in argument 'svm.L2.x' to be used by the",
                         " SVM classifier, is not in your census data.", sep = ""),
                   sep = ""))
        }
      }
    } else {
      # Check if svm.L2.x is NULL
      if (!is.null(svm.L2.x)) {
        warning(paste("The argument 'svm.L2.x', specifying the context-level",
                      " variables to be used by the SVM classifier, will be",
                      " ignored because 'svm' is set to FALSE.", sep = ""))
      }
    }
  } else {
    stop(paste("The logical argument 'svm', indicating whether the SVM",
               " classifier is to be used for predicting y,",
               " must be either TRUE or FALSE.", sep = ""))
  }

  # Check if mrp is logical
  if (is.logical(mrp)) {
    # Check if mrp is TRUE
    if (isTRUE(mrp)) {
      # Check if mrp.L2.x is NULL
      if (!is.null(mrp.L2.x)) {
        # Check if mrp.L2.x is a character vector
        if (!is.character(mrp.L2.x)) {
          stop(paste("The argument 'mrp.L2.x', specifying the context-level",
                     " variables to be used by the standard MRP classifier, must",
                     " be a character vector.", sep = ""))
        }

        # Check if mrp.L2.x is in survey data
        if (all(mrp.L2.x != "empty")){
          if (!all(mrp.L2.x %in% colnames(survey))) {
            stop(cat(paste("Context-level variable '",
                           mrp.L2.x[which(!(mrp.L2.x %in% colnames(survey)))],
                           "', specified in argument 'mrp.L2.x' to be used by the",
                           " standard MRP classifier, is not in your survey data.",
                           sep = ""), sep = ""))
          }
        }

        # Check if mrp.L2.x is in census data
        if (all(mrp.L2.x != "empty")){
          if (!all(mrp.L2.x %in% colnames(census))) {
            stop(cat(paste("Context-level variable '",
                           mrp.L2.x[which(!(mrp.L2.x %in% colnames(census)))],
                           "', specified in argument 'mrp.L2.x' to be used by the",
                           " standard MRP classifier, is not in your census data.",
                           sep = ""), sep = ""))
          }
        }
      }
    } else {
      # Check if mrp.L2.x is NULL
      if (!is.null(mrp.L2.x)) {
        warning(paste("The argument 'mrp.L2.x', specifying the context-level",
                      " variables to be used by the standard MRP classifier,",
                      " will be ignored because 'mrp' is set to FALSE.", sep = ""))
      }
    }
  } else {
    stop(paste("The logical argument 'mrp', indicating whether the standard",
               " MRP classifier is to be used for predicting y,",
               " must be either TRUE or FALSE.", sep = ""))
  }

  # Check if forward.select is logical
  if (!is.logical(forward.select)) {
    stop(paste("The logical argument 'forward.select', indicating whether to",
               " use forward selection instead of best subset selection,",
               " must be either TRUE or FALSE.", sep = ""))
  }
}


################################################################################
#                    Function to create EBMA hold-out fold                     #
################################################################################

#' Generates data fold to be used for EBMA tuning
#'
#' #' \code{ebma_folding()} generates a data fold that will not be used in
#' classifier tuning. It is data that is needed to determine the optimal
#' tolerance for EBMA.
#'
#' @param data The full survey data. A tibble.
#' @param L2.unit Geographic unit. A character scalar containing the column name
#'   of the geographic unit in \code{survey} and \code{census} at which outcomes
#'   should be aggregated.
#' @param ebma.size EBMA fold size. A number in the open unit interval
#'   indicating the proportion of respondents to be allocated to the EBMA fold.
#'   Default is \eqn{1/3}.

ebma_folding <- function(data, L2.unit, ebma.size) {
  # Add row number to data frame
  data <- data %>%
    dplyr::mutate(index = dplyr::row_number())

  # Split data by geographic unit into a list of data frames
  data_list <- data %>%
    dplyr::group_split(.data[[L2.unit]])

  # Sample one respondent per geographic unit
  one_per_unit <- lapply(data_list, function(x) {
    sample(x$index, size = 1, replace = FALSE)
  }) %>%
    unlist()

  # Create EBMA hold-out fold
  ebma_fold <- data %>%
    dplyr::filter(index %in% one_per_unit)

  data <- data %>%
    dplyr::filter(!(index %in% one_per_unit))

  remainder <- sample(data$index, size = ebma.size - length(one_per_unit),
                      replace = FALSE)

  ebma_remainder <- data %>%
    dplyr::filter(index %in% remainder)

  ebma_fold <- ebma_fold %>%
    dplyr::bind_rows(ebma_remainder)

  # Extract EBMA hold-out fold from survey sample
  cv_data <- data %>%
    dplyr::filter(!(index %in% ebma_fold$index))

  # Remove index variable
  ebma_fold <- ebma_fold %>%
    dplyr::select(-index)

  cv_data <- cv_data %>%
    dplyr::select(-index)

  # Function output
  out <- list(ebma_fold = ebma_fold,
              cv_data = cv_data)
  return(out)
}


################################################################################
#                         Function to create CV folds                          #
################################################################################

#' Generates folds for cross-validation
#'
#' \code{cv_folding} creates folds used in classifier training within the survey
#' data.
#'
#' @param data The survey data; must be a tibble.
#' @param L2.unit The column name of the factor variable identifying the
#'   context-level unit
#' @param k.folds An integer value indicating the number of folds to be
#'   generated.
#' @param cv.sampling Cross-validation sampling method. A character-valued
#'   scalar indicating whether cross-validation folds should be created by
#'   sampling individual respondents (\code{individuals}) or geographic units
#'   (\code{L2 units}). Default is \code{L2 units}. \emph{Note:} ignored if
#'   \code{folds} is provided, but must be specified otherwise.

cv_folding <- function(data, L2.unit, k.folds,
                       cv.sampling = c("individuals", "L2 units")) {

  if (cv.sampling == "individuals") {
    # Add row number to data frame
    data <- data %>%
      dplyr::mutate(index = row_number())

    # Randomize indices of individuals
    indices <- sample(data$index, size = length(data$index), replace = FALSE)

    # Define number of units per fold
    no_floor <- floor(length(indices) / k.folds)
    no_remaining <- length(indices) - no_floor * k.folds

    no_fold <- rep(no_floor, times = k.folds)

    if (no_remaining > 0) {
      no_fold[1:no_remaining] <- no_fold[1:no_remaining] + 1
    }

    # Split indices into folds
    fold_indices <- split(indices, rep(1:k.folds, times = no_fold))

    # Partition data according to indices
    out <- lapply(fold_indices, function(x) {
      data %>%
        dplyr::filter(index %in% x) %>%
        dplyr::select(-index)
    })
  } else {
    # Extract indices of geographic units
    indices <- data[[L2.unit]] %>%
      unique()

    # Randomize order of indices
    indices <- sample(indices, size = length(indices), replace = FALSE)

    # Define number of units per fold
    no_floor <- floor(length(indices) / k.folds)
    no_remaining <- length(indices) - no_floor * k.folds

    no_fold <- rep(no_floor, times = k.folds)

    if (no_remaining > 0) {
      no_fold[1:no_remaining] <- no_fold[1:no_remaining] + 1
    }

    # Split indices into folds
    fold_indices <- split(indices, rep(1:k.folds, times = no_fold))

    # Partition data according to indices
    out <- lapply(fold_indices, function(x) {
      data %>%
        dplyr::filter(.data[[L2.unit]] %in% x)
    })
  }

  # Function output
  return(out)
}


################################################################################
#           Function to create model list for best subset classifier           #
################################################################################

#' A list of models for the best subset selection.
#'
#' \code{model_list()} generates an exhaustive list of lme4 model formulas from
#' the individual level and context level variables as well as geographic unit
#' variables to be iterated over in best subset selection.
#'
#' @param y Outcome variable. A character vector containing the column names of
#'   the outcome variable.
#' @param L1.x Individual-level covariates. A character vector containing the
#'   column names of the individual-level variables in \code{survey} and
#'   \code{census} used to predict outcome \code{y}. Note that geographic unit
#'   is specified in argument \code{L2.unit}.
#' @param L2.x Context-level covariates. A character vector containing the
#'   column names of the context-level variables in \code{survey} and
#'   \code{census} used to predict outcome \code{y}.
#' @param L2.unit Geographic unit. A character scalar containing the column name
#'   of the geographic unit in \code{survey} and \code{census} at which outcomes
#'   should be aggregated.
#' @param L2.reg Geographic region. A character scalar containing the column
#'   name of the geographic region in \code{survey} and \code{census} by which
#'   geographic units are grouped (\code{L2.unit} must be nested within
#'   \code{L2.reg}). Default is \code{NULL}.

model_list <- function(y, L1.x, L2.x, L2.unit, L2.reg = NULL) {
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

  # Empty model
  empty_model <- list(as.formula(paste(y, " ~ ", all_re, sep = "")))

  # Remaining models
  L2_list <- lapply(seq_along(L2.x), function(x) {combn(L2.x, x)})
  L2_list <- lapply(L2_list, function(x) {
    apply(x, 2, function(c) {
      as.formula(paste(y, " ~ ", paste(c, collapse = " + "), " + ", all_re, sep = ""))
    })
  }) %>%
    unlist()

  # Combine models in list
  out <- c(empty_model, L2_list)

  # Function output
  return(out)
}


################################################################################
#               Function to create model list for PCA classifier               #
################################################################################

#' A list of models for the best subset selection with PCA.
#'
#' \code{model_list_pca()} generates an exhaustive list of lme4 model formulas
#' from the individual level and context level principal components as well as
#' geographic unit variables to be iterated over in best subset selection with
#' principal components.
#'
#' @param y Outcome variable. A character vector containing the column names of
#'   the outcome variable.
#' @param L1.x Individual-level covariates. A character vector containing the
#'   column names of the individual-level variables in \code{survey} and
#'   \code{census} used to predict outcome \code{y}. Note that geographic unit
#'   is specified in argument \code{L2.unit}.
#' @param L2.x Context-level covariates. A character vector containing the
#'   column names of the context-level variables in \code{survey} and
#'   \code{census} used to predict outcome \code{y}.
#' @param L2.unit Geographic unit. A character scalar containing the column name
#'   of the geographic unit in \code{survey} and \code{census} at which outcomes
#'   should be aggregated.
#' @param L2.reg Geographic region. A character scalar containing the column
#'   name of the geographic region in \code{survey} and \code{census} by which
#'   geographic units are grouped (\code{L2.unit} must be nested within
#'   \code{L2.reg}). Default is \code{NULL}.

model_list_pca <- function(y, L1.x, L2.x, L2.unit, L2.reg = NULL) {
  # Individual-level random effects
  L1_re <- paste(paste("(1 | ", L1.x, ")", sep = ""), collapse = " + ")

  # Geographic unit or Geographic unit-Geographic region random effects
  if (is.null(L2.reg)) {
    L2_re <- paste("(1 | ", L2.unit, ")", sep = "")
  } else {
    L2_re <- paste(paste("(1 | ", L2.reg, "/", L2.unit, ")", sep = ""),
                   collapse = " + ")
  }

  # Combine all random effects
  all_re <- paste(c(L1_re, L2_re), collapse = " + ")

  # Empty model
  empty_model <- list(as.formula(paste(y, " ~ ", all_re, sep = "")))

  # Remaining models
  L2_list <- lapply(seq_along(L2.x), function(x) {L2.x[1:x]})
  L2_list <- lapply(L2_list, function(x) {
    as.formula(paste(y, " ~ ", paste(x, collapse = " + "), " + ", all_re, sep = ""))
  })

  # Combine models in list
  out <- c(empty_model, L2_list)

  # Function output
  return(out)
}


################################################################################
#                           Prediction loss function                           #
################################################################################

#' Estimates loss value
#'
#' \code{loss_function()} estimates the loss based on a loss function.
#'
#' @param pred Predictions of outcome. A numeric vector of outcome predictions.
#' @param data.valid Test data set. A tibble of data that was not used for
#'   prediction.
#' @param loss.unit Loss function unit. A character-valued scalar indicating
#'   whether performance loss should be evaluated at the level of individual
#'   respondents (\code{individuals}) or geographic units (\code{L2 units}).
#'   Default is \code{individuals}.
#' @param loss.fun Loss function. A character-valued scalar indicating whether
#'   prediction loss should be measured by the mean squared error (\code{MSE})
#'   or the mean absolute error (\code{MAE}). Default is \code{MSE}.
#' @param y Outcome variable. A character vector containing the column names of
#'   the outcome variable.
#' @param L2.unit Geographic unit. A character scalar containing the column name
#'   of the geographic unit in \code{survey} and \code{census} at which outcomes
#'   should be aggregated.

loss_function <- function(pred, data.valid,
                          loss.unit = c("individuals", "L2 units"),
                          loss.fun = c("MSE", "MAE"),
                          y, L2.unit) {
  if (loss.unit == "individuals" & loss.fun == "MSE") {
    out <- mean((data.valid[[y]] - pred)^2)
  } else if (loss.unit == "individuals" & loss.fun == "MAE") {
    out <- mean(abs(data.valid[[y]] - pred))
  } else if (loss.unit == "L2 units" & loss.fun == "MSE") {
    data.valid <- data.valid %>%
      dplyr::mutate(pred = pred)

    out <- data.valid %>%
      dplyr::group_by_at(L2.unit) %>%
      dplyr::summarise_at(.vars = c(y, "pred"), mean) %>%
      dplyr::mutate(sqe = (.data[[y]] - pred)^2) %>%
      dplyr::pull(sqe)

    out <- mean(out)
  } else {
    data.valid <- data.valid %>%
      dplyr::mutate(pred = pred)

    out <- data.valid %>%
      dplyr::group_by_at(L2.unit) %>%
      dplyr::summarise_at(.vars = c(y, "pred"), mean) %>%
      dplyr::mutate(ae = abs(.data[[y]] - pred)) %>%
      dplyr::pull(ae)

    out <- mean(out)
  }

  # Function output
  return(out)
}


################################################################################
#                   Suppress cat in external package                           #
################################################################################

#' Suppress cat in external package
#'
#' \code{quiet()} suppresses cat output.
#'
#' @param x Input. It can be any kind.
#' @source The author is Hadley Wickham. We found the function here:
#' \url{http://r.789695.n4.nabble.com/Suppressing-output-e-g-from-cat-td859876.html}.

quiet <- function(x) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}

################################################################################
#                   Register Cores for multicore                               #
################################################################################

#' Register cores for multicore computing
#'
#' \code{multicore()} registers cores for parallel processing.
#'
#' @param cores Number of cores to be used. An integer. Default is \code{1}.
#' @param type Whether to start or end parallel processing. A character string.
#'   The possible values are \code{open}, \code{close}.
#' @param cl The registered cluster. Default is \code{NULL}

multicore <- function(cores = 1, type, cl = NULL) {

  # Start parallel processing
  if (type == "open"){
    # register clusters for windows
    if( Sys.info()["sysname"] == "Windows" ){
      cl <- parallel::makeCluster(cores)
      doParallel::registerDoParallel(cl)
    } else {
      cl <- parallel::makeForkCluster(cores)
      doParallel::registerDoParallel(cl)
    }
    return(cl)
  }

  # Stop parallel processing
  if (type == "close"){
    parallel::stopCluster(cl)
  }
}
