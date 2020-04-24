################################################################################
#            Function to check if arguments to auto_MrP() are valid            #
################################################################################

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
  if (!loss.fun %in% c("individuals", "L2 units")) {
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
        if (!all(mrp.L2.x %in% colnames(survey))) {
          stop(cat(paste("Context-level variable '",
                         mrp.L2.x[which(!(mrp.L2.x %in% colnames(survey)))],
                         "', specified in argument 'mrp.L2.x' to be used by the",
                         " standard MRP classifier, is not in your survey data.",
                         sep = ""), sep = ""))
        }

        # Check if mrp.L2.x is in census data
        if (!all(mrp.L2.x %in% colnames(census))) {
          stop(cat(paste("Context-level variable '",
                         mrp.L2.x[which(!(mrp.L2.x %in% colnames(census)))],
                         "', specified in argument 'mrp.L2.x' to be used by the",
                         " standard MRP classifier, is not in your census data.",
                         sep = ""), sep = ""))
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
