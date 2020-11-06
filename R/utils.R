################################################################################
#            Function to check if arguments to auto_MrP() are valid            #
################################################################################

#' Catches user input errors
#'
#' \code{error_checks()} checks for incorrect data entry in \code{autoMrP()}
#' call.
#'
#' @inheritParams auto_MrP

error_checks <- function(y, L1.x, L2.x, L2.unit, L2.reg, L2.x.scale, pcs,
                         folds, bin.proportion, bin.size, survey, census,
                         ebma.size, k.folds, cv.sampling, loss.unit, loss.fun,
                         best.subset, lasso, pca, gb, svm, mrp, forward.select,
                         best.subset.L2.x, lasso.L2.x, gb.L2.x, svm.L2.x,
                         mrp.L2.x, gb.L2.unit, gb.L2.reg, lasso.lambda,
                         lasso.n.iter, uncertainty, boot.iter, seed) {


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
      if (!(is.numeric(ebma.size) & ebma.size >= 0 & ebma.size < 1)) {
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
  if (!all(loss.unit %in% c("individuals", "L2 units"))) {
    stop(paste("The argument 'loss.unit', specifying the level at which to",
               " evaluate prediction performance, must be either",
               " 'individuals' or 'L2 units'.", sep = ""))
  }

  # Check if loss.fun is either "MSE" or "MAE"
  if (!all(loss.fun %in% c("MSE", "MAE", "cross-entropy", "f1", "msfe"))) {
    stop(paste("The argument 'loss.fun', specifying the loss function used",
               " to measure prediction performance, must be either",
               " 'MSE', 'MAE', 'cross-entropy', 'f1', or 'msfe'.", sep = ""))
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

      # Check if is provided but not a numeric vector
      if (!is.null(lasso.lambda)){
        if (!is.numeric(lasso.lambda)){
          stop("lasso.lambda must be 'NULL' or a non-negative numeric vector.")
        } else{
          # Check if lasso.lambda contains non-negative values
          if (!all(lasso.lambda > 0)){
            stop("lasso.lambda must not contain negative values")
          }
        }
      }

      # Check if lasso.n.iter is NULL
      if (!is.null(lasso.n.iter)) {
        if (!(dplyr::near(lasso.n.iter, as.integer(lasso.n.iter)) &
              length(lasso.n.iter) == 1)) {
          stop("lasso.n.iter specifies the Lasso grid size. It must be a non-negative integer valued scalar.")
        }
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
      # if (!isFALSE(gb.L2.unit)) {
      #   stop(paste("The argument 'gb.L2.unit', indicating whether 'L2.unit'",
      #              " should be included in the GB classifier, will be",
      #              " ignored because 'gb' is set to FALSE.", sep = ""))
      # }

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

  # Check if boot.iter corresponds to uncertainty, i.e. NULL if uncertainty = FALSE
  if (!uncertainty){
    if(!is.null(boot.iter)) {
      warning("boot.iter is ignored unless uncertainty = TRUE.")
    }
  }

  # Check if supplied seed is integer
  if (!is.null(seed)){
    if (isFALSE(dplyr::near(seed, as.integer(seed)))) {
      stop("Seed must be either NULL or an integer-valued scalar.")
    }
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

#' Estimates loss value.
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
                          loss.fun = c("MSE", "MAE", "cross-entropy"),
                          y, L2.unit) {

  ## Loss functions
  # MSE
  mse <- mean_squared_error(
    pred = pred, data.valid = data.valid,
    y = y, L2.unit = L2.unit)

  # MAE
  mae <- mean_absolute_error(
    pred = pred, data.valid = data.valid,
    y = y, L2.unit = L2.unit)

  # binary cross-entropy
  bce <- binary_cross_entropy(
    pred = pred, data.valid = data.valid,
    y = y, L2.unit = L2.unit)

  # f1 score
  f1 <- f1_score(
    pred = pred, data.valid = data.valid,
    L2.unit = L2.unit, y = y)

  # mean squared false error
  msfe <- mean_squared_false_error(
    pred = pred, data.valid = data.valid,
    y = y, L2.unit = L2.unit)

  # Combine loss functions
  score <- mse %>%
    dplyr::bind_rows(mae) %>%
    dplyr::bind_rows(bce) %>%
    dplyr::bind_rows(f1) %>%
    dplyr::bind_rows(msfe)

  # Filter score table by loss function and loss unit
  score <- score %>%
    dplyr::filter(measure %in% loss.fun) %>%
    dplyr::filter(level %in% loss.unit) %>%
    dplyr::group_by(measure) %>%
    dplyr::summarise(value = mean(value), .groups = "drop" )

  # Function output
  return(score)
}


###########################################################################
# mean squared error/ brier score -----------------------------------------
###########################################################################

#' Estimates the mean squared prediction error.
#'
#' \code{mean_squared_error()} estimates the mean squared error for the desired
#' loss unit.
#' @inheritParams loss_function

mean_squared_error <- function(pred, data.valid, y, L2.unit){

  # outcome
  out <- dplyr::tibble(
    measure = rep("MSE", 2),
    value = rep(NA, 2),
    level = c( "individuals", "L2 units")
  )

  # mse values
  values <- rep(NA, 2)

  # loss unit = "individual"
  values[1] <- mean((data.valid[[y]] - pred)^2)

  # loss unit = "L2 units"
  l2 <- data.valid %>%
    dplyr::mutate(pred = pred) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(sqe = (.data[[y]] - pred)^2 ) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(.dots = list(L2.unit)) %>%
    dplyr::summarise(mse = mean(sqe), .groups = "drop")

  values[2] <- mean(dplyr::pull(.data = l2, var = mse))

  out <- dplyr::mutate(out, value = values)

  return(out)
}


###########################################################################
# mean absolute error -----------------------------------------------------
###########################################################################

#' Estimates the mean absolute prediction error.
#'
#' \code{mean_absolute_error()} estimates the mean absolute error for the
#' desired loss unit.
#' @inheritParams loss_function

mean_absolute_error <- function(pred, data.valid, y, L2.unit){

  # outcome
  out <- dplyr::tibble(
    measure = rep("MAE", 2),
    value = rep(NA, 2),
    level = c( "individuals", "L2 units"))

  # mae values
  values <- rep(NA, 2)

  # loss unit = "individual"
  values[1] <- mean(abs(data.valid[[y]] - pred))

  # loss unit = "L2 units"
  l2 <- data.valid %>%
    dplyr::mutate(pred = pred) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(ae = abs(.data[[y]] - pred)) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(.dots = list(L2.unit)) %>%
    dplyr::summarise(mae = mean(ae), .groups = "drop")

  values[2] <- mean(dplyr::pull(.data = l2, var = mae))

  out <- dplyr::mutate(out, value = values)

  return(out = out)
}


###########################################################################
# binary cross-entropy ----------------------------------------------------
###########################################################################

#' Estimates the inverse binary cross-entropy, i.e. 0 is the best score and 1
#' the worst.
#'
#' \code{binary_cross_entropy()} estimates the inverse binary cross-entropy on
#' the individual and state-level.
#' @inheritParams loss_function

binary_cross_entropy <- function(pred, data.valid,
                                 loss.unit = c("individuals", "L2 units"),
                                 y, L2.unit){

  # outcome
  out <- dplyr::tibble(
    measure = rep("cross-entropy", 2),
    value = rep(NA, 2),
    level = c( "individuals", "L2 units")
  )

  # cross-entropy values
  values <- rep(NA, 2)

  # loss unit = "individual"
  values[1] <- (mean( data.valid[[y]] * log(pred) + (1 - data.valid[[y]]) * log(1 - pred)))*-1

  # loss unit = "L2 units"
  l2 <- data.valid %>%
    dplyr::mutate(pred = pred) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(ce = .data[[y]] * log(pred) + (1 - .data[[y]]) * log(1 - pred) ) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(.dots = list(L2.unit)) %>%
    dplyr::summarise(bce = mean(ce), .groups = "drop")

  values[2] <- mean(dplyr::pull(.data = l2, var = bce)) *-1

  out <- dplyr::mutate(out, value = values)
  return(out)
}


###########################################################################
# F1 score ----------------------------------------------------------------
###########################################################################
#' Estimates the inverse f1 score, i.e. 0 is the best score and 1 the worst.
#'
#' \code{f1_score()} estimates the inverse f1 scores on the individual and state
#' levels.
#' @inheritParams loss_function


f1_score <- function(pred, data.valid, y, L2.unit){

  ## individual level

  # true positives
  tp_ind <- data.valid %>%
    dplyr::mutate(pval = ifelse(test = pred > 0.5, yes = 1, no = 0)) %>%
    dplyr::select( !! rlang::sym(y), pval ) %>%
    dplyr::filter(pval == 1 & !! rlang::sym(y) == 1) %>%
    dplyr::summarise(tp = sum(pval)) %>%
    dplyr::pull(var = tp)

  # false positives
  fp_ind <- data.valid %>%
    dplyr::mutate(pval = ifelse(test = pred > 0.5, yes = 1, no = 0)) %>%
    dplyr::select( !! rlang::sym(y), pval ) %>%
    dplyr::filter(pval == 1 & !! rlang::sym(y) == 0 ) %>%
    dplyr::summarise(fp = sum(pval)) %>%
    dplyr::pull(var = fp)

  # false negatives
  fn_ind <- data.valid %>%
    dplyr::mutate(pval = ifelse(test = pred > 0.5, yes = 1, no = 0)) %>%
    dplyr::select( !! rlang::sym(y), pval ) %>%
    dplyr::filter(pval == 0 & !! rlang::sym(y) == 1 ) %>%
    dplyr::summarise(fn = sum(!! rlang::sym(y))) %>%
    dplyr::pull(var = fn)

  # f1 score
  f1 <- tp_ind / (tp_ind + 0.5 * (fp_ind + fn_ind) )

  # state-level f1 score
  state_out <- data.valid %>%
    # predicted values
    dplyr::mutate(pval = ifelse(test = pred > 0.5, yes = 1, no = 0)) %>%
    # select L2.unit, y, and predicted values
    dplyr::select( !! rlang::sym(L2.unit), !! rlang::sym(y), pval ) %>%
    # group by L2.unit
    dplyr::group_by( !! rlang::sym(L2.unit) ) %>%
    # nest data
    tidyr::nest() %>%
    # new column with state-level f1 values
    dplyr::mutate(
      f1 = purrr::map(data, function(x){
        # true positives
        tp <- x %>%
          dplyr::select( !! rlang::sym(y), pval ) %>%
          dplyr::filter(pval == 1 & !! rlang::sym(y) == 1) %>%
          dplyr::summarise(tp = sum(pval)) %>%
          dplyr::pull(var = tp)
        # false positives
        fp <- x %>%
          dplyr::select( !! rlang::sym(y), pval ) %>%
          dplyr::filter(pval == 1 & !! rlang::sym(y) == 0 ) %>%
          dplyr::summarise(fp = sum(pval)) %>%
          dplyr::pull(var = fp)
        # false negatives
        fn <- x %>%
          dplyr::select( !! rlang::sym(y), pval ) %>%
          dplyr::filter(pval == 0 & !! rlang::sym(y) == 1 ) %>%
          dplyr::summarise(fn = sum(!! rlang::sym(y))) %>%
          dplyr::pull(var = fn)
        # f1 score
        f1 <- tp / (tp + 0.5 * (fp + fn) ) })) %>%
    # unnest f1 values
    tidyr::unnest(f1) %>%
    dplyr::select( !! rlang::sym(L2.unit), f1 ) %>%
    dplyr::ungroup() %>%
    dplyr::summarise(f1 = mean(f1, na.rm = TRUE), .groups = "drop")

  # return
  out <- dplyr::tibble(
    measure = c("f1", "f1"),
    value = c(1 - f1, 1 - dplyr::pull(.data = state_out, var = f1)),
    level = c("individuals", "L2 units"))

  return(out)

}


###########################################################################
# Mean squared false error-------------------------------------------------
###########################################################################
#' Estimates the mean squared false error.
#'
#' \code{msfe()} estimates the inverse f1 scores on the individual and state
#' levels.
#' @inheritParams loss_function


mean_squared_false_error <- function(pred, data.valid, y, L2.unit){

  ## individual level
  msfe_l1 <- data.valid %>%
    dplyr::mutate(pval = ifelse(test = pred > 0.5, yes = 1, no = 0)) %>%
    dplyr::select( !! rlang::sym(y), pval ) %>%
    dplyr::group_by( !! rlang::sym(y) ) %>%
    dplyr::mutate(err = (!! rlang::sym(y) - pval) ) %>%
    dplyr::summarise(err_rates = mean(err), .groups = "drop") %>%
    dplyr::mutate(err_rates = err_rates^2) %>%
    dplyr::summarise( msfe = sum(err_rates)) %>%
    dplyr::pull(var = msfe)

  ## group level
  msfe_l2 <- data.valid %>%
    dplyr::mutate(pval = ifelse(test = pred > 0.5, yes = 1, no = 0)) %>%
    dplyr::select( !! rlang::sym(L2.unit), !! rlang::sym(y), pval ) %>%
    dplyr::group_by( !! rlang::sym(L2.unit) ) %>%
    tidyr::nest() %>%
    dplyr::mutate(msfe = purrr::map(data, function(x){
      msfe <- x %>%
        dplyr::group_by( !! rlang::sym(y) ) %>%
        dplyr::mutate(err = (!! rlang::sym(y) - pval) ) %>%
        dplyr::summarise(err_rates = mean(err), .groups = "drop") %>%
        dplyr::mutate(err_rates = err_rates^2) %>%
        dplyr::summarise( msfe = sum(err_rates)) %>%
        dplyr::pull(var = msfe)
    })) %>%
    tidyr::unnest(msfe) %>%
    dplyr::ungroup() %>%
    dplyr::summarise(msfe = mean(msfe), .groups = "drop") %>%
    dplyr::pull(var = msfe)

  # return
  out <- dplyr::tibble(
    measure = c("msfe", "msfe"),
    value = c(msfe_l1, msfe_l2),
    level = c("individuals", "L2 units"))

  return(out)

}

###########################################################################
# Loss score ranking ------------------------------------------------------
###########################################################################

#' Ranks tuning parameters according to loss functions
#'
#' \code{loss_score_ranking()} ranks tuning parameters according to the scores
#' received in multiple loss functions.
#'
#' @inheritParams loss_function
#' @param score A data set containing loss function names, the loss function
#'   values, and the tuning parameter values.

loss_score_ranking <- function(score, loss.fun){

  # tuning parameter names
  params <- names(score)[!names(score) %in% c("measure", "value")]

  ranking <- lapply(loss.fun, function(x){
    score %>%
      dplyr::filter(measure == x) %>%
      dplyr::arrange(value) %>%
      dplyr::mutate(rank = dplyr::row_number())
  })

  ranking <- dplyr::bind_rows(ranking) %>%
    dplyr::group_by(.dots = params) %>%
    dplyr::summarise(rank = sum(rank), .groups = "drop") %>%
    dplyr::arrange(rank)

  return(ranking)

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
#' \url{https://r.789695.n4.nabble.com/Suppressing-output-e-g-from-cat-td859876.html}.

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


################################################################################
#                 Predict function for glmmLasso                               #
################################################################################

#' Predicts on newdata from glmmLasso objects
#'
#' \code{glmmLasso()} predicts on newdata objects from a glmmLasso object.
#'
#' @inheritParams auto_MrP
#' @param m A \code{glmmLasso()} object.

predict_glmmLasso <- function(census, m, L1.x, lasso.L2.x, L2.unit, L2.reg) {

  # Fixed effects
  fixed_effects <- as.matrix(cbind(1, as.data.frame(census)[, lasso.L2.x])) %*% cbind(m$coefficients)

  # Individual-level random effects
  ind_ranef <- as.data.frame(census)[, L1.x]
  ind_ranef[] <- base::Map(paste, names(ind_ranef), ind_ranef, sep = '')
  ind_ranef <- cbind(apply(ind_ranef, 1, function(x){
    sum(m$ranef[which(names(m$ranef) %in% x)])
  }))

  # State random effects
  state_ranef <- cbind(paste(L2.unit, as.character(as.data.frame(census)[, L2.unit]), sep = ""))
  state_ranef <- cbind(apply(state_ranef, 1, function(x){
    m$ranef[names(m$ranef) == x]
  }))

  # Region random effect
  if(!is.null(L2.reg)){
    region_ranef <- cbind(paste(L2.reg, as.character(as.data.frame(census)[, L2.reg]), sep = ""))
    region_ranef <- cbind(apply(region_ranef, 1, function(x){
      m$ranef[names(m$ranef) == x]
    }))
  }

  # Predictions
  if(!is.null(L2.reg)){
    lasso_preds <- cbind(fixed_effects, ind_ranef, state_ranef, region_ranef)
  } else{
    lasso_preds <- cbind(fixed_effects, ind_ranef, state_ranef)
  }
  lasso_preds <- apply(lasso_preds, 1, sum)
  lasso_preds <- stats::pnorm(lasso_preds)

  return(lasso_preds)
}


################################################################################
#                 Plot method for autoMrP                                      #
################################################################################

#' A plot method for autoMrP objects. Plots unit-level preference estiamtes.
#'
#' \code{plot.autoMrP()} plots unit-level preference estimates and error bars.
#'
#' @param x An \code{autoMrP()} object.
#' @param algorithm The algorithm/classifier fo which preference estimates are
#'   desired. A character-valued scalar indicating either \code{ebma} or the
#'   classifier to be used. Allowed choices are: "ebma", "best_subset", "lasso",
#'   "pca", "gb", "svm", and "mrp". Default is \code{ebma}.
#' @param ci.lvl The level of the confidence intervals. A proportion. Default is
#'   \code{0.95}. Confidence intervals are based on bootstrapped estimates and
#'   will not be printed if bootstrapping was not carried out.
#' @param ... Additional arguments affecting the summary produced.
#' @export
#' @export plot.autoMrP

plot.autoMrP <- function(x, algorithm = "ebma", ci.lvl = 0.95, ...){

  # L2.unit identifier
  L2.unit <- names(x$classifiers)[1]


  # Error if requested algorithm was not fitted
  if(! algorithm %in% names(x$classifiers) & algorithm != "ebma" ){
    stop('The ', algorithm, ' classifier was not run. Re-run autoMrP() with the requested algorithm. Allowed choices are: "ebma", "best_subset", "lasso", "pca", "gb", "svm", and "mrp".')
  }

  # plot classifier if EBMA was not estimated
  if( "EBMA step skipped (only 1 classifier run)" %in% x$ebma ) {
    algorithm <- names(x$classifiers)[-1]
  }

  # plot data
  if(algorithm == "ebma"){
    plot_data <- x$ebma %>%
      dplyr::group_by(.dots = list(L2.unit)) %>%
      dplyr::summarise(median = stats::median(ebma, na.rm = TRUE),
                       lb = stats::quantile(x = ebma, p = (1 - ci.lvl)*.5, na.rm = TRUE),
                       ub = stats::quantile(x = ebma, p = ci.lvl + (1 - ci.lvl)*.5, na.rm = TRUE),
                       .groups = "drop") %>%
      dplyr::arrange(median) %>%
      dplyr::mutate(rank = dplyr::row_number()) %>%
      dplyr::mutate(rank = as.factor(rank))
  } else{
    plot_data <- x$classifiers %>%
      dplyr::group_by(.dots = list(L2.unit)) %>%
      dplyr::select(all_of(L2.unit), contains(algorithm)) %>%
      dplyr::summarise_all(.funs = list(median = ~ stats::quantile(x = ., probs = 0.5, na.rm = TRUE),
                                        lb = ~ stats::quantile(x = ., probs = (1 - ci.lvl) *.5, na.rm = TRUE),
                                        ub = ~ stats::quantile(x = ., probs = ci.lvl + (1 - ci.lvl) *.5, na.rm = TRUE))) %>%
      dplyr::arrange(median) %>%
      dplyr::mutate(rank = dplyr::row_number()) %>%
      dplyr::mutate(rank = as.factor(rank))
  }

  # y axis tick labels
  ylabs <- as.character(dplyr::pull(.data = plot_data, var = L2.unit))

  # plot (with/ without error bars)
  if(all(plot_data$median == plot_data$lb)){
    ggplot2::ggplot(data = plot_data, mapping = ggplot2::aes_string(x = "median", y = "rank", label = L2.unit)) +
      ggplot2::geom_point() +
      ggplot2::labs(x = "Estimates") +
      ggplot2::scale_y_discrete(breaks = rank, labels = ylabs, name = "States")
  } else{
    ggplot2::ggplot(data = plot_data, mapping = ggplot2::aes_string(x = "median", y = "rank", label = L2.unit)) +
      ggplot2::geom_point() +
      ggplot2::labs(x = "Estimates") +
      ggplot2::scale_y_discrete(breaks = rank, labels = ylabs, name = "States") +
      ggplot2::geom_errorbarh(mapping = ggplot2::aes(xmin = lb, xmax = ub))
  }
}


################################################################################
#                 Summary method for autoMrP                                   #
################################################################################

#' A summary method for autoMrP objects.
#'
#' \code{summary.autoMrP()} ...
#'
#' @param object An \code{autoMrP()} object for which a summary is desired.
#' @param ci.lvl The level of the confidence intervals. A proportion. Default is
#'   \code{0.95}. Confidence intervals are based on bootstrapped estimates and
#'   will not be printed if bootstrapping was not carried out.
#' @param digits The number of digits to be displayed. An integer scalar.
#'   Default is \code{4}.
#' @param format The table format. A character string passed to
#'   \code{\link[knitr]{kable}}. Default is \code{simple}.
#' @param classifiers Summarize a single classifier. A character string. Must be
#'   one of \code{best_subset}, \code{lasso}, \code{pca}, \code{gb}, \code{svm},
#'   or \code{mrp}. Default is \code{NULL}.
#' @param n Number of rows to be printed. An integer scalar. Default is
#'   \code{10}.
#' @param ... Additional arguments affecting the summary produced.
#' @export
#' @export summary.autoMrP

summary.autoMrP <- function(object, ci.lvl = 0.95, digits = 4, format = "simple",
                            classifiers = NULL, n = 10, ...){

  # weights
  if ( all(c("autoMrP", "weights") %in% class(object)) ){

    # error message if weights summary called without running multiple classifiers
    if (any(object == "EBMA step skipped (only 1 classifier run)")){
      stop("Weights are not reported if the EBMA step was skipped. Re-run autoMrP with multiple classifiers.")
    }

    # weights vector to tibble
    if( is.null(dim(object)) ){
      object <- dplyr::tibble(!!!object)
    }

    # summary statistics
    s_data <- object %>%
      tidyr::pivot_longer(
        cols = dplyr::everything(),
        names_to = "method",
        values_to = "estimates") %>%
      dplyr::group_by(method) %>%
      dplyr::summarise(
        min = base::min(estimates, na.rm = TRUE),
        quart1 = stats::quantile(x = estimates, probs = 0.25, na.rm = TRUE),
        median = stats::median(estimates, na.rm = TRUE),
        mean = base::mean(estimates, na.rm = TRUE),
        quart3 = stats::quantile(x = estimates, probs = 0.75, na.rm = TRUE),
        max = base::max(estimates, na.rm = TRUE),
        .groups = "drop") %>%
      dplyr::arrange(dplyr::desc(median))

    # weights with uncertainty
    if ( all(s_data$median != s_data$min) ){
      n <- ifelse(n <= nrow(s_data), yes = n, no = nrow(s_data) )
      cat( paste("\n", "# EBMA classifier weights:"), sep = "")
      # output table
      output_table(
        object = s_data[1:n, ],
        col.names = c(
          "Classifier",
          "Min.",
          "1st Qu.",
          "Median",
          "Mean",
          "3rd Qu.",
          "Max"),
        format = format,
        digits = digits)
      if (n < nrow(s_data)) cat( paste("... with", nrow(s_data)-n, " more rows", "\n", "\n"), sep = "")
    } else{
      s_data <- dplyr::select(.data = s_data, method, median)
      n <- ifelse(n <= nrow(s_data), yes = n, no = nrow(s_data) )
      cat( paste("\n", "# EBMA classifier weights:"), sep = "")
      output_table(
        object = s_data[1:n, ],
        col.names = c(
          "Classifier",
          "Weight"),
        format = format,
        digits = digits)
      if (n < nrow(s_data)) cat( paste("... with", nrow(s_data)-n, " more rows", "\n", "\n"), sep = "")
    }
  }

  # ensemble summary
  else if ( all(c("autoMrP", "ensemble") %in% class(object)) ) {

    # unit identifier
    L2.unit <- names(object)[1]

    # summary statistics
    s_data <- object %>%
      dplyr::group_by(.dots = list(L2.unit)) %>%
      dplyr::summarise(
        min = base::min(ebma, na.rm = TRUE),
        lb = stats::quantile(x = ebma, probs = (1 - ci.lvl)*.5, na.rm = TRUE),
        median = stats::quantile(x = ebma, probs = .5, na.rm = TRUE),
        ub = stats::quantile(x = ebma, probs = ci.lvl + (1 - ci.lvl)*.5, na.rm = TRUE),
        max = base::max(ebma, na.rm = TRUE),
        .groups = "drop"
      )

    # with or without uncertainty
    if ( all(s_data$median != s_data$lb) ){
      cat( paste("\n", "# EBMA estimates:"), sep = "")
      # output table
      output_table(
        object = s_data[1:n, ],
        col.names = c(
          L2.unit,
          "Min.",
          "Lower bound",
          "Median",
          "Upper bound",
          "Max"),
        format = format,
        digits = digits)
      if (n < nrow(s_data)) cat( paste("... with", nrow(s_data)-n, " more rows", "\n", "\n"), sep = "")

    } else{
      s_data <- dplyr::select(.data = s_data, one_of(L2.unit), median)
      n <- ifelse(n <= nrow(s_data), yes = n, no = nrow(s_data) )
      cat( paste("\n", "# EBMA estimates:"), sep = "")
      output_table(
        object = s_data[1:n, ],
        col.names = c(L2.unit, "Estimate"),
        format = format,
        digits = digits)
      if (n < nrow(s_data)) cat( paste("... with", nrow(s_data)-n, " more rows", "\n", "\n"), sep = "")
    }
  }

  # classifier summary
  else if ( all(c("autoMrP", "classifiers") %in% class(object)) ){

    # unit identifier
    L2.unit <- names(object)[1]

    # multiple classifiers
    if (base::is.null(classifiers)){

      # point estimates for all classifiers
      s_data <- object %>%
        dplyr::group_by(.dots = list(L2.unit)) %>%
        dplyr::summarise_all(.funs = median )

      # output table
      ests <- paste(names(object)[-1], collapse = ", ")
      n <- ifelse(n <= nrow(s_data), yes = n, no = nrow(s_data) )
      cat( paste("\n", "# estimates of classifiers: ", ests), sep = "")
      output_table(object = s_data[1:n, ],
                   col.names = names(s_data),
                   format = format,
                   digits = digits)
      if (n < nrow(s_data)) cat( paste("... with", nrow(s_data)-n, " more rows", "\n", "\n"), sep = "")
    } else{

      # summary statistics
      s_data <- object %>%
        dplyr::select(dplyr::one_of(L2.unit,classifiers)) %>%
        dplyr::group_by(.dots = list(L2.unit)) %>%
        dplyr::summarise_all(.funs = list(
          min = ~ base::min(x = ., na.rm = TRUE),
          lb = ~ stats::quantile(x = ., probs = (1 - ci.lvl)*.5, na.rm = TRUE),
          median = ~ stats::median(x = ., na.rm = TRUE),
          ub = ~ stats::quantile(x = ., probs = ci.lvl + (1 - ci.lvl)*.5, na.rm = TRUE),
          max = ~ base::max(x = ., na.rm = TRUE)
          ))

      # with or without uncertainty
      if( all(s_data$median != s_data$lb) ){
        n <- ifelse(n <= nrow(s_data), yes = n, no = nrow(s_data) )
        cat( paste("\n", "# estimates of", classifiers, "classifier"), sep = "")
        output_table(
          object = s_data[1:n, ],
          col.names = c(
            L2.unit,
            "Min.",
            "Lower bound",
            "Median",
            "Upper bound",
            "Max"),
          format = format,
          digits = digits)
        if (n < nrow(s_data)) cat( paste("... with", nrow(s_data)-n, " more rows", "\n", "\n"), sep = "")
      } else{
        s_data <- dplyr::select(.data = s_data, dplyr::one_of(L2.unit), "median")
        n <- ifelse(n <= nrow(s_data), yes = n, no = nrow(s_data) )
        cat( paste("\n", "# estimates of", classifiers, "classifier"), sep = "")
        output_table(
          object = s_data[1:n, ],
          col.names = c(L2.unit, "Estimate"),
          format = format,
          digits = digits)
        if (n < nrow(s_data)) cat( paste("... with", nrow(s_data)-n, " more rows", "\n", "\n"), sep = "")
      }
    }
  }

  # autoMrP list object
  else if ( all(c("autoMrP", "list") %in% class(object)) ){

    # unit identifier
    L2.unit <- names(object$classifiers)[1]

    # if EBMA was run
    if( !"EBMA step skipped (only 1 classifier run)" %in% object$ebma ){

      # Summarize EBMA or classifier specified in classifiers argument
      if( is.null(classifiers)) {
        s_data <- object$ebma
      } else{
        # check whether classifier was fitted
        if( !classifiers %in% names(object$classifiers) ){
          stop( classifiers, " was not fitted. Summary available for: ", paste(names(object$classifiers)[-1], collapse = ", "))
        }
        s_data <- object$classifiers %>%
          dplyr::select( one_of(L2.unit, classifiers) )
      }

      # summary statistics
      s_data <- s_data %>%
        dplyr::group_by(.dots = list(L2.unit)) %>%
        dplyr::summarise_all(.funs = list(
          min = ~ base::min(x = ., na.rm = TRUE),
          lb = ~ stats::quantile(x = ., probs = (1 - ci.lvl)*.5, na.rm = TRUE),
          median = ~ stats::median(x = ., na.rm = TRUE ),
          ub = ~ stats::quantile(x = ., probs = ci.lvl + (1 - ci.lvl)*.5, na.rm = TRUE),
          max = ~ base::max(x = ., na.rm = TRUE)
        ))

      # with or without uncertainty
      if( all(s_data$median != s_data$lb) ){
        n <- ifelse(n <= nrow(s_data), yes = n, no = nrow(s_data) )
        if( is.null(classifiers) ){
          cat( paste("\n", "# EBMA estimates:"), sep = "")
        } else{
          cat( paste("\n", "# ", classifiers, " estimates", sep = ""))
        }
        output_table(
          object = s_data[1:n, ],
          col.names = c( L2.unit, "Min.", "Lower bound", "Median", "Upper bound", "Max"),
        format = format,
        digits = digits)
        if (n < nrow(s_data)) cat( paste("... with", nrow(s_data)-n, " more rows", "\n", "\n"), sep = "")
      } else{
        s_data <- dplyr::select(.data = s_data, dplyr::one_of(L2.unit), median)
        n <- ifelse(n <= nrow(s_data), yes = n, no = nrow(s_data) )
        cat( paste("\n", "# EBMA estimates:"), sep = "")
        output_table(object = s_data[1:n, ], col.names = c(L2.unit, "Median"), format = format, digits = digits)
        if (n < nrow(s_data)) cat( paste("... with", nrow(s_data)-n, " more rows", "\n", "\n"), sep = "")
      }
    } else{

      # Summarize all classifiers or classifier specified in classifiers argument
      # or MrP model
      if( is.null(classifiers) ) {
        s_data <- object$classifiers
      } else{
        s_data <- object$classifiers %>%
          dplyr::select( one_of(L2.unit, classifiers) )
      }

      # summary statistics
      s_data <- s_data %>%
        dplyr::group_by(.dots = list(L2.unit)) %>%
        dplyr::summarise_all(.funs = list(
          min = ~ base::min(x = ., na.rm = TRUE),
          lb = ~ stats::quantile(x = ., probs = (1 - ci.lvl)*.5, na.rm = TRUE),
          median = ~ stats::median(x = ., na.rm = TRUE),
          ub = ~ stats::quantile(x = ., probs = ci.lvl + (1 - ci.lvl)*.5, na.rm = TRUE),
          max = ~ base::max(x = ., na.rm = TRUE)
        ))

      # with or without uncertainty
      comparison <- s_data %>%
        dplyr::select(one_of(grep(pattern = "median", x = names(s_data), value = "TRUE")[1],
                             grep(pattern = "lb", x = names(s_data), value = "TRUE")[1]))

        # summarize one classifier
        if ( sum(grepl(pattern = "best_subset|pca|lasso|gb|svm|mrp", x = names(s_data))) < 4 ){

          # without uncertainty
          if( all(comparison[,1] != comparison[,2]) ){

            n <- ifelse(n <= nrow(s_data), yes = n, no = nrow(s_data) )
            cat( paste("\n", "# ", names(object$classifiers)[2]," estimates:"), sep = "")
            output_table(
              object = s_data[1:n, ],
              col.names = c(
                L2.unit,
                "Min.",
                "Lower bound",
                "Median",
                "Upper bound",
                "Max"),
              format = format,
              digits = digits)
            if (n < nrow(s_data)) cat( paste("... with", nrow(s_data)-n, " more rows", "\n", "\n"), sep = "")
          } else{
            n <- ifelse(n <= nrow(s_data), yes = n, no = nrow(s_data) )
            cat( paste("\n", "# estimates of: ", paste(names(object$classifiers)[-1], collapse = ", ")), sep = "")
            s_data <- s_data %>%
              dplyr::select(one_of( L2.unit), contains("median"))
            output_table(
              object = s_data[1:n, ],
              col.names = names(s_data),
              format = format,
              digits = digits)
            if (n < nrow(s_data)) cat( paste("... with", nrow(s_data)-n, " more rows", "\n", "\n"), sep = "")
          }

      } else {
       # drop uncertainty columns
        if ( ncol(s_data) < 5 ){
          s_data <- dplyr::select(.data = s_data, dplyr::one_of(L2.unit), median)
          n <- ifelse(n <= nrow(s_data), yes = n, no = nrow(s_data) )
          cat( paste("\n", "# ", names(object$classifiers)[2]," estimates:"), sep = "")
          output_table(object = s_data[1:n, ], col.names = c(L2.unit, "Estimate"), format = format, digits = digits)
          if (n < nrow(s_data)) cat( paste("... with", nrow(s_data)-n, " more rows", "\n", "\n"), sep = "")
        } else{
          n <- ifelse(n <= nrow(s_data), yes = n, no = nrow(s_data) )
          s_data <- s_data %>%
            dplyr::select(one_of( L2.unit), contains("median"))
          output_table(
            object = s_data[1:n, ],
            col.names = names(s_data),
            format = format,
            digits = digits)
          if (n < nrow(s_data)) cat( paste("... with", nrow(s_data)-n, " more rows", "\n", "\n"), sep = "")
        }
      }
    }
  }
}

################################################################################
#                 Output table for summary                                     #
################################################################################

#' A table for the summary  function
#'
#' \code{output_table()} ...
#'
#' @inheritParams summary.autoMrP
#' @param col.names The column names of the table. A

output_table <- function(object, col.names, format, digits){

  # output table
  print( knitr::kable(x = object,
                      col.names = col.names,
                      format = format,
                      digits = digits))

}


################################################################################
#                 Equal spacing on the log scale                               #
################################################################################

#' @param min The minimum value of the sequence. A positive numeric scalar (min
#'   > 0).
#' @param max The maximum value of the sequence. a positive numeric scalar (max
#'   > 0).
#' @param n The length of the sequence. An integer valued scalar.

# Sequence that is equally spaced on the log scale
log_spaced <- function(min, max, n){
  return(base::exp( base::seq(from = base::log(min), to = base::log(max), length.out = n)))
}
