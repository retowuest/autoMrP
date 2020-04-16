################################################################################
#            Function to check if arguments to auto_MrP() are valid            #
################################################################################

error_checks <- function(y, L1.x) {
  # Check if y is a character scalar
  if (!(is.character(y) & length(y) == 1)) {
    stop(paste("Outcome '", y,
               "' must be a character scalar.", sep = ""))
  }

  # Check if y is in survey data
  if (!(y %in% colnames(survey))) {
    stop(paste("Outcome '", y,
               "' is not in your survey data.", sep = ""))
  }

  # Check if L1.x is a character vector


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
                          loss.measure = c("MSE", "MAE"),
                          y, L2.unit) {
  if (loss.unit == "individuals" & loss.measure == "MSE") {
    out <- mean((data.valid[[y]] - pred)^2)
  } else if (loss.unit == "individuals" & loss.measure == "MAE") {
    out <- mean(abs(data.valid[[y]] - pred))
  } else if (loss.unit == "L2 units" & loss.measure == "MSE") {
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
