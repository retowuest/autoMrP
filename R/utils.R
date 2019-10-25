# Function to create EBMA hold-out fold
ebma_folding <- function(data, geo.unit, n.ebma = NULL) {
  # Add row number to data frame
  data <- data %>%
    dplyr::mutate(index = row_number())

  # Split data by geographic unit into a list of data frames
  data_list <- data %>%
    dplyr::group_split(.data[[geo.unit]])

  # Sample one respondent per geographic unit
  one_per_unit <- lapply(data_list, function(x) sample(x$index, size = 1)) %>%
    unlist()

  # Create EBMA hold-out fold
  ebma_fold <- data %>%
    dplyr::filter(index %in% one_per_unit)

  data <- data %>%
    dplyr::filter(!(index %in% one_per_unit))

  remainder <- sample(data$index, size = n.ebma - length(one_per_unit))

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


# Function to create CV folds
cv_folding <- function(data, k.folds,
                       cv.sampling = c("respondents", "units")) {

  if (cv.sampling == "respondents") {
    # Add row number to data frame
    data <- data %>%
      dplyr::mutate(index = row_number())

    # Create folds
    fold_indices <- caret::createFolds(data$index, k = k.folds,
                                       list = TRUE, returnTrain = FALSE)

    out <- lapply(fold_indices, function(x) data %>%
                    dplyr::filter(index %in% x) %>%
                    dplyr::select(-index))
  } else {
    # Create folds
    fold_indices <- data[[geo.unit]] %>%
      unique() %>%
      caret::createFolds(k = k.folds, list = TRUE, returnTrain = FALSE)

    out <- lapply(fold_indices, function(x) data %>%
                    dplyr::filter(.data[[geo.unit]] %in% x))
  }

  # Function output
  return(out)
}


# Function to create model list for best subset classifier
model_list <- function(y, L1.x, L2.x) {
  # Individual-level random effects
  L1_re <- paste(paste("(1 | ", L1.x, ")", sep = ""), collapse = " + ")

  # Empty model
  all_models <- list(as.formula(paste(y, " ~ ", L1_re, sep = "")))

  # Remaining models
  L2_list <- lapply(seq_along(L2.x), function(x) combn(L2.x, x))
  L2_list <- lapply(L2_list, function(x) apply(x, 2, function(c)
    as.formula(paste(y, " ~ ", paste(c, collapse = " + "), " + ", L1_re, sep = "")))) %>%
    unlist()

  # Combine models in list
  out <- c(all_models, L2_list)

  # Function output
  return(out)
}


# Loss function
loss_function <- function(pred, data.valid,
                          unit = c("individual", "geo.unit"),
                          measure = c("mse", "mae"),
                          y, geo.unit) {
  if (unit == "individual" & measure == "mse") {
    out <- mean((data.valid[[y]] - pred)^2)
  } else if (unit == "individual" & measure == "mae") {
    out <- mean(abs(data.valid[[y]] - pred))
  } else if (unit == "geo.unit" & measure == "mse") {
    data.valid <- data.valid %>%
      dplyr::mutate(pred = pred)

    out <- data.valid %>%
      dplyr::group_by_at(geo.unit) %>%
      dplyr::summarise_at(.vars = c(y, "pred"), mean) %>%
      dplyr::mutate(sqe = (.data[[y]] - pred)^2) %>%
      dplyr::pull(sqe)

    out <- mean(out)
  } else {
    data.valid <- data.valid %>%
      dplyr::mutate(pred = pred)

    out <- data.valid %>%
      dplyr::group_by_at(geo.unit) %>%
      dplyr::summarise_at(.vars = c(y, "pred"), mean) %>%
      dplyr::mutate(ae = abs(.data[[y]] - pred)) %>%
      dplyr::pull(ae)

    out <- mean(out)
  }

  # Function output
  return(out)
}
