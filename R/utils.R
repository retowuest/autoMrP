# Function to create EBMA hold-out fold
ebma_folding <- function(data, geo.unit, n.ebma = NULL) {
  # Error checks
  if(is.null(n.ebma)) {
    n.ebma <- round(data / 4, digits = 0)
  } else if(!(is.numeric(n.ebma) & n.ebma == round(n.ebma, digits = 0))) {
      stop("n.ebma must be an integer number.")
    }

  # Add row number to data frame
  data <- data %>%
    dplyr::mutate(index = row_number())

  # Split data by geographic unit into a list of data frames
  data_list <- data %>%
    dplyr::group_split(get(geo.unit))

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

  # Remove index
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

