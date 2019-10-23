# Function to create EBMA hold-out fold
ebma_folding <- function(data, geo.unit) {
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
    dplyr::filter(index %in% one_per_unit) %>%
    dplyr::select(-index)

  # Extract EBMA hold-out fold from survey sample
  cv_data <- data %>%
    dplyr::filter(!(index %in% one_per_unit)) %>%
    dplyr::select(-index)

  # Function output
  out <- list(ebma_fold = ebma_fold,
              cv_data = cv_data)
  return(out)
}

