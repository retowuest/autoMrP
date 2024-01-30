deep_mrp_classifier <- function(y, form, data, verbose) {

  # Determine type of dependent variable
  if (
    data %>%
      dplyr::pull(!!y) %>%
      unique() %>%
      length() == 2
  ) {
    family <- "binomial"
  } else {
    family <- "linear"
  }

  # run vglmer model
  if (verbose) {
    out <- vglmer::vglmer(
      formula = form %>% as.formula(),
      data = data,
      family = family
    )
  } else {
    out <- suppressMessages(suppressWarnings(
      vglmer::vglmer(
        formula = form %>% as.formula(),
        data = data,
        family = family
      )
    ))
  }

  return(out)
}