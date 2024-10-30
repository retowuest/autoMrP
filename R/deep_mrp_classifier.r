#' Deep MrP classifier
#'
#' \code{deep_mrp_classifier} applies Deep MrP implemented in the \pkg{vglmer}
#' package to a data set.
#'
#' @inheritParams auto_MrP
#' @param form Model formula. A two-sided linear formula describing
#'   the model to be fit, with the outcome on the LHS and the covariates
#'   separated by + operators on the RHS.
#' @param data Data. A data.frame containing the data used to train the model.
#' @return A Deep MrP model. A \code{\link[vglmer]{vglmer}} object.

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
      formula = as.formula(form),
      data = data,
      family = family
    )
  } else {
    out <- suppressMessages(suppressWarnings(
      vglmer::vglmer(
        formula = as.formula(form),
        data = data,
        family = family
      )
    ))
  }
  return(out)
}