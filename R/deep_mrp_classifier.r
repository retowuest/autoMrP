deep_mrp_classifier <- function(form, data, verbose) {

  # run vglmer model
  if (verbose) {
    out <- vglmer::vglmer(
      formula = form %>% as.formula(),
      data = data,
      family = "binomial")
  } else {
    out <- suppressMessages(suppressWarnings(
      vglmer::vglmer(
        formula = form %>% as.formula(),
        data = data,
        family = "binomial")))
  }
  return(out)
}