apply_stack_weights <- function(ebma_out, stack_out, L2.unit, preds_all, y) {

  # initial binding of globals
  ebma_preds <- NULL

  # check whether stacking weights were calculated
  if (all(is.null(stack_out$stack_preds))) {

    # generate an object for individual level predicitions
    individual_level_predictions <- preds_all %>%
      dplyr::rename(!!rlang::sym(y) := y)

    # output object
    ebma_out <- list(
      ebma = ebma_out$ebma,
      classifiers = ebma_out$classifiers,
      weights = ebma_out$weights,
      stacking = "Stacking step skipped (only 1 classifier run)",
      stacking_weights = "Stacking step skipped (only 1 classifier run)",
      individual_level_predictions = individual_level_predictions
    )

  } else {

    # 1) glm stack
    glm_stack <- cbind(1, as.matrix(ebma_out$classifiers[, -1])) %*%
      stack_out$stack_weights$stack_glm

    # sigmoid transformation
    glm_stack <- 1 / (1 + exp(-glm_stack))

    # generate stack predictions
    stacked_preds <- tibble::tibble(
      !!rlang::sym(L2.unit) := dplyr::pull(
        .data = ebma_out$classifiers, var = 1
      ),
      stack_glm = as.numeric(glm_stack)
    )

    # 2) lme4 stack
    lme4_stack <- lapply(seq_len(nrow(ebma_out$classifiers)), function(.x) {

      # current post-stratified predictions
      c_pred <- ebma_out$classifiers %>%
        dplyr::slice(.x)

      # current state
      c_state <- dplyr::pull(c_pred, !!rlang::sym(L2.unit))

      # corresponding lme4 weights
      c_weights <- stack_out$stack_weights$stack_lme4 %>%
        dplyr::filter(!!rlang::sym(L2.unit) == c_state) %>%
        dplyr::select(-!!rlang::sym(L2.unit)) %>%
        as.matrix()

      # drop c_pred columns that are missing from c_weights
      c_pred <- c_pred[, colnames(c_pred) %in% colnames(c_weights)]

      # predictions without state and added intercept
      c_pred <- cbind(1, as.matrix(c_pred))

      # multiply predictions with weights
      c_pred <- c_pred %*% t(c_weights)

      # sigmoid transformation
      c_pred <- 1 / (1 + exp(-c_pred))

      # add state
      lme4_out <- tibble::tibble(
        !!rlang::sym(L2.unit) := c_state,
        stack_lme4 = as.numeric(c_pred)
      )

      return(lme4_out)

    }) %>%
      dplyr::bind_rows()

    # combine stack predictions
    stacked_preds <- stacked_preds %>%
      dplyr::left_join(y = lme4_stack, by = L2.unit)

    # 3) nlm stack
    nlm_stack <- as.matrix(ebma_out$classifiers[, -1]) %*%
      stack_out$stack_weights$stack_nlm

    # sigmoid transformation
    nlm_stack <- 1 / (1 + exp(-nlm_stack))

    # add nlm stack to stack predictions
    stacked_preds <- stacked_preds %>%
      dplyr::mutate(stack_nlm = as.numeric(nlm_stack))

    # 4) ornstein stack
    ornstein_stack <- as.matrix(ebma_out$classifiers[, -1]) %*%
      stack_out$stack_weights$stack_ornstein

    # add ornstein stack to stack predictions
    stacked_preds <- stacked_preds %>%
      dplyr::mutate(stack_ornstein = as.numeric(ornstein_stack))

    # 5) stack of stacks
    stack_of_stacks <- cbind(1, as.matrix(stacked_preds[, -1])) %*%
      stack_out$stack_weights$stack_of_stacks

    # sigmoid transformation
    stack_of_stacks <- 1 / (1 + exp(-stack_of_stacks))

    # add stack of stacks to stack predictions
    stacked_preds <- stacked_preds %>%
      dplyr::mutate(stack_of_stacks = as.numeric(stack_of_stacks))

    # 6) stack of stacks with ebma
    stack_of_stacks_ebma <- as.matrix(
      cbind(1, stacked_preds[, "stack_of_stacks"], ebma_out$ebma[, "ebma"])
    ) %*%
      stack_out$stack_weights$stack_of_stacks_with_ebma

    # sigmoid transformation
    stack_of_stacks_ebma <- 1 / (1 + exp(-stack_of_stacks_ebma))

    # add stack of stacks with ebma to stack predictions
    stacked_preds <- stacked_preds %>%
      dplyr::mutate(stack_of_stacks_ebma = as.numeric(stack_of_stacks_ebma))

    # generate an object for individual level predicitions
    individual_level_predictions <- preds_all %>%
      dplyr::mutate(
        ebma = ebma_out$individual_level_predictions %>%
          dplyr::pull(var = ebma_preds)
      ) %>%
      dplyr::bind_cols(
        stack_out$stack_preds %>%
          dplyr::select(-1)
      ) %>%
      dplyr::rename(!!rlang::sym(y) := y)

    # combine everything
    ebma_out <- list(
      ebma = ebma_out$ebma,
      classifiers = ebma_out$classifiers,
      weights = ebma_out$weights,
      stacking = stacked_preds,
      stacking_weights = stack_out$stack_weights,
      individual_level_predictions = individual_level_predictions
    )
  }

  return(ebma_out)
}