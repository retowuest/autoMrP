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

    # 1) non-negative least squares stack
    nnls_stack <- as.matrix(ebma_out$classifiers[, -1]) %*%
      stack_out$stack_weights$stack_nnls

    # generate stack predictions
    stacked_preds <- tibble::tibble(
      !!rlang::sym(L2.unit) := dplyr::pull(
        .data = ebma_out$classifiers, var = 1
      ),
      stack_nnls = as.numeric(nnls_stack)
    )

    # 2) optim with constraints stack
    optim_stack <- as.matrix(ebma_out$classifiers[, -1]) %*%
      stack_out$stack_weights$stack_optim

    # add optim stack to stack predictions
    stacked_preds <- stacked_preds %>%
      dplyr::mutate(stack_optim = as.numeric(optim_stack))

    # 3) quadratic programming stack
    qp_stack <- as.matrix(ebma_out$classifiers[, -1]) %*%
      stack_out$stack_weights$stack_qp

    # add qp stack to stack predictions
    stacked_preds <- stacked_preds %>%
      dplyr::mutate(stack_qp = as.numeric(qp_stack))

    # 4) ornstein stack
    ornstein_stack <- as.matrix(ebma_out$classifiers[, -1]) %*%
      stack_out$stack_weights$stack_ornstein

    # add ornstein stack to stack predictions
    stacked_preds <- stacked_preds %>%
      dplyr::mutate(stack_ornstein = as.numeric(ornstein_stack))

    # 5) stack of stacks
    stack_of_stacks <- as.matrix(stacked_preds[, -1]) %*%
      stack_out$stack_weights$stack_of_stacks

    # add stack of stacks to stack predictions
    stacked_preds <- stacked_preds %>%
      dplyr::mutate(stack_of_stacks = as.numeric(stack_of_stacks))

    # 6) stack of stacks with ebma
    stack_of_stacks_ebma <- as.matrix(
      cbind(stacked_preds[, "stack_of_stacks"], ebma_out$ebma[, "ebma"])
    ) %*%
      stack_out$stack_weights$stack_of_stacks_ebma

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
          dplyr::select(-id, -y, -dplyr::all_of(L2.unit))
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