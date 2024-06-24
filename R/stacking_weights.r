stacking_weights <- function(preds, ebma_out, L2.unit) {

  # initial binding of globals
  stack_lme4 <- stack_nlm <- y <- stack_of_stacks <- ebma_preds <- NULL

  # check whether more than one classifier was tuned
  if (ncol(preds) <= 3) {
    message("Skipping stacking, only one classifier was tuned.")
    # return the stacked predictions and weights
    return(list(stack_preds = NULL, stack_weights = NULL))
  } else {

    # true labels
    true_labels <- preds$y

    # grouping variable
    group_var <- preds %>% dplyr::pull(!!rlang::sym(L2.unit))

    # matrix of predictions
    stack_predictions <- as.matrix(preds[, -c(1, 2)])

    # is the dpendent variable binary or not?
    binary_dv <- if (length(unique(true_labels)) == 2) {
      TRUE
    } else {
      FALSE
    }

    # stacking for binary dependent variables
    if (binary_dv) {

      #--------------------------------------------------
      # regression for binary response
      #--------------------------------------------------
      # meta formula
      meta_formula <- as.formula(
        paste("y ~ ", paste(colnames(stack_predictions), collapse = " + "))
      )

      # logistic regression model
      logit_m <- glm(meta_formula, data = preds, family = binomial)

      # coefficients of the logistic model
      logit_coefs <- coef(logit_m)

      # weights for base models
      stacking_weights <- logit_coefs  # Exclude the intercept
      stacking_weights[is.na(stacking_weights)] <- 0

      # predicions using logistic regression weights
      final_prediction <- cbind(1L, stack_predictions) %*% stacking_weights

      # apply sigmoid
      final_prediction <- 1 / (1 + exp(-final_prediction))

      # stacked predictions output
      stack_out <- tibble::tibble(
        stack_glm = as.numeric(final_prediction)
      )

      # stacking weights container
      stack_weights <- list(
        stack_glm = stacking_weights
      )

      #----------------------------------------------------------------
      # random intercepts model for binary response
      #----------------------------------------------------------------

      # data for glmer
      stack_data <- tibble::tibble(
        true_labels = true_labels,
        !!rlang::sym(L2.unit) := group_var
      ) %>%
        dplyr::bind_cols(stack_predictions)

      # random intercepts formula
      glmer_form <- sprintf(
        "true_labels ~ %s + (1 | %s)",
        paste0(colnames(stack_predictions), collapse = " + "), L2.unit
      )

      # random intercepts model
      glmer_m <- lme4::glmer(
        formula = glmer_form, data = stack_data, family = binomial
      )

      # Predict and apply sigmoid for binary classification
      stack_out <- stack_out %>%
        dplyr::mutate(
          stack_lme4 = stats::predict(
            glmer_m, newdata = as.data.frame(stack_data), re.form = NULL,
          )
        ) %>%
        dplyr::mutate(
          stack_lme4 = 1 / (1 + exp(-stack_lme4))
        )

      # extract fixed effects (common weights)
      fixed_effects <- lme4::fixef(glmer_m)

      # extract random effects (group-specific deviations)
      random_effects <- lme4::ranef(glmer_m)[[L2.unit]]

      # compute group-specific intercepts
      group_intercepts <- random_effects + fixed_effects["(Intercept)"]

      # Group-specific weights
      group_specific_weights <- tibble::tibble(
        !!rlang::sym(L2.unit) := levels(stack_data[[L2.unit]])
      )
      group_specific_weights <- lapply(
        seq_len(nrow(group_specific_weights)), function(x) {
          # current group
          c_group <- dplyr::slice(group_specific_weights, x)
          # group-specific weights
          c_weights <- c(
            intercept = group_intercepts[
              c_group[[L2.unit]], "(Intercept)"
            ], fixed_effects[-1]
          )
          # combine group and weights
          c_group <- dplyr::bind_cols(c_group, t(c_weights))
          return(c_group)
        }
      ) %>%
        dplyr::bind_rows()

      # add group-specific weights to the stack_weights list
      stack_weights$stack_lme4 <- group_specific_weights

      #----------------------------------------------------------------
      # opitmization for binary response
      #----------------------------------------------------------------

      # objective function to minimize (log loss)
      objective <- function(weights) {
        weighted_preds <- stack_predictions %*% weights
        # apply sigmoid
        weighted_preds <- 1 / (1 + exp(-weighted_preds))
        # log loss
        -mean(
          true_labels * log(weighted_preds) + (1 - true_labels) *
            log(1 - weighted_preds)
        )
      }

      # initial guess for weights - same weights for all models
      initial_weights <- rep(
        x = 1 / (ncol(stack_predictions) - 1),
        times = ncol(stack_predictions - 1)
      )

      # minimize the objective function
      result <- stats:::nlm(objective, initial_weights)

      # opptimal weights
      optimal_weights <- result$estimate
      names(optimal_weights) <- colnames(stack_predictions)

      # stacked prediction
      stack_out <- stack_out %>%
        dplyr::mutate(
          stack_nlm = as.numeric(stack_predictions %*% optimal_weights)
        ) %>%
        dplyr::mutate(
          stack_nlm = 1 / (1 + exp(-stack_nlm))
        )

      # add optimal weights to the stack_weights list
      stack_weights$stack_nlm <- optimal_weights

      #----------------------------------------------------------------
      # Ornstein hillclimbing for binary response
      #----------------------------------------------------------------

      # initial weights for hill climbing
      weights_ornstein <- rep(x = 1, times = ncol(stack_predictions))

      # predictions based on hillclimbing weights
      stack_ornstein <- stack_predictions %*% weights_ornstein /
        sum(weights_ornstein)

      # compute log loss
      c_log_loss <- -mean(
        true_labels * log(stack_ornstein) + (1 - true_labels) *
          log(1 - stack_ornstein)
      )

      # container for best log loss
      best_log_loss <- c_log_loss

      # while loop control
      keep_going <- TRUE

      # hill climbing until number of iterations is reached or no improvement
      # through addition of models
      while (keep_going) {

        # keep track of which model would be best 'greedy' addition to ensemble
        best_addition <- 0L

        # iterate over all models
        for (j in seq_along(weights_ornstein)) {

          # increase model weight by 1
          weights_ornstein[j] <- weights_ornstein[j] + 1L

          # generate predictions with updated weights
          preds_update <- stack_predictions %*% weights_ornstein /
            sum(weights_ornstein)

          # compute log loss of updated predictions
          updated_loglos <- -mean(
            true_labels * log(preds_update) + (1 - true_labels) *
              log(1 - preds_update)
          )

          # check whether updated log loss is better than current best
          if (updated_loglos < best_log_loss) {
            best_log_loss <- updated_loglos
            # model addition that best improves log loss
            best_addition <- j
          }

          # reset weights
          weights_ornstein[j] <- weights_ornstein[j] - 1L
        } # end of loop over all models

        # if we found an improvement (and we're below the cutoff number of
        # iterations), then add the best model
        if (sum(weights_ornstein) < 1000 && best_log_loss < c_log_loss) {

          weights_ornstein[best_addition] <- weights_ornstein[best_addition] + 1L
          c_log_loss <- best_log_loss

        } else {
          # stop the while loop if no model addition improves log loss
          keep_going <- FALSE
        }
      } # end of while loop

      # normalize the weights
      weights_ornstein <- weights_ornstein %>%
        magrittr::divide_by(sum(weights_ornstein))
      names(weights_ornstein) <- colnames(stack_predictions)

      # add optimal weights to the stack_weights list
      stack_weights$stack_ornstein <- weights_ornstein

      # predictions based on hillclimbing weights
      stack_out <- stack_out %>%
        dplyr::mutate(
          stack_ornstein = as.numeric(stack_predictions %*% weights_ornstein)
        )

      #----------------------------------------------------------------
      # Stack of staacks for binary response without EBMA
      #----------------------------------------------------------------

      # add the true labels to stacked predictions
      stack_out <- dplyr::mutate(.data = stack_out, y = true_labels) %>%
        dplyr::relocate(y, .before = 1)

      # use a logistic regression to combine the stacked predictions
      m_blend <- glm(y ~ ., data = stack_out, family = binomial)

      # add model predictions
      stack_out <- stack_out %>%
        dplyr::mutate(
          stack_of_stacks = predict(
            m_blend, newdata = stack_out, type = "response"
          )
        )

      # add optimal weights to the stack_weights list
      stack_weights$stack_of_stacks <- coef(m_blend)

      #----------------------------------------------------------------
      # Stack of staacks for binary response with EBMA
      #----------------------------------------------------------------

      # add the true labels to stacked predictions
      stack_with_ebma <- stack_out %>%
        dplyr::mutate(
          ebma_preds = ebma_out$individual_level_predictions$ebma_preds
        ) %>%
        dplyr::select(y, stack_of_stacks, ebma_preds)

      # use a logistic regression to combine the stacked predictions
      m_blend <- glm(y ~ ., data = stack_with_ebma, family = binomial)

      # add model predictions
      stack_out <- stack_out %>%
        dplyr::mutate(
          stack_of_stacks_with_ebma = predict(
            m_blend, newdata = stack_with_ebma, type = "response"
          )
        )

      # add optimal weights to the stack_weights list
      stack_weights$stack_of_stacks_with_ebma <- coef(m_blend)

    } # end if binary_dv

    # stacking for continuous dependent variables
    if (binary_dv == FALSE) {
      stop("Stacking for continuous dependent variables not implemented yet.")
    }

    # return the stacked predictions and weights
    return(list(stack_preds = stack_out, stack_weights = stack_weights))
  }
}