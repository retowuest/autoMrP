stacking_weights <- function(preds) {

  #--------------------------------------------------
  # regression for binary response
  #--------------------------------------------------
  # meta formula
  meta_formula <- as.formula(
    paste("y ~ ", paste(names(preds)[-1], collapse = " + "))
  )

  # logistic regression model
  logit_m <- glm(meta_formula, data = preds, family = binomial)

  # coefficients of the logistic model
  logit_coefs <- coef(logit_m)

  # weights for base models
  stacking_weights <- logit_coefs[-1]  # Exclude the intercept
  stacking_weights[is.na(stacking_weights)] <- 0

  # predicions using logistic regression weights
  final_prediction <- as.matrix(preds[, -1]) %*% stacking_weights

  # apply sigmoid
  final_prediction <- 1 / (1 + exp(-final_prediction))

  preds$stacking <- final_prediction

  #----------------------------------------------------------------
  # opitmization for binary response
  #----------------------------------------------------------------
stack_predictions <- as.matrix(stack_predictions)
true_labels <- preds$y
  
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

  # initial guess for weights
  initial_weights <- rep(1 / (ncol(stack_predictions) - 1), ncol(stack_predictions - 1))

  # minimize the objective function
  result <- stats:::nlm(objective, initial_weights)

  # Optimal weights
  optimal_weights <- result$estimate
  optimal_weights

}