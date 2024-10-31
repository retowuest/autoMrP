stacking_weights <- function(preds, ebma_out, L2.unit, k.folds, cores) {

  # initial binding of globals
  stack_lme4 <- stack_nlm <- y <- stack_of_stacks <- ebma_preds <- k <- NULL
  id <- NULL

  # check whether more than one classifier was tuned
  if (ncol(preds) <= 3) {
    message("Skipping stacking, only one classifier was tuned.")
    # return the stacked predictions and weights
    return(list(stack_preds = NULL, stack_weights = NULL))
  } else {

    # log loss
    log_loss_fun <- function(true_labels, predictions) {
      l_loss <- -mean(
        true_labels * log(predictions) + (1 - true_labels) *
          log(1 - predictions)
      )
      return(l_loss)
    }

    # Objective function to minimize (log loss) with non-negative least squares
    # adds a small number to avoid log(0)
    objective <- function(weights) {
      weighted_preds <- stack_predictions %*% weights
      # Apply sigmoid
      weighted_preds <- 1 / (1 + exp(-weighted_preds))
      # Log loss
      -mean(
        true_labels * log(weighted_preds + 1e-15) +
          (1 - true_labels) * log(1 - weighted_preds + 1e-15)
      )
    }

    # add IDs
    preds <- preds %>%
      dplyr::mutate(id = seq_len(nrow(.)))
    ebma_out$individual_level_predictions <-
      ebma_out$individual_level_predictions %>%
      dplyr::mutate(id = seq_len(nrow(.)))

    # generate folds for stacking
    cv_folds <- cv_folding(
      data = preds,
      L2.unit = L2.unit,
      k.folds = k.folds,
      cv.sampling = "L2 units"
    )

    # is the dpendent variable binary or not?
    binary_dv <- if (length(unique(preds$y)) == 2) {
      TRUE
    } else {
      FALSE
    }

    # stacking for binary dependent variables
    if (binary_dv) {

      #--------------------------------------------------
      # cross-validation stacking
      #--------------------------------------------------

      # Register cores
      cl <- multicore(cores = cores, type = "open", cl = NULL)

      # loop over k folds
      k_weights <- foreach::foreach(
        k = seq_along(cv_folds), .packages = "autoMrP"
      ) %dorng% {

        `%>%` <- magrittr::`%>%`

        # generate weights based on training data
        data_train <- dplyr::bind_rows(cv_folds[-k])
        stack_predictions <- data_train %>%
          dplyr::select(-y, -!!rlang::sym(L2.unit), -id) %>%
          as.matrix()

        # evaluate based on test data
        data_test <- cv_folds[[k]] %>%
          dplyr::select(-y, -!!rlang::sym(L2.unit), -id) %>%
          as.matrix()

        #--------------------------------------------------
        # non-negative least squares
        #--------------------------------------------------
        ols_coefs <- nnls::nnls(
          as.matrix(stack_predictions), data_train$y
        )$x
        names(ols_coefs) <- colnames(stack_predictions)

        # normalize the weights
        ols_coefs <- ols_coefs / sum(ols_coefs)

        # weights for base models (replace NAs with 0)
        stacking_weights <- ols_coefs
        stacking_weights[is.na(ols_coefs)] <- 0

        # test data log loss
        test_preds <- data_test %*% stacking_weights
        test_error <- log_loss_fun(cv_folds[[k]]$y, test_preds)

        # store results
        nnls_out <- list(
          weights = stacking_weights,
          test_error = test_error
        )

        #--------------------------------------------------
        # optim
        #--------------------------------------------------

        # true labels
        true_labels <- data_train$y

        # Initial guess for weights
        initial_weights <- rep(
          1 / ncol(stack_predictions), ncol(stack_predictions)
        )

        # Define bounds for weights
        lower_bounds <- rep(0, ncol(stack_predictions))
        upper_bounds <- rep(Inf, ncol(stack_predictions))

        # Run the optimization
        result <- nloptr::nloptr(
          x0 = initial_weights,
          eval_f = objective,
          lb = lower_bounds,
          ub = upper_bounds,
          opts = list(
            algorithm = "NLOPT_LN_COBYLA",
            xtol_rel = 1.0e-8,
            maxeval = 1000
          )
        )

        # Optimal weights
        stacking_weights <- result$solution
        names(stacking_weights) <- colnames(stack_predictions)

        # normalize the weights
        stacking_weights <- stacking_weights / sum(stacking_weights)

        # test data log loss
        test_preds <- data_test %*% stacking_weights
        test_error <- log_loss_fun(cv_folds[[k]]$y, test_preds)

        # store results
        optim_out <- list(
          weights = stacking_weights,
          test_error = test_error
        )

        #--------------------------------------------------
        # quadratic programming
        #--------------------------------------------------
        d_m <- t(stack_predictions) %*% stack_predictions
        d_v <- t(stack_predictions) %*% data_train$y

        # regularize
        epsilon <- 1e-8
        d_m <- d_m + epsilon * diag(ncol(d_m))

        # constraints (non-negative weights)
        a_m <- rbind(
          rep(1, ncol(stack_predictions)), diag(ncol(stack_predictions))
        )
        b_v <- c(1, rep(0, ncol(stack_predictions)))

        # solve the quadratic programming problem
        qp_solution <- quadprog::solve.QP(
          Dmat = d_m,
          dvec = d_v,
          Amat = t(a_m),
          bvec = b_v,
          meq = 1 # Equality constraint for sum(weights) = 1
        )

        # get weights
        stacking_weights <- qp_solution$solution
        names(stacking_weights) <- colnames(stack_predictions)

        # evaluate based on test data
        test_preds <- data_test %*% stacking_weights
        test_error <- log_loss_fun(cv_folds[[k]]$y, test_preds)

        # store results
        qp_out <- list(
          weights = stacking_weights,
          test_error = test_error
        )

        #--------------------------------------------------
        # Ornstein hillclimbing
        #--------------------------------------------------

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

          # keep track of which model would be best 'greedy' addition to
          # ensemble
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

            weights_ornstein[best_addition] <- weights_ornstein[best_addition] +
              1L
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

        # evaluate based on test data
        test_preds <- data_test %*% weights_ornstein
        test_error <- log_loss_fun(cv_folds[[k]]$y, test_preds)

        # store results
        ornstein_out <- list(
          weights = weights_ornstein,
          test_error = test_error
        )

        return(
          list(
            nnls = nnls_out,
            optim = optim_out,
            qp = qp_out,
            ornstein = ornstein_out
          )
        )
      }

      # De-register cluster
      multicore(cores = cores, type = "close", cl = cl)

      # all data
      stack_predictions <- cv_folds %>%
        dplyr::bind_rows() %>%
        dplyr::select(-y, -!!rlang::sym(L2.unit), -id) %>%
        as.matrix()

      # non-negative least squares stacking
      nnls_out <- do.call(rbind, k_weights)[, "nnls"]
      # weights
      nnls_w <- do.call(rbind, nnls_out)[, "weights"] %>%
        dplyr::bind_rows()
      # test error
      nnls_e <- do.call(rbind, nnls_out)[, "test_error"] %>%
        unlist()
      nnls_e <- 1 / nnls_e
      nnls_e <- nnls_e / sum(nnls_e)
      # genearte a weighted mean of stacking_weights based on test error
      nnls_w <- apply(nnls_w, 2, function(x) {
        weighted.mean(x = x, w = nnls_e)
      })

      # individual level predictions
      final_prediction <- stack_predictions %*% nnls_w

      # stacked predictions output
      stack_out <- tibble::tibble(
        id = cv_folds %>%
          dplyr::bind_rows() %>%
          dplyr::pull(id),
        y = cv_folds %>%
          dplyr::bind_rows() %>%
          dplyr::pull(y),
        !!rlang::sym(L2.unit) := cv_folds %>%
          dplyr::bind_rows() %>%
          dplyr::pull(!!rlang::sym(L2.unit)),
        stack_nnls = as.numeric(final_prediction)
      )

      # stacking weights container
      stack_weights <- list(
        stack_nnls = nnls_w
      )

      # optim stacking weights
      optim_out <- do.call(rbind, k_weights)[, "optim"]
      # weights
      optim_w <- do.call(rbind, optim_out)[, "weights"] %>%
        dplyr::bind_rows()
      # test error
      optim_e <- do.call(rbind, optim_out)[, "test_error"] %>%
        unlist()
      optim_e <- 1 / optim_e
      optim_e <- optim_e / sum(optim_e)
      # genearte a weighted mean of stacking_weights based on test error
      optim_w <- apply(optim_w, 2, function(x) {
        weighted.mean(x = x, w = optim_e)
      })

      # individual level predictions
      final_prediction <- stack_predictions %*% optim_w

      # stacked predictions output
      stack_out <- stack_out %>%
        dplyr::mutate(
          stack_optim = as.numeric(final_prediction)
        )

      # stacking weights container
      stack_weights$stack_optim <- optim_w

      # quadratic programming stacking
      qp_out <- do.call(rbind, k_weights)[, "qp"]
      # weights
      qp_w <- do.call(rbind, qp_out)[, "weights"] %>%
        dplyr::bind_rows()
      # test error
      qp_e <- do.call(rbind, qp_out)[, "test_error"] %>%
        unlist()
      qp_e <- 1 / qp_e
      qp_e <- qp_e / sum(qp_e)
      # genearte a weighted mean of stacking_weights based on test error
      qp_w <- apply(qp_w, 2, function(x) {
        weighted.mean(x = x, w = qp_e)
      })

      # individual level predictions
      final_prediction <- stack_predictions %*% qp_w

      # stacked predictions output
      stack_out <- stack_out %>%
        dplyr::mutate(
          qp_optim = as.numeric(final_prediction)
        )

      # stacking weights container
      stack_weights$stack_qp <- qp_w

      # Orstein stacking
      ornstein_out <- do.call(rbind, k_weights)[, "ornstein"]
      # weights
      ornstein_w <- do.call(rbind, ornstein_out)[, "weights"] %>%
        dplyr::bind_rows()
      # test error
      ornstein_e <- do.call(rbind, ornstein_out)[, "test_error"] %>% unlist()
      ornstein_e <- 1 / ornstein_e
      ornstein_e <- ornstein_e / sum(ornstein_e)
      # genearte a weighted mean of stacking_weights based on test error
      ornstein_w <- apply(ornstein_w, 2, function(x) {
        weighted.mean(x = x, w = ornstein_e)
      })

      # individual level predictions
      final_prediction <- stack_predictions %*% ornstein_w

      # stacked predictions output
      stack_out <- stack_out %>%
        dplyr::mutate(
          stack_ornstein = as.numeric(final_prediction)
        )

      # stacking weights container
      stack_weights$stack_ornstein <- ornstein_w

      # sort the stacked predictions by id
      stack_out <- stack_out %>%
        dplyr::arrange(id)

      #----------------------------------------------------------------
      # 2nd round of stacking: stack of stacking algorithms
      #----------------------------------------------------------------

      # generate folds
      cv_folds <- cv_folding(
        data = stack_out,
        L2.unit = L2.unit,
        k.folds = k.folds,
        cv.sampling = "L2 units"
      )

      # Register cores
      cl <- multicore(cores = cores, type = "open", cl = NULL)

      # loop over k folds
      k_weights <- foreach::foreach(
        k = seq_along(cv_folds), .packages = "autoMrP"
      ) %dorng% {

        `%>%` <- magrittr::`%>%`

        # generate weights based on training data
        data_train <- dplyr::bind_rows(cv_folds[-k])
        stack_predictions <- data_train %>%
          dplyr::select(-y, -!!rlang::sym(L2.unit), -id) %>%
          as.matrix()

        # evaluate based on test data
        data_test <- cv_folds[[k]] %>%
          dplyr::select(-y, -!!rlang::sym(L2.unit), -id) %>%
          as.matrix()

        true_labels <- data_train$y

        #--------------------------------------------------
        # Stacks without EBMA via Ornstein hillclimbing
        #--------------------------------------------------

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

          # keep track of which model would be best 'greedy' addition to
          # ensemble
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

            weights_ornstein[best_addition] <- weights_ornstein[best_addition] +
              1L
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

        # evaluate based on test data
        test_preds <- data_test %*% weights_ornstein
        test_error <- log_loss_fun(cv_folds[[k]]$y, test_preds)

        # store results
        return(
          list(
            weights = weights_ornstein,
            test_error = test_error
          )
        )
      }

      # De-register cluster
      multicore(cores = cores, type = "close", cl = cl)

      # all data
      stack_predictions <- cv_folds %>%
        dplyr::bind_rows() %>%
        dplyr::select(-y, -!!rlang::sym(L2.unit), -id) %>%
        as.matrix()

      # 2nd round of stacking output
      # weights
      stack2_w <- do.call(rbind, k_weights)[, "weights"] %>%
        dplyr::bind_rows()
      # test error
      stack2_e <- do.call(rbind, k_weights)[, "test_error"] %>%
        unlist()
      stack2_e <- 1 / stack2_e
      stack2_e <- stack2_e / sum(stack2_e)
      # genearte a weighted mean of stacking_weights based on test error
      stack2_w <- apply(stack2_w, 2, function(x) {
        weighted.mean(x = x, w = stack2_e)
      })

      # individual level predictions
      final_prediction <- stack_predictions %*% stack2_w

      # stck of stacks output
      stack_of_stacks <- tibble::tibble(
        id = cv_folds %>%
          dplyr::bind_rows() %>%
          dplyr::pull(id),
        stack_of_stacks = as.numeric(final_prediction)
      ) %>%
        dplyr::arrange(id)

      # merge predictions into stack_out
      stack_out <- stack_out %>%
        dplyr::left_join(stack_of_stacks, by = "id")

      # stacking weights container
      stack_weights$stack_of_stacks <- stack2_w

      #----------------------------------------------------------------
      # 3rd round of stacking: stack of stacks with ebma
      #----------------------------------------------------------------

      # data for third round of stacking
      stack_data3 <- stack_out %>%
        dplyr::select(id, y, !!rlang::sym(L2.unit), stack_of_stacks)

      # ebma predictions
      ebma_preds <- ebma_out$individual_level_predictions %>%
        dplyr::select(id, ebma_preds)

      # join ebma predictions
      stack_data3 <- stack_data3 %>%
        dplyr::left_join(ebma_preds, by = "id")

      # generate folds
      cv_folds <- cv_folding(
        data = stack_data3,
        L2.unit = L2.unit,
        k.folds = k.folds,
        cv.sampling = "L2 units"
      )

      # Register cores
      cl <- multicore(cores = cores, type = "open", cl = NULL)

      # loop over k folds
      k_weights <- foreach::foreach(
        k = seq_along(cv_folds), .packages = "autoMrP"
      ) %dorng% {

        `%>%` <- magrittr::`%>%`

        # generate weights based on training data
        data_train <- dplyr::bind_rows(cv_folds[-k])
        stack_predictions <- data_train %>%
          dplyr::select(-y, -!!rlang::sym(L2.unit), -id) %>%
          as.matrix()

        # evaluate based on test data
        data_test <- cv_folds[[k]] %>%
          dplyr::select(-y, -!!rlang::sym(L2.unit), -id) %>%
          as.matrix()

        true_labels <- data_train$y

        #--------------------------------------------------
        # Stacks without EBMA via Ornstein hillclimbing
        #--------------------------------------------------

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

          # keep track of which model would be best 'greedy' addition to
          # ensemble
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

            weights_ornstein[best_addition] <- weights_ornstein[best_addition] +
              1L
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

        # evaluate based on test data
        test_preds <- data_test %*% weights_ornstein
        test_error <- log_loss_fun(cv_folds[[k]]$y, test_preds)

        # store results
        return(
          list(
            weights = weights_ornstein,
            test_error = test_error
          )
        )
      }

      # De-register cluster
      multicore(cores = cores, type = "close", cl = cl)

      # all data
      stack_predictions <- cv_folds %>%
        dplyr::bind_rows() %>%
        dplyr::select(-y, -!!rlang::sym(L2.unit), -id) %>%
        as.matrix()

      # 3rd round of stacking output
      # weights
      stack3_w <- do.call(rbind, k_weights)[, "weights"] %>%
        dplyr::bind_rows()
      # test error
      stack3_e <- do.call(rbind, k_weights)[, "test_error"] %>%
        unlist()
      stack3_e <- 1 / stack3_e
      stack3_e <- stack3_e / sum(stack3_e)
      # genearte a weighted mean of stacking_weights based on test error
      stack3_w <- apply(stack3_w, 2, function(x) {
        weighted.mean(x = x, w = stack3_e)
      })

      # individual level predictions
      final_prediction <- stack_predictions %*% stack3_w

      # stck of stacks output
      stack_of_stacks_ebma <- tibble::tibble(
        id = cv_folds %>%
          dplyr::bind_rows() %>%
          dplyr::pull(id),
        stack_of_stacks_ebma = as.numeric(final_prediction)
      ) %>%
        dplyr::arrange(id)

      # merge predictions into stack_out
      stack_out <- stack_out %>%
        dplyr::left_join(stack_of_stacks_ebma, by = "id")


      # stacking weights container
      stack_weights$stack_of_stacks_ebma <- stack3_w

    } # end if binary_dv

    # stacking for continuous dependent variables
    if (binary_dv == FALSE) {
      stop("Stacking for continuous dependent variables not implemented yet.")
    }

    # return the stacked predictions and weights
    return(list(stack_preds = stack_out, stack_weights = stack_weights))
  }
}