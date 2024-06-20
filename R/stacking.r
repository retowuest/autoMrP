# classifier predictions often correlate highly, therefore we use
# lasso to stack them
stackweights <- function(s_data, L2.unit) {

  # model formula
  form_stack <- as.formula(
    paste(
      # DV
      colnames(s_data)[1], "~",
      # ignore intercept
      "-1 +",
      # IVs
      paste(
        colnames(s_data)[2:length(colnames(s_data))], collapse = " + "
      )
    )
  )

  # stacking via a simple glm
  m_stack <- glm(
    formula = form_stack, data = s_data, family = binomial(link = probit)
  )

  preds <- predict(
    m_stack, newdata = census, type = "response"
  )


  # # add best subset prediction to cenus
  # test <- census %>%
  #   dplyr::mutate(
  #     stacking = stats::predict(
  #       object = m_stack,
  #       newdata = .,
  #       allow.new.levels = TRUE,
  #       type = "response"
  #     )
  #   )
    
    
  #    %>%
  #   dplyr::group_by(!! rlang::sym(L2.unit)) %>%
  #   dplyr::summarize(staticking = mean(stacking))
    
    
    #  %>%
    # dplyr::summarize(stacking = stats::weighted.mean(
    #   x = stacking, w = prop
    # ), .groups = "keep") %>%
    # dplyr::ungroup()




}