rolling_validate_regression <- function(data, formula, model_fitting_config) {
  rv <- model_fitting_config$rv
  n <- nrow(data)

  r <- rv$window_size
  h <- ifelse(is.null(rv$horizon), 1, rv$horizon)
  step <- ifelse(is.null(rv$validation_step), 1, rv$validation_step)

  y_name <- as.character(formula[[2]])

  valid_values <- seq(r + h, n, by = step)
  sq_errors <- numeric(length(valid_values))

  for (j in seq_along(valid_values)) {
    i <- valid_values[j]

    train_index <- 1:(i - h)
    valid_index <- i

    fit <- lm(formula, data = data[train_index, , drop = FALSE])
    pred <- predict(fit, newdata = data[valid_index, , drop = FALSE])

    y_val <- data[[y_name]][valid_index]

    sq_errors[j] <- (y_val - pred)^2
  }

  list(
    mean_rv_mse = mean(sq_errors),
    n_steps = length(sq_errors),
    w = r
  )
}

rolling_validate_polynomial <- function(data, response, predictor, model_fitting_config) {
  formula <- as.formula(sprintf(
    "`%s` ~ poly(`%s`, degree = %d, raw = TRUE)",
    response,
    predictor,
    as.integer(model_fitting_config$model$poly_degree)
  ))
  rolling_validate_regression(data, formula, model_fitting_config)
}
