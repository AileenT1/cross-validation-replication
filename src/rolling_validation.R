rolling_validate_regression <- function(data, formula, model_fitting_config) {
  rv <- model_fitting_config$rv
  n <- nrow(data)

  if (rv$window_type == "fixed") {
    w <- rv$window_size
  } else {
    w <- floor(rv$train_frac * n)
  }

  h <- rv$horizon
  step <- ifelse(is.null(rv$validation_step), 1L, rv$validation_step)

  y_name <- as.character(formula[[2L]])

  if (rv$window_type == "expanding") {
    t_values <- seq(w, n - h, by = step)
    sq_errors <- numeric(length(t_values))

    for (j in seq_along(t_values)) {
      t <- t_values[j]
      train_index <- 1:t
      valid_index <- t + h

      fit <- stats::lm(formula, data = data[train_index, , drop = FALSE])
      pred <- stats::predict(fit, newdata = data[valid_index, , drop = FALSE])
      y_val <- data[[y_name]][valid_index]

      sq_errors[j] <- (y_val - pred)^2
    }

  } else {
    t_values <- seq(1L, n - w - h + 1L, by = step)
    sq_errors <- numeric(length(t_values))

    for (j in seq_along(t_values)) {
      t <- t_values[j]
      train_index <- t:(t + w - 1L)
      valid_index <- t + w + h - 1L

      fit <- stats::lm(formula, data = data[train_index, , drop = FALSE])
      pred <- stats::predict(fit, newdata = data[valid_index, , drop = FALSE])
      y_val <- data[[y_name]][valid_index]

      sq_errors[j] <- (y_val - pred)^2
    }
  }

  list(
    mean_rv_mse = mean(sq_errors),
    n_steps = length(sq_errors),
    w = w
  )
}

rolling_validate_polynomial <- function(data, response, predictor, model_fitting_config) {
  formula <- stats::as.formula(sprintf(
    "`%s` ~ poly(`%s`, degree = %d, raw = TRUE)",
    response,
    predictor,
    as.integer(model_fitting_config$model$poly_degree)
  ))
  rolling_validate_regression(data, formula, model_fitting_config)
}
