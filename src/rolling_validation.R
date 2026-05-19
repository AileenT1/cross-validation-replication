rolling_validate_regression <- function(data, formula, rv_config) {
  n <- nrow(data)
  w <- floor(rv_config$train_frac * n)
  y_name <- as.character(formula[[2L]])

  if (rv_config$window_type == "expanding") {
    n_steps <- n - w
    sq_errors <- numeric(n_steps)
    j <- 1L
    for (t in w:(n - 1L)) {
      train_index <- 1:t
      valid_index <- t + 1L
      train_data <- data[train_index, , drop = FALSE]
      val_row <- data[valid_index, , drop = FALSE]
      fit <- stats::lm(formula, data = train_data)
      pred <- stats::predict(fit, newdata = val_row)
      y_val <- val_row[[y_name]]
      sq_errors[j] <- (y_val - pred)^2
      j <- j + 1L
    }
  } else {
    n_steps <- n - w
    sq_errors <- numeric(n_steps)
    j <- 1L
    for (t in 1L:(n - w)) {
      train_index <- t:(t + w - 1L)
      valid_index <- t + w
      train_data <- data[train_index, , drop = FALSE]
      val_row <- data[valid_index, , drop = FALSE]
      fit <- stats::lm(formula, data = train_data)
      pred <- stats::predict(fit, newdata = val_row)
      y_val <- val_row[[y_name]]
      sq_errors[j] <- (y_val - pred)^2
      j <- j + 1L
    }
  }

  list(
    mean_rv_mse = mean(sq_errors),
    n_steps = n_steps,
    w = w
  )
}

rolling_validate_polynomial <- function(data, response, predictor, rv_config) {
  formula <- stats::as.formula(sprintf(
    "`%s` ~ poly(`%s`, degree = %d, raw = TRUE)",
    response,
    predictor,
    as.integer(rv_config$poly_degree)
  ))
  rolling_validate_regression(data, formula, rv_config)
}
