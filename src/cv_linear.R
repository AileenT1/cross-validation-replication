kfold_cv_linear <- function(data, formula, K, fold_id = NULL) {
  n <- nrow(data)
  stopifnot(K >= 2L, n >= K)
  if (is.null(fold_id)) {
    fold_id <- sample(rep(seq_len(K), length.out = n))
  } else {
    stopifnot(length(fold_id) == n)
  }

  lhs <- formula[[2L]]
  if (!is.symbol(lhs)) {
    stop("Left-hand side of `formula` must be a single column name (e.g. y ~ x), not a transform.", call. = FALSE)
  }
  y_name <- as.character(lhs)

  mse_by_fold <- numeric(K)

  for (k in seq_len(K)) {
    train <- fold_id != k
    val <- fold_id == k

    train_data <- data[train, , drop = FALSE]
    val_data <- data[val, , drop = FALSE]

    fit <- stats::lm(formula, data = train_data)
    pred <- stats::predict(fit, newdata = val_data)
    y_val <- val_data[[y_name]]
    mse_by_fold[k] <- mean((y_val - pred)^2)
  }

  list(
    mean_cv_mse = mean(mse_by_fold),
    mse_by_fold = mse_by_fold,
    fold_id = fold_id
  )
}
