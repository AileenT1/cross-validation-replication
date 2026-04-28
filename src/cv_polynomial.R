if (!exists("kfold_cv_linear", mode = "function", inherits = TRUE)) {
  source("src/cv_linear.R")
}

kfold_cv_polynomial <- function(
    data,
    response,
    predictor,
    degree,
    K,
    raw = TRUE,
    fold_id = NULL) {
  formula <- stats::as.formula(sprintf(
    "`%s` ~ poly(`%s`, degree = %d, raw = %s)",
    response,
    predictor,
    as.integer(degree),
    if (isTRUE(raw)) "TRUE" else "FALSE"
  ))

  kfold_cv_linear(data, formula, K = K, fold_id = fold_id)
}
