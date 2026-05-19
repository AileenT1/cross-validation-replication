# Set.seed() before run_rolling_selection_error_rate_sim() for reproducibility.
here::i_am("src/rolling_error_rate.R")
source(here::here("src", "rolling_validation.R"))

run_rolling_selection_error_rate_sim <- function(rv_config) {
  sample_sizes <- rv_config$sample_sizes
  mc_runs <- rv_config$mc_runs

  by_n_list <- vector("list", length(sample_sizes))
  for (i in seq_along(sample_sizes)) {
    n_now <- sample_sizes[[i]]
    run_results <- lapply(seq_len(mc_runs), function(replication_id) {
      x <- rnorm(n_now)
      eps <- rnorm(n_now, sd = rv_config$sigma)
      if (rv_config$truth_type == "linear") {
        y <- rv_config$beta0 + rv_config$beta1 * x + eps
      } else {
        y <- rv_config$beta0 + rv_config$beta1 * x + rv_config$beta2 * x^2 + rv_config$beta3 * x^3 + eps
      }
      dat <- data.frame(y = y, x = x)
      dat <- dat[order(dat$x, dat$y), , drop = FALSE]

      linear_rv <- rolling_validate_regression(
        data = dat,
        formula = stats::as.formula("y ~ x"),
        rv_config = rv_config
      )
      poly_rv <- rolling_validate_polynomial(
        data = dat,
        response = "y",
        predictor = "x",
        rv_config = rv_config
      )

      mean_rv_linear <- linear_rv$mean_rv_mse
      mean_rv_poly3 <- poly_rv$mean_rv_mse
      selected_model <- if (mean_rv_linear <= mean_rv_poly3) "linear" else "poly3"
      wrong <- if (rv_config$truth_type == "linear") {
        selected_model == "poly3"
      } else {
        selected_model == "linear"
      }

      data.frame(
        n = n_now,
        replication = replication_id,
        wrong = wrong,
        stringsAsFactors = FALSE
      )
    })
    df <- do.call(rbind, run_results)
    wrong_rate <- mean(df$wrong)
    mc_se <- sqrt(wrong_rate * (1 - wrong_rate) / mc_runs)
    by_n_list[[i]] <- data.frame(
      n = n_now,
      error_rate = wrong_rate,
      mc_se = mc_se,
      stringsAsFactors = FALSE
    )
    if (i %% 5 == 0 || i == length(sample_sizes)) {
      message("Finished sample size index ", i, " / ", length(sample_sizes), " (n = ", n_now, ")")
    }
  }
  do.call(rbind, by_n_list)
}
