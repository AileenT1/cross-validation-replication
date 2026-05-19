here::i_am("src/cv_error_rate.R")
source(here::here("src", "cv_linear.R"))
source(here::here("src", "cv_polynomial.R"))

run_cv_selection_error_rate_sim <- function(sim_config) {
  sample_sizes <- sim_config$sample_sizes
  mc_runs <- sim_config$mc_runs

  by_n_list <- vector("list", length(sample_sizes))
  for (i in seq_along(sample_sizes)) {
    n_now <- sample_sizes[[i]]
    run_results <- lapply(seq_len(mc_runs), function(replication_id) {
      x <- rnorm(n_now)
      eps <- rnorm(n_now, sd = sim_config$sigma)
      if (sim_config$truth_type == "linear") {
        y <- sim_config$beta0 + sim_config$beta1 * x + eps
      } else {
        y <- sim_config$beta0 + sim_config$beta1 * x + sim_config$beta2 * x^2 + sim_config$beta3 * x^3 + eps
      }
      dat <- data.frame(y = y, x = x)
      fold_id <- sample(rep(seq_len(sim_config$K), length.out = n_now))

      linear_cv <- kfold_cv_linear(
        data = dat,
        formula = stats::as.formula("y ~ x"),
        K = sim_config$K,
        fold_id = fold_id
      )
      poly_cv <- kfold_cv_polynomial(
        data = dat,
        response = "y",
        predictor = "x",
        degree = sim_config$poly_degree,
        K = sim_config$K,
        raw = TRUE,
        fold_id = fold_id
      )

      mean_cv_linear <- linear_cv$mean_cv_mse
      mean_cv_poly3 <- poly_cv$mean_cv_mse
      selected_model <- if (mean_cv_linear <= mean_cv_poly3) "linear" else "poly3"
      wrong <- if (sim_config$truth_type == "linear") {
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
    if (i %% 25 == 0 || i == length(sample_sizes)) {
      message("Finished sample size index ", i, " / ", length(sample_sizes), " (n = ", n_now, ")")
    }
  }
  do.call(rbind, by_n_list)
}
