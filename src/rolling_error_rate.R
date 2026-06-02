here::i_am("src/rolling_error_rate.R")
source(here::here("src", "rolling_validation.R"))
source(here::here("src", "simulate_regression_data.R"))

run_rolling_selection_error_rate_sim <- function(simulation_config, model_fitting_config) {
  sample_sizes <- simulation_config$sample_sizes
  mc_runs <- simulation_config$mc_runs
  truth <- simulation_config$truth
  candidate_models <- model_fitting_config$candidate_models

  by_n_list <- vector("list", length(sample_sizes))
  for (i in seq_along(sample_sizes)) {
    n_now <- sample_sizes[[i]]
    run_results <- lapply(seq_len(mc_runs), function(replication_id) {
      dat <- simulate_regression_data(n_now, simulation_config)
      dat <- dat[order(dat$x, dat$y), , drop = FALSE]

      linear_rv <- rolling_validate_regression(
        data = dat,
        formula = candidate_models$linear,
        model_fitting_config = model_fitting_config
      )
      poly_rv <- rolling_validate_regression(
        data = dat,
        formula = candidate_models$polynomial,
        model_fitting_config = model_fitting_config
      )

      mean_linear <- linear_rv$mean_rv_mse
      mean_poly3 <- poly_rv$mean_rv_mse
      selected_model <- if (mean_linear <= mean_poly3) "linear" else "poly3"
      wrong <- if (truth$type == "linear") {
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
    if (i %% 5L == 0L || i == length(sample_sizes)) {
      message("Finished sample size index ", i, " / ", length(sample_sizes), " (n = ", n_now, ")")
    }
  }
  do.call(rbind, by_n_list)
}
