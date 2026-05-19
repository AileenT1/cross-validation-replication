library(dplyr)
library(ggplot2)

here::i_am("code/03_expanding_linear.R")
source(here::here("src", "rolling_error_rate.R"))

#Main simulation configuration
linear_truth_rv_config <- list(
  window_type = "expanding",
  train_frac = 0.5,
  horizon = 1,
  mc_runs = 100,
  sample_sizes = seq(50, 1000, by = 50),
  beta0 = 1,
  beta1 = 2,
  sigma = 4,
  poly_degree = 3,
  truth_type = "linear",
  truth = "linear: y = beta0 + beta1*x + eps",
  selection_rule = "Choose model with lower mean RV MSE; tie -> linear"
)

results_dir <- here::here("data", "final")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
rds_path <- file.path(results_dir, "03_expanding_linear_results.rds")

expanding_linear_results <- if (file.exists(rds_path)) readRDS(rds_path) else NULL
if (!is.null(expanding_linear_results) &&
    all(c("rv_config", "error_rate_by_n") %in% names(expanding_linear_results))) {
  linear_truth_rv_config <- expanding_linear_results$rv_config
  error_rate_by_n <- expanding_linear_results$error_rate_by_n
  message("Loaded results from ", rds_path)
} else {
  set.seed(2026)
  error_rate_by_n <- run_rolling_selection_error_rate_sim(linear_truth_rv_config)
  saveRDS(
    list(rv_config = linear_truth_rv_config, error_rate_by_n = error_rate_by_n),
    rds_path
  )
  message("Saved results to ", rds_path)
}

#Graphing and tables
expanding_linear_settings_tbl <- data.frame(
  setting = c(
    "seed", "window type", "train_frac (w = floor(train_frac * n))",
    "horizon", "Monte Carlo replications per n", "Sample sizes (grid)",
    "beta0", "beta1", "sigma", "poly_degree", "truth", "selection rule"
  ),
  value = c(
    2026,
    linear_truth_rv_config$window_type,
    linear_truth_rv_config$train_frac,
    linear_truth_rv_config$horizon,
    linear_truth_rv_config$mc_runs,
    "50, 100, ..., 1000 (step 50)",
    linear_truth_rv_config$beta0,
    linear_truth_rv_config$beta1,
    linear_truth_rv_config$sigma,
    linear_truth_rv_config$poly_degree,
    linear_truth_rv_config$truth,
    linear_truth_rv_config$selection_rule
  )
)

expanding_linear_error_rate_summary_tbl <- error_rate_by_n |> arrange(n)
expanding_linear_error_rate_plot <- ggplot(error_rate_by_n, aes(x = n, y = error_rate)) +
  geom_line(linewidth = 0.35) +
  labs(
    x = "Sample size",
    y = "RV selection error rate",
    title = "RV selection error rate (expanding window, linear truth)"
  ) +
  scale_x_continuous(
    limits = c(0, 1000),
    breaks = seq(0, 1000, by = 50)
  ) +
  scale_y_continuous(
    limits = c(0, 0.3),
    expand = expansion(mult = c(0, 0.02))
  ) +
  theme_bw()
