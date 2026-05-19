library(dplyr)
library(ggplot2)

here::i_am("code/03_sliding_poly.R")
source(here::here("src", "rolling_error_rate.R"))

#Main simulation configuration
poly_truth_rv_config <- list(
  window_type = "sliding",
  train_frac = 0.5,
  horizon = 1,
  mc_runs = 1000,
  sample_sizes = seq(50, 1000, by = 50),
  beta0 = 1,
  beta1 = 2,
  beta2 = 1,
  beta3 = 1,
  sigma = 4,
  poly_degree = 3,
  truth_type = "polynomial",
  truth = "polynomial: y = beta0 + beta1*x + beta2*x^2 + beta3*x^3 + eps",
  selection_rule = "Choose model with lower mean RV MSE; tie -> linear"
)

results_dir <- here::here("data", "final")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
rds_path <- file.path(results_dir, "03_sliding_poly_results.rds")

sliding_poly_results <- if (file.exists(rds_path)) readRDS(rds_path) else NULL
if (!is.null(sliding_poly_results) &&
    all(c("rv_config", "error_rate_by_n") %in% names(sliding_poly_results))) {
  poly_truth_rv_config <- sliding_poly_results$rv_config
  error_rate_by_n <- sliding_poly_results$error_rate_by_n
  message("Loaded results from ", rds_path)
} else {
  set.seed(2026)
  error_rate_by_n <- run_rolling_selection_error_rate_sim(poly_truth_rv_config)
  saveRDS(
    list(rv_config = poly_truth_rv_config, error_rate_by_n = error_rate_by_n),
    rds_path
  )
  message("Saved results to ", rds_path)
}

#Graphing and tables
sliding_poly_settings_tbl <- data.frame(
  setting = c(
    "seed", "window type", "train_frac (w = floor(train_frac * n))",
    "horizon", "Monte Carlo replications per n", "Sample sizes (grid)",
    "beta0", "beta1", "beta2", "beta3", "sigma", "poly_degree", "truth", "selection rule"
  ),
  value = c(
    2026,
    poly_truth_rv_config$window_type,
    poly_truth_rv_config$train_frac,
    poly_truth_rv_config$horizon,
    poly_truth_rv_config$mc_runs,
    "50, 100, ..., 1000 (step 50)",
    poly_truth_rv_config$beta0,
    poly_truth_rv_config$beta1,
    poly_truth_rv_config$beta2,
    poly_truth_rv_config$beta3,
    poly_truth_rv_config$sigma,
    poly_truth_rv_config$poly_degree,
    poly_truth_rv_config$truth,
    poly_truth_rv_config$selection_rule
  )
)

sliding_poly_error_rate_summary_tbl <- error_rate_by_n |> arrange(n)
sliding_poly_error_rate_plot <- ggplot(error_rate_by_n, aes(x = n, y = error_rate)) +
  geom_line(linewidth = 0.35) +
  labs(
    x = "Sample size",
    y = "RV selection error rate",
    title = "RV selection error rate (sliding window, polynomial truth)"
  ) +
  scale_x_continuous(
    limits = c(0, 1000),
    breaks = seq(0, 1000, by = 50)
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2),
    expand = expansion(mult = c(0, 0.02))
  ) +
  theme_bw()
