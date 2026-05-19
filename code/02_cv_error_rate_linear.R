library(dplyr)
library(ggplot2)

here::i_am("code/02_cv_error_rate_linear.R")
source(here::here("src", "cv_error_rate.R"))

#Main simulation configuration
linear_truth_sim_config <- list(
  K = 10,
  mc_runs = 1000,
  sample_sizes = seq(50, 1000, by = 50),
  beta0 = 1,
  beta1 = 2,
  sigma = 4,
  poly_degree = 3,
  truth_type = "linear",
  truth = "linear: y = beta0 + beta1*x + eps",
  selection_rule = "Choose poly3 if mean CV MSE strictly lower than linear; else linear"
)

results_dir <- here::here("data", "final")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
rds_path <- file.path(results_dir, "02_cv_error_rate_results.rds")

cv_error_rate_linear_results <- if (file.exists(rds_path)) readRDS(rds_path) else NULL
if (!is.null(cv_error_rate_linear_results) &&
    all(c("linear_truth_config", "linear_truth_error_rate_by_n") %in% names(cv_error_rate_linear_results))) {
  linear_truth_sim_config <- cv_error_rate_linear_results$linear_truth_config
  linear_truth_error_rate_by_n <- cv_error_rate_linear_results$linear_truth_error_rate_by_n
  message("Loaded results from ", rds_path)
} else {
  set.seed(2026)
  linear_truth_error_rate_by_n <- run_cv_selection_error_rate_sim(linear_truth_sim_config)
  cv_error_rate_linear_results <- list(
    linear_truth_config = linear_truth_sim_config,
    linear_truth_error_rate_by_n = linear_truth_error_rate_by_n
  )
  saveRDS(cv_error_rate_linear_results, rds_path)
  message("Saved results to ", rds_path)
}

#Graphing and tables
settings_tbl <- data.frame(
  setting = c(
    "seed", "K", "Monte Carlo replications per n", "Sample sizes (grid)",
    "beta0", "beta1", "sigma", "poly_degree", "truth", "selection rule"
  ),
  value = c(
    2026,
    linear_truth_sim_config$K,
    linear_truth_sim_config$mc_runs,
    "50, 100, ..., 10000 (step 50)",
    linear_truth_sim_config$beta0,
    linear_truth_sim_config$beta1,
    linear_truth_sim_config$sigma,
    linear_truth_sim_config$poly_degree,
    linear_truth_sim_config$truth,
    linear_truth_sim_config$selection_rule
  )
)

error_rate_summary_tbl <- linear_truth_error_rate_by_n |> arrange(n)
error_rate_plot <- ggplot(linear_truth_error_rate_by_n, aes(x = n, y = error_rate)) +
  geom_line(linewidth = 0.35) +
  labs(
    x = "Sample size",
    y = "CV selection error rate",
    title = "Rate CV picks polynomial under linear truth"
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
