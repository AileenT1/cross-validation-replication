library(dplyr)
library(ggplot2)

here::i_am("code/02_cv_error_rate_poly.R")
source(here::here("src", "cv_error_rate.R"))

#Main simulation configuration
poly_truth_sim_config <- list(
  K = 10,
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
  selection_rule = "Choose poly3 if mean CV MSE strictly lower than linear; else linear"
)

results_dir <- here::here("data", "final")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
rds_path <- file.path(results_dir, "02_cv_error_rate_poly_results.rds")

cv_error_rate_poly_results <- if (file.exists(rds_path)) readRDS(rds_path) else NULL
if (!is.null(cv_error_rate_poly_results) &&
    all(c("poly_truth_config", "poly_truth_error_rate_by_n") %in% names(cv_error_rate_poly_results))) {
  poly_truth_sim_config <- cv_error_rate_poly_results$poly_truth_config
  poly_truth_error_rate_by_n <- cv_error_rate_poly_results$poly_truth_error_rate_by_n
  message("Loaded results from ", rds_path)
} else {
  set.seed(2026)
  poly_truth_error_rate_by_n <- run_cv_selection_error_rate_sim(poly_truth_sim_config)
  cv_error_rate_poly_results <- list(
    poly_truth_config = poly_truth_sim_config,
    poly_truth_error_rate_by_n = poly_truth_error_rate_by_n
  )
  saveRDS(cv_error_rate_poly_results, rds_path)
  message("Saved results to ", rds_path)
}

#Graphing and tables
settings_tbl <- data.frame(
  setting = c(
    "seed", "K", "Monte Carlo replications per n", "Sample sizes (grid)",
    "beta0", "beta1", "beta2", "beta3", "sigma", "poly_degree", "truth", "selection rule"
  ),
  value = c(
    2026,
    poly_truth_sim_config$K,
    poly_truth_sim_config$mc_runs,
    "50, 100, ..., 1000 (step 50)",
    poly_truth_sim_config$beta0,
    poly_truth_sim_config$beta1,
    poly_truth_sim_config$beta2,
    poly_truth_sim_config$beta3,
    poly_truth_sim_config$sigma,
    poly_truth_sim_config$poly_degree,
    poly_truth_sim_config$truth,
    poly_truth_sim_config$selection_rule
  )
)

error_rate_summary_tbl <- poly_truth_error_rate_by_n |> arrange(n)
error_rate_plot <- ggplot(poly_truth_error_rate_by_n, aes(x = n, y = error_rate)) +
  geom_line(linewidth = 0.35) +
  labs(
    x = "Sample size",
    y = "CV selection error rate",
    title = "Rate CV picks linear under polynomial truth"
  ) +
  scale_x_continuous(
    limits = c(0, 1000),
    breaks = seq(0, 1000, by = 50)
  ) +
  scale_y_continuous(
    limits = c(0, 0.2),
    expand = expansion(mult = c(0, 0.02))
  ) +
  theme_bw()
