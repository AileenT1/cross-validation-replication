library(dplyr)
library(ggplot2)

here::i_am("code/04_rv_fixed_poly.R")
source(here::here("src", "rolling_error_rate.R"))

#code to run the simulation
simulation_config <- list(
  mc_runs = 500,
  sample_sizes = seq(100, 1000, by = 100),
  truth = list(
    type = "polynomial",
    beta = c(beta0 = 1, beta1 = 2, beta2 = 1, beta3 = 1)
  ),
  noise = list(
    type = "gaussian",
    sigma = 4
  )
)
model_fitting_config <- list(
  candidate_models = list(
    linear = as.formula("y ~ x"),
    polynomial = as.formula("y ~ poly(x, 3, raw = TRUE)")
  ),
  rv = list(
    window_type = "fixed",
    window_size = 50,
    horizon = 1,
    validation_step = 10,
    selection_rule = "minimum average validation error"
  )
)

results_dir <- here::here("data", "final")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
rds_path <- file.path(results_dir, "04_rv_fixed_poly_results.rds")

fixed_poly_results <- if (file.exists(rds_path)) readRDS(rds_path) else NULL
if (is.null(fixed_poly_results) ||
    !all(c("simulation_config", "model_fitting_config", "error_rate_by_n") %in% names(fixed_poly_results))) {
  set.seed(2026)
  fixed_poly_error_rate_by_n <- run_rolling_selection_error_rate_sim(
    simulation_config,
    model_fitting_config
  )
  saveRDS(
    list(
      simulation_config = simulation_config,
      model_fitting_config = model_fitting_config,
      error_rate_by_n = fixed_poly_error_rate_by_n
    ),
    rds_path
  )
  message("Saved results to ", rds_path)
} else {
  simulation_config <- fixed_poly_results$simulation_config
  model_fitting_config <- fixed_poly_results$model_fitting_config
  fixed_poly_error_rate_by_n <- fixed_poly_results$error_rate_by_n
  message("Loaded results from ", rds_path)
}

#tables&graphs
fixed_poly_simulation_settings_tbl <- data.frame(
  setting = c(
    "seed",
    "Monte Carlo replications per n",
    "Sample sizes (grid)",
    "truth type",
    "beta0",
    "beta1",
    "beta2",
    "beta3",
    "noise type",
    "sigma"
  ),
  value = c(
    2026,
    simulation_config$mc_runs,
    paste0(
      min(simulation_config$sample_sizes), ", ",
      min(simulation_config$sample_sizes) + 50, ", ..., ",
      max(simulation_config$sample_sizes),
      " (step 50)"
    ),
    simulation_config$truth$type,
    simulation_config$truth$beta["beta0"],
    simulation_config$truth$beta["beta1"],
    simulation_config$truth$beta["beta2"],
    simulation_config$truth$beta["beta3"],
    simulation_config$noise$type,
    simulation_config$noise$sigma
  )
)

fixed_poly_error_rate_summary_tbl <- fixed_poly_error_rate_by_n |> arrange(n)
fixed_poly_error_rate_plot <- ggplot(fixed_poly_error_rate_by_n, aes(x = n, y = error_rate)) +
  geom_line(linewidth = 0.35) +
  labs(
    x = "Sample size",
    y = "RV selection error rate",
    title = "RV selection error rate (fixed window w = 50, polynomial truth)"
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
