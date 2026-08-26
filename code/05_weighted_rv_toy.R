library(dplyr)
library(ggplot2)
library(kableExtra)

here::i_am("code/05_weighted_rv_toy.R")
source(here::here("src", "weighted_rv_toy.R"))

# simulation settings (broad experiment)
simulation_config <- list(
  mc_runs = 1000,
  sample_sizes = c(100, 200, 500, 1000, 2000, 5000, 10000),
  alphas = c(0, -1, -0.5, 0.5, 1)
)

# near-zero alpha experiment (same MC and n grid)
near_zero_config <- list(
  mc_runs = simulation_config$mc_runs,
  sample_sizes = simulation_config$sample_sizes,
  alphas = c(-0.5, -0.25, -0.1, 0, 0.1, 0.25, 0.5)
)

results_dir <- here::here("data", "final")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
rds_path <- file.path(results_dir, "05_weighted_rv_toy_results.rds")
rds_near_path <- file.path(results_dir, "05_weighted_rv_toy_near_results.rds")

# broad power-weight grid simulation
weighted_rv_results <- if (file.exists(rds_path)) readRDS(rds_path) else NULL

if (is.null(weighted_rv_results) ||
    !all(c("simulation_config", "results_by_n_alpha") %in% names(weighted_rv_results))) {
  message("Running power-weight grid (alpha = 0 control first in alphas)...")
  results_by_n_alpha <- run_weight_grid(
    sample_sizes = simulation_config$sample_sizes,
    alphas = simulation_config$alphas,
    mc_runs = simulation_config$mc_runs
  )

  saveRDS(
    list(
      simulation_config = simulation_config,
      results_by_n_alpha = results_by_n_alpha
    ),
    rds_path
  )
  message("Saved results to ", rds_path)
} else {
  simulation_config <- weighted_rv_results$simulation_config
  results_by_n_alpha <- weighted_rv_results$results_by_n_alpha
  message("Loaded results from ", rds_path)
}

# near-zero alpha grid simulation
near_results <- if (file.exists(rds_near_path)) readRDS(rds_near_path) else NULL

if (is.null(near_results) ||
    !all(c("near_zero_config", "results_near_by_n_alpha") %in% names(near_results))) {
  message("Running near-zero alpha grid...")
  results_near_by_n_alpha <- run_weight_grid(
    sample_sizes = near_zero_config$sample_sizes,
    alphas = near_zero_config$alphas,
    mc_runs = near_zero_config$mc_runs
  )
  saveRDS(
    list(
      near_zero_config = near_zero_config,
      results_near_by_n_alpha = results_near_by_n_alpha
    ),
    rds_near_path
  )
  message("Saved near-zero results to ", rds_near_path)
} else {
  near_zero_config <- near_results$near_zero_config
  results_near_by_n_alpha <- near_results$results_near_by_n_alpha
  message("Loaded near-zero results from ", rds_near_path)
}

# ---- tables ----
weighted_rv_simulation_settings_tbl <- data.frame(
  setting = c(
    "seed",
    "Monte Carlo replications",
    "Sample sizes (grid)",
    "weight family",
    "broad alphas (0 = control first)",
    "near-zero alphas",
    "incorrect selection"
  ),
  value = c(
    2026,
    simulation_config$mc_runs,
    paste(simulation_config$sample_sizes, collapse = ", "),
    "w_i = i^alpha",
    paste(simulation_config$alphas, collapse = ", "),
    paste(near_zero_config$alphas, collapse = ", "),
    "Xi_n < 0"
  )
)

control_results_tbl <- results_by_n_alpha |>
  filter(alpha == 0) |>
  arrange(n) |>
  transmute(
    n,
    p_neg,
    mean_xi,
    theory_mean,
    log_n = log(n),
    var_xi,
    R
  )

weighted_rv_results_tbl <- results_by_n_alpha |>
  arrange(alpha, n)

near_zero_results_tbl <- results_near_by_n_alpha |>
  arrange(alpha, n)

# simulation interpretation tables
weighted_rv_summary_tbl <- data.frame(
  alpha = c(-1, -0.5, 0, 0.5, 1),
  weight = c("i^-1", "i^-0.5", "1", "sqrt(i)", "i"),
  expected_signal_order = c(
    "bounded O(1) (derivation)",
    "bounded O(1) (derivation)",
    "~ log(n) (derivation)",
    "~ 2 sqrt(n) (derivation)",
    "~ n (derivation)"
  ),
  error_probability_behavior = c(
    "approximately flat near 0.08 (simulation)",
    "decreases then levels near 0.01 (simulation)",
    "decreases toward 0 (simulation)",
    "noisy and remains positive (simulation)",
    "remains around 0.1 (simulation)"
  ),
  R_behavior = c(
    "no clear tendency toward 0 (simulation)",
    "no clear tendency toward 0 (simulation)",
    "decreases with n (simulation)",
    "no clear tendency toward 0 (simulation)",
    "no clear tendency toward 0 (simulation)"
  ),
  interpretation = c(
    "appears bounded away from 0 (simulation interpretation)",
    "decreases initially, then appears to level off (simulation interpretation)",
    "clear evidence toward 0 (simulation interpretation)",
    "unclear; no clear convergence to 0 (simulation interpretation)",
    "appears bounded away from 0 (simulation interpretation)"
  ),
  stringsAsFactors = FALSE
)

# Near-zero alpha summary (hand-written; derivation + simulation reading)
weighted_rv_near_summary_tbl <- data.frame(
  alpha = c(-0.5, -0.25, -0.1, 0, 0.1, 0.25, 0.5),
  weight = c("i^-0.5", "i^-0.25", "i^-0.1", "1", "i^0.1", "i^0.25", "sqrt(i)"),
  expected_signal_order = c(
    "bounded O(1) (derivation)",
    "bounded O(1) (derivation)",
    "bounded O(1) (derivation)",
    "~ log(n) (derivation)",
    "~ c n^0.1 (derivation)",
    "~ c n^0.25 (derivation)",
    "~ 2 sqrt(n) (derivation)"
  ),
  error_probability_behavior = c(
    "decreases then levels near 0.02 (simulation)",
    "decreases toward ~0 (simulation)",
    "decreases to ~0 (simulation)",
    "decreases toward 0 (simulation)",
    "decreases toward ~0 (simulation)",
    "decreases then stays noisy near 0.02 (simulation)",
    "noisy and remains positive (simulation)"
  ),
  R_behavior = c(
    "slow decrease; stays large (simulation)",
    "decreases with n (simulation)",
    "decreases with n (simulation)",
    "decreases with n (simulation)",
    "decreases with n (simulation)",
    "mild decrease; noisy (simulation)",
    "no clear tendency toward 0 (simulation)"
  ),
  interpretation = c(
    "levels off away from 0 (simulation interpretation)",
    "looks similar to control at these n (simulation interpretation)",
    "looks similar to control (simulation interpretation)",
    "clear evidence toward 0 (simulation interpretation)",
    "looks similar to control (simulation interpretation)",
    "weaker / unclear vs alpha = 0 (simulation interpretation)",
    "no clear convergence to 0 (simulation interpretation)"
  ),
  stringsAsFactors = FALSE
)

# plots (broad)
results_by_n_alpha$alpha_lab <- factor(
  results_by_n_alpha$alpha,
  levels = simulation_config$alphas
)

n_breaks <- simulation_config$sample_sizes

make_weighted_rv_plot <- function(df, y_aes, ylab, title, log_y = FALSE) {
  p <- ggplot(df, aes(x = n, y = .data[[y_aes]], color = alpha_lab, group = alpha_lab)) +
    geom_line(linewidth = 0.4) +
    geom_point(size = 1.2) +
    scale_x_log10(breaks = n_breaks) +
    labs(
      x = "Sample size n (log scale)",
      y = ylab,
      color = "alpha",
      title = title
    ) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  if (isTRUE(log_y)) {
    p <- p + scale_y_log10()
  } else {
    p <- p + expand_limits(y = 0)
  }
  p
}

weighted_rv_pneg_plot <- make_weighted_rv_plot(
  results_by_n_alpha, "p_neg",
  "P(Xi_n < 0)",
  "Incorrect-selection probability vs n (broad alphas)"
)

weighted_rv_mean_plot <- ggplot(results_by_n_alpha, aes(x = n, color = alpha_lab, group = alpha_lab)) +
  geom_line(aes(y = mean_xi), linewidth = 0.4) +
  geom_point(aes(y = mean_xi), size = 1.2) +
  geom_line(aes(y = theory_mean), linetype = "dashed", linewidth = 0.35) +
  scale_x_log10(breaks = n_breaks) +
  expand_limits(y = 0) +
  labs(
    x = "Sample size n (log scale)",
    y = "Mean of Xi_n",
    color = "alpha",
    title = "Empirical mean (solid) and theoretical mean (dashed)"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

weighted_rv_var_plot <- make_weighted_rv_plot(
  results_by_n_alpha, "var_xi",
  "Var(Xi_n) (log scale)",
  "Empirical variance vs n",
  log_y = TRUE
)

weighted_rv_R_plot <- make_weighted_rv_plot(
  results_by_n_alpha, "R",
  "R_n = Var / (mean)^2",
  "Chebyshev ratio R_n vs n"
)

control_pneg_plot <- ggplot(
  control_results_tbl,
  aes(x = n, y = p_neg)
) +
  geom_line(linewidth = 0.4) +
  geom_point(size = 1.5) +
  scale_x_log10(breaks = n_breaks) +
  expand_limits(y = 0) +
  labs(
    x = "Sample size n (log scale)",
    y = "P(Xi_n < 0)",
    title = "Control (alpha = 0): incorrect-selection probability"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# plots (near-zero)
results_near_by_n_alpha$alpha_lab <- factor(
  results_near_by_n_alpha$alpha,
  levels = near_zero_config$alphas
)

weighted_rv_near_pneg_plot <- ggplot(
  results_near_by_n_alpha,
  aes(x = n, y = p_neg, color = alpha_lab, group = alpha_lab)
) +
  geom_line(linewidth = 0.4) +
  geom_point(size = 1.2) +
  scale_x_log10(breaks = near_zero_config$sample_sizes) +
  expand_limits(y = 0) +
  labs(
    x = "Sample size n (log scale)",
    y = "P(Xi_n < 0)",
    color = "alpha",
    title = "Near-zero alpha experiment: incorrect-selection probability"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

near_zero_settings_tbl <- data.frame(
  setting = c(
    "seed",
    "Monte Carlo replications",
    "Sample sizes (grid)",
    "alphas (near zero)"
  ),
  value = c(
    2026,
    near_zero_config$mc_runs,
    paste(near_zero_config$sample_sizes, collapse = ", "),
    paste(near_zero_config$alphas, collapse = ", ")
  )
)
