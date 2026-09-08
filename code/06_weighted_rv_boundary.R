library(dplyr)
library(ggplot2)
library(kableExtra)

here::i_am("code/06_weighted_rv_boundary.R")
source(here::here("src", "06_weighted_rv_boundary.R"))

simulation_config <- list(
  seed = 2026L,
  mc_runs = 5000L,
  batch_size = 250L,
  sample_sizes = c(100L, 200L, 500L, 1000L, 2000L, 5000L, 10000L),
  alphas = c(-1, -0.5, -0.1, 0, 0.1, 0.25, 0.5, 1)
)

theory_config <- list(
  sample_sizes = sort(unique(c(
    simulation_config$sample_sizes,
    as.integer(round(exp(seq(log(50), log(1000000), length.out = 120L))))
  ))),
  alphas = simulation_config$alphas,
  negative_limit_alphas = c(-1, -0.5, -0.25, -0.1, -0.05, -0.025),
  positive_limit_alphas = seq(0.01, 1, by = 0.01),
  limit_truncation = 1000000L
)

results_dir <- here::here("data", "final")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
rds_path <- file.path(results_dir, "06_weighted_rv_boundary_results.rds")

required_result_names <- c(
  "simulation_config", "theory_config", "theory_results",
  "mc_results", "limit_results"
)
analysis_results <- if (file.exists(rds_path)) readRDS(rds_path) else NULL
cache_is_current <- !is.null(analysis_results) &&
  all(required_result_names %in% names(analysis_results)) &&
  identical(analysis_results$simulation_config, simulation_config) &&
  identical(analysis_results$theory_config, theory_config)

if (!cache_is_current) {
  validate_weighted_rv_identities()

  theory_results <- run_weighted_rv_theory_grid(
    sample_sizes = theory_config$sample_sizes,
    alphas = theory_config$alphas
  )
  mc_results <- run_weighted_rv_mc(
    sample_sizes = simulation_config$sample_sizes,
    alphas = simulation_config$alphas,
    mc_runs = simulation_config$mc_runs,
    seed = simulation_config$seed,
    batch_size = simulation_config$batch_size
  )

  limit_results <- run_weighted_rv_limit_grid(
    negative_alphas = theory_config$negative_limit_alphas,
    positive_alphas = theory_config$positive_limit_alphas,
    truncation = theory_config$limit_truncation
  )
  coarse_negative_limits <- run_weighted_rv_limit_grid(
    negative_alphas = theory_config$negative_limit_alphas,
    positive_alphas = numeric(0),
    truncation = theory_config$limit_truncation %/% 2L
  ) |>
    filter(alpha < 0) |>
    select(alpha, limit_R_coarse = limit_R)
  limit_results <- limit_results |>
    left_join(coarse_negative_limits, by = "alpha") |>
    mutate(limit_abs_change = abs(limit_R - limit_R_coarse))

  saveRDS(
    list(
      simulation_config = simulation_config,
      theory_config = theory_config,
      theory_results = theory_results,
      mc_results = mc_results,
      limit_results = limit_results
    ),
    rds_path
  )
  message("Saved results to ", rds_path)
} else {
  simulation_config <- analysis_results$simulation_config
  theory_config <- analysis_results$theory_config
  theory_results <- analysis_results$theory_results
  mc_results <- analysis_results$mc_results
  limit_results <- analysis_results$limit_results
  validate_weighted_rv_identities()
  message("Loaded results from ", rds_path)
}

stopifnot(
  nrow(mc_results) == length(simulation_config$sample_sizes) * length(simulation_config$alphas),
  all(is.finite(mc_results$emp_mean)),
  all(is.finite(mc_results$emp_var)),
  all(mc_results$emp_var >= 0),
  all(is.finite(theory_results$theory_R)),
  all(theory_results$theory_var >= 0)
)

alpha_labels <- setNames(
  paste0("alpha = ", simulation_config$alphas),
  as.character(simulation_config$alphas)
)
theory_results <- theory_results |>
  mutate(alpha_lab = factor(alpha_labels[as.character(alpha)], levels = alpha_labels))
mc_results <- mc_results |>
  mutate(alpha_lab = factor(alpha_labels[as.character(alpha)], levels = alpha_labels))

mc_theory_results <- mc_results |>
  left_join(
    theory_results |>
      select(alpha, n, theory_mean, theory_var, theory_R),
    by = c("alpha", "n")
  )

regime_summary_tbl <- data.frame(
  regime = c("$\\alpha<0$", "$\\alpha=0$", "$\\alpha>0$"),
  mean = c(
    "$\\to M_\\infty(\\alpha)\\in(0,\\infty)$",
    "$\\sim\\log n$",
    "$\\sim n^\\alpha/\\alpha$"
  ),
  variance = c(
    "$\\to V_\\infty(\\alpha)\\in(0,\\infty)$",
    "$O(1)$",
    "$\\sim 2n^{2\\alpha}/(1+\\alpha)$"
  ),
  ratio = c(
    "$\\to R_\\infty(\\alpha)>0$",
    "$\\to0$",
    "$\\to2\\alpha^2/(1+\\alpha)$"
  ),
  stringsAsFactors = FALSE
)

simulation_settings_tbl <- data.frame(
  setting = c("Seed", "Monte Carlo replications", "Batch size", "Sample sizes", "Power exponents"),
  value = c(
    simulation_config$seed,
    simulation_config$mc_runs,
    simulation_config$batch_size,
    paste(simulation_config$sample_sizes, collapse = ", "),
    paste(simulation_config$alphas, collapse = ", ")
  ),
  stringsAsFactors = FALSE
)

limit_check_tbl <- limit_results |>
  filter(alpha < 0) |>
  select(alpha, limit_R, limit_R_coarse, limit_abs_change)

weighted_rv_R_plot <- ggplot(
  theory_results,
  aes(x = n, y = theory_R, color = alpha_lab, group = alpha_lab)
) +
  geom_line(linewidth = 0.55) +
  geom_point(
    data = mc_theory_results,
    aes(x = n, y = emp_R, color = alpha_lab),
    inherit.aes = FALSE,
    size = 1.25,
    alpha = 0.75
  ) +
  scale_x_log10() +
  labs(
    x = "Sample size n (log scale)",
    y = expression(R[n](alpha)),
    color = NULL
  ) +
  theme_bw() +
  theme(legend.position = "bottom", legend.text = element_text(size = 7))

mean_comparison <- mc_theory_results |>
  transmute(alpha_lab, n, empirical = emp_mean, theoretical = theory_mean)
variance_comparison <- mc_theory_results |>
  transmute(alpha_lab, n, empirical = emp_var, theoretical = theory_var)
comparison_long <- rbind(
  data.frame(
    alpha_lab = mean_comparison$alpha_lab,
    n = mean_comparison$n,
    empirical = mean_comparison$empirical,
    theoretical = mean_comparison$theoretical,
    quantity = "Mean"
  ),
  data.frame(
    alpha_lab = variance_comparison$alpha_lab,
    n = variance_comparison$n,
    empirical = variance_comparison$empirical,
    theoretical = variance_comparison$theoretical,
    quantity = "Variance"
  )
)

weighted_rv_moment_plot <- ggplot(comparison_long, aes(x = n, color = alpha_lab)) +
  geom_line(aes(y = theoretical, group = alpha_lab), linewidth = 0.45) +
  geom_point(aes(y = empirical), size = 1.1, alpha = 0.75) +
  facet_wrap(~quantity, scales = "free_y", ncol = 1) +
  scale_x_log10(breaks = simulation_config$sample_sizes) +
  scale_y_log10() +
  labs(
    x = "Sample size n (log scale)",
    y = "Value (log scale)",
    color = NULL
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    legend.text = element_text(size = 7)
  )

weighted_rv_pneg_plot <- ggplot(
  mc_results,
  aes(x = n, y = p_neg, color = alpha_lab, group = alpha_lab)
) +
  geom_errorbar(
    aes(
      ymin = pmax(0, p_neg - 1.96 * p_neg_mc_se),
      ymax = pmin(1, p_neg + 1.96 * p_neg_mc_se)
    ),
    width = 0,
    linewidth = 0.25,
    alpha = 0.65
  ) +
  geom_line(linewidth = 0.45) +
  geom_point(size = 1.2) +
  scale_x_log10(breaks = simulation_config$sample_sizes) +
  expand_limits(y = 0) +
  labs(
    x = "Sample size n (log scale)",
    y = expression(hat(P)(Xi[n]^((alpha)) < 0)),
    color = NULL
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    legend.text = element_text(size = 7)
  )

floor_plot_data <- limit_results |>
  mutate(source = ifelse(alpha < 0, "Numerical infinite-series limit", "Closed-form limit"))

weighted_rv_floor_plot <- ggplot(floor_plot_data, aes(x = alpha, y = limit_R)) +
  geom_line(
    data = filter(floor_plot_data, alpha >= 0),
    color = "black",
    linewidth = 0.6
  ) +
  geom_point(aes(shape = source), size = 1.7) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.35) +
  scale_shape_manual(values = c("Closed-form limit" = 16, "Numerical infinite-series limit" = 1)) +
  labs(
    x = expression(alpha),
    y = expression(R[infinity](alpha)),
    shape = NULL
  ) +
  theme_bw() +
  theme(legend.position = "bottom")
