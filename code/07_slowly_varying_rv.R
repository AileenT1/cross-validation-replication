library(dplyr)
library(ggplot2)
library(kableExtra)
here::i_am("code/07_slowly_varying_rv.R")
source(here::here("src", "07_slowly_varying_rv.R"))

registry07 <- sv_registry()
config07 <- list(
  version = 1L, seed = 2026L, mc_runs = 5000L, batch_size = 250L,
  mc_n = c(100L, 200L, 500L, 1000L, 2000L, 5000L, 10000L, 20000L, 50000L, 100000L),
  theory_n = NULL, limit_n = c(100000L, 1000000L, 10000000L)
)
config07$theory_n <- sort(unique(c(config07$mc_n, config07$limit_n,
  as.integer(round(exp(seq(log(50), log(1e7), length.out = 160L)))))))
cache07_path <- here::here("data", "final", "07_slowly_varying_rv_results.rds")
results07 <- if (file.exists(cache07_path)) readRDS(cache07_path) else NULL
cache07_valid <- !is.null(results07) &&
  isTRUE(tryCatch({
    sv_validate_cache(results07, config07, registry07)
    TRUE
  }, error = function(e) FALSE))

if (!cache07_valid) {
  validation07 <- sv_validate()
  message("Calculating exact trajectories through n = ", max(config07$theory_n))
  exact07 <- sv_theory(config07$theory_n, registry07)
  mc07 <- sv_mc(config07$mc_n, registry07, config07$mc_runs,
                 config07$batch_size, config07$seed)
  limits07 <- sv_limits(config07$limit_n, registry07)
  results07 <- list(config = config07, registry = registry07, exact = exact07,
                    mc = mc07, limits = limits07, validation = validation07,
                    session = capture.output(sessionInfo()))
  sv_validate_cache(results07, config07, registry07)
  dir.create(dirname(cache07_path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(results07, cache07_path)
  message("Saved Analysis 07 results.")
} else {
  message("Using validated Analysis 07 cache; no Monte Carlo rerun.")
}

exact07 <- results07$exact |> left_join(registry07, by = "weight_id")
mc07 <- results07$mc |>
  left_join(results07$exact, by = c("weight_id", "n")) |>
  left_join(registry07, by = "weight_id") |>
  mutate(ratio_exact_mean_mc = variance_mc / mean^2)
limits07 <- results07$limits
main_labels07 <- registry07$label[registry07$main]
exact07$label <- factor(exact07$label, levels = registry07$label)
mc07$label <- factor(mc07$label, levels = registry07$label)
style07 <- theme_bw(base_size = 10) +
  theme(legend.position = "bottom", strip.text = element_text(size = 9),
        panel.grid.minor = element_blank(), plot.caption = element_text(hjust = 0))

ratio_plot07 <- ggplot(filter(exact07, main), aes(n, ratio)) +
  geom_line(linewidth = 0.5, colour = "#2166AC") +
  geom_point(data = filter(mc07, main), aes(y = ratio_exact_mean_mc),
             size = 1.1, colour = "#B35806") +
  facet_wrap(~label, ncol = 2, scales = "free_y", drop = TRUE) +
  scale_x_log10(breaks = c(1e2, 1e4, 1e6)) +
  labs(x = "Sample size n (log scale)", y = "Chebyshev ratio",
       caption = "Blue: exact ratio. Orange: empirical variance / exact mean squared.\nPanel-specific y scales; the ratio is not an error probability.") + style07

probability_plot07 <- ggplot(filter(mc07, main), aes(n, p_wrong)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.07,
                colour = "#636363") +
  geom_line(linewidth = 0.4, colour = "#2166AC") +
  geom_point(size = 1, colour = "#2166AC") +
  facet_wrap(~label, ncol = 2, drop = TRUE) +
  scale_x_log10(breaks = c(1e2, 1e3, 1e4, 1e5)) +
  labs(x = "Sample size n (log scale)", y = "Estimated P(score < 0)",
       caption = "Pointwise 95% Wilson intervals from 5,000 replications.\nPaths are shared across weights and sample sizes; plotted estimates are dependent.") +
  style07

moments_data07 <- bind_rows(
  filter(exact07, main, n %in% config07$mc_n) |>
    transmute(weight_id, n, label, moment = "Mean", exact_value = mean),
  filter(exact07, main, n %in% config07$mc_n) |>
    transmute(weight_id, n, label, moment = "Variance", exact_value = variance))
moments_mc07 <- bind_rows(
  filter(mc07, main) |> transmute(weight_id, n, label, moment = "Mean", empirical = mean_mc),
  filter(mc07, main) |> transmute(weight_id, n, label, moment = "Variance", empirical = variance_mc))
moment_plot07 <- ggplot(moments_data07, aes(n, exact_value, colour = label)) +
  geom_line(linewidth = 0.55) +
  geom_point(data = moments_mc07, aes(y = empirical), size = 1.2) +
  facet_wrap(~moment, scales = "free_y") + scale_x_log10() + scale_y_log10() +
  labs(x = "Sample size n (log scale)", y = "Moment (log scale)", colour = "Weight",
       caption = "Lines: exact moments. Points: Monte Carlo moments.") +
  style07 + guides(colour = guide_legend(ncol = 2))

family_plot07 <- ggplot(filter(exact07, family == "log_power"),
                        aes(n, ratio, colour = factor(p))) +
  geom_line(linewidth = 0.65) + scale_x_log10() + scale_y_log10() +
  labs(x = "Sample size n (log scale)", y = "Exact ratio (log scale)", colour = "p",
       caption = "The ratio tends to zero for p >= -1; it has a positive limit for p < -1.\nSlow finite-sample change near p = -1 does not establish a plateau.") + style07

reference_plot07 <- ggplot(filter(exact07, weight_id %in%
  c("logp_1", "logp_2", "loglog", "power001")), aes(n, ratio, colour = label)) +
  geom_line(linewidth = 0.65) +
  geom_hline(yintercept = 2 * 0.01^2 / 1.01, linetype = "dashed", colour = "#555555") +
  scale_x_log10() + scale_y_log10() +
  labs(x = "Sample size n (log scale)", y = "Exact ratio (log scale)", colour = "Weight",
       caption = "Dashed: the positive-power limit 2(0.01)^2 / 1.01, not a finite-n approximation.\nThe log-based ratios tend to zero, despite being larger over the displayed range.") +
  style07 + guides(colour = guide_legend(ncol = 2))

endpoint_table07 <- filter(mc07, main, n == max(config07$mc_n)) |>
  transmute(Weight = as.character(label), Mean = mean, Variance = variance,
            Ratio = ratio, Error = p_wrong, MCSE = mcse, Lower = lower, Upper = upper)
limit_table07 <- limits07 |> filter(p < -1) |>
  transmute(p, N = n, M_est = mean_limit_est, V_est = variance_limit_est,
            R_est = ratio_limit_est)
