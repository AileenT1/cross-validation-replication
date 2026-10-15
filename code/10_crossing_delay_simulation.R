library(ggplot2)
here::i_am("code/10_crossing_delay_simulation.R")
source(here::here("src", "10_crossing_delay.R"), local = TRUE)

config10 <- list(
  alphas = c(0.5, 1, 2, 4),
  inverse_theta = c(5L, 10L, 20L, 40L, 80L, 160L),
  factor_n = sort(unique(as.integer(round(exp(seq(log(20), log(10000), length.out = 100L)))))),
  mc_alpha = 1,
  sample_sizes = c(200L, 500L, 1000L, 2000L),
  c_values = c(0.5, 1, 1.5, 1.8, 2.2, 2.5, 3),
  replications = 5000L, bootstraps = 400L, seed = 20261014L
)

cache10 <- here::here("data", "final", "10_crossing_delay_results.rds")
signature10 <- unname(tools::md5sum(c(
  here::here("src", "10_crossing_delay.R"),
  here::here("code", "10_crossing_delay_simulation.R"))))
results10 <- if (file.exists(cache10)) readRDS(cache10) else NULL
cache_current10 <- !is.null(results10) && identical(results10$config, config10) &&
  identical(results10$signature, signature10)
if (!cache_current10) {
  message("Running Analysis 10: exact calculations and 5,000 Gaussian Monte Carlo paths.")
  start10 <- proc.time()[["elapsed"]]
  validation10 <- ten_validate()
  deterministic10 <- ten_deterministic(config10)
  mc10 <- ten_mc(config10)
  results10 <- list(config = config10, signature = signature10,
    crossing = deterministic10$crossing, factor = deterministic10$factor,
    mc = mc10, validation_error = validation10,
    elapsed_seconds = proc.time()[["elapsed"]] - start10)
  dir.create(dirname(cache10), recursive = TRUE, showWarnings = FALSE)
  saveRDS(results10, cache10)
} else message("Using cached Analysis 10 simulation results (matching code and configuration).")
stopifnot(nrow(results10$crossing) == 24L, nrow(results10$mc) == 28L,
          results10$validation_error < 1e-10)

style10 <- theme_bw(base_size = 10) + theme(
  legend.position = "bottom", legend.title = element_blank(),
  panel.grid.minor = element_blank(), strip.text = element_text(size = 9))
label10 <- function(alpha) paste0("alpha = ", alpha, "; lambda = ", 1 + 1 / alpha)
crossing10 <- results10$crossing
crossing10$label <- factor(label10(crossing10$alpha), levels = label10(config10$alphas))
factor10 <- results10$factor
factor10$label <- factor(label10(factor10$alpha), levels = label10(config10$alphas))
mc10 <- results10$mc
mc10$panel <- factor(paste("n =", mc10$n), levels = paste("n =", config10$sample_sizes))
mc10$side <- ifelse(mc10$c_value < 2, "before", "after")

crossings_plot10 <- ggplot(crossing10, aes(-log10(theta), ratio)) +
  geom_hline(aes(yintercept = lambda), linetype = 2, colour = "#B35806") +
  geom_line(colour = "#2166AC", linewidth = 0.6) +
  geom_point(colour = "#2166AC", size = 1.4) +
  scale_x_continuous(breaks = -log10(c(0.2, 0.05, 0.0125)),
    labels = c("0.2", "0.05", "0.0125")) +
  facet_wrap(~label, ncol = 2, scales = "free_y") +
  labs(x = "|theta| (decreasing)", y = "Expected RV / population crossing") + style10

factor_plot10 <- ggplot(factor10, aes(n, factor)) +
  geom_hline(aes(yintercept = target), linetype = 2, colour = "#B35806") +
  geom_line(colour = "#2166AC", linewidth = 0.6) + scale_x_log10() +
  facet_wrap(~label, ncol = 2, scales = "free_y") +
  labs(x = "Training size n", y = "Square-weight factor") + style10

curve10 <- expand.grid(c_value = c(seq(0.5, 1.8, length.out = 150L),
                                  seq(2.2, 3, length.out = 150L)),
                       panel = levels(mc10$panel))
curve10$panel <- factor(curve10$panel, levels = levels(mc10$panel))
curve10$side <- ifelse(curve10$c_value < 2, "before", "after")
curve10$ratio_limit <- (4 + 16 * curve10$c_value / 3) / (2 - curve10$c_value)^2
ratio_plot10 <- ggplot(mc10, aes(c_value)) +
  geom_vline(xintercept = 2, linetype = 3, colour = "grey50") +
  geom_line(data = curve10, aes(y = ratio_limit, group = side, colour = "Asymptotic limit"), linewidth = 0.55) +
  geom_line(aes(y = ratio_exact, group = side, colour = "Exact finite-n"), linetype = 2, linewidth = 0.5) +
  geom_point(aes(y = ratio_exact, colour = "Exact finite-n"), shape = 15, size = 1.4) +
  geom_errorbar(aes(ymin = low, ymax = high, colour = "Monte Carlo (95% bootstrap)"),
                width = 0.035, linewidth = 0.35) +
  geom_point(aes(y = ratio_mc, colour = "Monte Carlo (95% bootstrap)"), size = 1.7) +
  scale_colour_manual(values = c("Asymptotic limit" = "#222222",
    "Exact finite-n" = "#2166AC", "Monte Carlo (95% bootstrap)" = "#B35806")) +
  scale_y_log10() + scale_x_continuous(breaks = c(0.5, 1, 1.5, 2, 2.5, 3)) +
  facet_wrap(~panel, ncol = 2) +
  labs(x = "c = n * theta^2 (weighted crossing at c = 2)", y = "Chebyshev ratio (log scale)") +
  style10 + theme(legend.text = element_text(size = 8))

last_crossing10 <- crossing10[crossing10$inverse == max(config10$inverse_theta), ]
last_factor10 <- factor10[factor10$n == max(config10$factor_n), ]
last_mc10 <- mc10[mc10$n == max(config10$sample_sizes), ]
near10 <- last_mc10[abs(last_mc10$c_value - 2) < 0.3, ]
away10 <- last_mc10[abs(last_mc10$c_value - 2) >= 0.5, ]
finite_error10 <- max(abs(last_mc10$ratio_exact / last_mc10$ratio_limit - 1))
mc_error10 <- max(abs(away10$ratio_mc / away10$ratio_exact - 1))
mean_error10 <- max(near10$mean_se / abs(near10$mean_exact))
