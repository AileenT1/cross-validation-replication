library(dplyr)
library(ggplot2)
library(knitr)
library(kableExtra)

here::i_am("code/09_crossing_point_simulation.R")
source(here::here("src", "09_crossing_point.R"), local = TRUE)

config09 <- list(
  theta = c(0.2, 0.1, 0.05, 0.025),
  weights = c("unweighted", "log", "log_squared", "loglog", "geometric"),
  geometric_base = 1.05,
  n_max = 40000L,
  c_values = c(0.75, 1, 1.25, 2),
  mc_theta = 0.1,
  mc_replications = 2000L,
  mc_batch = 200L,
  seed = 2026L
)

validate09 <- function() {
  n <- 5L
  theta <- 0.2
  v <- nine_weight(n, "log")
  tails <- rev(cumsum(rev(v / seq_len(n)^2)))
  b <- c(tails[-1L], 0) - v / seq_len(n)
  mat <- matrix(0, n + 1L, n + 1L)
  diag(mat)[seq_len(n)] <- tails
  for (l in 2:(n + 1L)) {
    mat[seq_len(l - 1L), l] <- b[l - 1L]
    mat[l, seq_len(l - 1L)] <- b[l - 1L]
  }
  linear <- c(0, -2 * theta * v)
  xi <- c(-0.6, 0.2, 0.3, -0.4, 0.5, -0.1)
  y <- theta + xi
  score_direct <- sum(v * ((y[-1L] - cumsum(y)[seq_len(n)] /
                             seq_len(n))^2 - y[-1L]^2))
  score_matrix <- as.numeric(t(xi) %*% mat %*% xi) +
    sum(linear * xi) - theta^2 * sum(v)
  expected_matrix <- sum(diag(mat)) - theta^2 * sum(v)
  variance_matrix <- 2 * sum(mat^2) + sum(linear^2)
  moment <- nine_variance(n, theta, "log")
  discrepancy <- max(
    abs(score_direct - score_matrix),
    abs(expected_matrix / sum(v) - moment$mean),
    abs(variance_matrix / sum(v)^2 - moment$total)
  )
  stopifnot(is.finite(discrepancy), discrepancy < 1e-12)
  discrepancy
}

compute09 <- function(config) {
  validation_error <- validate09()
  crossings <- list()
  paths <- list()
  index <- 0L
  for (theta in config$theta) {
    pop <- nine_population_crossing(theta)
    for (kind in config$weights) {
      index <- index + 1L
      path <- nine_exact_path(theta, kind, config$n_max,
                              config$geometric_base)
      rv <- nine_crossing(path)
      stopifnot(!is.na(rv), rv >= pop, path$delta_expected[rv] < 0,
                rv == 1L || path$delta_expected[rv - 1L] >= 0)
      crossings[[index]] <- data.frame(
        theta = theta, weight = kind, n_pop = pop, n_rv = rv,
        difference = rv - pop, crossing_ratio = rv / pop,
        improves_over_unweighted = NA
      )
      keep <- sort(unique(c(
        round(exp(seq(log(2), log(config$n_max), length.out = 240))),
        round(config$c_values * pop), pop, rv, rv - 1L
      )))
      keep <- keep[keep >= 1L & keep <= config$n_max]
      paths[[index]] <- path[keep, , drop = FALSE]
      paths[[index]]$n_pop <- pop
      paths[[index]]$n_over_pop <- paths[[index]]$n / pop
      paths[[index]]$delta_scaled <- paths[[index]]$delta_expected / theta^2
      paths[[index]]$population_scaled <- paths[[index]]$delta_pop / theta^2
    }
  }
  crossing <- do.call(rbind, crossings)
  for (theta in config$theta) {
    baseline <- crossing$difference[
      crossing$theta == theta & crossing$weight == "unweighted"
    ]
    crossing$improves_over_unweighted[crossing$theta == theta] <-
      crossing$difference[crossing$theta == theta] < baseline
  }
  theory_path <- do.call(rbind, paths)
  local <- list()
  index <- 0L
  for (theta in config$theta) {
    pop <- nine_population_crossing(theta)
    for (c_value in config$c_values) {
      n <- max(1L, as.integer(round(c_value * pop)))
      for (kind in config$weights) {
        index <- index + 1L
        moment <- nine_variance(n, theta, kind, config$geometric_base)
        local[[index]] <- data.frame(
          theta = theta, c_target = c_value, n = n, n_over_pop = n / pop,
          weight = kind, delta_pop = 1 / n - theta^2,
          delta_expected = moment$mean,
          variance_quadratic = moment$quadratic,
          variance_linear = moment$linear,
          variance_total = moment$total,
          chebyshev_ratio = if (abs(moment$mean) < 1e-14) NA_real_
                            else moment$total / moment$mean^2,
          effective_n = moment$effective_n
        )
      }
    }
  }
  local <- do.call(rbind, local)
  mc_theta <- config$mc_theta
  mc_cross <- crossing[crossing$theta == mc_theta, ]
  mc_pop <- nine_population_crossing(mc_theta)
  mc_max <- as.integer(ceiling(1.15 * max(mc_cross$n_rv)))
  mc_n <- sort(unique(c(
    round(exp(seq(log(5), log(mc_max), length.out = 26))),
    round(config$c_values * mc_pop), mc_pop,
    mc_cross$n_rv, mc_cross$n_rv - 1L, mc_cross$n_rv + 1L
  )))
  mc_n <- mc_n[mc_n >= 1L & mc_n <= mc_max]
  mc <- nine_mc(mc_theta, config$weights, mc_n,
                config$mc_replications, config$mc_batch,
                config$seed, config$geometric_base)
  mc$mean_exact <- vapply(seq_len(nrow(mc)), function(j)
    nine_exact_path(mc_theta, mc$weight[j], mc$n[j],
                    config$geometric_base)$delta_expected[mc$n[j]],
    numeric(1))
  stopifnot(all(is.finite(mc$mean_mc)), all(is.finite(mc$variance_mc)),
            all(mc$variance_mc >= 0), all(mc$p_model0 >= 0 & mc$p_model0 <= 1),
            all(is.finite(local$variance_total)), all(local$variance_total > 0))
  list(config = config, validation_error = validation_error,
       crossing = crossing, theory_path = theory_path,
       local = local, mc = mc)
}

cache09 <- here::here("data", "final", "09_crossing_point_results.rds")
if (file.exists(cache09)) {
  results09 <- readRDS(cache09)
  stopifnot(identical(results09$config, config09),
            all(c("crossing", "theory_path", "local", "mc",
                  "validation_error") %in% names(results09)))
} else {
  results09 <- compute09(config09)
  saveRDS(results09, cache09)
}
stopifnot(nrow(results09$crossing) ==
            length(config09$theta) * length(config09$weights),
          nrow(results09$local) ==
            length(config09$theta) * length(config09$weights) *
            length(config09$c_values),
          results09$validation_error < 1e-12)

weight_labels09 <- c(
  unweighted = "Unweighted: 1",
  log = "Log: log(i+1)",
  log_squared = "Squared log: [log(i+1)]^2",
  loglog = "Log-log: log log(i+e)",
  geometric = "Geometric: 1.05^i"
)
label09 <- function(x) factor(weight_labels09[x], levels = weight_labels09)

crossing09 <- results09$crossing |>
  mutate(label = label09(weight))
path09 <- results09$theory_path |>
  mutate(label = label09(weight))
local09 <- results09$local |>
  mutate(label = label09(weight))
mc09 <- results09$mc |>
  mutate(label = label09(weight), n_over_pop =
           n / nine_population_crossing(config09$mc_theta))

style09 <- theme_bw(base_size = 10) +
  theme(legend.position = "bottom", strip.text = element_text(size = 8),
        panel.grid.minor = element_blank())

comparison09 <- path09 |>
  filter(theta == config09$mc_theta, n_over_pop >= 0.5,
         n_over_pop <= 9) |>
  ggplot(aes(n_over_pop)) +
  geom_hline(yintercept = 0, colour = "grey55") +
  geom_line(aes(y = population_scaled, colour = "Population risk"),
            linewidth = 0.5) +
  geom_line(aes(y = delta_scaled, colour = "Expected weighted RV"),
            linewidth = 0.6) +
  geom_vline(aes(xintercept = n_rv / n_pop), data = crossing09 |>
               filter(theta == config09$mc_theta),
             colour = "#B35806", linetype = 3, linewidth = 0.35) +
  geom_vline(xintercept = 1, colour = "grey40", linetype = 2,
             linewidth = 0.35) +
  facet_wrap(~label, ncol = 3, scales = "free_y") +
  scale_colour_manual(values = c("Population risk" = "#2166AC",
                                 "Expected weighted RV" = "#B35806")) +
  labs(x = "Training size / population crossing",
       y = "Loss difference / theta^2", colour = NULL) + style09

crossings_plot09 <- crossing09 |>
  ggplot(aes(1 / theta^2, crossing_ratio, colour = label)) +
  geom_hline(yintercept = 1, linetype = 2, colour = "grey50") +
  geom_line(linewidth = 0.6) + geom_point(size = 1.5) +
  scale_x_log10() +
  labs(x = "1 / theta^2 (log scale)",
       y = "Expected RV crossing / population crossing",
       colour = "Weight") + style09

tradeoff09 <- local09 |>
  filter(c_target %in% c(0.75, 1.25, 2)) |>
  ggplot(aes(1 / theta^2, chebyshev_ratio, colour = label)) +
  geom_line(linewidth = 0.55) + geom_point(size = 1.35) +
  facet_wrap(~c_target, nrow = 1,
             labeller = labeller(c_target = function(x) paste("c =", x))) +
  scale_x_log10() + scale_y_log10() +
  labs(x = "1 / theta^2 (log scale)",
       y = "Exact variance / expected RV difference squared",
       colour = "Weight") + style09

mc_mean_plot09 <- mc09 |>
  ggplot(aes(n_over_pop)) +
  geom_hline(yintercept = 0, colour = "grey60") +
  geom_line(aes(y = mean_exact), colour = "#2166AC", linewidth = 0.5) +
  geom_point(aes(y = mean_mc), colour = "#B35806", size = 1.2) +
  facet_wrap(~label, ncol = 3, scales = "free_y") +
  scale_x_log10() +
  labs(x = "Training size / population crossing (log scale)",
       y = "Normalized expected RV difference") + style09

mc_select_plot09 <- mc09 |>
  ggplot(aes(n_over_pop, p_model0)) +
  geom_line(colour = "#2166AC", linewidth = 0.5) +
  geom_point(colour = "#2166AC", size = 1.2) +
  geom_errorbar(aes(ymin = pmax(0, p_model0 - 1.96 * se_model0),
                    ymax = pmin(1, p_model0 + 1.96 * se_model0)),
                width = 0.04, colour = "#636363") +
  facet_wrap(~label, ncol = 3, scales = "free_y") +
  scale_x_log10() +
  labs(x = "Training size / population crossing (log scale)",
       y = "Monte Carlo frequency selecting Model 0") + style09
