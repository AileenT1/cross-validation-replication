here::i_am("code/02_cv_error_rate_poly.R")

report_only <- isTRUE(getOption("cv_error_rate_02_poly_report_only", FALSE))

if (report_only) {
  suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
  })
  rds_path <- here::here("data", "final", "02_cv_error_rate_poly_results.rds")
  if (!file.exists(rds_path)) {
    stop(
      "Missing results file: data/final/02_cv_error_rate_poly_results.rds. ",
      "Run code/02_cv_error_rate_poly.R first (e.g. Rscript code/02_cv_error_rate_poly.R)."
    )
  }
  output <- readRDS(rds_path)
  params <- output$parameters
  by_n <- output$by_n
  settings_tbl <- data.frame(
    setting = c(
      "seed", "K", "Monte Carlo replications per n", "Sample sizes (grid)",
      "beta0", "beta1", "beta2", "beta3", "sigma", "poly_degree", "truth", "selection rule"
    ),
    value = c(
      params$seed,
      params$K,
      params$mc_runs,
      params$sample_sizes_label,
      params$beta0,
      params$beta1,
      params$beta2,
      params$beta3,
      params$sigma,
      params$poly_degree,
      params$truth,
      params$selection_rule
    )
  )
  error_rate_summary_tbl <- by_n |>
    arrange(n)
  error_rate_plot <- ggplot2::ggplot(by_n, ggplot2::aes(x = n, y = error_rate)) +
    ggplot2::geom_line(linewidth = 0.35) +
    ggplot2::geom_ribbon(
      ggplot2::aes(
        ymin = pmax(0, error_rate - mc_se),
        ymax = pmin(0.3, error_rate + mc_se)
      ),
      alpha = 0.15,
      linewidth = 0
    ) +
    ggplot2::labs(
      x = "Sample size",
      y = "CV selection error rate",
      title = "Rate CV picks linear under polynomial truth"
    ) +
    ggplot2::scale_x_continuous(
      limits = c(0, 1000),
      breaks = seq(0, 1000, by = 50)
    ) +
    ggplot2::scale_y_continuous(limits = c(0, 0.2), expand = ggplot2::expansion(mult = c(0, 0.02))) +
    ggplot2::theme_bw()
} else {
  if (exists("kfold_cv_linear", mode = "function", inherits = TRUE)) {
    rm(kfold_cv_linear, envir = .GlobalEnv)
  }
  if (exists("kfold_cv_polynomial", mode = "function", inherits = TRUE)) {
    rm(kfold_cv_polynomial, envir = .GlobalEnv)
  }

  source(here::here("src", "cv_linear.R"), local = .GlobalEnv)
  source(here::here("src", "cv_polynomial.R"), local = .GlobalEnv)

  set.seed(2026)

  K <- 10L
  mc_runs <- 1000L
  sample_sizes <- seq(50L, 1000L, by = 50L)
  beta0 <- 1
  beta1 <- 2
  beta2 <- 1
  beta3 <- 1
  sigma <- 4
  poly_degree <- 3L

  run_one_replication <- function(replication_id, n_now) {
    x <- rnorm(n_now)
    eps <- rnorm(n_now, sd = sigma)
    y <- beta0 + beta1 * x + beta2 * x^2 + beta3 * x^3 + eps
    dat <- data.frame(y = y, x = x)

    fold_id <- sample(rep(seq_len(K), length.out = n_now))

    linear_cv <- kfold_cv_linear(
      data = dat,
      formula = stats::as.formula("y ~ x"),
      K = K,
      fold_id = fold_id
    )

    poly_cv <- kfold_cv_polynomial(
      data = dat,
      response = "y",
      predictor = "x",
      degree = poly_degree,
      K = K,
      raw = TRUE,
      fold_id = fold_id
    )

    mean_cv_linear <- linear_cv$mean_cv_mse
    mean_cv_poly3 <- poly_cv$mean_cv_mse
    selected_model <- if (mean_cv_linear <= mean_cv_poly3) "linear" else "poly3"
    wrong <- selected_model == "linear"

    data.frame(
      n = n_now,
      replication = replication_id,
      wrong = wrong,
      stringsAsFactors = FALSE
    )
  }

  by_n_list <- vector("list", length(sample_sizes))
  for (i in seq_along(sample_sizes)) {
    n_now <- sample_sizes[[i]]
    run_results <- lapply(seq_len(mc_runs), function(replication_id) {
      run_one_replication(replication_id = replication_id, n_now = n_now)
    })
    df <- do.call(rbind, run_results)
    wrong_rate <- mean(df$wrong)
    mc_se <- sqrt(wrong_rate * (1 - wrong_rate) / mc_runs)
    by_n_list[[i]] <- data.frame(
      n = n_now,
      error_rate = wrong_rate,
      mc_se = mc_se,
      stringsAsFactors = FALSE
    )
    if (i %% 25L == 0L || i == length(sample_sizes)) {
      message("Finished sample size index ", i, " / ", length(sample_sizes), " (n = ", n_now, ")")
    }
  }
  by_n <- do.call(rbind, by_n_list)

  results_dir <- here::here("data", "final")
  if (!dir.exists(results_dir)) {
    dir.create(results_dir, recursive = TRUE)
  }

  output <- list(
    parameters = list(
      seed = 2026L,
      K = K,
      mc_runs = mc_runs,
      sample_sizes_label = "50, 100, ..., 1000 (step 50)",
      beta0 = beta0,
      beta1 = beta1,
      beta2 = beta2,
      beta3 = beta3,
      sigma = sigma,
      poly_degree = poly_degree,
      truth = "polynomial: y = beta0 + beta1*x + beta2*x^2 + beta3*x^3 + eps",
      selection_rule = "Choose poly3 if mean CV MSE strictly lower than linear; else linear"
    ),
    by_n = by_n
  )

  saveRDS(output, file = file.path(results_dir, "02_cv_error_rate_poly_results.rds"))

  message("Saved results to ", file.path(results_dir, "02_cv_error_rate_poly_results.rds"))
}
