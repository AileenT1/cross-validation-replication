if (exists("kfold_cv_linear", mode = "function", inherits = TRUE)) {
  rm(kfold_cv_linear, envir = .GlobalEnv)
}
if (exists("kfold_cv_polynomial", mode = "function", inherits = TRUE)) {
  rm(kfold_cv_polynomial, envir = .GlobalEnv)
}

source("src/cv_linear.R", local = .GlobalEnv)
source("src/cv_polynomial.R", local = .GlobalEnv)

set.seed(2026)

K <- 10L
mc_runs <- 1000L
sample_sizes <- c(10000L, 30000L)
beta0 <- 1
beta1 <- 2
sigma <- 4
poly_degree <- 3L

run_one_replication <- function(replication_id, n_now) {
  x <- rnorm(n_now)
  eps <- rnorm(n_now, sd = sigma)
  y <- beta0 + beta1 * x + eps
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

  data.frame(
    n = n_now,
    replication = replication_id,
    mean_cv_linear = mean_cv_linear,
    mean_cv_poly3 = mean_cv_poly3,
    cv_gap_poly_minus_linear = mean_cv_poly3 - mean_cv_linear,
    selected_model = selected_model,
    stringsAsFactors = FALSE
  )
}

all_results <- vector("list", length(sample_sizes))
for (i in seq_along(sample_sizes)) {
  n_now <- sample_sizes[[i]]
  run_results <- lapply(seq_len(mc_runs), function(replication_id) {
    run_one_replication(replication_id = replication_id, n_now = n_now)
  })
  all_results[[i]] <- do.call(rbind, run_results)
}
results <- do.call(rbind, all_results)

selection_summary <- aggregate(
  replication ~ n + selected_model,
  data = results,
  FUN = length
)
names(selection_summary)[names(selection_summary) == "replication"] <- "count"
selection_summary <- selection_summary[order(selection_summary$n, selection_summary$selected_model), ]
selection_summary$selection_rate <- selection_summary$count / ave(selection_summary$count, selection_summary$n, FUN = sum)

gap_summary <- aggregate(
  cbind(mean_cv_linear, mean_cv_poly3, cv_gap_poly_minus_linear) ~ n,
  data = results,
  FUN = mean
)
gap_sd <- aggregate(
  cv_gap_poly_minus_linear ~ n,
  data = results,
  FUN = sd
)
names(gap_sd)[names(gap_sd) == "cv_gap_poly_minus_linear"] <- "sd_gap_poly_minus_linear"
gap_summary <- merge(gap_summary, gap_sd, by = "n", sort = TRUE)
names(gap_summary)[names(gap_summary) == "cv_gap_poly_minus_linear"] <- "mean_gap_poly_minus_linear"

results_dir <- "report/results"
if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}

output <- list(
  parameters = list(
    seed = 2026L,
    K = K,
    mc_runs = mc_runs,
    sample_sizes = sample_sizes,
    beta0 = beta0,
    beta1 = beta1,
    sigma = sigma,
    poly_degree = poly_degree
  ),
  results = results,
  selection_summary = selection_summary,
  gap_summary = gap_summary
)

saveRDS(output, file = file.path(results_dir, "01_cv_linear_truth_results.rds"))

message("Saved results to report/results/01_cv_linear_truth_results.rds")
