# Exact finite-n calculations and actual Gaussian-path Monte Carlo for Analysis 10.
# Validation time i=2,...,n+1; training size i-1. No image files are read or written.

ten_moments <- function(n, alpha, theta) {
  training <- seq_len(n)
  weight <- (training + 1)^alpha
  tail <- rev(cumsum(rev(weight / training^2)))
  off <- c(tail[-1L], 0) - weight / training
  mean <- sum(weight / training) - theta^2 * sum(weight)
  variance <- 2 * sum(tail^2) + 4 * sum(training * off^2) +
    4 * theta^2 * sum(weight^2)
  list(mean = mean, variance = variance, ratio = variance / mean^2)
}

ten_validate <- function() {
  n <- 6L
  alpha <- 1
  theta <- 0.23
  xi <- c(-0.6, 0.2, 0.3, -0.4, 0.5, -0.1, 0.8)
  training <- seq_len(n)
  weight <- (training + 1)^alpha
  y <- theta + xi
  direct <- sum(weight * ((y[-1] - cumsum(y)[training] / training)^2 - y[-1]^2))
  past <- cumsum(xi)[training] / training
  expanded <- sum(weight * (past^2 - 2 * xi[-1] * past - 2 * theta * xi[-1] - theta^2))
  tail <- rev(cumsum(rev(weight / training^2)))
  off <- c(tail[-1L], 0) - weight / training
  mat <- matrix(0, n + 1L, n + 1L)
  diag(mat)[training] <- tail
  for (j in 2:(n + 1L)) mat[seq_len(j - 1L), j] <- mat[j, seq_len(j - 1L)] <- off[j - 1L]
  linear <- c(0, -2 * theta * weight)
  matrix_score <- as.numeric(t(xi) %*% mat %*% xi) + sum(linear * xi) - theta^2 * sum(weight)
  moment <- ten_moments(n, alpha, theta)
  error <- max(abs(c(direct, expanded) - matrix_score),
               abs(moment$mean - (sum(diag(mat)) - theta^2 * sum(weight))),
               abs(moment$variance - (2 * sum(mat^2) + sum(linear^2))))
  stopifnot(is.finite(error), error < 1e-10)
  error
}

ten_deterministic <- function(config) {
  crossing <- factor <- list()
  j <- k <- 0L
  for (alpha in config$alphas) {
    lambda <- 1 + 1 / alpha
    n_max <- ceiling(2 * lambda * max(config$inverse_theta)^2)
    training <- seq_len(n_max)
    weight <- (training + 1)^alpha
    total <- cumsum(weight)
    reciprocal <- cumsum(weight / training)
    for (inverse in config$inverse_theta) {
      theta <- 1 / inverse
      # Exact integer population tie on this reciprocal-integer theta grid.
      n_pop <- inverse^2 + 1L
      delta <- reciprocal - theta^2 * total
      n_rv <- which(delta < 0)[1L]
      stopifnot(!is.na(n_rv), n_rv >= n_pop, delta[n_rv] < 0,
                n_rv == 1L || delta[n_rv - 1L] >= 0)
      j <- j + 1L
      crossing[[j]] <- data.frame(alpha, lambda, theta, inverse, n_pop, n_rv,
                                  ratio = n_rv / n_pop)
    }
    n <- config$factor_n
    k <- k + 1L
    factor[[k]] <- data.frame(alpha, lambda, n,
      factor = n * cumsum(weight^2)[n] / total[n]^2,
      target = (alpha + 1)^2 / (2 * alpha + 1))
  }
  list(crossing = do.call(rbind, crossing), factor = do.call(rbind, factor))
}

ten_mc <- function(config) {
  set.seed(config$seed)
  runs <- config$replications
  past_sum <- rnorm(runs)
  quadratic <- weighted_noise <- numeric(runs)
  snapshots <- list()
  # Each vector element is an independent Gaussian path. Shared paths pair
  # comparisons across c and n; they do not increase replication counts.
  for (n in seq_len(max(config$sample_sizes))) {
    fresh <- rnorm(runs)
    weight <- (n + 1)^config$mc_alpha
    past_mean <- past_sum / n
    quadratic <- quadratic + weight * (past_mean^2 - 2 * fresh * past_mean)
    weighted_noise <- weighted_noise + weight * fresh
    if (n %in% config$sample_sizes)
      snapshots[[as.character(n)]] <- list(quadratic = quadratic, noise = weighted_noise)
    past_sum <- past_sum + fresh
  }
  set.seed(config$seed + 1L)
  bootstrap <- matrix(sample.int(runs, runs * config$bootstraps, replace = TRUE), nrow = runs)
  rows <- list()
  j <- 0L
  lambda <- 1 + 1 / config$mc_alpha
  for (n in config$sample_sizes) {
    path <- snapshots[[as.character(n)]]
    total <- sum((seq_len(n) + 1)^config$mc_alpha)
    for (c_value in config$c_values) {
      theta <- sqrt(c_value / n)
      # This is the realized squared-loss score, with the identity checked
      # independently against direct predictions in ten_validate().
      score <- path$quadratic - 2 * theta * path$noise - theta^2 * total
      mean_mc <- mean(score)
      variance_mc <- var(score)
      exact <- ten_moments(n, config$mc_alpha, theta)
      sampled <- matrix(score[bootstrap], nrow = runs)
      bootstrap_ratio <- apply(sampled, 2L, var) / colMeans(sampled)^2
      interval <- quantile(bootstrap_ratio, c(0.025, 0.975), names = FALSE)
      j <- j + 1L
      rows[[j]] <- data.frame(n, c_value, theta, mean_mc, variance_mc,
        ratio_mc = variance_mc / mean_mc^2,
        mean_exact = exact$mean, variance_exact = exact$variance,
        ratio_exact = exact$ratio,
        ratio_limit = (2 * lambda / (lambda - 1) +
          4 * c_value * lambda^2 / (lambda^2 - 1)) / (lambda - c_value)^2,
        low = interval[1L], high = interval[2L], mean_se = sqrt(variance_mc / runs))
    }
  }
  do.call(rbind, rows)
}
