here::i_am("src/weighted_rv_toy.R")

# Power weight: w_i = i^alpha
make_power_weight_fun <- function(alpha) {
  force(alpha)
  function(i, n) {
    i^alpha
  }
}

# Unweighted per-step loss difference between Model 0 and Model 1
unweighted_increments <- function(y) {
  n <- length(y)
  cumsum_y <- cumsum(y)
  k <- seq_len(n - 1) # training size
  mu_hat <- cumsum_y[k] / k
  y_val <- y[k + 1]
  -2 * y_val * mu_hat + mu_hat^2
}

# Theoretical mean: sum_{i=2}^n w(i,n) / (i-1)
# Theoretical mean of weighted RV
theory_mean_weighted_rv <- function(n, weight_fun) {
  i <- seq.int(2, n)
  sum(weight_fun(i, n) / (i - 1))
}

# Simulate weighted RV once, for a given sample size and weight function, returns one weighted RV score
simulate_weighted_rv_once <- function(n, weight_fun) {
  y <- rnorm(n)
  inc <- unweighted_increments(y)
  i <- seq.int(2, n)
  w <- weight_fun(i, n)
  sum(w * inc)
}

estimate_weighted_rv <- function(n, mc_runs, weight_fun) {
  xi <- replicate(mc_runs, simulate_weighted_rv_once(n, weight_fun))
  mean_xi <- mean(xi)
  var_xi <- var(xi)
  theory_mean <- theory_mean_weighted_rv(n, weight_fun)
  list(
    n = n,
    mc_runs = mc_runs,
    p_neg = mean(xi < 0),  # incorrect selection rate
    mean_xi = mean_xi,
    var_xi = var_xi,
    theory_mean = theory_mean,
    R = var_xi / (mean_xi^2)
  )
}

# Efficient grid: for each MC path of length n_max, reuse increments across n and alpha
run_weight_grid <- function(sample_sizes, alphas, mc_runs) {
  sample_sizes <- sort(as.integer(sample_sizes))
  n_max <- max(sample_sizes)
  alphas <- as.numeric(alphas)

  set.seed(2026)

  # Store Xi for each (alpha, n, replication)
  # Build tidy rows
  rows <- vector("list", length(alphas) * length(sample_sizes))
  row_id <- 0

  # Precompute weights for each alpha and each n (depends only on i and n, not on Y)
  weight_cache <- lapply(alphas, function(a) {
    wf <- make_power_weight_fun(a)
    lapply(sample_sizes, function(n) {
      i <- seq.int(2, n)
      wf(i, n)
    })
  })

  # Precompute theoretical mean of weighted RV for each alpha and each n
  theory_cache <- lapply(seq_along(alphas), function(a_idx) {
    sapply(seq_along(sample_sizes), function(n_idx) {
      n <- sample_sizes[[n_idx]]
      w <- weight_cache[[a_idx]][[n_idx]]
      i <- seq.int(2, n)
      sum(w / (i - 1))
    })
  })

  # Accumulators: for each (alpha, n) store vector of Xi of length mc_runs
  xi_store <- lapply(alphas, function(a) {
    lapply(sample_sizes, function(n) numeric(mc_runs))
  })

  for (rep in seq_len(mc_runs)) {
    y <- rnorm(n_max)
    inc_full <- unweighted_increments(y) # length n_max - 1, for i=2..n_max

    for (a_idx in seq_along(alphas)) {
      for (n_idx in seq_along(sample_sizes)) {
        n <- sample_sizes[[n_idx]]
        # increments for i=2..n are first (n-1) entries
        inc <- inc_full[seq_len(n - 1)]
        w <- weight_cache[[a_idx]][[n_idx]]
        xi_store[[a_idx]][[n_idx]][rep] <- sum(w * inc)
      }
    }

    if (rep %% 100L == 0L || rep == mc_runs) {
      message("Finished MC replication ", rep, " / ", mc_runs)
    }
  }

  for (a_idx in seq_along(alphas)) {
    for (n_idx in seq_along(sample_sizes)) {
      xi <- xi_store[[a_idx]][[n_idx]]
      mean_xi <- mean(xi)
      var_xi <- var(xi)
      theory_mean <- theory_cache[[a_idx]][[n_idx]]
      row_id <- row_id + 1L
      rows[[row_id]] <- data.frame(
        alpha = alphas[[a_idx]],
        n = sample_sizes[[n_idx]],
        mc_runs = mc_runs,
        p_neg = mean(xi < 0),
        mean_xi = mean_xi,
        var_xi = var_xi,
        theory_mean = theory_mean,
        R = if (abs(mean_xi) < .Machine$double.eps) NA_real_ else var_xi / (mean_xi^2),
        stringsAsFactors = FALSE
      )
    }
  }

  do.call(rbind, rows)
}
