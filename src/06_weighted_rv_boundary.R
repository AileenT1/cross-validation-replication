here::i_am("src/06_weighted_rv_boundary.R")
source(here::here("src", "weighted_rv_toy.R"))

# Exact coefficients in the Gaussian quadratic-form representation.
weighted_rv_coefficients <- function(n, alpha) {
  n <- as.integer(n)
  stopifnot(length(n) == 1L, n >= 2L, length(alpha) == 1L, is.finite(alpha))

  i <- seq.int(2L, n)
  weights <- i^alpha
  terms <- weights / (i - 1)^2
  reverse_sums <- rev(cumsum(rev(terms)))

  a <- reverse_sums
  future_sums <- c(reverse_sums[-1L], 0)
  b <- future_sums - weights / (i - 1)

  list(a = a, b = b)
}

weighted_rv_exact_theory <- function(n, alpha) {
  n <- as.integer(n)
  coefficients <- weighted_rv_coefficients(n, alpha)
  i <- seq.int(2L, n)

  theory_mean <- sum(i^alpha / (i - 1))
  diag_var <- 2 * sum(coefficients$a^2)
  offdiag_var <- 4 * sum((i - 1) * coefficients$b^2)
  theory_var <- diag_var + offdiag_var

  data.frame(
    alpha = alpha,
    n = n,
    theory_mean = theory_mean,
    diag_var = diag_var,
    offdiag_var = offdiag_var,
    theory_var = theory_var,
    theory_R = theory_var / theory_mean^2,
    asymptotic_R = if (alpha > 0) 2 * alpha^2 / (1 + alpha) else if (alpha == 0) 0 else NA_real_,
    stringsAsFactors = FALSE
  )
}

run_weighted_rv_theory_grid <- function(sample_sizes, alphas) {
  rows <- vector("list", length(sample_sizes) * length(alphas))
  row_id <- 0L
  for (alpha in alphas) {
    for (n in sample_sizes) {
      row_id <- row_id + 1L
      rows[[row_id]] <- weighted_rv_exact_theory(n, alpha)
    }
  }
  do.call(rbind, rows)
}

# Batched simulation avoids retaining all Gaussian paths in memory. Within a
# path, every requested n and alpha uses the same nested observations.
run_weighted_rv_mc <- function(
    sample_sizes,
    alphas,
    mc_runs,
    seed = 2026,
    batch_size = 250L) {
  sample_sizes <- sort(unique(as.integer(sample_sizes)))
  alphas <- as.numeric(alphas)
  mc_runs <- as.integer(mc_runs)
  batch_size <- as.integer(batch_size)
  stopifnot(min(sample_sizes) >= 2L, mc_runs >= 2L, batch_size >= 1L)

  n_max <- max(sample_sizes)
  n_count <- length(sample_sizes)
  alpha_count <- length(alphas)
  xi_store <- array(NA_real_, dim = c(mc_runs, n_count, alpha_count))
  target_lookup <- match(seq_len(n_max), sample_sizes, nomatch = 0L)
  alpha_weights <- outer(seq_len(n_max), alphas, `^`)

  set.seed(seed)
  first_rep <- 1L
  while (first_rep <= mc_runs) {
    last_rep <- min(first_rep + batch_size - 1L, mc_runs)
    current_size <- last_rep - first_rep + 1L
    y <- matrix(rnorm(current_size * n_max), nrow = current_size)
    past_sum <- y[, 1L]
    xi <- matrix(0, nrow = current_size, ncol = alpha_count)

    for (i in seq.int(2L, n_max)) {
      mu_hat <- past_sum / (i - 1)
      increment <- -2 * y[, i] * mu_hat + mu_hat^2
      xi <- xi + increment * rep(alpha_weights[i, ], each = current_size)

      target_id <- target_lookup[i]
      if (target_id > 0L) {
        xi_store[first_rep:last_rep, target_id, ] <- xi
      }
      past_sum <- past_sum + y[, i]
    }
    message("Finished Monte Carlo replications ", first_rep, "--", last_rep, " / ", mc_runs)
    first_rep <- last_rep + 1L
  }

  rows <- vector("list", n_count * alpha_count)
  row_id <- 0L
  for (alpha_id in seq_along(alphas)) {
    for (n_id in seq_along(sample_sizes)) {
      xi_values <- xi_store[, n_id, alpha_id]
      exact <- weighted_rv_exact_theory(sample_sizes[n_id], alphas[alpha_id])
      p_neg <- mean(xi_values < 0)
      emp_var <- stats::var(xi_values)
      row_id <- row_id + 1L
      rows[[row_id]] <- data.frame(
        alpha = alphas[alpha_id],
        n = sample_sizes[n_id],
        mc_runs = mc_runs,
        p_neg = p_neg,
        p_neg_mc_se = sqrt(p_neg * (1 - p_neg) / mc_runs),
        emp_mean = mean(xi_values),
        emp_var = emp_var,
        emp_R = emp_var / exact$theory_mean^2,
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

# Midpoint Euler--Maclaurin tail approximations for the two power series used
# below. Six expansion terms make the truncation error negligible on the
# grids used in Analysis 06.
weighted_mean_tail <- function(alpha, truncation, expansion_terms = 6L) {
  x <- truncation + 0.5
  k <- 0:expansion_terms
  sum(x^(alpha - k) / (k - alpha))
}

weighted_square_tail <- function(alpha, truncation, expansion_terms = 6L) {
  x <- truncation + 0.5
  k <- 0:expansion_terms
  sum((k + 1) * x^(alpha - 1 - k) / (1 + k - alpha))
}

weighted_rv_negative_limit <- function(alpha, truncation = 1000000L) {
  stopifnot(length(alpha) == 1L, alpha < 0, truncation >= 1000L)
  truncation <- as.integer(truncation)
  i <- seq.int(2L, truncation)
  terms <- i^alpha / (i - 1)^2
  square_tail <- weighted_square_tail(alpha, truncation)
  reverse_sums <- rev(cumsum(rev(terms))) + square_tail

  # a_{infinity,t}, t=1,...,N, including the first omitted-series tail.
  a <- c(reverse_sums, square_tail)
  l <- seq.int(2L, truncation)
  b <- a[l] - l^alpha / (l - 1)

  mean_limit <- sum(i^alpha / (i - 1)) + weighted_mean_tail(alpha, truncation)
  diagonal_tail <- 2 * truncation^(2 * alpha - 1) /
    ((1 - alpha)^2 * (1 - 2 * alpha))
  offdiagonal_tail <- -2 * alpha * truncation^(2 * alpha) / (1 - alpha)^2
  variance_limit <- 2 * sum(a^2) + diagonal_tail +
    4 * sum((l - 1) * b^2) + offdiagonal_tail

  data.frame(
    alpha = alpha,
    truncation = truncation,
    mean_limit = mean_limit,
    variance_limit = variance_limit,
    limit_R = variance_limit / mean_limit^2,
    stringsAsFactors = FALSE
  )
}

run_weighted_rv_limit_grid <- function(
    negative_alphas,
    positive_alphas,
    truncation = 1000000L) {
  negative_rows <- lapply(negative_alphas, weighted_rv_negative_limit, truncation = truncation)
  negative_results <- if (length(negative_rows)) do.call(rbind, negative_rows) else NULL
  nonnegative_results <- data.frame(
    alpha = c(0, positive_alphas),
    truncation = NA_integer_,
    mean_limit = NA_real_,
    variance_limit = NA_real_,
    limit_R = c(0, 2 * positive_alphas^2 / (1 + positive_alphas)),
    stringsAsFactors = FALSE
  )
  rbind(negative_results, nonnegative_results)
}

weighted_rv_quadratic_matrix <- function(n, alpha) {
  coefficients <- weighted_rv_coefficients(n, alpha)
  A <- matrix(0, nrow = n, ncol = n)
  diag(A)[seq_len(n - 1L)] <- coefficients$a
  for (l in seq.int(2L, n)) {
    A[seq_len(l - 1L), l] <- coefficients$b[l - 1L]
    A[l, seq_len(l - 1L)] <- coefficients$b[l - 1L]
  }
  A
}

validate_weighted_rv_identities <- function() {
  for (alpha in c(-0.5, 0, 0.5, 1)) {
    for (n in c(4L, 9L, 25L)) {
      exact <- weighted_rv_exact_theory(n, alpha)
      coefficients <- weighted_rv_coefficients(n, alpha)
      A <- weighted_rv_quadratic_matrix(n, alpha)
      stopifnot(
        isTRUE(all.equal(sum(coefficients$a), exact$theory_mean, tolerance = 1e-12)),
        isTRUE(all.equal(sum(diag(A)), exact$theory_mean, tolerance = 1e-12)),
        isTRUE(all.equal(2 * sum(A^2), exact$theory_var, tolerance = 1e-12))
      )
    }
  }
  invisible(TRUE)
}
