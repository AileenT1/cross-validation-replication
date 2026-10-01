# Exact calculations and simulation for the new Analysis 09 model.
# Model 0: sample mean; Model 1: the fixed-zero predictor.

nine_weight <- function(n, kind, geometric_base = 1.05) {
  stopifnot(length(n) == 1L, n >= 1L, geometric_base > 1)
  i <- seq_len(n) + 1
  raw <- switch(
    kind,
    unweighted = rep(1, n),
    log = log(i + 1),
    log_squared = log(i + 1)^2,
    loglog = log(log(i + exp(1))),
    geometric = exp((seq_len(n) - n) * log(geometric_base)),
    stop("Unknown weight kind: ", kind)
  )
  raw / max(raw)
}

nine_population_crossing <- function(theta) {
  stopifnot(is.finite(theta), abs(theta) > 0, abs(theta) < 1)
  floor(1 / theta^2 + 1e-8) + 1L
}

nine_exact_path <- function(theta, kind, n_max, geometric_base = 1.05) {
  stopifnot(is.finite(theta), abs(theta) > 0, n_max >= 2L)
  m <- seq_len(n_max)
  if (kind == "geometric") {
    discount <- 1 / geometric_base
    a <- w <- w2 <- numeric(n_max)
    for (j in m) {
      a[j] <- (if (j > 1L) discount * a[j - 1L] else 0) + 1 / j
      w[j] <- (if (j > 1L) discount * w[j - 1L] else 0) + 1
      w2[j] <- (if (j > 1L) discount^2 * w2[j - 1L] else 0) + 1
    }
  } else {
    v <- nine_weight(n_max, kind, geometric_base)
    a <- cumsum(v / m)
    w <- cumsum(v)
    w2 <- cumsum(v^2)
  }
  data.frame(
    n = m, theta = theta, weight = kind,
    delta_pop = 1 / m - theta^2,
    average_inverse_training = a / w,
    delta_expected = a / w - theta^2,
    effective_n = w^2 / w2,
    linear_variance = 4 * theta^2 * w2 / w^2
  )
}

nine_crossing <- function(path) {
  hit <- which(path$delta_expected < 0)
  if (length(hit)) path$n[hit[1L]] else NA_integer_
}

nine_variance <- function(n, theta, kind, geometric_base = 1.05) {
  stopifnot(n >= 1L, is.finite(theta))
  v <- nine_weight(n, kind, geometric_base)
  m <- seq_len(n)
  tails <- rev(cumsum(rev(v / m^2)))
  b <- c(tails[-1L], 0) - v / m
  quadratic <- 2 * sum(tails^2) + 4 * sum(m * b^2)
  linear <- 4 * theta^2 * sum(v^2)
  w <- sum(v)
  list(
    quadratic = quadratic / w^2,
    linear = linear / w^2,
    total = (quadratic + linear) / w^2,
    mean = sum(v / m) / w - theta^2,
    effective_n = w^2 / sum(v^2)
  )
}

nine_mc <- function(theta, kinds, n_values, replications = 2000L,
                    batch_size = 200L, seed = 2026L,
                    geometric_base = 1.05) {
  stopifnot(length(n_values) > 0L, all(n_values >= 1L),
            replications > 1L, batch_size > 0L)
  n_values <- sort(unique(as.integer(n_values)))
  n_max <- max(n_values)
  set.seed(seed)
  key <- expand.grid(weight = kinds, n = n_values,
                     KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  scores <- matrix(NA_real_, nrow = replications, ncol = nrow(key))
  columns <- split(seq_len(nrow(key)), key$weight)
  fixed_weights <- lapply(setNames(kinds[kinds != "geometric"],
                                   kinds[kinds != "geometric"]),
                          function(kind) nine_weight(n_max, kind, geometric_base))
  for (start in seq.int(1L, replications, by = batch_size)) {
    finish <- min(start + batch_size - 1L, replications)
    for (r in start:finish) {
      xi <- rnorm(n_max + 1L)
      m <- seq_len(n_max)
      past <- cumsum(xi)[m]
      increment <- (past / m)^2 -
        2 * xi[m + 1L] * past / m -
        2 * theta * xi[m + 1L] - theta^2
      for (kind in kinds) {
        if (kind == "geometric") {
          discount <- 1 / geometric_base
          numer <- denom <- numeric(n_max)
          for (j in m) {
            numer[j] <- (if (j > 1L) discount * numer[j - 1L] else 0) +
              increment[j]
            denom[j] <- (if (j > 1L) discount * denom[j - 1L] else 0) + 1
          }
          values <- numer[n_values] / denom[n_values]
        } else {
          v <- fixed_weights[[kind]]
          values <- cumsum(v * increment)[n_values] / cumsum(v)[n_values]
        }
        scores[r, columns[[kind]]] <- values
      }
    }
  }
  key$mean_mc <- colMeans(scores)
  key$variance_mc <- apply(scores, 2, var)
  key$p_model0 <- colMeans(scores < 0)
  key$se_model0 <- sqrt(key$p_model0 * (1 - key$p_model0) / replications)
  key$replications <- replications
  key$theta <- theta
  key
}
