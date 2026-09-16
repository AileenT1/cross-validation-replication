here::i_am("src/07_slowly_varying_rv.R")
source(here::here("src", "weighted_rv_toy.R"))

# All weights are deterministic functions of validation time, not final n.
sv_registry <- function() {
  p <- c(0, 1, 2, -1, -2, -1.25, -0.75, -0.5, 0.5)
  data.frame(
    weight_id = c(paste0("logp_", p), "poly2", "loglog", "power001"),
    label = c("1", "log(i+1)", "log(i+1)^2", "1/log(i+1)",
              "1/log(i+1)^2", "log(i+1)^(-1.25)", "log(i+1)^(-0.75)",
              "log(i+1)^(-0.5)", "log(i+1)^0.5",
              "1 + log(i+1) + log(i+1)^2", "log(log(i+e))", "i^0.01"),
    family = c(rep("log_power", length(p)), "polynomial", "loglog", "power"),
    p = c(p, NA_real_, NA_real_, NA_real_),
    main = c(rep(TRUE, 5), rep(FALSE, 4), TRUE, TRUE, FALSE)
  )
}

sv_weight <- function(i, spec) {
  switch(spec$family,
    log_power = log1p(i)^spec$p,
    polynomial = 1 + log1p(i) + log1p(i)^2,
    loglog = log(log(i + exp(1))),
    power = i^0.01,
    stop("Unknown weight family.")
  )
}

sv_coefficients <- function(w) {
  stopifnot(length(w) >= 1L, all(is.finite(w)), all(w > 0))
  m <- seq_along(w)
  a <- rev(cumsum(rev(w / m^2)))
  b <- c(a[-1L], 0) - w / m
  list(a = a, b = b)
}

sv_moments <- function(w) {
  z <- sv_coefficients(w)
  mean <- sum(w / seq_along(w))
  variance <- 2 * sum(z$a^2) + 4 * sum(seq_along(w) * z$b^2)
  c(mean = mean, variance = variance, ratio = variance / mean^2)
}

# O(max(n)) time and O(chunk_size) working storage per weight.
# From the loss-derived quadratic matrix:
# V_n - V_(n-1) = 2 w_n/(n-1)^2 * ((2n-1)w_n - 2 sum_(i<n)w_i).
sv_theory <- function(sample_sizes, registry = sv_registry(), chunk_size = 100000L) {
  ns <- sort(unique(as.integer(sample_sizes)))
  stopifnot(length(ns) > 0, min(ns) >= 2L, chunk_size >= 1L)
  rows <- lapply(seq_len(nrow(registry)), function(j) {
    spec <- registry[j, ]
    mean0 <- variance0 <- weight0 <- 0
    ans <- matrix(NA_real_, length(ns), 3L)
    for (lo in seq.int(2L, max(ns), by = chunk_size)) {
      i <- seq.int(lo, min(max(ns), lo + chunk_size - 1L))
      w <- sv_weight(i, spec)
      previous <- weight0 + c(0, head(cumsum(w), -1L))
      means <- mean0 + cumsum(w / (i - 1))
      variances <- variance0 + cumsum(2 * w / (i - 1)^2 *
        ((2 * i - 1) * w - 2 * previous))
      use <- which(ns >= lo & ns <= tail(i, 1))
      if (length(use)) {
        pos <- ns[use] - lo + 1L
        ans[use, ] <- cbind(means[pos], variances[pos],
                            variances[pos] / means[pos]^2)
      }
      mean0 <- tail(means, 1)
      variance0 <- tail(variances, 1)
      weight0 <- weight0 + sum(w)
    }
    data.frame(weight_id = spec$weight_id, n = ns,
               mean = ans[, 1], variance = ans[, 2], ratio = ans[, 3])
  })
  do.call(rbind, rows)
}

# Shared Gaussian paths across all weights and n; independent across replicates.
# Matrix products are restricted to disjoint n intervals (no repeated prefixes).
sv_mc <- function(sample_sizes, registry = sv_registry(), mc_runs = 5000L,
                  batch_size = 250L, seed = 2026L) {
  ns <- sort(unique(as.integer(sample_sizes)))
  stopifnot(min(ns) >= 2L, mc_runs >= 2L, batch_size >= 1L)
  set.seed(seed)
  nw <- nrow(registry)
  store <- array(NA_real_, c(mc_runs, length(ns), nw))
  for (start in seq.int(1L, mc_runs, by = batch_size)) {
    ids <- start:min(mc_runs, start + batch_size - 1L)
    b <- length(ids)
    running_sum <- rnorm(b)
    scores <- matrix(0, b, nw)
    left <- 2L
    for (k in seq_along(ns)) {
      # Small time blocks bound memory even for large MC horizons.
      for (lo in seq.int(left, ns[k], by = 2000L)) {
        i <- seq.int(lo, min(ns[k], lo + 1999L))
        y <- matrix(rnorm(b * length(i)), nrow = length(i), ncol = b)
        cs <- apply(y, 2L, cumsum)
        dim(cs) <- dim(y)
        past <- sweep(cs - y, 2L, running_sum, "+")
        mu <- past / (i - 1)
        increments <- mu^2 - 2 * y * mu
        weights <- vapply(seq_len(nw), function(j) sv_weight(i, registry[j, ]),
                          numeric(length(i)))
        dim(weights) <- c(length(i), nw)
        scores <- scores + crossprod(increments, weights)
        running_sum <- running_sum + cs[nrow(cs), ]
      }
      store[ids, k, ] <- scores
      left <- ns[k] + 1L
    }
    message("Monte Carlo: ", max(ids), "/", mc_runs)
  }
  rows <- list()
  for (j in seq_len(nw)) for (k in seq_along(ns)) {
    x <- store[, k, j]
    count <- sum(x < 0)
    phat <- count / mc_runs
    z <- qnorm(0.975)
    center <- (phat + z^2 / (2 * mc_runs)) / (1 + z^2 / mc_runs)
    half <- z * sqrt(phat * (1 - phat) / mc_runs +
                     z^2 / (4 * mc_runs^2)) / (1 + z^2 / mc_runs)
    rows[[length(rows) + 1L]] <- data.frame(
      weight_id = registry$weight_id[j], n = ns[k], mc_runs = mc_runs,
      mean_mc = mean(x), variance_mc = var(x), ratio_mc = var(x) / mean(x)^2,
      wrong_count = count, p_wrong = phat,
      mcse = sqrt(phat * (1 - phat) / mc_runs),
      lower = max(0, center - half), upper = min(1, center + half))
  }
  do.call(rbind, rows)
}

# Negative-p limits: the mean tail is essential (converges very slowly near -1).
# Integral-tail corrections for the variance recurrence are asymptotic estimates,
# not rigorous error bounds. Multiple truncations are reported for transparency.
sv_limits <- function(truncations, registry = sv_registry()) {
  neg <- registry[registry$family == "log_power" & !is.na(registry$p) &
                    registry$p < 0, ]
  exact <- sv_theory(truncations, neg)
  out <- merge(exact, neg[c("weight_id", "p")], by = "weight_id")
  L <- log(out$n)
  # V infinity - V_N ~ -2 (log N)^(2p), from summing the recurrence.
  out$variance_limit_est <- out$variance - 2 * L^(2 * out$p)
  out$mean_limit_est <- ifelse(out$p < -1,
                              out$mean - L^(out$p + 1) / (out$p + 1), NA_real_)
  out$ratio_limit_est <- ifelse(out$p < -1,
                               out$variance_limit_est / out$mean_limit_est^2,
                               0)
  out
}

sv_validate <- function() {
  registry <- sv_registry()
  max_matrix <- max_moment <- max_score <- max_recurrence <- 0
  set.seed(17)
  for (n in c(2L, 3L, 7L, 19L)) for (j in seq_len(nrow(registry))) {
    w <- sv_weight(2:n, registry[j, ])
    # Independent construction: (e_i - past-average)(...)' - e_i e_i'.
    Q <- matrix(0, n, n)
    for (i in 2:n) {
      e <- numeric(n)
      e[i] <- 1
      prediction <- numeric(n)
      prediction[seq_len(i - 1L)] <- 1 / (i - 1)
      Q <- Q + w[i - 1L] * (tcrossprod(e - prediction) - tcrossprod(e))
    }
    cf <- sv_coefficients(w)
    A <- diag(c(cf$a, 0))
    for (l in 2:n) A[seq_len(l - 1), l] <- A[l, seq_len(l - 1)] <- cf$b[l - 1]
    moments <- sv_moments(w)
    max_matrix <- max(max_matrix, abs(A - Q))
    max_moment <- max(max_moment, abs(moments[1] - sum(diag(Q))),
                      abs(moments[2] - 2 * sum(Q^2)))
    y <- rnorm(n)
    direct <- sum(w * ((y[-1] - cumsum(y)[1:(n - 1)] / (1:(n - 1)))^2 - y[-1]^2))
    max_score <- max(max_score, abs(direct - sum(w * unweighted_increments(y))),
                     abs(direct - drop(crossprod(y, Q %*% y))))
    rec <- sv_theory(n, registry[j, ], chunk_size = 3L)
    max_recurrence <- max(max_recurrence, abs(unlist(rec[c("mean", "variance", "ratio")]) - moments))
  }
  errors <- c(matrix = max_matrix, moments = max_moment,
              scores = max_score, recurrence = max_recurrence)
  stopifnot(all(errors < 1e-9))
  base <- sv_theory(c(100L, 10000L), registry[registry$weight_id == "logp_0", ])
  stopifnot(max(abs(base$variance - vapply(base$n,
    function(n) 6 * sum(1 / seq_len(n - 1)^2), numeric(1)))) < 1e-9)
  small1 <- suppressMessages(sv_mc(c(5L, 10L), registry, 20L, 7L, 2026L))
  small2 <- suppressMessages(sv_mc(c(5L, 10L), registry, 20L, 7L, 2026L))
  stopifnot(identical(small1, small2))
  data.frame(test = names(errors), max_absolute_error = unname(errors))
}

sv_validate_cache <- function(x, config, registry) {
  stopifnot(identical(x$config, config), identical(x$registry, registry),
            is.data.frame(x$validation), is.data.frame(x$limits))
  for (name in c("exact", "mc")) {
    d <- x[[name]]
    sizes <- if (name == "exact") config$theory_n else config$mc_n
    required <- if (name == "exact") c("weight_id", "n", "mean", "variance", "ratio") else
      c("weight_id", "n", "mc_runs", "mean_mc", "variance_mc", "ratio_mc",
        "wrong_count", "p_wrong", "mcse", "lower", "upper")
    stopifnot(is.data.frame(d), nrow(d) == nrow(registry) * length(sizes),
              all(required %in% names(d)),
              !anyDuplicated(d[c("weight_id", "n")]),
              setequal(d$n, sizes), setequal(d$weight_id, registry$weight_id))
    numeric_columns <- vapply(d, is.numeric, logical(1))
    stopifnot(all(is.finite(as.matrix(d[numeric_columns]))))
  }
  stopifnot(all(x$exact$mean > 0), all(x$exact$variance > 0),
            all(x$exact$ratio > 0), all(x$mc$variance_mc >= 0),
            all(x$mc$p_wrong >= 0 & x$mc$p_wrong <= 1),
            all(x$mc$lower <= x$mc$p_wrong & x$mc$upper >= x$mc$p_wrong),
            all(x$mc$mc_runs == config$mc_runs),
            all(x$mc$wrong_count == round(x$mc$p_wrong * config$mc_runs)),
            all(x$mc$wrong_count >= 0 & x$mc$wrong_count <= config$mc_runs),
            all(x$mc$lower >= 0 & x$mc$upper <= 1))
  lim <- x$limits
  neg <- registry[registry$family == "log_power" & !is.na(registry$p) & registry$p < 0, ]
  stopifnot(nrow(lim) == nrow(neg) * length(config$limit_n),
            !anyDuplicated(lim[c("weight_id", "n")]),
            setequal(lim$n, config$limit_n),
            setequal(lim$weight_id, neg$weight_id),
            all(is.finite(lim$variance_limit_est)), all(lim$variance_limit_est > 0),
            all(is.finite(lim$ratio_limit_est)), all(lim$ratio_limit_est >= 0),
            all(is.na(lim$mean_limit_est) == (lim$p >= -1)),
            all(is.finite(lim$mean_limit_est[lim$p < -1])),
            all(lim$mean_limit_est[lim$p < -1] > 0),
            nrow(x$validation) == 4L,
            all(is.finite(x$validation$max_absolute_error)),
            all(x$validation$max_absolute_error < 1e-9))
  invisible(TRUE)
}
