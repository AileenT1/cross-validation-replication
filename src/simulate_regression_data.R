here::i_am("src/simulate_regression_data.R")

simulate_regression_data <- function(n, simulation_config) {
  truth <- simulation_config$truth
  beta <- truth$beta
  noise <- simulation_config$noise

  x <- rnorm(n)
  eps <- rnorm(n, sd = noise$sigma)

  if (truth$type == "linear") {
    y <- beta["beta0"] +
      beta["beta1"] * x +
      eps
  } else if (truth$type == "polynomial") {
    y <- beta["beta0"] +
      beta["beta1"] * x +
      beta["beta2"] * x^2 +
      beta["beta3"] * x^3 +
      eps
  } else {
    stop("Unknown truth type.")
  }

  data.frame(y = y, x = x)
}
