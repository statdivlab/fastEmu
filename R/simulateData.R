#' @importFrom stats rnorm rbinom rpois rnbinom
# based on function to simulate data for simulations in Clausen & Willis (2024) in radEmu
simulateData <- function(X,
                         B,
                         distn,
                         zinb_size = NULL,
                         zinb_zero_prop = NULL,
                         mean_count_before_ZI,
                         scale_ZINB_mean = TRUE) {

  n <- nrow(X)
  J <- ncol(B)
  log_means <- do.call(cbind,
                       lapply(1:J,
                              function(j) X %*% B[, j, drop = FALSE]))

  row_means <- rowSums(exp(log_means)) / J

  z <- sapply(row_means, function(x) log(mean_count_before_ZI) - log(x) + stats::rnorm(1))
  Y <- matrix(0, ncol = J, nrow = n)

  for (i in 1:n) {
    log_means[i,] <- log_means[i,] + z[i]
  }

  for (i in 1:n) {
    accepted <- FALSE
    while (!accepted) {
      for (j in 1:J) {
        if (distn == "Poisson") {
          Y[i,j] <- stats::rpois(1, lambda = exp(log_means[i, j]))
        }
        if(distn == "ZINB"){
          if (scale_ZINB_mean) {
            Y[i,j] <- stats::rnbinom(1, mu = exp(log_means[i, j]) / (1 - zinb_zero_prop),
                                     size = zinb_size / ((1 - zinb_zero_prop)^2)) *
              (1 - stats::rbinom(1, 1, prob = zinb_zero_prop))
          } else {
            Y[i,j] <- stats::rnbinom(1,mu = exp(log_means[i, j]), size = zinb_size) *
              (1 - stats::rbinom(1, 1, prob = zinb_zero_prop))
          }
        }
        if (sum(Y[i,]) > 0) {
          accepted <- TRUE
        }
      }
    }
  }
  return(Y)
}
