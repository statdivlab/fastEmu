test_that("same mean from either distribution", {
  # simulate data
  n <- 2
  J <- 100
  b0 <- rnorm(J)
  b1 <- rnorm(J)
  b2 <- rnorm(J)
  B <- rbind(b0, b1, b2)
  X <- cbind(1, rep(c(0, 1), each = n/2), rnorm(n))
  sim_Y <- data.frame(Y1_P = rep(NA, 1000),
                      Y1_Z = rep(NA, 1000))
  for (i in 1:1000) {
    Y_P <- simulateData(X = X,
                        B = B,
                        distn = "Poisson",
                        zinb_size = 5,
                        zinb_zero_prop = 0.6,
                        mean_count_before_ZI = 50)
    sim_Y$Y1_P[i] <- Y_P[1]
    Y_Z <- simulateData(X = X,
                        B = B,
                        distn = "ZINB",
                        zinb_size = 5,
                        zinb_zero_prop = 0.6,
                        mean_count_before_ZI = 50)
    sim_Y$Y1_Z[i] <- Y_Z[1]
  }
  means <- colMeans(sim_Y)
  expect_true(base::all.equal(as.vector(means[1]), as.vector(means[2]), tol = 0.25))
})

