n <- 50
J <- 10
b0 <- rnorm(J)
b1 <- rnorm(J)
X <- cbind(1, rep(0:1, each = n / 2))
Y <- simulateData(X, rbind(b0, b1), distn = "Poisson", mean_count_before_ZI = 50)
rownames(X) <- paste0("sample", 1:nrow(X))
rownames(Y) <- rownames(X)

test_that("radEmu and fastEmu give same results when using full model and all categories in constraint", {

  radEmu_res <- suppressWarnings(emuFit(Y, X, test_kj = data.frame(k = 2, j = 6), match_row_names = FALSE))
  fastEmu_res <- suppressWarnings(fastEmuTest(constraint_cats = 1:10, Y = Y, X = X,
                             test_kj = data.frame(k = 2, j = 6), model = "full", match_row_names = FALSE))

  expect_true(all.equal(radEmu_res$coef, fastEmu_res$coef[, -10]))
})

test_that("fastEmu runs with and without estimation, and results don't change", {

  colnames(Y) <- paste0("taxon", 1:10)
  suppressWarnings({res_est <- fastEmuTest(constraint_cats = 1:5, Y = Y, X = X,
                            test_kj = data.frame(k = 2, j = 6), estimate_full_model = TRUE,
                            match_row_names = FALSE)})
  suppressWarnings({res_no_est <- fastEmuTest(constraint_cats = 1:5, Y = Y, X = X,
                            test_kj = data.frame(k = 2, j = 6), estimate_full_model = FALSE,
                            match_row_names = FALSE)})
  expect_true(res_est$coef$pval[6] == res_no_est$coef$pval)

  # also check we're returning what we want

  # all constraint categories are in "included_categories"
  expect_true(sum(paste0("taxon", 1:5) %in% unlist(res_est$included_categories)) == 5)
  # categories in coef are in the same order as they are in Y
  expect_true(all.equal(colnames(Y), res_est$coef$category))

})

test_that("fastEmu can work with additional arguments", {

  res <- suppressWarnings(fastEmuTest(constraint_cats = 1:5, Y = Y, X = X, test_kj = data.frame(k = 2, j = 6),
                            penalize = FALSE, return_nullB = TRUE, match_row_names = FALSE))

  # check that penalty is FALSE
  expect_false(res$penalized)

  # check that nullB is returned
  expect_true("null_B" %in% names(res))

})

test_that("fastEmu still works with deprecated aggregated model", {

  res_agg <- suppressWarnings(fastEmuTest(constraint_cats = 1:5, Y = Y, X = X, test_kj = data.frame(k = 2, j = 6),
                            model = "agg", match_row_names = FALSE))

  res_drop <- suppressWarnings(fastEmuTest(constraint_cats = 1:5, Y = Y, X = X, test_kj = data.frame(k = 2, j = 6)))
  expect_true(all.equal(res_agg$coef$pval[6], res_drop$coef$pval[6], tol = 0.05, match_row_names = FALSE))

})

test_that("fastEmu controls Type I error rate when it should", {

  skip("skipping this in automatic tests because it is slow")

  n <- 40
  J <- 10
  set.seed(1569)
  nsim <- 100
  ps <- rep(NA, nsim)
  bs <- simulateBs(J = J, constraint_cats = 1:5, test_j = 6, constraint_mag = 1,
                   other_mag = 5, under_null = TRUE)
  X = cbind(1, rep(0:1, each = n /2))

  for (i in 1:nsim) {

    dat <- simulateData(X = X,
                        B = rbind(bs$b0, bs$b1), distn = "Poisson", mean_count_before_ZI = 50)
    emu_res <- suppressWarnings(fastEmuTest(constraint_cats = 1:5, Y = dat, X = X, test_kj = data.frame(k = 2, j = 6),
                           estimate_full_model = FALSE))
    ps[i] <- emu_res$coef$pval
  }

  expect_true(mean(ps <= 0.05) < 0.05 + 1.96 * sqrt(0.05 * 0.95 / nsim))

})

test_that("fastEmu has power that increases with sample size and signal magnitude", {

  skip("skipping this in automatic tests because it is slow")

  n <- 20
  J <- 10
  set.seed(1569)
  nsim <- 100
  ps <- rep(NA, nsim)
  bs <- simulateBs(J = J, constraint_cats = 1:5, test_j = 6, constraint_mag = 1,
                   other_mag = 5, under_null = FALSE, alt_val = 0.1)
  X = cbind(1, rep(0:1, each = n /2))

  for (i in 1:nsim) {
    print(i)
    dat <- simulateData(X = X,
                        B = rbind(bs$b0, bs$b1), distn = "Poisson", mean_count_before_ZI = 50)
    emu_res <- suppressWarnings(fastEmuTest(constraint_cats = 1:5, Y = dat, X = X, test_kj = data.frame(k = 2, j = 6),
                           estimate_full_model = FALSE))
    ps[i] <- emu_res$coef$pval
  }

  n <- 50
  set.seed(1570)
  nsim <- 100
  new_ps <- rep(NA, nsim)
  bs <- simulateBs(J = J, constraint_cats = 1:5, test_j = 6, constraint_mag = 1,
                   other_mag = 5, under_null = FALSE, alt_val = 1)
  X = cbind(1, rep(0:1, each = n /2))

  for (i in 1:nsim) {
    print(i)
    dat <- simulateData(X = X,
                        B = rbind(bs$b0, bs$b1), distn = "Poisson", mean_count_before_ZI = 50)
    emu_res <- suppressWarnings(fastEmuTest(constraint_cats = 1:5, Y = dat, X = X, test_kj = data.frame(k = 2, j = 6),
                           estimate_full_model = FALSE))
    new_ps[i] <- emu_res$coef$pval
  }

  expect_true(mean(ps <= 0.05) < mean(new_ps <= 0.05))

})
