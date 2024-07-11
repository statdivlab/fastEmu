n <- 50
J <- 10
b0 <- rnorm(J)
b1 <- rnorm(J)
X <- cbind(1, rep(0:1, each = n / 2))
Y <- simulateData(X, rbind(b0, b1), distn = "Poisson", mean_count_before_ZI = 50)

test_that("radEmu and fastEmu give same results when using full model and all categories in constraint", {

  radEmu_res <- emuFit(Y, X, test_kj = data.frame(k = 2, j = 6))
  fastEmu_res <- fastEmuTest(constraint_cats = 1:10, Y = Y, X = X,
                             test_kj = data.frame(k = 2, j = 6), model = "full")

  expect_true(all.equal(radEmu_res$coef, fastEmu_res$coef[, -10]))
})

test_that("fastEmu runs with and without estimation, and results don't change", {

  colnames(Y) <- paste0("taxon", 1:10)
  expect_silent({res_est <- fastEmuTest(constraint_cats = 1:5, Y = Y, X = X,
                            test_kj = data.frame(k = 2, j = 6), estimate_full_model = TRUE)})
  expect_silent({res_no_est <- fastEmuTest(constraint_cats = 1:5, Y = Y, X = X,
                            test_kj = data.frame(k = 2, j = 6), estimate_full_model = FALSE)})
  expect_true(res_est$coef$pval[6] == res_no_est$coef$pval)

  # also check we're returning what we want

  # all constraint categories are in "included_categories"
  expect_true(sum(paste0("taxon", 1:5) %in% unlist(res_est$included_categories)) == 5)
  # categories in coef are in the same order as they are in Y
  expect_true(all.equal(colnames(Y), res_est$coef$category))

})

test_that("fastEmu can work with additional arguments", {

  expect_silent({res <- fastEmuTest(constraint_cats = 1:5, Y = Y, X = X, test_kj = data.frame(k = 2, j = 6),
                            penalize = FALSE, return_nullB = TRUE)})

  # check that penalty is FALSE
  expect_false(res$penalized)

  # check that nullB is returned
  expect_true("null_B" %in% names(res))

})

test_that("fastEmu still works with deprecated aggregated model", {

  expect_silent({res_agg <- fastEmuTest(constraint_cats = 1:5, Y = Y, X = X, test_kj = data.frame(k = 2, j = 6),
                            model = "agg")})

  res_drop <- fastEmuTest(constraint_cats = 1:5, Y = Y, X = X, test_kj = data.frame(k = 2, j = 6))
  expect_true(all.equal(res_agg$coef$pval[6], res_drop$coef$pval[6], tol = 0.001))

})
