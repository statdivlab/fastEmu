n <- 50
J <- 10
b0 <- rnorm(J)
b1 <- rnorm(J)
X <- cbind(1, rep(0:1, each = n / 2))
Y <- simulateData(X, rbind(b0, b1), distn = "Poisson", mean_count_before_ZI = 50)
rownames(X) <- paste0("sample", 1:nrow(X))
rownames(Y) <- rownames(X)
colnames(Y) <- paste0("taxon", 1:ncol(Y))

test_that("fastEmuFit() gives nearly the same fit as fastEmuTest() for given reference set", {
  test_res <- suppressWarnings(fastEmuTest(constraint_cats = 1:5, Y = Y, X = X,
                                           test_kj = data.frame(k = 2, j = 6)))
  fit_res <- fastEmuFit(reference_set = 1:5, Y = Y, X = X, test_kj = data.frame(k = 2, j = 6))
  # check that estimates are approximately equal
  expect_true(all.equal(test_res$coef$estimate[6], fit_res$coef$estimate[6], 0.001))
  # check that score test p-values are approximately equal
  expect_true(all.equal(test_res$coef$pval[6], fit_res$coef$pval[6], 0.001))
})

test_that("fastEmuFit() works with data-driven reference set", {
  dd_ref_res <- fastEmuFit(Y = Y, X = X, run_score_tests = FALSE, reference_set_size = 5)
  rad_res <- radEmu::emuFit(Y = Y, X = X, run_score_tests = FALSE)
  # check that reference set is as expected
  rad_ref <- data.frame(ind = 1:ncol(Y), abs_est = abs(rad_res$coef$estimate))
  rad_ord <- order(rad_ref$abs_est)
  ref <- rad_ref$ind[rad_ord[1:5]]
  expect_equal(ref, dd_ref_res$reference_set)
})

test_that("when fastEmuFit() is run, warning from match_row_names = TRUE is only run once", {
  X1 <- X
  rownames(X1) <- NULL
  fit_message <- capture_messages(fastEmuFit(X = X1, Y = Y, test_kj = data.frame(k = 2, j = 6), reference_set_size = 5))
  expect_true(sum(fit_message == "Row names are missing from the covariate matrix X. We will assume the rows are in the same order as in the response matrix Y. You are responsible for ensuring the order of your observations is the same in both matrices.\n")
              == 1)
})

test_that("fastEmuFit() gives appropriate errors/warnings/messages related to reference set covariate, reference set number, and reference set vs constraint function mismatch", {
  # errors and messages related to reference set covariate
  expect_error(fastEmuFit(X = X, Y = Y, reference_set_covariate = "covariateA"))
  expect_error(fastEmuFit(X = X, Y = Y, reference_set_covariate = 3))
  expect_message(fit <- fastEmuFit(X = X, Y = Y, reference_set_covariate = 1, run_score_tests = FALSE, reference_set_size = 1))
  # message related to reference set size
  expect_message(fit <- fastEmuFit(X = X, Y = Y, run_score_tests = FALSE, reference_set_size = 11))
  # errors related to reference set
  expect_error(fastEmuFit(X = X, Y = Y, reference_set = 1:11))
  expect_error(fastEmuFit(X = X, Y = Y, reference_set = c("OTU1", "OTU2")))
  expect_error(fastEmuFit(X = X, Y = Y, reference_set = 1:2, constraint_fn = 3))
  expect_error(fastEmuFit(X = X, Y = Y, reference_set = 2, constraint_fn = 3))
  # check that correctly giving reference set as indices or column names provides the same results
  fit1 <- fastEmuFit(X = X, Y = Y, run_score_tests = FALSE, reference_set = 1:2)
  fit2 <- fastEmuFit(X = X, Y = Y, run_score_tests = FALSE, reference_set = c("taxon2", "taxon1"))
  expect_true(all.equal(fit1$coef, fit2$coef))
})

test_that("estimation is the same whether no fitted object, fastEmuFit fitted object, or radEmu fitted object", {
  rad_res <- radEmu::emuFit(Y = Y, X = X, run_score_tests = FALSE)
  fast_res <- fastEmuFit(Y = Y, X = X, reference_set_size = 5, run_score_tests = FALSE)
  fast_res_with_rad <- fastEmuFit(Y = Y, X = X, fitted_model = rad_res, run_score_tests = FALSE, reference_set_size = 5)
  fast_res_with_fast <- suppressMessages(fastEmuFit(Y = Y, X = X, fitted_model = fast_res, run_score_tests = FALSE, reference_set_size = 5))
  expect_true(all.equal(fast_res$coef, fast_res_with_rad$coef))
  expect_true(all.equal(fast_res$coef, fast_res_with_fast$coef))
  expect_true(all.equal(fast_res$reference_set, fast_res_with_rad$reference_set))
  expect_true(all.equal(fast_res$reference_set, fast_res_with_fast$reference_set))

  fast_res <- fastEmuFit(Y = Y, X = X, reference_set = 1:5, run_score_tests = FALSE)
  fast_res_with_rad <- fastEmuFit(Y = Y, X = X, fitted_model = rad_res, run_score_tests = FALSE, reference_set = 1:5)
  fast_res_with_fast <- suppressMessages(fastEmuFit(Y = Y, X = X, fitted_model = fast_res, run_score_tests = FALSE, reference_set = 1:5))
  expect_true(all.equal(fast_res$coef, fast_res_with_rad$coef))
  expect_true(all.equal(fast_res$coef, fast_res_with_fast$coef))
})

test_that("estimation is the same is use default constraints or provide them with constraint arguments", {
  constraint_fn <- (function(x) radEmu:::pseudohuber_center(x, d = 0.1))
  constraint_grad_fn <- (function(x) radEmu:::dpseudohuber_center_dx(x, d = 0.1))

  default_res <- fastEmuFit(Y = Y, X = X, reference_set_size = 5, run_score_tests = FALSE)
  set_res <- fastEmuFit(Y = Y, X = X, reference_set_size = 5, run_score_tests = FALSE,
                        constraint_fn = constraint_fn, constraint_grad_fn = constraint_grad_fn,
                        constraint_param = NA)
  expect_true(all.equal(default_res$coef, set_res$coef))
})

test_that("estimation removes arguments that we don't want (and throws warning)", {
  expect_warning(fast_res <- fastEmuFit(Y = Y, X = X, reference_set_size = 5, test_kj = data.frame(k = 2, j = 6),
                                        return_both_score_pvals = TRUE, use_fullmodel_cov = TRUE))
  expect_false("score_stat_full_info" %in% names(fast_res$coef))
})

test_that("multiple score tests can be run and results are given in the right order, and if null_B is provided that it is used", {
  res1 <- fastEmuFit(Y = Y, X = X, reference_set = 1:5, test_kj = data.frame(k = 2, j = 3), return_nullB = TRUE)
  res2 <- fastEmuFit(Y = Y, X = X, reference_set = 1:5, test_kj = data.frame(k = 2, j = c(10, 3)))
  expect_true(sum(!(is.na(res2$coef$pval))) == 2)
  expect_true(sum(!(is.na(res1$coef$pval))) == 1)
  expect_true(res1$coef$pval[3] == res2$coef$pval[3])

  res3 <- fastEmuFit(Y = Y, X = X, reference_set = 1:5, test_kj = data.frame(k = 2, j = 3),
                     B_null_list = list(res1$null_B[[1]]))
  expect_true(all.equal(res1$coef, res3$coef, 0.001))
})

test_that("with large number of taxa, fastEmu results with data-driven reference set results are similar to radEmu results with default constraint", {

  skip("skipping this in automatic tests because with large J, radEmu robust score tests are slow.")

  n <- 50
  J <- 500
  b0 <- rnorm(J)
  b1 <- rnorm(J)
  Y <- simulateData(X, rbind(b0, b1), distn = "Poisson", mean_count_before_ZI = 50)
  rownames(X) <- paste0("sample", 1:nrow(X))
  rownames(Y) <- rownames(X)
  colnames(Y) <- paste0("taxon", 1:ncol(Y))

  dd_ref_res <- fastEmuFit(Y = Y, X = X, test_kj = data.frame(k = 2, j = 1:5), verbose = TRUE, tolerance = 0.001, B_null_tol = 0.01)
  rad_res <- radEmu::emuFit(Y = Y, X = X, test_kj = data.frame(k = 2, j = 1:5), verbose = TRUE, tolerance = 0.001, B_null_tol = 0.01)
  expect_true(dd_ref_res$constraint_diff[2] < 0.05)
  # there may be some larger differences, but this is helpful to look at
  # all.equal(dd_ref_res$coef$pval[1:5], rad_res$coef$pval[1:5], 0.25))
})

test_that("fastEmuFit controls Type I error rate, with given reference set", {

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
    emu_res <- fastEmuFit(reference_set = 1:5, Y = dat, X = X, test_kj = data.frame(k = 2, j = 6),
                                           match_row_names = FALSE)
    ps[i] <- emu_res$coef$pval[6]
  }

  expect_true(mean(ps <= 0.05) < 0.05 + 1.96 * sqrt(0.05 * 0.95 / nsim))

})

test_that("fastEmuFit controls Type I error rate, with data-driven reference set", {

  skip("skipping this in automatic tests because it is slow")

  n <- 40
  J <- 100
  set.seed(1569)
  nsim <- 100
  ps <- rep(NA, nsim)
  bs <- simulateBs(J = J, constraint_cats = 1:5, test_j = 6, constraint_mag = 5,
                   other_mag = 5, under_null = TRUE)
  X = cbind(1, rep(0:1, each = n /2))

  for (i in 1:nsim) {

    dat <- simulateData(X = X,
                        B = rbind(bs$b0, bs$b1), distn = "Poisson", mean_count_before_ZI = 50)
    emu_res <- fastEmuFit(Y = dat, X = X, test_kj = data.frame(k = 2, j = 6),
                                           match_row_names = FALSE)
    ps[i] <- emu_res$coef$pval[6]
  }

  expect_true(mean(ps <= 0.05) < 0.05 + 1.96 * sqrt(0.05 * 0.95 / nsim))

})

test_that("fastEmuFit has power that increases with sample size and signal magnitude", {

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


