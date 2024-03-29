test_that("test-fastEmuTest works", {

  library(fastEmu)

  # simulate data
  n <- 100
  J <- 200
  b0 <- rnorm(J)
  b1 <- rnorm(J)
  X <- cbind(1, rep(c(0, 1), each = n/2),
             rnorm(n))
  Y <- radEmu:::simulate_data(n = n,
                              J = J,
                              b0 = b0,
                              b1 = b1,
                              distn = "Poisson",
                              zinb_size = 5,
                              zinb_zero_prop = 0.6,
                              mean_count_before_ZI = 50)
  colnames(Y) <- paste0("category", 1:ncol(Y))

  # test 100th category
  # randomly select 25 categories to use as constraint
  cats <- sample((1:ncol(Y))[-100], 25)

  # run full model, cat categories as constaint, test 100th
  res_full <- fastEmuTest(model = "full",
                          constraint_cats = cats,
                          Y = Y,
                          X = X,
                          test_kj = data.frame(k = 2:3, j = 100),
                          tau = 2,
                          B_null_tol = 0.005,
                          tolerance = 0.005,
                          constraint_tol = 0.001,
                          return_wald_p = TRUE,
                          use_both_cov = FALSE,
                          use_fullmodel_info = TRUE,
                          return_both_score_pvals = TRUE)
  expect_false(is.na(res_full$p_vals$score_pval_null_info[1]))

  # run dropped model
  # run dropped model, cat categories as constraint, test 100th
  res_drop <- fastEmuTest(model = "drop",
                         constraint_cats = cats,
                         Y = Y,
                         X = X,
                         test_kj = data.frame(k = 2, j = 100),
                         tau = 2,
                         B_null_tol = 0.005,
                         tolerance = 0.005,
                         constraint_tol = 0.001,
                         return_wald_p = TRUE,
                         use_both_cov = FALSE,
                         use_fullmodel_info = TRUE,
                         return_both_score_pvals = TRUE)
  expect_false(is.na(res_drop$p_vals$score_pval_null_info))

  # run aggregated model
  # run aggregated model, cat categories as constaint, test 100th
  res_agg <- fastEmuTest(model = "agg",
                          constraint_cats = cats,
                          Y = Y,
                          X = X,
                          test_kj = data.frame(k = 2, j = 100),
                          tau = 2,
                          B_null_tol = 0.005,
                          tolerance = 0.005,
                          constraint_tol = 0.001,
                          return_wald_p = TRUE,
                          use_both_cov = FALSE,
                          use_fullmodel_info = TRUE,
                          return_both_score_pvals = TRUE)
  expect_false(is.na(res_agg$p_vals$score_pval_null_info))

  # the correct number of categories are included for each model
  expect_true(nrow(res_full$coef) == 2*ncol(Y) &
                nrow(res_drop$coef) == 2*(1 + length(cats)) &
                nrow(res_agg$coef) == 2*(2 + length(cats)))

  # we can add other categories with the categories_in_model argument
  cats_in_model <- rep(2, ncol(Y))
  cats_in_model[cats] <- 0
  cats_in_model[100] <- 1
  cats_to_add <- sample((1:ncol(Y))[-c(100, cats)], 25)
  cats_in_model[cats_to_add] <- 1
  # run dropped model, cat categories as constaint, test 100th, add in more categories
  res_drop_more <- fastEmuTest(model = "drop",
                               constraint_cats = cats,
                               categories_in_model = cats_in_model,
                               Y = Y,
                               X = X,
                               test_kj = data.frame(k = 2, j = 100),
                               tau = 2,
                               B_null_tol = 0.005,
                               tolerance = 0.005,
                               constraint_tol = 0.001,
                               return_wald_p = TRUE,
                               use_both_cov = FALSE,
                               use_fullmodel_info = TRUE,
                               return_both_score_pvals = TRUE)
  expect_true(nrow(res_drop_more$coef) == 2*(1 + length(cats) + length(cats_to_add)))
  expect_false(is.na(res_drop_more$p_vals$score_pval_null_info))

})
