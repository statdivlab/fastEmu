#' Fit radEmu model after collapsing data into smaller joint model, runs score tests
#'
#' @param model either "full" to run score tests with the full joint model, "drop" to run score tests with
#' a conditional joint model, or "agg" to run score tests with a joint marginal model
#' @param constraint_cats a vector of category indices for categories that should be involved in the constraint
#' and retained in any smaller joint model
#' @param Y an n x J matrix or dataframe of nonnegative observations, or a phyloseq object containing an otu table and sample data.
#' @param X an n x p matrix or dataframe of covariates (optional)
#' @param formula a one-sided formula specifying the form of the mean model to be fit
#' @param data an n x p data frame containing variables given in \code{formula}
#' @param cluster a numeric vector giving cluster membership for each row of Y to
#' be used in computing GEE test statistics. Default is NULL, in which case rows of
#' Y are treated as independent.
#' @param penalize logical: should Firth penalty be used in fitting model? Default is TRUE.
#' @param B starting value of coefficient matrix (p x J). If not provided,
#' B will be initiated as a zero matrix.
#' @param fitted_model a fitted model produced by a separate call to emuFit; to
#' be provided if score tests are to be run without refitting the full unrestricted model.
#' Default is NULL.
#' @param refit logical: if B or fitted_model is provided, should full model be fit (TRUE) or
#' should fitting step be skipped (FALSE), e.g., if score tests are to be run on an already
#' fitted model. Default is TRUE.
#' @param test_kj a data frame whose rows give coordinates (in category j and
#' covariate k) of elements of B to construct hypothesis tests for. If \code{test_kj}
#' is not provided, all elements of B save the intercept row will be tested.
#' @param alpha nominal type 1 error level to be used to construct confidence intervals. Default is 0.05
#' (corresponding to 95% confidence intervals)
#' @param return_wald_p logical: return p-values from Wald tests? Default is FALSE.
#' @param compute_cis logical: compute and return Wald CIs? Default is TRUE.
#' @param run_score_tests logical: perform robust score testing? Default is TRUE.
#' @param use_fullmodel_info logical: TODO? Default is FALSE.
#' @param use_fullmodel_cov logical: use information matrix and empirical score covariance
#' computed for full model fit? Defaults to FALSE, in which case these quantities are
#' recomputed for each null model fit for score testing.
#' @param use_both_cov logical: should score tests be run using information and
#' empirical score covariance evaluated both under the null and full models?
#' Used in simulations
#' @param verbose provide updates as model is being fitted? Defaults to TRUE.
#' @param tolerance tolerance for stopping criterion in full model fitting; once
#' no element of B is updated by more than this value in a single step, we exit
#' optimization. Defaults to 1e-3.
#' @param rho_init numeric: value at which to initiate rho parameter in augmented Lagrangian
#' algorithm. Default is 1.
#' @param tau numeric: value to scale rho by in each iteration of augmented Lagrangian
#' algorithm that does not move estimate toward zero sufficiently. Default is 2.
#' @param kappa numeric: value between 0 and 1 that determines the cutoff on the ratio
#' of current distance from feasibility over distance in last iteration triggering
#' scaling of rho. If this ratio is above kappa, rho is scaled by tau to encourage
#' estimate to move toward feasibility.
#' @param constraint_tol numeric: constraint tolerance for fits under null hypotheses
#' (tested element of B must be equal to constraint function to within this tolerance for
#' a fit to be accepted as a solution to constrained optimization problem). Default is 1e-5.
#' @param maxit maximum number of outer iterations of augmented lagrangian algorithm to perform before
#' exiting optimization. Default is 1000.
#' @param inner_maxit maximum number of coordinate descent passes through columns of B to make within each
#' outer iteration of augmented lagrangian algorithm before exiting inner loop
#' @param max_step maximum stepsize; update directions computed during optimization
#' will be rescaled if a step in any parameter exceeds this value. Defaults to 0.5.
#' @param B_null_tol numeric: convergence tolerance for null model fits for score testing (if max of absolute difference in
#' B across outer iterations is below this threshold, we declare convergence).
#' Default is 0.01.
#' @param inner_tol numeric: convergence tolerance for augmented Lagrangian subproblems
#' within null model fitting. Default is 1.
#' @param ntries numeric: how many times should optimization be tried in null
#' models where at least one optimization attempt fails? Default is 4.
#' @param c1 numeric: parameter for Armijo line search. Default is 1e-4.
#' @param trackB logical: should values of B be recorded across optimization
#' iterations and be returned? Primarily used for debugging. Default is FALSE.
#' @param return_both_score_pvals logical: should score p-values be returned using both
#' information matrix computed from full model fit and from null model fits? Default is
#' FALSE. This parameter is used for simulations - in any applied analysis, type of
#' p-value to be used should be chosen before conducting tests.
#'
#' @return A list containing elements 'coef', 'B', 'penalized', 'Y_augmented',
#' 'I', and 'Dy'.  Parameter estimates by
#' covariate and outcome category (e.g., taxon for microbiome data), as well as
#' optionally confidence intervals and p-values, are contained in 'coef'. 'B'
#' contains parameter estimates in matrix format (rows indexing covariates and
#' columns indexing outcome category / taxon). 'penalized' is equal to TRUE
#' if Firth penalty is used in estimation (default) and FALSE otherwise. 'I' and
#' 'Dy' contain an information matrix and empirical score covariance matrix
#' computed under the full model.
#'
#' @importFrom stats cov median model.matrix optim pchisq qnorm weighted.mean
#' @import Matrix
#' @import MASS
#' @import radEmu
#'
#' @export
#'
fastEmuTest <- function(model = "full",
                        constraint_cats,
                        Y,
                        X = NULL,
                        formula = NULL,
                        data = NULL,
                        cluster = NULL,
                        penalize = TRUE,
                        B = NULL,
                        fitted_model = NULL,
                        refit = TRUE,
                        test_kj,
                        alpha = 0.05,
                        return_wald_p = FALSE,
                        compute_cis = TRUE,
                        run_score_tests = TRUE,
                        use_fullmodel_info = FALSE,
                        use_fullmodel_cov = FALSE,
                        use_both_cov = FALSE,
                        #constraint_fn = pseudohuber_center,
                        #constraint_grad_fn = dpseudohuber_center_dx,
                        #constraint_param = 0.1,
                        verbose = FALSE,
                        tolerance = 1e-4,
                        B_null_tol = 1e-3,
                        rho_init = 1,
                        inner_tol = 1,
                        ntries = 4,
                        tau = 2,
                        kappa = 0.8,
                        constraint_tol = 1e-5,
                        c1 = 1e-4,
                        maxit = 1000,
                        inner_maxit = 25,
                        max_step = 1,
                        trackB = FALSE,
                        return_both_score_pvals = FALSE) {

  if (sum(rowSums(Y) == 0) > 0) {
    stop("There is at least one sample with no counts in any category. Please remove samples that have no counts in any category.")
  }

  # currently only implemented for a single category
  if (length(unique(test_kj$j)) > 1) {
    stop("Currently testing is not implemented for more than one category at once")
  }
  # # currently only implemented for testing j not in constraint set
  # if (test_kj$j[1] %in% constraint_cats) {
  #   stop("Currently testing is not implemented for categories included in the constraint.")
  # }

  # add category names if they aren't already there
  if (is.null(colnames(Y))) {
    colnames(Y) <- paste0("category", 1:ncol(Y))
  }
  test_cat <- colnames(Y)[test_kj$j[1]]

  # set constraints
  # reorder categories in Y so that constraint categories are in the beginning
  ind_keep <- c(constraint_cats, test_kj$j[1])
  ind_keep <- unique(ind_keep)
  new_order <- c(ind_keep, (1:ncol(Y))[-ind_keep])
  # set constraint as single category constraint if there is one value in constraint_cats
  if (length(constraint_cats) == 1) {
    constraint_fn_sub <- (function(x) x[1])
    constraint_grad_fn_sub <- function(x) {
      grad <- c(1, rep(0, length(x) - 1))
    }
  # set constraint as pseudo-Huber over constraint categories if there are multiple
  } else {
    constraint_fn_sub <- (function(x) {
      radEmu:::pseudohuber_center(x[1:length(constraint_cats)], d = .1)
      })
    constraint_grad_fn_sub <- (function(x) {
      grad <- rep(0, length(x))
      grad[1:length(constraint_cats)] <-
        radEmu:::dpseudohuber_center_dx(x[1:length(constraint_cats)], d = .1)
      return(grad)
    })
  }

  added_cats <- NULL

  # check for separation in categorical data, cases in which the baseline category and
  # category k both only have zeros
  if (is.null(X)) {
    if (is.null(formula) | is.null(data)) {
      stop("If design matrix X not provided, both formula and data containing\ncovariates in formula must be provided.")
    }
    X <- model.matrix(formula, data)
  }
  if ("data.frame" %in% class(X)) {
    X <- as.matrix(X)
    if (!is.numeric(X)) {
      stop("X is a data frame that cannot be coerced to a numeric matrix. Please fix and try again.")
    }
  }


  # update model based on model chosen and categories to use
  if (model == "full") {
    mod_Y <- Y[, new_order]
  } else if (model == "drop") {
    mod_Y <- Y[, ind_keep]
    if (sum(rowSums(mod_Y) == 0) > 0) {
      message("There is at least one observation with no counts across the categories included in this model. Additional categories will be added to the model")
      count_df <- data.frame(id = (1:ncol(Y))[-ind_keep],
                             counts = colSums(Y)[-ind_keep])
      count_df$rank <- rank(-count_df$counts, ties.method = "first")
      added_cats <- c()
      curr_rank <- 1
      while (sum(rowSums(mod_Y) == 0) > 0) {
        cat <- count_df$id[count_df$rank == curr_rank]
        ind_keep <- c(ind_keep, cat)
        mod_Y <- Y[, ind_keep]
        added_cats <- c(added_cats, cat)
        curr_rank <- curr_rank + 1
      }
    }
  } else if (model == "agg") {
    other_cat_sum <- rowSums(Y[, -ind_keep])
    mod_Y <- cbind(Y[, ind_keep], other_cat_sum)
  } else {
    stop("Please use 'full', 'drop', or 'agg' for method argument.")
  }

  emuObj <- emuFit(Y = mod_Y,
                   X = X,
                   formula = formula,
                   data = data,
                   cluster = cluster,
                   penalize = penalize,
                   B = B,
                   fitted_model = fitted_model,
                   refit = refit,
                   test_kj = data.frame(k = test_kj$k,
                                        j = which(new_order == test_kj$j[1])),
                   alpha = alpha,
                   return_wald_p = return_wald_p,
                   compute_cis = compute_cis,
                   run_score_tests = run_score_tests,
                   use_fullmodel_info = use_fullmodel_info,
                   use_fullmodel_cov = use_fullmodel_cov,
                   use_both_cov = use_both_cov,
                   constraint_fn = constraint_fn_sub,
                   constraint_grad_fn = constraint_grad_fn_sub,
                   constraint_param = NA,
                   verbose = verbose,
                   tolerance = tolerance,
                   B_null_tol = B_null_tol,
                   rho_init = rho_init,
                   inner_tol = inner_tol,
                   ntries = ntries,
                   tau = tau,
                   kappa = kappa,
                   constraint_tol = constraint_tol,
                   c1 = c1,
                   maxit = maxit,
                   inner_maxit = inner_maxit,
                   max_step = max_step,
                   trackB = trackB,
                   return_both_score_pvals = return_both_score_pvals)

  res <- emuObj
  res$model <- model
  pval_cols <- which(stringr::str_detect(names(emuObj$coef), "pval") |
                       stringr::str_detect(names(emuObj$coef), "_p"))
  covs <- unique(emuObj$coef$covariate)
  res$p_vals <- emuObj$coef[emuObj$coef$category == test_cat &
                              emuObj$coef$covariate %in% covs[test_kj$k - 1], pval_cols,
                            drop = FALSE]
  if (!is.null(added_cats)) {
    res$added_cats <- added_cats
  }
  return(res)
}
