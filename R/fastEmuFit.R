#' A fast approximation to \code{emuFit} from \code{radEmu}.
#'
#' @param reference_set The reference set to use in the identifiability constraint.
#' The user can input a reference set as a vector of numbers that represent indices
#' for columns of the \code{Y} matrix, or names that correspond with column names of
#' the \code{Y} matrix. If a reference set is not provided, by default, this is set
#' to \code{data_driven}, and \code{fastEmuFit} will identify a reference set of typical
#' taxa of size \code{reference_set_size}. If \code{data_driven_ss} or
#' \code{data_driven_thin}, a data-driven reference set will be determined using sample
#' splitting or Poisson thinning respectively.
#' @param reference_set_size The size of the reference set if it is data-driven, default
#' is set to \code{50}. We recommend a reference set of size 30-100 for the best balance
#' of computational efficiency and estimation precision.
#' @param reference_set_covariate If the reference set is data-driven, which covariates should
#' it be chosen relative to. By default, this will be all covariates in the model (ignoring
#' the intercept). However, if a model includes a main covariate that will be tested and several
#' precision variables, we recommend choosing the reference set with respect to the main
#' covariate of interest. This argument should be a vector of numbers that correspond to
#' column indices in the \code{X} design matrix. If you don't know which columns correspond to
#' each covariate in your design matrix, run the function \code{radEmu::make_design_matrix()} to
#' see the design matrix for your model.
#' @param Y an n x J matrix or dataframe of nonnegative observations, or a \code{phyloseq}
#' or \code{TreeSummarizedExperiment} object containing an otu table and sample data.
#' @param X an n x p matrix or dataframe of covariates (optional, either include \code{X}
#' or \code{formula} and \code{data})
#' @param formula a one-sided formula specifying the form of the mean model to be fit
#' @param data an n x p data frame containing variables given in \code{formula}
#' @param test_kj a data frame whose rows give coordinates (in category j and
#' covariate k) of elements of B to construct hypothesis tests for. If \code{test_kj}
#' is not provided, all elements of B save the intercept row will be tested.
#' @param cluster a numeric vector giving cluster membership for each row of Y to
#' be used in computing GEE test statistics. Default is NULL, in which case rows of
#' Y are treated as independent.
#' @param penalize logical: should Firth penalty be used in fitting model? Default is TRUE.
#' @param B starting value of coefficient matrix (p x J). If not provided,
#' B will be initiated as a zero matrix.
#' @param fitted_model a fitted model produced by a call to fastEmu::fastEmuFit or radEmu::emuFit;
#' to be provided if score tests are to be run without refitting the full unrestricted model.
#' Default is NULL.
#' @param refit logical: if B or fitted_model is provided, should full model be fit (TRUE) or
#' should fitting step be skipped (FALSE), e.g., if score tests are to be run on an already
#' fitted model. Default is TRUE.
#' @param return_wald_p logical: return p-values from Wald tests? Default is FALSE. These can only be
#' returned if \code{estimate_full_model} is TRUE.
#' @param compute_cis logical: compute and return Wald CIs? Default is TRUE. These can only be
#' returned if \code{estimate_full_model} is TRUE.
#' @param verbose provide updates as model is being fitted and score tests are run? Defaults to FALSE.
#' @param ... Additional arguments to radEmu:::emuFit. See possible arguments with \code{?radEmu::emuFit}.
#'
#' @return A list that includes all elements of an \code{emuFit} object from \code{radEmu::emuFit()}, as
#' well as additional elements. See the documentation in \code{?radEmu::emuFit} for a full description of the
#' elements in an \code{emuFit} object.The \code{emuFit} object includes the matrix \code{coef}, which provides
#' estimates for all parameters and score statistics and p-values for all parameters that were tested.
#' The returned object also includes \code{reference_set} and \code{reference_set_names}, which give the
#' indices of the reference set in terms of columns of the \code{Y} matrix and category names respectively,
#' of the categories (taxa) that were used as a reference set of "typical taxa" for the identifiability
#' constraint. Other elements of the list correspond to score tests. \code{included_categories} gives the
#' set of categories used for the reduced model for each score test, \code{score_test_hyperparams} provides
#' the hyperparameters related to estimation under the null hypothesis for each score test. If \code{return_null_B}
#' or \code{return_score_components} were set to \code{TRUE}, then \code{null_B} or \code{score_components}
#' will also be returned, which respectively give the estimated B values under the null hypothesis and the
#' components of the robust score test that are run, for each score test.
#'
#' @importFrom stats model.matrix
#' @import radEmu
#'
#' @export
#'
fastEmuFit <- function(reference_set = "data_driven",
                       reference_set_size = 50,
                       reference_set_covariate = NULL,
                       Y,
                       X = NULL,
                       formula = NULL,
                       data = NULL,
                       test_kj = NULL,
                       cluster = NULL,
                       penalize = TRUE,
                       B = NULL,
                       fitted_model = NULL,
                       refit = TRUE,
                       return_wald_p = FALSE,
                       compute_cis = TRUE,
                       verbose = FALSE,
                       ...) {

  # record call
  call <- match.call(expand.dots = FALSE)

  # record extra arguments
  extra_args <- list(...)

  # if unsupported arguments included, warn user
  unsupported <- c("use_both_cov", "use_fullmodel_info",
                   "use_fullmodel_cov", "return_both_score_pvals")
  if (sum(names(extra_args) %in% unsupported) > 0) {
    warning("You have used one of the following arguments: `use_both_cov`, `use_fullmodel_info`, `use_fullmodel_cov`, or `return_both_score_pvals`. These arguments are not supported in `fastEmu` and will be ignored.")
  }

  # start by running radEmu checks
  if (!("assay_name" %in% names(extra_args))) {
    assay_name <- NULL
  } else {
    assay_name <- extra_args$assay_name
  }
  if (!("B_null_list" %in% names(extra_args))) {
    B_null_list <- NULL
  } else {
    B_null_list <- extra_args$B_null_list
  }
  if (!("match_row_names" %in% names(extra_args))) {
    match_row_names <- TRUE
  } else {
    match_row_names <- extra_args$match_row_names
  }
  if (!("remove_zero_comparison_pvals" %in% names(extra_args))) {
    remove_zero_comparison_pvals <- 0.01
  } else {
    remove_zero_comparison_pvals <- extra_args$remove_zero_comparison_pvals
  }
  if (!("unobserved_taxon_error" %in% names(extra_args))) {
    unobserved_taxon_error <- TRUE
  } else {
    unobserved_taxon_error <- extra_args$unobserved_taxon_error
  }
  if (!("constraint_fn" %in% names(extra_args))) {
    constraint_fn <- radEmu:::pseudohuber_center
  } else {
    constraint_fn <- extra_args$constraint_fn
  }
  if (!("constraint_grad_fn" %in% names(extra_args))) {
    constraint_grad_fn <- radEmu:::dpseudohuber_center_dx
  } else {
    constraint_grad_fn <- extra_args$constraint_grad_fn
  }
  if (!("constraint_param" %in% names(extra_args))) {
    constraint_param <- 0.1
  } else {
    constraint_param <- extra_args$constraint_param
  }
  # run score tests
  if ("run_score_tests" %in% names(extra_args)) {
    run_score_tests <- extra_args$run_score_tests
  } else {
    run_score_tests <- TRUE
  }
  check_results <- radEmu:::emuFit_check(Y = Y,
                                         X = X,
                                         formula = formula,
                                         data = data,
                                         assay_name = assay_name,
                                         cluster = cluster,
                                         B_null_list = B_null_list,
                                         test_kj = test_kj,
                                         match_row_names = match_row_names,
                                         verbose = verbose,
                                         remove_zero_comparison_pvals = remove_zero_comparison_pvals,
                                         unobserved_taxon_error = unobserved_taxon_error,
                                         constraint_fn = constraint_fn,
                                         constraint_grad_fn = constraint_grad_fn,
                                         constraint_param = constraint_param,
                                         run_score_tests = run_score_tests)
  Y <- check_results$Y
  X <- check_results$X
  cluster <- check_results$cluster
  B_null_list <- check_results$B_null_list
  n <- nrow(X)
  J <- ncol(Y)
  p <- ncol(X)
  constraint_fn <- check_results$constraint_fn
  constraint_grad_fn <- check_results$constraint_grad_fn
  constraint_param <- check_results$constraint_param

  # check for valid reference set covariate (if data-driven)
  if (length(reference_set) == 1 && reference_set == "data_driven") {
    if (is.null(reference_set_covariate)) {
      reference_set_covariate <- 2:p
    } else {
      if (!is.numeric(reference_set_covariate)) {
        stop("Make sure that `reference_set_covariate` is a numeric value or vector of numeric values. If you don't know which column of the design matrix corresponds to the covariate(s) you want to use, use `radEmu::make_design_matrix()` to view your design matrix.")
      }
      if (1 %in% reference_set_covariate) {
        message("You've included the intercept in your set of covariates to use to determine the data-driven reference set. While we don't recommend this by default, we'll let you do this and hope that you have a good reason!")
      }
      if (sum(!(reference_set_covariate %in% 1:p)) > 0) {
        stop("You've included one or more values in `reference_set_covariate` that is larger than the number of columns of the design matrix. Please try again, and if you don't know which column of the design matrix corresponds to the covariate(s) you want to use, use `radEmu::make_design_matrix()` to view your design matrix.")
      }
    }
  }

  # check for valid reference set
  if (length(reference_set) == 1 && reference_set %in% c("data_driven_thin", "data_driven_ss")) {
    stop("Sarah hasn't implemented data-driven reference sets with sample splitting or thinning yet. If you need this, please open a Github issue.")
  } else if (length(reference_set) == 1 && reference_set == "data_driven") {
    reference_set <- "data_driven"
  } else {
    if (is.numeric(reference_set)) {
      if (sum(!(reference_set %in% 1:ncol(Y))) > 0) {
        stop("At least one index provided in the reference set is larger than the number of taxa included in `Y`. Please double-check your reference set.")
      }
      if (!(is.null(colnames(Y)))) {
        reference_set_names <- colnames(Y)[reference_set]
      } else {
        reference_set_names <- NULL
      }
    } else if (is.character(reference_set)) {
      if (sum(!(reference_set %in% colnames(Y))) > 0) {
        stop("You have provided a character vector for the reference set, but at least one character string does not match a column name of `Y`. Please double-check your reference set.")
      }
      reference_set_names <- reference_set
      reference_set <- which(colnames(Y) %in% reference_set)
    }
  }

  # check on constraint
  if (!("constraint_param" %in% names(extra_args))) {
    constraint_param <- 0.1
  } else {
    constraint_param <- extra_args$constraint_param
  }
  if (!("constraint_fn" %in% names(extra_args))) {
    constraint_fn <- (function(x) radEmu:::pseudohuber_center(x, d = constraint_param))
  } else {
    constraint_fn <- extra_args$constraint_fn
  }
  if (!("constraint_grad_fn" %in% names(extra_args))) {
    constraint_grad_fn <- (function(x) radEmu:::dpseudohuber_center_dx(x, d = constraint_param))
  } else {
    constraint_grad_fn <- extra_args$constraint_grad_fn
  }
  if (length(reference_set) == 1 && reference_set == "data-driven" & is.numeric(constraint_fn)) {
    stop("You've set a data-driven reference set, but chosen a constraint function that consists of a single reference category. Please rerun this function, and set `reference_set` to be the index of the taxa that you would like to the be the reference category, or remove the `constraint_fn` argument.")
  } else if (is.numeric(constraint_fn) & is.numeric(reference_set)) {
    if (length(reference_set) > 1 || reference_set != constraint_fn) {
      stop("You've set a single reference taxon with the `constraint_fn` argument but this doesn't match the argument `reference_set`. Please remove one of these arguments.")
    }
  }

  # if reference set is data_driven, determine this reference set
  if (length(reference_set) == 1 && reference_set == "data_driven") {
    if (verbose %in% c(TRUE, "development")) {
      message("Running estimation algorithm to determine reference set.")
    }
    if (!(is.null(fitted_model)) & "reference_set" %in% names(fitted_model)) {
      message("You have set `reference_set = data_driven` but provided a `fitted_model` that already includes a `reference_set`. This reference set will be used. If you would like the reference set to be recalculated, set `fitted_model$reference_set <- NULL` and rerun this function.")
      reference_set <- fitted_model$reference_set
      reference_set_names <- fitted_model$reference_set_names
      if ("constraint_diff" %in% names(fitted_model)) {
        constraint_diff <- fitted_model$constraint_diff
      } else {
        constraint_diff <- NULL
      }
    } else {
      # remove arguments from ... that we don't want for estimation for reference set or for re-estimation for confidence intervals
      est_args_rm <- which(names(extra_args) %in%
                             c("run_score_tests", "use_both_cov", "use_fullmodel_info",
                               "use_fullmodel_cov", "return_both_score_pvals", "match_row_names"))
      est_args <- list(Y = Y,
                       X = X,
                       cluster = cluster,
                       penalize = penalize,
                       B = B,
                       fitted_model = fitted_model,
                       refit = refit,
                       test_kj = test_kj,
                       return_wald_p = FALSE,
                       compute_cis = FALSE,
                       verbose = ifelse(verbose == "development", "development", FALSE),
                       run_score_tests = FALSE,
                       use_both_cov = FALSE,
                       use_fullmodel_info = FALSE,
                       use_fullmodel_cov = FALSE,
                       return_both_score_pvals = FALSE,
                       match_row_names = FALSE)
      if (length(est_args_rm) > 0) {
        modified_extra_args <- extra_args[-est_args_rm]
      } else {
        modified_extra_args <- extra_args
      }
      est_args <- c(est_args, modified_extra_args)
      if (is.null(fitted_model) & is.null(B)) {

        # estimate with radEmu::emuFit
        fitted_model <- do.call(emuFit, est_args)
      }

      # choose reference set
      if (is.null(fitted_model)) {
        fitted_model$B <- B
      }
      ref_set_res <- chooseRefSet(fitted_model = fitted_model,
                                  reference_set_size = reference_set_size,
                                  reference_set_covariate = reference_set_covariate,
                                  constraint_fn = constraint_fn)
      reference_set <- ref_set_res$reference_set
      reference_set_names <- colnames(Y)[reference_set]
      constraint_diff <- ref_set_res$constraint_diff
      new_B <- ref_set_res$new_B

      # rerun emuFit (using new B) to update coef data frame
      # update constraint functions based on reference set
      constraint_fn_est <- (function(x) {
        constraint_fn(x[reference_set])
      })
      constraint_grad_fn_est <- (function(x) {
        grad <- rep(0, length(x))
        grad[reference_set] <-
          constraint_grad_fn(x[reference_set])
        return(grad)
      })
      re_est_args <- est_args
      re_est_args$constraint_fn <- constraint_fn_est
      re_est_args$constraint_param <- NA
      re_est_args$constraint_grad_fn <- constraint_grad_fn_est
      re_est_args$fitted_model <- NULL
      re_est_args$B <- new_B
      re_est_args$refit <- FALSE
      re_est_args$compute_cis <- compute_cis
      re_est_args$return_wald_p <- return_wald_p
      if (verbose == "development") {
        message("Optionally computing confidence intervals and Wald p-values")
      }
      fitted_model <- do.call(emuFit, re_est_args)
    }

  } else {
    constraint_diff <- NULL

    # update constraint functions based on reference set
    constraint_fn_est <- (function(x) {
      constraint_fn(x[reference_set])
    })
    constraint_grad_fn_est <- (function(x) {
      grad <- rep(0, length(x))
      grad[reference_set] <-
        constraint_grad_fn(x[reference_set])
      return(grad)
    })

    # check if fitted model (with reference set) exists, otherwise fit model
    if (is.null(fitted_model) & is.null(B)) {
      if (verbose %in% c(TRUE, "development")) {
        message("Running estimation algorithm.")
      }

      # remove arguments from ... that we don't want for estimation for reference set
      est_args_rm <- which(names(extra_args) %in%
                             c("run_score_tests", "use_both_cov", "use_fullmodel_info",
                               "use_fullmodel_cov", "return_both_score_pvals", "match_row_names",
                               "constraint_fn", "constraint_grad_fn", "constraint_param"))
      est_args <- list(Y = Y,
                       X = X,
                       cluster = cluster,
                       penalize = penalize,
                       B = B,
                       fitted_model = fitted_model,
                       refit = refit,
                       test_kj = test_kj,
                       verbose = ifelse(verbose == "development", "development", FALSE),
                       run_score_tests = FALSE,
                       constraint_fn = constraint_fn_est,
                       constraint_grad_fn = constraint_grad_fn_est,
                       constraint_param = NA,
                       use_both_cov = FALSE,
                       use_fullmodel_info = FALSE,
                       use_fullmodel_cov = FALSE,
                       return_both_score_pvals = FALSE,
                       match_row_names = FALSE)
      if (length(est_args_rm) > 0) {
        modified_extra_args <- extra_args[-est_args_rm]
      } else {
        modified_extra_args <- extra_args
      }
      est_args <- c(est_args, modified_extra_args)
      fitted_model <- do.call(emuFit, est_args)
    } else {
      if (is.null(fitted_model)) {
        fitted_model$B <- B
      }
      if (!("reference_set" %in% names(fitted_model))) {
        # check whether constraint has already been set
        if (all.equal(constraint_fn(fitted_model$B[1, reference_set]), 0, 0.0001) == TRUE) {
          new_B <- fitted_model$B
        } else {
          new_B <- fitted_model$B
          p <- nrow(fitted_model$B)
          for (k in 1:p) {
            new_B[k, ] <- new_B[k, ] - constraint_fn(new_B[k, reference_set])
          }
        }

        # remove arguments from ... that we don't want for estimation for reference set
        est_args_rm <- which(names(extra_args) %in%
                               c("run_score_tests", "use_both_cov", "use_fullmodel_info",
                                 "use_fullmodel_cov", "return_both_score_pvals", "match_row_names",
                                 "constraint_fn", "constraint_grad_fn", "constraint_param"))
        est_args <- list(Y = Y,
                         X = X,
                         cluster = cluster,
                         penalize = penalize,
                         B = new_B,
                         fitted_model = NULL,
                         test_kj = test_kj,
                         refit = FALSE,
                         verbose = ifelse(verbose == "development", "development", FALSE),
                         run_score_tests = FALSE,
                         constraint_fn = constraint_fn_est,
                         constraint_grad_fn = constraint_grad_fn_est,
                         constraint_param = NA,
                         use_both_cov = FALSE,
                         use_fullmodel_info = FALSE,
                         use_fullmodel_cov = FALSE,
                         return_both_score_pvals = FALSE,
                         match_row_names = FALSE)
        if (length(est_args_rm) > 0) {
          modified_extra_args <- extra_args[-est_args_rm]
        } else {
          modified_extra_args <- extra_args
        }
        est_args <- c(est_args, modified_extra_args)
        fitted_model <- do.call(emuFit, est_args)
      }
    }
  }
  result <- fitted_model
  result$reference_set <- reference_set
  result$reference_set_names <- reference_set_names
  if (!(is.null(constraint_diff))) {
    result$constraint_diff <- constraint_diff
  }


  if (run_score_tests) {
    if (is.null(test_kj)) {
      test_kj <- data.frame(expand.grid(k = 2:p, j = 1:J))
    }

    n_test <- nrow(test_kj)
    included_categories <- vector(mode = "list", length = n_test)
    score_test_hyperparams <- vector(mode = "list", length = n_test)
    null_B <- vector(mode = "list", length = n_test)
    score_pieces <- vector(mode = "list", length = n_test)

    # check if list of starting values for B under the null has been provided
    if ("B_null_list" %in% names(extra_args)) {
      B_null_list <- extra_args$B_null_list
      if (length(B_null_list) != nrow(test_kj)) {
        warning("Length of 'B_null_list' is different than the number of tests specified in 'test_kj'. Ignoring object 'B_null_list'.")
        B_null_list <- NULL
      }
    } else {
      B_null_list <- NULL
    }

    # set constraint function and gradient over reference set
    if (length(reference_set) == 1) {
      constraint_fn_inf <- (function(x) x[1])
      constraint_grad_fn_inf <- function(x) {
        grad <- c(1, rep(0, length(x) - 1))
      }
    } else {
      constraint_fn_inf <- (function(x) {
        constraint_fn(x[1:length(reference_set)])
      })
      constraint_grad_fn_inf <- (function(x) {
        grad <- rep(0, length(x))
        grad[1:length(reference_set)] <-
          constraint_grad_fn(x[1:length(reference_set)])
        return(grad)
      })
    }

    # run score tests
    for (i_test in 1:n_test) {
      if (verbose %in% c(TRUE, "development")) {
        message(paste0("Running score test ", i_test, " of ", n_test, "."))
        start <- proc.time()
      }

      # reorder categories in Y so that reference set categories are in the beginning
      ind_keep <- c(reference_set, test_kj$j[i_test])
      ind_keep <- unique(ind_keep)
      if (length(ind_keep) == 1) {
        stop("If you only have one category in the reference set, this category cannot also be tested.")
      }
      new_order <- c(ind_keep, (1:ncol(Y))[-ind_keep])
      # add in additional categories if needed (so that each observation has at least one
      # category with non-zero counts)
      added_cats <- NULL
      mod_Y <- Y[, ind_keep]
      if (sum(rowSums(mod_Y) == 0) > 0) {
        if (verbose == "development") {
          message("There is at least one observation with no counts across the categories included in this model. Additional categories will be added to the model")
        }
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

      j_ind <- which(new_order == test_kj$j[i_test])

      # use estimated B to start
      upd_result <- result
      upd_result$B <- upd_result$B[, ind_keep]
      if (penalize) {
        upd_result$Y_augmented <- upd_result$Y_augmented[, ind_keep]
      }

      # remove arguments we don't want for score tests
      inf_args_rm <- which(names(extra_args) %in%
                             c("constraint_fn", "constraint_grad_fn", "constraint_param",
                               "B_null_list", "use_fullmodel_info", "use_fullmodel_cov",
                               "use_both_cov", "return_both_score_pvals", "match_row_names"))
      inf_args <- list(Y = mod_Y,
                       X = X,
                       formula = formula,
                       data = data,
                       cluster = cluster,
                       penalize = penalize,
                       fitted_model = upd_result,
                       refit = FALSE,
                       test_kj = data.frame(k = test_kj$k[i_test],
                                            j = j_ind),
                       return_wald_p = FALSE,
                       compute_cis = FALSE,
                       constraint_fn = constraint_fn_inf,
                       constraint_grad_fn = constraint_grad_fn_inf,
                       constraint_param = NA,
                       use_both_cov = FALSE,
                       use_fullmodel_info = FALSE,
                       use_fullmodel_cov = FALSE,
                       return_both_score_pvals = FALSE,
                       match_row_names = FALSE,
                       verbose = ifelse(verbose == "development", "development", FALSE))
      if (!(is.null(B_null_list))) {
        inf_args$B_null_list <- list(B_null_list[[i_test]])
      }
      if (length(inf_args_rm) > 0) {
        modified_extra_args <- extra_args[-inf_args_rm]
      } else {
        modified_extra_args <- extra_args
      }
      inf_args <- c(inf_args, modified_extra_args)
      # run score tests
      emuObj <- do.call(emuFit, inf_args)
      cols_rm <- which(names(emuObj$coef) %in% c("category_num", "estimate", "se", "lower", "upper", "zero_comparison"))
      row_data <- emuObj$coef[j_ind, -cols_rm]
      ind <- which(result$coef$covariate == row_data$covariate &
                       result$coef$category == row_data$category)
      result$coef[ind, names(row_data)] <- row_data
      included_categories[[i_test]] <- emuObj$coef$category
      score_test_hyperparams[[i_test]] <- emuObj$score_test_hyperparams
      if ("null_B" %in% names(emuObj)) {
        null_B[[i_test]] <- emuObj$null_B[[1]]
      }
      if ("score_components" %in% names(emuObj)) {
        score_pieces[[i_test]] <- emuObj$score_components[[1]]
      }
      if (verbose %in% c(TRUE, "development")) {
        end <- proc.time() - start
        sec <- round(end[3])
        if (sec <= 300) {
          time <- paste0(sec, " seconds")
        } else if (sec <= 18000) {
          min <- round(sec / 60)
          time <- paste0(min, " minutes")
        } else {
          hour <- round(sec / (60^2))
          time <- paste0(hour, " hours")
        }
        message(paste("Score test ", i_test, " of ", n_test,
                      " has completed in approximately ", time, ".",sep = ""))

      }
    }
    # add additional information from running score tests
    result$included_categories <- included_categories
    result$score_test_hyperparams <- score_test_hyperparams
    if (!is.null(null_B[[1]])) {
      result$null_B <- null_B
    }
    if (!is.null(score_pieces[[1]])) {
      result$score_components <- score_pieces
    }
  }

  # add column to coef data frame for reference set
  result$coef$reference_set_category <- result$coef$category %in% reference_set_names

  # return object
  return(result)
}
