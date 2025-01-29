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
#' @return Sarah UPDATE!!! A list containing elements 'coef', 'included_categories', and 'score_test_hyperparams'.
#' 'coef' is a matrix with score test statistics and p-values for all parameters tested. 'Included_categories'
#' is a list of the categories used in the reduced model for each score test that is fit.
#' 'score_test_hyperparameters' contains parameters and hyperparameters related to estimation under the
#' null, including whether or not the algorithm converged, which can be helpful for debugging.
#' The elements in 'Included_categories' and 'score_test_hyperparams' correspond to the tests described
#' in \code{test_kj} and are in the same order.
#' If \code{esttimate_full_model} is TRUE, then 'coef' will also include the estimates and optionally
#' confidence intervals and Wald test p-values for all parameters, estimated from the full model.
#' Additionally, if \code{estimate_full_model} is TRUE, the return list will also include 'B', 'penalized', 'Y_augmented',
#' 'z_hat', 'I', and 'Dy'. See the documentation in \code{?radEmu::emuFit} for a description of these
#' elements.
#'
#' @importFrom stats model.matrix
#' @import radEmu
#'
#' @export
#'
fastEmuFit <- function(reference_set = "data_driven",
                       reference_set_size = 50,
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

  # Record call
  call <- match.call(expand.dots = FALSE)

  # start by running radEmu checks
  extra_args <- list(...)
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
                                         unobserved_taxon_error = unobserved_taxon_error)
  Y <- check_results$Y
  X <- check_results$X
  cluster <- check_results$cluster
  B_null_list <- check_results$B_null_list
  n <- nrow(X)
  J <- ncol(Y)
  p <- ncol(X)

  # check for valid reference set
  if (reference_set %in% c("data_driven_thin", "data_driven_ss")) {
    stop("Sarah hasn't implemented data-driven reference sets with sample splitting or thinning yet. If you need this, please open a Github issue.")
  } else if (reference_set == "data_driven") {
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
  if (!("constraint_fn" %in% names(extra_args))) {
    if (!("constraint_param" %in% names(extra_args))) {
      constraint_param <- 0.1
    } else {
      constraint_param <- extra_args$constraint_param
    }
    constraint_fn <- (function(x) radEmu:::pseudohuber_center(x, d = constraint_param))
  } else {
    constraint_fn <- extra_args$constraint_fn
  }
  if (reference_set == "data-driven" & is.numeric(constraint_fn)) {
    stop("You've set a data-driven reference set, but chosen a constraint function that consists of a single reference category. Please rerun this function, and set `reference_set` to be the index of the taxa that you would like to the be the reference category, or remove the `constraint_fn` argument.")
  } else if (is.numeric(constraint_fn) & is.numeric(reference_set)) {
    if (reference_set != constraint_fn) {
      stop("You've set a single reference taxon with the `constraint_fn` argument but this doesn't match the argument `reference_set`. Please remove one of these arguments.")
    }
  }

  # right now can't do data-driven reference set for p > 2, add this functionality in later!
  if (reference_set == "data_driven" & p > 2) {
    stop("Sarah hasn't implemented data-driven reference sets for multiple covariates or covariate levels yet. Sarah please implement this!!!")
  }



  # if reference set is data_driven, determine this reference set
  if (reference_set == "data_driven") {
    if (verbose %in% c(TRUE, "development")) {
      message("Running estimation algorithm to determine reference set.")
    }
    if (!(is.null(fitted_model)) & "reference_set" %in% names(fitted_model)) {
      message("You have set `reference_set = data_driven` but provided a `fitted_model` that already includes a `reference_set`. This reference set will be used. If you would like the reference set to be recalculated, set `fitted_model$reference_set <- NULL` and rerun this function.")
      reference_set <- fitted_model$reference_set
    } else {
      if (is.null(fitted_model) & is.null(B)) {
        # estimate with radEmu::emuFit

        # remove arguments from ... that we don't want for estimation for reference set
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
                         return_wald_p = FALSE,
                         compute_cis = FALSE,
                         verbose = (verbose == "development"),
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
        fitted_model <- do.call(emuFit, est_args)
      }

      # choose reference set
      if (is.null(fitted_model)) {
        fitted_model$B <- B
      }
      ref_set_res <- chooseRefSet(fitted_model = fitted_model,
                                  reference_set_size = reference_set_size,
                                  constraint_fn = constraint_fn)
      reference_set <- ref_set_res$reference_set
      reference_set_names <- colnames(Y)[reference_set]
      constraint_diff <- ref_set_res$constraint_diff
      new_B <- ref_set_res$new_B

      # rerun emuFit (using new B) to update coef data frame
      re_est_args <- est_args
      re_est_args$fitted_model <- NULL
      re_est_args$B <- new_B
      re_est_args$refit <- FALSE
      re_est_args$compute_cis <- compute_cis
      re_est_args$return_wald_p <- return_wald_p
      fitted_model <- do.call(emuFit, re_est_args)
    }
  } else {
    # check if fitted model (with reference set) exists, otherwise fit model
    if (is.null(fitted_model) & is.null(B)) {
      if (verbose %in% c(TRUE, "development")) {
        message("Running estimation algorithm.")
      }
      # remove arguments from ... that we don't want for estimation for reference set
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
                       verbose = (verbose == "development"),
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
      fitted_model <- do.call(emuFit, est_args)
    } else {
      if (is.null(fitted_model)) {
        fitted_model$B <- B
      }
      if (!("reference_set" %in% names(fitted_model))) {
        # check whether constraint has already been set
        if (all.equal(constraint_fn(fitted_model$B[1, reference_set]), 0, 0.0001)) {
          new_B <- fitted_model$B
        } else {
          new_B <- fitted_model$B
          p <- nrow(fitted_model$B)
          for (k in 1:p) {
            new_B[k, ] <- new_B[k, ] - constraint_fn(new_B[k, ])
          }
        }

        # remove arguments from ... that we don't want for estimation for reference set
        est_args_rm <- which(names(extra_args) %in%
                               c("run_score_tests", "use_both_cov", "use_fullmodel_info",
                                 "use_fullmodel_cov", "return_both_score_pvals", "match_row_names"))
        est_args <- list(Y = Y,
                         X = X,
                         cluster = cluster,
                         penalize = penalize,
                         B = new_B,
                         fitted_model = NULL,
                         refit = FALSE,
                         verbose = (verbose == "development"),
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

  # run score tests
  if ("run_score_tests" %in% names(extra_args)) {
    run_score_tests <- extra_args$run_score_tests
  } else {
    run_score_tests <- TRUE
  }

  if (run_score_tests) {
    if (is.null(test_kj)) {
      test_kj <- data.frame(expand.grid(k = 2:p, j = 1:J))
    }

    n_test <- nrow(test_kj)
    score_res <- NULL
    included_categories <- vector(mode = "list", length = n_test)
    score_test_hyperparams <- vector(mode = "list", length = n_test)
    null_B <- vector(mode = "list", length = n_test)
    score_pieces <- vector(mode = "list", length = n_test)
  }

  # next - figure out how to update constraint if not pseudo-Huber

  # set constraint_fn to index if reference set size is 1

  # -------------- OLD ----------------------

  if (run_score) {
    # prepare to run score tests

    # get X matrix
    if (is.null(X)) {
      if (is.null(formula) | is.null(data)) {
        stop("If design matrix X not provided, both formula and data containing
covariates in formula must be provided.")
      }
      X <- model.matrix(formula, data)
    }
    if ("data.frame" %in% class(X)) {
      X <- as.matrix(X)
      if (!is.numeric(X)) {
        stop("X is a data frame that cannot be coerced to a numeric matrix. Please fix and try again.")
      }
    }
    p <- ncol(X)

    if (is.null(test_kj)) {
      test_kj <- data.frame(expand.grid(k = 2:p, j = 1:J))
    }

    n_test <- nrow(test_kj)
    score_res <- NULL
    included_categories <- vector(mode = "list", length = n_test)
    score_test_hyperparams <- vector(mode = "list", length = n_test)
    null_B <- vector(mode = "list", length = n_test)
    score_pieces <- vector(mode = "list", length = n_test)

    if (is.null(model)) {
      model = "drop"
    }

    # set constraints for inference
    # for full model, don't do reordering
    if (model == "full") {
      if (estimate_full_model) {
        constraint_fn_inf <- constraint_fn_est
        constraint_grad_fn_inf <- constraint_grad_fn_est
      } else {
        if (length(constraint_cats) == 1) {
          constraint_fn_inf <- (function(x) x[constraint_cats])
          constraint_grad_fn_inf <- function(x) {
            grad <- rep(0, length(x))
            grad[constraint_cats] <- 1
          }
          # set constraint as pseudo-Huber over constraint categories if there are multiple
        } else {
          constraint_fn_inf <- (function(x) {
            radEmu:::pseudohuber_center(x[constraint_cats], d = .1)
          })
          constraint_grad_fn_inf <- (function(x) {
            grad <- rep(0, length(x))
            grad[constraint_cats] <-
              radEmu:::dpseudohuber_center_dx(x[constraint_cats], d = .1)
            return(grad)
          })
        }
      }
      # for reduced model, reorder so that constraint categories come first
    } else {
      # set constraint as single category constraint if there is one value in constraint_cats
      if (length(constraint_cats) == 1) {
        constraint_fn_inf <- (function(x) x[1])
        constraint_grad_fn_inf <- function(x) {
          grad <- c(1, rep(0, length(x) - 1))
        }
        # set constraint as pseudo-Huber over constraint categories if there are multiple
      } else {
        constraint_fn_inf <- (function(x) {
          radEmu:::pseudohuber_center(x[1:length(constraint_cats)], d = .1)
        })
        constraint_grad_fn_inf <- (function(x) {
          grad <- rep(0, length(x))
          grad[1:length(constraint_cats)] <-
            radEmu:::dpseudohuber_center_dx(x[1:length(constraint_cats)], d = .1)
          return(grad)
        })
      }
    }

    # run score tests
    if (verbose) {
      message("Running score tests")
    }
    for (i_test in 1:n_test) {
      if (verbose) {
        message(paste0("Running score test ", i_test))
      }

      # update model based on model chosen and categories to use
      if (model == "full") {
        mod_Y <- Y
        j_ind <- test_kj$j[i_test]
      } else {
        # reorder categories in Y so that constraint categories are in the beginning
        ind_keep <- c(constraint_cats, test_kj$j[i_test])
        ind_keep <- unique(ind_keep)
        if (length(ind_keep) == 1) {
          stop("If you only have one constraint category, this category cannot also be tested.")
        }
        new_order <- c(ind_keep, (1:ncol(Y))[-ind_keep])

        added_cats <- NULL

        if (model == "drop") {
          mod_Y <- Y[, ind_keep]
          # add in additional categories if needed
          if (sum(rowSums(mod_Y) == 0) > 0) {
            if (verbose) {
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
        } else if (model == "agg") {
          other_cat_sum <- rowSums(Y[, -ind_keep])
          mod_Y <- cbind(Y[, ind_keep], other_cat_sum)
        } else {
          stop("Please use 'full', 'drop', or 'agg' for method argument.")
        }
        j_ind <- which(new_order == test_kj$j[i_test])
      }

      # if we've estimated full model or have B or fitted model as inputs, use these as start
      upd_B <- NULL
      if (estimate_full_model) {
        upd_B <- emu_est$B
      } else if (!is.null(fitted_model)) {
        upd_B <- fitted_model$B
      } else if (!is.null(B)) {
        upd_B <- B
      }

      if (!is.null(upd_B)) {
        if (model != "full") {
          upd_B <- upd_B[, ind_keep]
        }
      }

      if (estimate_full_model & penalize) {
        upd_fitted_model <- emu_est
        upd_fitted_model$B <- upd_B
        if (model != "full") {
          upd_fitted_model$Y_augmented <- upd_fitted_model$Y_augmented[, ind_keep]
        }
      } else if (!is.null(fitted_model) & penalize) {
        upd_fitted_model <- fitted_model
        upd_fitted_model$B <- upd_B
        if (model != "full") {
          upd_fitted_model$Y_augmented <- upd_fitted_model$Y_augmented[, ind_keep]
        }
      } else {
        upd_fitted_model <- NULL
      }

      # remove arguments we don't want for score tests
      inf_args_rm <- which(names(extra_args) %in%
                             c("constraint_fn", "constraint_grad_fn", "constraint_param"))
      inf_args <- list(Y = mod_Y,
                       X = X,
                       formula = formula,
                       data = data,
                       cluster = cluster,
                       penalize = penalize,
                       B = upd_B,
                       fitted_model = upd_fitted_model,
                       refit = refit,
                       test_kj = data.frame(k = test_kj$k[i_test],
                                            j = j_ind),
                       return_wald_p = FALSE,
                       compute_cis = FALSE,
                       constraint_fn = constraint_fn_inf,
                       constraint_grad_fn = constraint_grad_fn_inf,
                       constraint_param = NA,
                       verbose = verbose)
      if (length(inf_args_rm) > 0) {
        modified_extra_args <- extra_args[-inf_args_rm]
      } else {
        modified_extra_args <- extra_args
      }
      inf_args <- c(inf_args, modified_extra_args)
      # run score tests
      emuObj <- do.call(emuFit, inf_args)
      cols_rm <- which(names(emuObj$coef) %in% c("category_num", "estimate", "se", "lower", "upper"))
      row_data <- emuObj$coef[j_ind, -cols_rm]
      score_res <- rbind(score_res, row_data)
      if (estimate_full_model) {
        ind <- which(emu_est$coef$covariate == row_data$covariate &
                       emu_est$coef$category == row_data$category)
        emu_est$coef[ind, names(row_data)] <- row_data
      }
      included_categories[[i_test]] <- emuObj$coef$category
      score_test_hyperparams[[i_test]] <- emuObj$score_test_hyperparams
      if ("null_B" %in% names(emuObj)) {
        null_B[[i_test]] <- emuObj$null_B[[1]]
      }
      if ("score_components" %in% names(emuObj)) {
        score_pieces[[i_test]] <- emuObj$score_components[[1]]
      }
    }
  }

  if (estimate_full_model) {

    # add column in coef about constraint categories
    emu_est$coef$constraint_cat <- emu_est$coef$category %in% colnames(Y)[constraint_cats]

    res <- emu_est
  } else {
    # reorder coef
    cov_order <- unique(score_res$covariate)
    cat_order <- colnames(Y)
    reordered <- order(sapply(score_res$covariate, function(x) {which(cov_order == x)}),
                       sapply(score_res$category, function(x) {which(cat_order == x)}))
    score_res <- score_res[reordered, ]

    # add column in coef about constraint categories
    score_res$constraint_cat <- score_res$category %in% colnames(Y)[constraint_cats]

    res <- list()
    res$coef <- score_res
  }
  if (run_score) {
    res$included_categories <- included_categories
    res$score_test_hyperparams <- score_test_hyperparams
    if (!is.null(null_B[[1]])) {
      res$null_B <- null_B
    }
    if (!is.null(score_pieces[[1]])) {
      res$score_components <- score_pieces
    }
  }

  return(res)
}
