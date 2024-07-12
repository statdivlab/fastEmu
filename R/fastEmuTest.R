#' Run robust score tests using a reduced version of the radEmu model
#'
#' @param constraint_cats a vector of category indices for categories that should be involved in the constraint
#' and retained in any smaller joint model
#' @param estimate_full_model If TRUE, this function will estimate parameters from full model (recommended)
#' and run score tests with the reduced model. If FALSE, score tests will be run with the reduced model but
#' parameter estimates will not be produced. Set to FALSE if you are running score tests in parallel
#' and don't want to estimate parameters for the full model every time. Default is TRUE.
#' @param Y an n x J matrix or dataframe of nonnegative observations, or a phyloseq object containing an otu table and sample data.
#' @param X an n x p matrix or dataframe of covariates (optional, either include \code{X} or \code{formula} and \code{data})
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
#' @param fitted_model a fitted model produced by a call to radEmu::emuFit; to
#' be provided if score tests are to be run without refitting the full unrestricted model.
#' Default is NULL.
#' @param refit logical: if B or fitted_model is provided, should full model be fit (TRUE) or
#' should fitting step be skipped (FALSE), e.g., if score tests are to be run on an already
#' fitted model. Default is TRUE.
#' @param return_wald_p logical: return p-values from Wald tests? Default is FALSE. These can only be
#' returned if \code{estimate_full_model} is TRUE.
#' @param compute_cis logical: compute and return Wald CIs? Default is TRUE. These can only be
#' returned if \code{estimate_full_model} is TRUE.
#' @param verbose provide updates as model is being fitted? Defaults to TRUE.
#' @param ... Additional arguments to radEmu:::emuFit. See possible arguments with \code{?radEmu::emuFit}.
#' @param model deprecated argument, default is "drop" to run the reduced model that drops
#' all categories not in the constraint or being tested for each score test, "full" to run tests
#' with the full joint model (equivalent to radEmu), or "agg" to run a different version of the reduced
#' model in which all categories not in the constraint or being tested are aggregated together
#'
#' @return A list containing elements 'coef', 'included_categories', and 'score_test_hyperparams'.
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
fastEmuTest <- function(constraint_cats,
                        estimate_full_model = TRUE,
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
                        ...,
                        model = NULL) {


  run_score <- TRUE
  extra_args <- list(...)
  if ("run_score_tests" %in% names(extra_args)) {
    if (!extra_args$run_score_tests) {
      run_score <- FALSE
    }
  }
  if (run_score == FALSE & estimate_full_model == FALSE) {
    stop("You are not running estimation or score tests. Please choose at least one of these
         functionalities to use fastEmu. Do this by setting 'estimate_full_model' to TRUE or
         'run_score_tests' to TRUE.")
  }

  if ("constraint_fn" %in% names(extra_args) | "constraint_grad_fn" %in% names(extra_args) |
      "constraint_param" %in% names(extra_args)) {
    warning("You've included one of the parameters 'constraint_fn', 'constraint_grad_fn', or 'constraint_param'. These will be ignored.
            By default, fastEmu uses the pseudo-median over the categories specified in 'constraint_cats' as a constraint.")
  }

  # check if Y is a phyloseq object
  if ("phyloseq" %in% class(Y)) {
    if (requireNamespace("phyloseq", quietly = TRUE)) {
      if (is.null(formula)) {
        stop("If Y is a `phyloseq` object, make sure to include the formula argument.")
      } else {
        data <- data.frame(phyloseq::sample_data(Y))
        X <- model.matrix(formula, data)
        taxa_are_rows <- Y@otu_table@taxa_are_rows
        Y <- as.matrix(phyloseq::otu_table(Y))
        if (taxa_are_rows) {
          Y <- t(Y)
        }
      }
    } else {
      stop("You are trying to use a `phyloseq` data object or `phyloseq` helper function without having the `phyloseq` package installed. Please either install the package or use a standard data frame.")
    }
  } else if ("data.frame" %in% class(Y)) {
    Y <- as.matrix(Y)
    if (!is.numeric(Y)) {
      stop("Y is a data frame that cannot be coerced to a numeric matrix. Please fix and try again.")
    }
  }

  if (sum(rowSums(Y) == 0) > 0) {
    stop("There is at least one sample with no counts in any category. Please remove samples that have no counts in any category.")
  }

  # add category names if they aren't already there
  if (is.null(colnames(Y))) {
    colnames(Y) <- paste0("category_", 1:ncol(Y))
  }

  # estimate using full radEmu model
  if (estimate_full_model) {
    if (verbose) {
      message("Estimating parameters from full model.")
    }

    # set constraints for estimation
    if (length(constraint_cats) == 1) {
      constraint_fn_est <- (function(x) x[constraint_cats])
      constraint_grad_fn_est <- function(x) {
        grad <- rep(0, length(x))
        grad[constraint_cats] <- 1
      }
    # set constraint as pseudo-Huber over constraint categories if there are multiple
    } else {
      constraint_fn_est <- (function(x) {
        radEmu:::pseudohuber_center(x[constraint_cats], d = .1)
      })
      constraint_grad_fn_est <- (function(x) {
        grad <- rep(0, length(x))
        grad[constraint_cats] <-
          radEmu:::dpseudohuber_center_dx(x[constraint_cats], d = .1)
        return(grad)
      })
    }

    # remove arguments from ... that we don't want for estimation
    est_args_rm <- which(names(extra_args) %in%
                           c("constraint_fn", "constraint_grad_fn", "constraint_param",
                             "run_score_tests", "use_both_cov", "use_fullmodel_info",
                             "use_fullmodel_cov", "return_both_score_pvals"))
    est_args <- list(Y = Y,
                     X = X,
                     formula = formula,
                     data = data,
                     cluster = cluster,
                     penalize = penalize,
                     B = B,
                     fitted_model = fitted_model,
                     refit = refit,
                     return_wald_p = return_wald_p,
                     compute_cis = compute_cis,
                     constraint_fn = constraint_fn_est,
                     constraint_grad_fn = constraint_grad_fn_est,
                     constraint_param = NA,
                     verbose = verbose,
                     run_score_tests = FALSE,
                     use_both_cov = FALSE,
                     use_fullmodel_info = FALSE,
                     use_fullmodel_cov = FALSE,
                     return_both_score_pvals = FALSE)
    if (length(est_args_rm) > 0) {
      modified_extra_args <- extra_args[-est_args_rm]
    } else {
      modified_extra_args <- extra_args
    }
    est_args <- c(est_args, modified_extra_args)
    emu_est <- do.call(emuFit, est_args)
  }

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
                       fitted_model = NULL,
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
    res$score_rest_hyperparams <- score_test_hyperparams
    if (!is.null(null_B[[1]])) {
      res$null_B <- null_B
    }
  }

  return(res)
}
