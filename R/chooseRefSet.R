#' A function to choose a data-driven reference set.
#'
#' @param fitted_model The output from \code{radEmu::emuFit} when applied to data with a
#' constraint over all taxa.
#' @param k The index of the column in the design matrix for which to choose the reference set.
#' @param reference_set_size The size of the reference set if it is data-driven, default
#' is set to \code{50}. We recommend a reference set of size 30-100 for the best balance
#' of computational efficiency and estimation precision.
#' @param constraint_fn The constraint function (by default the smoothed median).
#'
#' @importFrom utils head
#'
#' @return A list including the set of taxa of size \code{reference_set_size} with the smallest
#' absolute L2 norm of estimated log fold-differences for the specified covariate(s), relative to
#' the smoothed median log fold-difference over all taxa. Also a vector of length p (the number of
#' columns in the design matrix) of differences between the constraint function over log fold-differences
#' across all taxa and the constraint function log fold-difference over the reference set. Also
#' the estimated B matrix, with each row shifted according to the constraint over the new
#' reference set.
#'
chooseRefSet <- function(fitted_model,
                         k,
                         reference_set_size,
                         constraint_fn) {

  B_df <- data.frame(index = 1:ncol(fitted_model$B),
                    l2_lfd = sqrt(fitted_model$B[k, ]^2))
  B_ord <- order(B_df$l2_lfd)
  B_df <- B_df[B_ord, ]
  if (reference_set_size > nrow(B_df)) {
    reference_set_size <- nrow(B_df)
    message("Reference set size is larger than the number of taxa. All taxa will be used as the reference set, but consider decreasing reference set size or using `radEmu` instead of `fastEmu`.")
  }
  ref_set <- head(B_df$index, reference_set_size)
  p <- nrow(fitted_model$B)
  new_B <- fitted_model$B
  constraint_diff <- constraint_fn(fitted_model$B[k, ref_set])
  new_B[k, ] <- new_B[k, ] - constraint_diff[k]

  return(list(reference_set = ref_set, constraint_diff = constraint_diff, new_B = new_B))
}
