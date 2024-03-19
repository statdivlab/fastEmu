# function to simulate data under the null
# adapted from `get_sim_bs` in `radEmu`
simulateBs <- function(J, constraint_cats, test_j,
                       constraint_mag, other_mag,
                       under_null = TRUE, alt_val = NULL) {
  # set b0 as in radEmu
  evens <- ((1:J)%%2 == 0)
  b0 <- numeric(J)
  b0[!evens] <- seq(-3, 3, length.out = sum(!evens))
  b0[evens] <- seq(3, -3, length.out = sum(evens))
  # set b1 separately for constraint categories, test category, other categories
  b1 <- rep(0, J)
  b1[constraint_cats] <- constraint_mag *
    sinh(seq(-10, 10, length.out = length(constraint_cats)))/sinh(10)
  b1[-c(constraint_cats, test_j)] <- other_mag *
    sinh(seq(-10, 10, length.out = (J - 1 - length(constraint_cats))))/sinh(10)
  if (!under_null) {
    if (is.null(alt_val)) {
      stop("If simulating under alternative, provide value for alternative B1")
    } else {
      b1[test_j] <- alt_val
    }
  }
  return(list(b0 = b0, b1 = b1))
}
