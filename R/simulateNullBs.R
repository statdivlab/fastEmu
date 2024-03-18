# function to simulate data under the null
# adapted from `get_sim_bs` in `radEmu`
simulateNullBs <- function(J, constraint_cats, test_j,
                           constraint_mag, other_mag) {
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

  return(list(b0 = b0, b1 = b1))
}
