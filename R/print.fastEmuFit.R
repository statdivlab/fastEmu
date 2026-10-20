#' Print function
#'
#' @param x Object of class \code{fastEmuFit}
#' @param n The number of coefficient estimates to be printed (ordered by largest absolute value to smallest)
#' @param ... No optional arguments are accepted at this time.
#'
#' @return \code{NULL}. Displays printed model summary.
#'
#' @method print fastEmuFit
#'
#' @export
print.fastEmuFit <- function(x,
                             n = 20,
                             ...) {
  radEmu:::print.emuFit(x, n = n, ...)
}
