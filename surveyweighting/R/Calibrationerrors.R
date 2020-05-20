#' Calibration errors
#'
#' @description Calculate maximum absolute error and maximum relative error
#'
#' @param Xs matrix of calibration variables.
#' @param d vector of initial weights.
#' @param total vector of population totals.
#' @param g vector of g-weights.
#'
#' @return Calculate maximum absolute error and maximum relative error.
#'
#' @seealso \code{\link{calib2}}
#'
#' @examples
#' \donttest{
#' library(sampling)
#'Xs <- cbind(c(1,1,1,1,1,0,0,0,0,0),c(0,0,0,0,0,1,1,1,1,1),c(1,2,3,4,5,6,7,8,9,10))
# inclusion probabilities
# pik <- rep(0.2, times = 10)
# total <- c(24, 26, 280)
# the g-weights
# g <- calib(Xs, d = 1/pik, total, method = "raking")
# Calibrationerrors(Xs, d=1/pik, total, g)
#'
#' }
#'
#' @export Calibrationerrors
#'
#' @import matrixcalc
#' @import data.table



Calibrationerrors <- function(Xs, d, total, g) {
    if (is.null(g)) stop("the g-weight vector is null")
    if (!is.matrix(Xs)) Xs = as.matrix(Xs)
    tr <- crossprod(Xs, g * d)
    matrica <- data.table(Max_absolute_error = max(abs(tr - total)),
                          Max_relative_error = max(abs(tr - total) / total))
    return(matrica)
}
