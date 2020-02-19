#' Calibration errors
#'
#' Computes the g-weights of the calibration estimator.
#'
#' @param Xs_old matrix of calibration variables.
#' @param d vector of initial weights.
#' @param total_old vector of population totals.
#' @param q vector of positive values accounting for heteroscedasticity; the variation of the g-weights is reduced for small values of q.
#' @param method calibration method (linear, raking, logit, truncated).
#' @param bounds vector of bounds for the g-weights used in the truncated and logit methods; 'low' is the smallest value and 'upp' is the largest value.
#' @param description if description=TRUE, summary of initial and final weights are printed, and their boxplots and histograms are drawn; by default, its value is FALSE.
#' @param max_iter maximum number of iterations in the Newton's method.
#' @param EPS1 tolerance in the algorithm.(default is eps1=1e-06)
#' @param cut_off A numeric value for the pair-wise absolute correlation cutoff
#' @param verbose A boolean for printing the details of function findCorrelations.
#' @param relative_error if relative_error=TRUE, checks quality of the relative error. If relative_error = FALSE, checks quality of the absolute error.
#' @param constant if constant=TRUE, its remove constant variables in the matrix of calibration.
#' @param digit_level integer which shows how many digits rounds the results in the result table.
#' @param print_indvar if print_indvar=TRUE, Independent variables are printed; by default, its value is FALSE.
#'
#' @return Calculate maximum absolute error and maximum relative error.
#'
#' @seealso \code{\link{Calibrationerrors}}
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
# g2 <- calib2(Xs, d = 1/pik, total, method = "raking")
# Calibrationerrors(Xs, d=1/pik, total, g)
#'
#' }
#'
#' @export calib2
#'
#' @import MASS
#' @import sampling
#' @import quantreg
#' @import matrixcalc
#' @import data.table
#' @importFrom caret  findCorrelation
#' @importFrom caret  findLinearCombos
#' @importFrom graphics boxplot
#' @importFrom graphics hist
#' @importFrom graphics par
#' @importFrom stats cor
#'
calib2 <-  function(Xs_old, d, total_old, q = rep(1, length(d)),
                    method = c("linear", "raking", "truncated", "logit"),
                    bounds = c(low = 0, upp = 10), description = FALSE,
                    max_iter = 500, EPS1 = 1e-06, cut_off = 1,
                    constant = TRUE, verbose = FALSE, relative_error = FALSE,
                    digit_level = 6, print_indvar = FALSE) {


  if (any(is.na(Xs_old)) | any(is.na(d)) | any(is.na(total_old)) |
      any(is.na(q)))
    stop("the input should not contain NAs")

  if (!(is.matrix(Xs_old) | is.array(Xs_old)))
    Xs_old <- as.matrix(Xs_old)
  if (is.matrix(Xs_old))
    if (length(total_old) != ncol(Xs_old))
      stop("Xs_old and total_old have different dimensions")

  if (is.vector(Xs_old) & length(total_old) != 1)
    stop("Xs_old and total_old have different dimensions")
  if (any(is.infinite(q)))
    stop("there are Inf values in the q vector")
  if (missing(method))
    stop("Specify a method")
  if (!(method %in% c("linear", "raking", "logit", "truncated")))
    stop("the specified method is not in the list")
  if (method %in% c("linear", "raking") & !missing(bounds))
    stop("for the linear and raking the bounds are not allowed")

  EPS = .Machine$double.eps
  n = length(d)

  total_old2 = total_old
  Xs_old2 = Xs_old

  var2 = nrow(Xs_old) / (nrow(Xs_old) - 1) * (colMeans(Xs_old * Xs_old) - colMeans(Xs_old)^2)

  if (!constant) {
    Xs_old2 = as.matrix(Xs_old[,var2 != 0])
    total_old2 = total_old[var2 != 0]
  } else {total_old2 = total_old}

  v3 = findLinearCombos(Xs_old2)$remove
  Xs_old3 = as.matrix(Xs_old2)
  total_old3 = total_old2

  if (!is.null(v3)) {
    Xs_old3 = as.matrix(Xs_old2[,-v3])
    total_old3 = total_old2[-v3]
  }

  Xs = as.matrix(Xs_old3)
  total = total_old3
  v4 = integer(0)

  if (ncol(Xs) > 1) {
    var2 = nrow(Xs) / (nrow(Xs) - 1) * (colMeans(Xs * Xs) - colMeans(Xs)^2)
    Xs3cor = cor(Xs[, var2 != 0])
    Xs3cor[Xs3cor > 1] = 1
    Xs3cor[Xs3cor < -1] = -1
    v4 = findCorrelation(Xs3cor, cutoff = cut_off, verbose = verbose )
  }

  if (1 - !(length(v4))) {
    Xs = as.matrix(Xs_old3[, -v4])
    total = total_old3[-v4] }

  if (description) {
    cat("Matrix Xs column count: ")
    print(ncol(Xs))
    cat("Input matrix Xs column count: ")
    print(ncol(Xs_old))
    cat("\n")
  }

  if (print_indvar) {
    cat("Independent variables:\n")
    print(colnames(Xs))
    cat("\n")}

  lambda = as.matrix(rep(0, n))
  lambda1 = ginv(t(Xs * d * q) %*% Xs, tol = EPS) %*% (total -
                                                         as.vector(t(d) %*% Xs))
  if (relative_error | method == "linear" | max(abs(lambda1)) < EPS) {
    g <- calib(Xs = Xs, d = d, total = total,
               q = q, method = method,
               bounds = bounds, description = FALSE,
               max_iter = max_iter)
  } else if (!relative_error & method == "truncated") {
    if (!missing(bounds)) {
      if (bounds[2] <= 1 || bounds[1] >= 1 || bounds[1] >
          bounds[2])
        stop("The conditions low<1<upp are not fulfilled")
    }
    else stop("Specify the bounds")
    g = 1 + q * as.vector(Xs %*% lambda1)
    list = 1:length(g)
    l = 0
    g1 = g
    Xs1 = Xs
    d1 = d
    t2 = total
    list1 = 1:length(g)
    q1 = q
    while (l < max_iter) {
      l = l + 1
      if (any(g < bounds[1]) | any(g > bounds[2])) {
        g[g < bounds[1]] = bounds[1]
        g[g > bounds[2]] = bounds[2]
        list = (1:length(g))[g > bounds[1] & g < bounds[2]]
        if (length(list) != 0) {
          g1 = g[list]
          t2 = total - as.vector(t(g[-list] * d[-list]) %*%
                                   Xs[-list, ])
          Xs1 = Xs[list, ]
          d1 = d[list]
          q1 = q[list]
          list1 = list
        }
      }
      t1 = as.vector(t(d1) %*% Xs1)
      lambda1 = ginv(t(Xs1 * d1 * q1) %*% Xs1, tol = EPS) %*%
        (t2 - t1)
      if (length(list1) > 1)
        g1 = 1 + q1 * as.vector(Xs1 %*% lambda1)
      else if (length(list1) == 1) {
        g1 = 1 + q1 * as.vector(as.vector(Xs1) %*% as.vector(lambda1))
      }
      g[list1] = g1
      tr = crossprod(Xs, g * d)
      if (max(abs(tr - total)) < EPS1 & all(g >=
                                            bounds[1] & g <= bounds[2]))
        break
    }
    if (l == max_iter) {
      cat("No convergence in", max_iter, "iterations with the given bounds. \n")
      cat("The bounds for the g-weights are:", min(g),
          " and ", max(g), "\n")
      cat(" and the g-weights are given by g\n")
    }
  } else if (!relative_error & method == "raking") {
    lambda = as.matrix(rep(0, ncol(Xs)))
    w1 = as.vector(d * exp(Xs %*% lambda * q))
    for (l in 1:max_iter) {
      phi = t(Xs) %*% w1 - total
      T1 = t(Xs * w1)
      phiprim = T1 %*% Xs
      lambda = lambda - ginv(phiprim, tol = EPS) %*% phi
      w1 = as.vector(d * exp(Xs %*% lambda * q))
      if (any(is.na(w1)) | any(is.infinite(w1))) {
        warning("No convergence")
        g = NULL
        break
      }
      tr = crossprod(Xs, w1)
      if (max(abs(tr - total)) < EPS1)
        break
    }
    if (l == max_iter) {
      warning("No convergence")
      g = NULL
    }
    else g = w1/d
  } else if (!relative_error & method == "logit") {
    if (bounds[2] <= 1 || bounds[1] >= 1 || bounds[1] > bounds[2])
      stop("The conditions low < 1 < upp are not fulfilled")
    A = (bounds[2] - bounds[1])/((1 - bounds[1]) * (bounds[2] -
                                                      1))
    u = rep(1, length(d))
    F = (bounds[1] * (bounds[2] - 1) + bounds[2] * (1 - bounds[1]) *
           u)/(bounds[2] - 1 + (1 - bounds[1]) * u)
    w1 = as.vector(d * F)
    T = t(Xs * w1)
    phiprim = ginv(T %*% Xs, tol = EPS)
    g = F
    tr = crossprod(Xs, w1)
    if (max(abs(tr - total)) < EPS1 | any(g < bounds[1]) | any(g > bounds[2])) {
      lambda1 = rep(0, ncol(Xs))
      list = 1:length(g)
      t2 = total
      Xs1 = Xs
      d1 = d
      g1 = g
      q1 = q
      list1 = 1:length(g)
      for (l in 1:max_iter) {
        if (any(g < bounds[1]) | any(g > bounds[2])) {
          g[g < bounds[1]] = bounds[1]
          g[g > bounds[2]] = bounds[2]
          list = (1:length(g))[g > bounds[1] & g < bounds[2]]
          if (length(list) != 0) {
            g1 = g[list]
            t2 = total - as.vector(t(g[-list] * d[-list]) %*%
                                     Xs[-list, ])
            Xs1 = Xs[list, ]
            d1 = d[list]
            q1 = q[list]
            list1 = list
          }
        }
        t1 = as.vector(t(d1) %*% Xs1)
        phi = t(Xs1) %*% as.vector(d1 * g1) - t1
        T = t(Xs1 * as.vector(d1 * g1))
        phiprime = T %*% Xs1
        lambda1 = lambda1 - ginv(phiprime, tol = EPS) %*%
          (as.vector(phi) - t2 + t1)
        u = exp(A * (Xs1 %*% lambda1 * q1))
        F = g1 = (bounds[1] * (bounds[2] - 1) + bounds[2] *
                    (1 - bounds[1]) * u)/(bounds[2] - 1 + (1 -
                                                             bounds[1]) * u)
        if (any(is.na(g1))) {
          warning("no convergence")
          g1 = g = NULL
          break
        }
        g[list1] = g1
        tr = crossprod(Xs, g * d)
        if (max(abs(tr - total)) < EPS1 & all(g >=
                                              bounds[1] & g <= bounds[2]))
          break
      }
      if (l == max_iter) {
        cat("no convergence in", max_iter, "iterations with the given bounds. \n")
        cat("the bounds for the g-weights are:", min(g),
            " and ", max(g), "\n")
        cat(" and the g-weights are given by g\n")
        g = NULL
      }
    }
  }

  result <- NULL
  if (!is.null(g)) {
    result <- data.table(Xs_d = crossprod(Xs_old, d),
                         Xs_d_g = crossprod(Xs_old, d * g),
                         real = total_old)
    result$change <- round(result$real - result$Xs_d_g, digit_level)
    result$diff <- round(result$real - result$Xs_d_g, digit_level)
    names(result)[1:2] <- list("Xs*d", "Xs*d*g")
  }
  if (description && !is.null(g)) {
    par(mfrow = c(3, 2), pty = "s")
    hist(g)
    boxplot(g, main = "Boxplot of g")
    hist(d)
    boxplot(d, main = "Boxplot of d")
    hist(g * d)
    boxplot(g * d, main = "Boxplot of w=g*d")

    if (method %in% c("truncated", "raking", "logit"))
                       cat("number of iterations ", l, "\n")

    cat("summary - initial weigths d\n")
    print(summary(d))
    cat("\n")
    cat("summary - calibrated weigths g \n")
    print(summary(as.vector(g)))
    cat("\n")
    cat("summary - final weigths w=g*d\n")
    print(summary(as.vector(g * d)))
    cat("\n")
    cat("\n")

    print(result)

  }

  drop_var <- toString(colnames(Xs_old)[!(colnames(Xs_old) %in% colnames(Xs))])

  calibrationerrors <- Calibrationerrors(Xs = Xs_old, d = d, total = total_old, g = g)

  return(list(g = g, drop_var = drop_var,
              calib_Xnames = colnames(Xs),
              result = result,
              calibrationerrors =  calibrationerrors))
}


