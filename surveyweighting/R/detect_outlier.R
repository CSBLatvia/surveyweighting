#' Calculated Outliers for single stage surveys
#'
#' @description Detected outliers for single stage surveys.
#'
#' @param id (\code{character}) - unit id variable name.
#' @param vars (\code{character}) - variable name.
#' @param group (\code{character}) - group variable name.
#' @param b (\code{numeric}) - the user-dependent attitude to distances between succeeding values.
#' @param a (\code{numeric}) - the threshold.
#' @param nd_log (\code{numeric}) - correction for calculating average density of consecutive values
#' @param nd_log2 (\code{numeric}) - second correction for calculating average density of consecutive values
#' @param percthreshold (\code{numeric}) - percentage of the threshold;
#' @param dataset data.table;
#'
#' @return A list with two objects
#' \itemize{
#'    \item datao - data.table with outliers
#'    \item count - count of the outliers.
#' }
#'
#' @references
#' Mark Last, Abraham Kandelm, Automated Detection of Outliers in Real-World Data, \url{http://www.ise.bgu.ac.il/faculty/mlast/papers/outliers2.pdf}
#' @return dataset - survey data with outliers.
#'
#' @seealso \code{\link{Calibrationerrors}}
#'
#' @examples
#' \donttest{
#' library(data.table)
#' }
#'
#' @export detect_outlier
#'
#' @import data.table


detect_outlier <- function(id, vars, group, dataset, b = 1,
                           a = 0.1, nd_log = 2, nd_log2 = 4,
                           percthreshold = 0.95){

  N <- casen <-  d <- grps <- lag_1 <- lag_m1 <- NULL
  ll <- m <- mu <- nd <- nv <- nv_sum <- outl <- NULL
  outl_v <- powh <- uu <- vv <- z2 <- NULL

  dataset <- data.table(dataset)
  dataset[, grps := .GRP, by = c(vars, group)]
  data_w <- dataset[!is.na(get(vars)), c(id, group, "grps", vars), with = FALSE]

  # Calculate indicator vars (vars) values frequencies in group
  temp2_vars <- data_w[, .N, by = c(vars, group, "grps")]

  setkeyv(temp2_vars, c(group, vars))
  temp2_vars[, d := sum(N), by = group]
  temp2_vars[, nd := .N, by = group]

  temp2_vars[, m := round(max(1, log(nd/log(nd_log)))), by = group]


  #  "a" value
  # with "a"  you can adjust the number of outlier
  # if there are too many bounces, the value of "a" is reduced
  # If a larger number of bounces is allowed, the value of "a" can be increased

  temp2_vars[, casen := 1:.N, by = group]
  temp2_vars[, lag_1 := shift(.SD[, vars, with = FALSE], n = 1, fill = NA, type = "lag"), by = group]
  temp2_vars[, lag_m1 := shift(.SD[, vars, with = FALSE], n = m + 1, fill = NA, type = "lag"), by = group]

  temp2_vars[nd > 3, powh := ((b * m * (get(vars) - lag_1)) / (N * (lag_1 - lag_m1)))]
  temp2_vars[nd > 3, mu := 2/(1 + exp(powh))]

  temp2_vars[casen < m + 2, mu := NA]
  temp2_vars[, ll := round(max(1, log(nd/log(nd_log2)))), by = group]
  temp2_vars[, nv := N * get(vars)]
  temp2_vars[, nv_sum := sum(nv), by = group]

  # if group is only 2 or 3 values, then outliers is deteted by other algorithm
  # must calculated variable "z2"
  temp2_vars[nd < 4 & nd > 1, z2 := nv_sum * percthreshold]
  temp2_vars[nd >= 4 | nd <= 1, z2 := NA]

  temp2_vars[, outl := 0]
  temp2_vars[casen > nd - ll & mu < a, outl := 1]
  temp2_vars[z2 < get(vars), outl := 1]

  temp2_vars[outl == 1, uu := min(get(vars)), by = c(group)]
  temp2_vars[is.na(uu), uu := 0]
  temp2_vars[, vv := max(uu, na.rm = TRUE), by = c(group)]
  temp2_vars[get(vars) >= vv & vv > 0, outl_v := 1]
  #temp2[apgr_isd < vv | vv == 0, outl := 0]
  temp2_vars[is.na(outl_v), outl_v := 0]

  outliers <- table(temp2_vars$outl_v, useNA = "ifany")

  temp2_vars[, c("N", "d", "nd", "m", "lag_1", "lag_m1", "powh", "mu", "casen",
                 "ll", "nv", "nv_sum", "z2", "outl", "uu", "vv") := NULL]

  data_w <- merge(dataset, temp2_vars, by = c(vars, group, "grps"), all.x = TRUE)
  data_w <- data_w[, c(id, "outl_v"), with = FALSE]

  setnames(data_w, "outl_v", paste0("outl_", vars))
  return(list(datao = data_w[], count = outliers))
}
