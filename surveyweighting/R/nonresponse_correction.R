#' Nonresponse correction for single stage surveys
#'
#' @description Calcutate nonresponse correction factor for single stage surveys
#'
#' @param in_system (\code{character}) - indicator that tell us, is data in system.
#' @param sample (\code{character}) - sampling indicator.
#' @param frame (\code{character}) - frame indicator variable name.
#' @param big_units (\code{character}) - big unit indicator variable name.
#' @param response_ind (\code{character}) - response indicator variable name,  that indicates which respondents are weighted in the weight calculation process.
#' @param summary_w (\code{character}) summary weight indicator variable name, that indicates which respondents have a weight greater than zero. Weight will be used to calculate summaries.
#' @param design_weight (\code{character}) - design weight variable name.
#' @param pop_size (\code{character}) - population size in strata variable name.
#' @param strata (\code{character}) - strata variable name;
#' @param dataset data.table;
#'
#' @return dataset - survey data with weight (koef) and summary weight (koef_nodot).
#'
#' @seealso \code{\link{Calibrationerrors}}
#'
#' @examples
#' library(data.table)
#' \donttest{
#'
#' }
#'
#' @export nonresponse_correction
#'
#' @import data.table
#'

nonresponse_correction <- function(in_system = "in_isdavs",
                                     sample = "izlase",
                                     frame = "frame",
                                     big_units = "lielu",
                                     response_ind = "ievad",
                                     summary_w = "kops_sv",
                                     design_weight = "diz_sv",
                                     pop_size = NULL,
                                     strata = "strata",
                                     dataset){

  ..pop_size <- ..design_weight <- koef <- ievad <- ievad_sum <- koef_nodot <- NULL

  skaits_stratas <- dataset[, lapply(.SD, sum, na.rm = TRUE), keyby = c(strata),
                            .SDcols = c(response_ind, summary_w, frame,
                                        big_units, in_system, sample)]

  setnames(skaits_stratas, c(response_ind, summary_w, frame,
                             big_units, in_system, sample),
                    paste0(c(response_ind, summary_w, frame,
                             big_units, in_system, sample), "_sum"))

  dataset <- merge(dataset, skaits_stratas, by = c(strata), all = TRUE)

  if (!is.null(dataset$koef)) dataset[, koef := NULL]
  dataset[, koef := 0]
  if (is.null(pop_size)) {
    dataset[get(response_ind) == 1, koef := get(response_ind) * get(..design_weight) * get(paste0(sample, "_sum")) / get(paste0(response_ind, "_sum"))]
  } else {
    dataset[get(response_ind) == 1 & get(big_units) == 1, koef := 1]
    dataset[get(response_ind) == 1 & get(big_units) == 0, koef := get(response_ind) * (get(..pop_size) - get(paste0(big_units, "_sum"))) / (get(paste0(response_ind, "_sum")) - get(paste0(big_units, "_sum")))]
  }


  dataset[, koef_nodot := get(summary_w) * koef]

  dataset[]
}


