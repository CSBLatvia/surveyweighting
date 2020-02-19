#' Nonresponse correction for single stage surveys
#'
#' Caclutate nonresponse correction factor for single stage surveys
#'
#' @param in_system (\code{character}) - ISDAVS sistēma, kurā veikts apsekojums
#' @param sample (\code{character}) - sampling indicator;
#' @param frame (\code{character}) - frame indicator variable name.
#' @param response_ind (\code{character}) - response indicator variable name,  that indicates which respondents are weighted in the weight calculation process.
#' @param summary_w (\code{character}) summary weight indicator variable name, that indicates which respondents have a weight greater than zero. Weight will be used to calculate summaries.
#' @param design_weight teksts (\code{character}) - design weight variable name.
#' @param strata (\code{character}) - strata variable name;
#' @param dataset data.table;
#'
#' @return dataset - survey data with weight (koef) and summary weight (koef_nodot).
#'
#' @seealso \code{\link{Calibrationerrors}}
#'
#' @examples
#' \donttest{
#' library(data.table)
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
                                     response_ind = "ievad",
                                     summary_w = "kops_sv",
                                     design_weight = "diz_sv",
                                     strata = "strata",
                                     dataset){

  koef <- ievad <- ievad_sum <- koef_nodot <- kops_sv <- NULL

  skaits_stratas <- dataset[, lapply(.SD, sum, na.rm = TRUE), keyby = strata,
                            .SDcols = c(response_ind, summary_w, frame, in_system, sample)]

  setnames(skaits_stratas, c(response_ind, summary_w, frame, in_system, sample),
                    paste0(c(response_ind, summary_w, frame, in_system, sample), "_sum"))

  dataset <- merge(dataset, skaits_stratas, by = c(strata), all = TRUE)

  dataset[, koef := 0]
  dataset[get(response_ind) == 1, koef := get(response_ind) * get(design_weight) * get(paste0(sample, "_sum")) / get(paste0(response_ind, "_sum"))]
  dataset[, koef_nodot := get(summary_w) * koef]

  dataset[]
}



