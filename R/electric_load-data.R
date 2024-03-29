#' Electricity forecasting data set
#'
#'
#' Electricity forecasting data set provided by EDF R&D.
#' It contains weekly measurements of the total electricity consumption in France from 1996 to 2009,
#' together with several covariates, including temperature, industrial production indices (source: INSEE) and calendar information.
#'
#' @docType data
#'
#' @usage data(electric_load)
#'
#' @keywords datasets
#'
#'
#'
#' @examples
#'  data(electric_load)
#'  # a few graphs to display the data
#'  attach(electric_load)
#'  plot(Load, type = 'l')
#'  plot(Temp, Load, pch = 16, cex = 0.5)
#'  plot(NumWeek, Load, pch = 16, cex = 0.5)
#'  plot(Load, Load1, pch = 16, cex = 0.5)
#'  acf(Load, lag.max = 20)
#'  detach(electric_load)
"electric_load" 
