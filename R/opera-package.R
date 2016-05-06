#' Online Prediction by ExpeRt Aggregation
#' 
#' The package \code{opera} performs, for regression-oriented time-series,
#' predictions by combining a finite set of forecasts provided by the user.
#' More formally, it considers a sequence of observations \code{Y} (such as
#' electricity consumption, or any bounded time series) to be predicted step
#' by step. At each time instance \code{t}, a finite set of experts
#' (basicly some based forecasters) provide predictions \code{x} of the next
#' observation in \code{y}. This package proposes several adaptive and robust
#' methods to combine the expert forecasts based on their past performance.
#' 
#' 
#' @name opera-package
#' @aliases opera-package opera
#' @docType package
#' @author Pierre Gaillard <pierre@@gaillard.me>
#' @references Prediction, Learning, and Games. N. Cesa-Bianchi and G. Lugosi.
#' \cr
#' 
#' Forecasting the electricity consumption by aggregating specialized experts;
#' a review of sequential aggregation of specialized experts, with an
#' application to Slovakian an French contry-wide one-day-ahead (half-)hourly
#' predictions, Machine Learning, in press, 2012. Marie Devaine, Pierre
#' Gaillard, Yannig Goude, and Gilles Stoltz
#' \cr
#'  
#' Contributions to online robust aggregation: work on the approximation error and on 
#' probabilistic forecasting. Pierre Gaillard. PhD Thesis, University Paris-Sud, 2015.
#' \cr
#'
#' 
#' @keywords package
#' @template example

NULL 
