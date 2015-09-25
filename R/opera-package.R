#' Online Prediction by ExpeRts Aggregation
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
#' \tabular{ll}{ Package: \tab opera\cr Type: \tab Package\cr Version: \tab
#' 1.0\cr Date: \tab 2014-02\cr License: \tab Property of EDF R&D and CNRS\cr }
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
#' @keywords package
#' @examples
#' 
#'library('opera')  # load the package
#'set.seed(1)
#'
#'# Example: find the best one week ahead forecasting strategy (weekly data)
#'# packages
#'library(mgcv)
#'library(caret)
#'
#'# import data
#'data(electric_load)
#'idx_data_test <- 680:nrow(electric_load)
#'data_train <- electric_load[-idx_data_test, ]
#'data_test <- electric_load[idx_data_test, ]
#'
#'# Medium term model to remove trend and seasonality (using generalized additive model)
#'detrend.fit <- gam(Load ~ s(Time,k=3) + s(NumWeek) + s(Temp) + s(IPI), data = data_train)
#'electric_load$Trend <- c(predict(detrend.fit), predict(detrend.fit,newdata = data_test))
#'electric_load$Load.detrend <- electric_load$Load - electric_load$Trend
#'
#'# a few graphs to display the data
#'attach(data_train)
#'plot(Load, type = 'l')
#'plot(Temp, Load, pch = 16, cex = 0.5)
#'plot(NumWeek, Load, pch = 16, cex = 0.5)
#'plot(Load, Load1, pch = 16, cex = 0.5)
#'acf(Load, lag.max = 20)
#'detach(data_train)
#'
#'# Build the expert forecasts 
#'# ##########################
#'
#'# A generalized additive model
#'gam.fit <- gam(Load ~ s(IPI) + s(Temp) + s(Time, k=3) + 
#'                s(Load1) + as.factor(NumWeek), data = data_train)
#'gam.forecast <- predict(gam.fit, newdata = data_test)
#'rmse(gam.forecast,data_test$Load)
#'
#'# An online autoregressive model on the residuals of the medium term model
#'ar.forecast <- numeric(length(idx_data_test))
#'for (i in seq(idx_data_test)) {
#'  ar.fit <- ar(electric_load$Load.detrend[1:(idx_data_test[i] - 1)])
#'  ar.forecast[i] <- as.numeric(predict(ar.fit)$pred) + electric_load$Trend[idx_data_test[i]]
#'}
#'rmse(ar.forecast,data_test$Load) 
#'
#'# A GBM
#'gbm0.fit <- train(Load ~ IPI + IPI_CVS + Temp + Temp1 + Time + Load1 + NumWeek, 
#'                  data = data_train, method = 'gbm')
#'gbm.forecast <- predict(gbm0.fit, newdata = data_test)
#'rmse(gbm.forecast,data_test$Load) 
#'
#'
#'# Aggregation of experts
#'###########################
#'
#'X <- cbind(gam.forecast, ar.forecast, gbm.forecast)
#'colnames(X) <- c('gam', 'ar', 'gbm')
#'
#'Y <- data_test$Load
#'
#'matplot(cbind(Y, X), type = 'l', col = 1:6, ylab = 'Weekly load', xlab = 'Week')
#'
#'
#'# How good are the expert? Look at the oracles
#'oracle.convex <- oracle(Y = Y, experts = X, loss.type = 'square', model = 'convex')
#'plot(oracle.convex)
#'oracle.convex
#'
#'# Is a single expert the best over time ? Are there breaks ?
#'oracle.shift <- oracle(Y = Y, experts = X, loss.type = 'percentage', model = 'shifting')
#'plot(oracle.shift)
#'oracle.shift
#'
#'# Online aggregation of the experts with MLpol
#'#############################################
#'
#'# Initialize the aggregation rule
#'m0.MLpol <- mixture(model = 'MLpol', loss.type = 'square')
#'
#'# Perform online prediction using MLpol There are 3 equivalent possibilities 1)
#'# start with an empty model and update the model sequentially
#'m1.MLpol <- m0.MLpol
#'for (i in 1:length(Y)) {
#'  m1.MLpol <- predict(m1.MLpol, newexperts = X[i, ], newY = Y[i])
#'}
#'
#'# 2) perform online prediction directly from the empty model
#'m2.MLpol <- predict(m0.MLpol, newexpert = X, newY = Y, online = TRUE)
#'
#'# 3) perform the online aggregation directly
#'m3.MLpol <- mixture(Y = Y, experts = X, model = 'MLpol', loss.type = 'square')
#'
#'# These predictions are equivalent:
#'identical(m1.MLpol, m2.MLpol)  # TRUE
#'identical(m1.MLpol, m3.MLpol)  # TRUE
#'
#'# Display the results
#'summary(m3.MLpol)
#'plot(m3.MLpol) 

NULL 
