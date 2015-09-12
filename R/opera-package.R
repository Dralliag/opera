

#' Online Prediction by ExpeRts Aggregation
#' 
#' The package \code{opera} performs, for regression-oriented time-series,
#' predictions by combining a finite set of forecasts provided by the user.
#' More formally, it considers a sequence of observations \code{y} (such as
#' electricity consumption, or any bounded time serie) to be predicted instance
#' after instance. At each time instance \code{t}, a finite set of experts
#' (basicly some based forecasters) provide predictions \code{x} of the next
#' observation in \code{y}. This package proposes several adaptive and robust
#' methods to combine the experts' forecasts based on their past performance.
#' 
#' \tabular{ll}{ Package: \tab opera\cr Type: \tab Package\cr Version: \tab
#' 1.0\cr Date: \tab 2014-02\cr License: \tab Property of EDF R&D and CNRS\cr }
#' 
#' @name opera-package
#' @aliases opera-package opera
#' @docType package
#' @author Pierre Gaillard <pierre@gaillard.me>
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
#' library('opera')              # load the package
#' set.seed(1)                   
#' 
#' # -----------------------------------------------------------------
#' #              EASY IID DATA WITHOUT BREAKS
#' # -----------------------------------------------------------------
#' 
#' T = 100                       # number of instances
#' t = 1:T                       # instances
#' Y = cos(5*2*pi*t / T)         # sequence to be predicted
#' 
#' X1 = Y + 0.1*rnorm(T)         # first expert (with small average error)
#' X2 = Y + 0.3*rnorm(T)         # second expert
#' awake1 = rep(c(rep(1,9),0),T/10) # the first expert is not always available
#' awake2 = rep(1,T)             # the second expert is always available
#' 
#' X = cbind(X1,X2)              # matrix of experts
#' awake = cbind(awake1,awake2)  # activation matrix
#' 
#' matplot(X, type='l', col=2:3) # plot experts' predictions
#' lines(Y)                      # plot observations
#' 
#' # Performance of the experts
#' cat('Expert 1, rmse :', rmse(X1,Y,awake=awake1), '\n')
#' cat('Expert 2, rmse :', rmse(X2,Y,awake=awake2), '\n')
#' 
#' # Performance of taking expert 1 if available, expert 2 otherwise
#' X3 = X1 * awake[,1] + X2 * (1-awake[,1])
#' cat("Best sequence of experts in hindsight, rmse :", rmse(X3,Y), '\n\n')
#' 
#' 
#' # EWA with fixed learning rate
#' mod = mixture(y=Y, experts=X, 
#'                aggregationRule=list(name="EWA", eta=1, loss.type='square', loss.gradient=FALSE), 
#'                awake=awake)
#' # plot weights assigned to both experts (when an expert is not available its weight is 0)
#' matplot(mod$weights, type='l', main='EWA with fixed learning rate', col=2:3) 
#' cat('EWA mod, rmse :', rmse(mod$prediction,Y), '\n')
#' 
#' # EWA algorithm with gradient loss function
#' mod = mixture(y=Y, experts=X, 
#'                aggregationRule=list(name="EWA", eta=1, loss.type='square', loss.gradient=TRUE), 
#'                awake=awake)
#' matplot(mod$weights, type='l', main='EWA with gradient losses', col=2:3) 
#' cat('EWA mod with gradient losses, rmse :', rmse(mod$prediction,Y), '\n')
#' 
#' # EWA algorithm with automatic calibration of the learning parameter
#' mod = mixture(y=Y, experts=X, aggregationRule="EWA", awake=awake)
#' matplot(mod$weights, type='l', main = 'Automatic EWA', col=2:3) 
#' cat('EWA mod with automatic tuning, rmse :', rmse(mod$prediction,Y), '\n')
#' 
#' # MLpol aggregation rule
#' mod = mixture(y=Y, experts=X, aggregationRule="MLpol", awake=awake)
#' mod$prediction = apply(mod$weights*X, 1, sum)
#' matplot(mod$weights, type='l', main = 'MLpol mod', col=2:3, ylim = c(0,1))
#' cat('MLpol mod, rmse :', rmse(mod$prediction,Y), '\n')
#' 
#' 
#' # -----------------------------------------------------------------
#' #                 TIME-SERIES WITH BREAKS
#' # -----------------------------------------------------------------
#' 
#' # We now assume that there is a break in the time series and 
#' # that experts are swaped after alpha*T instances
#' alpha = 1/2
#' X[floor(alpha*T):T,] = X[floor(alpha*T):T,2:1]
#' awake[floor(alpha*T):T,] = awake[floor(alpha*T):T,2:1]
#' 
#' # Performances of the experts
#' cat('Expert 1, rmse :', rmse(X1,Y,awake=awake1), '\n')
#' cat('Expert 2, rmse :', rmse(X2,Y,awake=awake2), '\n')
#' cat("Best sequence of experts in hindsight, rmse :", rmse(X3,Y), '\n\n')
#' 
#' 
#' # EWA with fixed learning rate
#' mod = mixture(y=Y, experts=X, 
#'                aggregationRule=list(name="EWA", eta=1, loss.type='square', loss.gradient=FALSE), 
#'                awake=awake) 
#' # plot weights assigned to both experts (when an expert is not available its weight is 0)
#' matplot(mod$weights, type='l', main='EWA with fixed learning rate', col=2:3) 
#' cat('EWA mod, rmse :', rmse(mod$prediction,Y), '\n')
#' 
#' 
#' # Fixed-share with automatic tuning of learning rate
#' mod = mixture(y=Y, experts=X, aggregationRule="FS", awake=awake)
#' # plot weights assigned to both experts (when an expert is not available its weight is 0)
#' matplot(mod$weights, type='l', main='Fixed-share with automatic tuning', col=2:3) 
#' cat('Fixed-share mod, rmse :', rmse(mod$prediction,Y), '\n')
#' 
#' # MLpol mod
#' mod = mixture(y=Y, experts=X, aggregationRule="MLpol", awake=awake)
#' matplot(mod$weights, type='l', main = 'MLpol mod', col=2:3)
#' cat('MLpol mod, rmse :', rmse(mod$prediction,Y), '\n')
#' 
#' 
#' 
NULL



