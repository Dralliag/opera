\name{opera-package}
\alias{opera-package}
\alias{opera}
\docType{package}
\title{
Online Prediction by ExpeRts Aggregation
}
\description{
The package \code{opera} performs, for regression-oriented time-series, predictions by combining a finite set of forecasts provided by the user. More formally, it considers a sequence of observations \code{y} (such as electricity consumption, or any bounded time serie) to be predicted instance after instance. At each time instance \code{t}, a finite set of experts (basicly some based forecasters) provide predictions \code{x} of the next observation in \code{y}. This package proposes several adaptive and robust methods to combine the experts' forecasts based on their past performance.
}
\details{
\tabular{ll}{
Package: \tab opera\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-02\cr
License: \tab Property of EDF R&D and CNRS\cr
}
}
\author{
  Pierre Gaillard <pierre-p.gaillard@edf.fr>
}
\references{
Prediction, Learning, and Games. N. Cesa-Bianchi and G. Lugosi. \cr

Forecasting the electricity consumption by aggregating specialized experts; a review of sequential aggregation of specialized experts, with an application to Slovakian an French contry-wide one-day-ahead (half-)hourly predictions, Machine Learning, in press, 2012. Marie Devaine, Pierre Gaillard, Yannig Goude, and Gilles Stoltz 
}

\keyword{ package }

\examples{
library('opera')              # load the package
set.seed(1)                   

# -----------------------------------------------------------------
#              EASY IID DATA WITHOUT BREAKS
# -----------------------------------------------------------------

T = 100                       # number of instances
t = 1:T                       # instances
Y = cos(5*2*pi*t / T)         # sequence to be predicted

X1 = Y + 0.1*rnorm(T)         # first expert (with small average error)
X2 = Y + 0.3*rnorm(T)         # second expert
awake1 = rep(c(rep(1,9),0),T/10) # the first expert is not always available
awake2 = rep(1,T)             # the second expert is always available

X = cbind(X1,X2)              # matrix of experts
awake = cbind(awake1,awake2)  # activation matrix

matplot(X, type='l', col=2:3) # plot experts' predictions
lines(Y)                      # plot observations

# Performance of the experts
cat('Expert 1, rmse :', rmse(X1,Y,awake=awake1), '\n')
cat('Expert 2, rmse :', rmse(X2,Y,awake=awake2), '\n')

# Performance of taking expert 1 if available, expert 2 otherwise
X3 = X1 * awake[,1] + X2 * (1-awake[,1])
cat("Best sequence of experts in hindsight, rmse :", rmse(X3,Y), '\n\n')


# EWA with fixed learning rate
mixture = ewa(y=Y, experts=X, eta=1, awake=awake, loss.type='squareloss', loss.gradient=FALSE) 
# plot weights assigned to both experts (when an expert is not available its weight is 0)
matplot(mixture$weights, type='l', main='EWA with fixed learning rate', col=2:3) 
cat('EWA mixture, rmse :', rmse(mixture$prediction,Y), '\n')

# EWA algorithm with gradient loss function
mixture = ewa(y=Y, experts=X, eta=1, awake=awake, loss.type='squareloss', loss.gradient=TRUE) 
matplot(mixture$weights, type='l', main='EWA with gradient losses', col=2:3) 
cat('EWA mixture with gradient losses, rmse :', rmse(mixture$prediction,Y), '\n')

# EWA algorithm with automatic calibration of the learning parameter
mixture = ewaCalib(y=Y, experts=X, awake=awake, loss.type='squareloss', loss.gradient=TRUE) 
matplot(mixture$weights, type='l', main = 'Automatic EWA', col=2:3) 
cat('EWA mixture with automatic tuning, rmse :', rmse(mixture$prediction,Y), '\n')

# MLpol aggregation rule
mixture = MLpol(y=Y, experts=X, awake=awake, loss.gradient=TRUE)
mixture$prediction = apply(mixture$weights*X, 1, sum)
matplot(mixture$weights, type='l', main = 'MLpol mixture', col=2:3, ylim = c(0,1))
cat('MLpol mixture, rmse :', rmse(mixture$prediction,Y), '\n')


# -----------------------------------------------------------------
#                 TIME-SERIES WITH BREAKS
# -----------------------------------------------------------------

# We now assume that there is a break in the time series and that experts are swaped after alpha*T instances
alpha = 1/2
X[floor(alpha*T):T,] = X[floor(alpha*T):T,2:1]
awake[floor(alpha*T):T,] = awake[floor(alpha*T):T,2:1]

# Performances of the experts
cat('Expert 1, rmse :', rmse(X1,Y,awake=awake1), '\n')
cat('Expert 2, rmse :', rmse(X2,Y,awake=awake2), '\n')
cat("Best sequence of experts in hindsight, rmse :", rmse(X3,Y), '\n\n')


# EWA with fixed learning rate
mixture = ewa(y=Y, experts=X, eta=1, awake=awake, loss.type='squareloss', loss.gradient=FALSE) 
# plot weights assigned to both experts (when an expert is not available its weight is 0)
matplot(mixture$weights, type='l', main='EWA with fixed learning rate', col=2:3) 
cat('EWA mixture, rmse :', rmse(mixture$prediction,Y), '\n')


# Fixed-share with automatic tuning of learning rate
mixture = fixedshareCalib(y=Y, experts=X, awake=awake, loss.type='squareloss', loss.gradient=FALSE) 
# plot weights assigned to both experts (when an expert is not available its weight is 0)
matplot(mixture$weights, type='l', main='Fixed-share with automatic tuning', col=2:3) 
cat('Fixed-share mixture, rmse :', rmse(mixture$prediction,Y), '\n')

# MLpol mixture
mixture = MLpol(y=Y, experts=X, awake=awake, loss.gradient=TRUE)
matplot(mixture$weights, type='l', main = 'MLpol mixture', col=2:3)
cat('MLpol mixture, rmse :', rmse(mixture$prediction,Y), '\n')


}
