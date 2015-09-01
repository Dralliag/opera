#' Compute an aggregation rule
#' 
#' The function \code{mixture} performs an
#' aggregation rule chosen by the user. It considers a sequence \code{y} of observations to be predicted
#' sequentially with the help of experts advices \code{x}.  The forms at each
#' instance \code{t} a prediction by assigning weight to the experts advices
#' and combining them.
#' 
#' 
#' @param y  A vector containing the observations
#' to be predicted.
#' @param experts A matrix containing the experts
#' forecasts. Each column corresponds to the predictions proposed by an expert
#' to predict \code{Y}.  It has as many columns as there are experts.
#' @param aggregationRule Either a character string specifying the aggregation rule to use or a list with a component \code{name} specifying the aggregation rule and any additional parameter needed. 

#' Currently available aggregation rules are:
#' \describe{
#'    \item{"EWA"}{Exponentially weighted average aggregation rule. A positive learning rate \strong{eta} can be chosen by the user. The
#' bigger it is the faster the aggregation rule will learn from observations
#' and experts performances. However too hight values lead to unstable weight
#' vectors and thus unstable predictions. If it is not specified, the learning rate is calibrated online. }
#'    \item{"FS"}{Fixed-share aggregation rule. As for \code{ewa}, a learning rate \strong{eta} can be chosen by the user or calibrated online. The main difference with \code{ewa} aggregation
#' rule rely in the mixing rate \strong{alpha}\eqn{\in [0,1]} wich considers at
#' each instance a small probability \code{alpha} to have a rupture in the
#' sequence and that the best expert may change. Fixed-share aggregation rule
#' can thus compete with the best sequence of experts that can change a few
#' times (see \code{\link{bestShifts}}), while \code{ewa} can only
#' compete with the best fixed expert. The mixing rate is either chosen by the user either calibrated online.}
#'    \item{"Ridge"}{Ridge regression. It minimizes at
#' each instance a penalized criterion.  It forms at each instance linear
#' combination of the experts' forecasts and can assign negative weights that
#' not necessarily sum to one.  It is useful if the experts are biased or
#' correlated. It cannot be used with specialized experts. A positive regularization coefficient \strong{lambda} can either be chosen by the user or calibrated online. }
#'    \item{"MLpol"}{Polynomial Potential aggregation rule
#' with different learning rates for each expert.  The learning rates are
#' calibrated using theoretical values. There are similar aggregation rules like "BOA" (Bernstein online Aggregation see [Wintenberger, 2014] "MLewa", and "MLprod" (see [Gaillard, Erven, and Stoltz, 2014])}
#'    \item{"pinball"}{ It performs a mixing aggregation rule for quantile
#' regression.  At each instance, it forms the mixture by performing a convex
#' minimisation. It chooses the mixture that minimizes among the past a
#' penalized criterion based on cumulated pinball loss.}
#'    \item{"gamMixture"}{#'  Fits a general additive model (GAM) on the data
#' to form weight vectors for the experts that can depend on exogeneous data.
#' The process is however not currently stable. We advice not using it yet. Use
#' only at most one exogeneous variable.}
#' }
#' Possible optional additional parameters are:
#' \describe{
#'    \item{loss.type}{(not possible for "ridge", "pinball", and "gamMixture") a string specifying
#' the loss function considered to evaluate the performance.  It can be
#' "squareloss", "mae", "mape", or "pinballloss". See \code{\link{loss}} for
#' more details. If "pinballloss" is chosen, an additional parameter \code{tau} in \code{[0,1]} is required. }   
#'    \item{loss.gradient}{A boolean. If
#' TRUE (default) the aggregation rule will not be directly applied to the loss
#' function at hand but to a gradient version of it.  The aggregation rule is
#' then similar to gradient descent aggregation rule.}
#' }
#' @param w0 A vector containing the prior weights of the experts. (not possible for "MLpol", "BOA", "MLewa", and "MLprod")
#' @param awake A matrix specifying the
#' activation coefficients of the experts. Its entries lie in \code{[0,1]}.
#' Needed if some experts are specialists and do not always form and suggest
#' prediction.  If the expert number \code{k} at instance \code{t} does not
#' form any prediction of observation \code{Y_t}, we can put
#' \code{awake[t,k]=0} so that the mixture does not consider expert \code{k} in
#' the mixture to predict \code{Y_t}.
#' @param href A number in \code{[1,period]}
#' specifying the instant in the day when the aggregation rule can update its
#' weights.  It should lie in the interval \code{c(1,period)}.
#' @param period The number of instants in
#' each day.
#' @param delay A positive number that indicates the number of instants before
#' the mixture has access to the true observations in \code{y}.  If \code{delay
#' > 0}, \code{y.ETR} can not be \code{NULL}.
#' @param y.ETR A vector containing real time estimations of \code{y} to be
#' used at each instant \code{t} instead of \code{y_t} to predict instants
#' \code{t+1,...,t+delay}.

#' @return  \item{weights }{ A matrix of dimension \code{c(T,N)}, with
#' \code{T} the number of instances to be predicted and \code{N} the number of
#' experts.  Each row contains the convex combination to form the predictions }
#' \item{prediction }{ A vector of length \code{T} that contains the
#' predictions outputted by the aggregation rule.  } \item{cumulativeLoss }{ The
#' cumulated loss suffered by the aggregation rule.  } \item{regret }{ An array
#' that contains the cumulated regret suffered by the aggregation rule against
#' each expert.  }
#' @author Pierre Gaillard <pierre-p.gaillard@@edf.fr>
#' @keywords ~kwd1 ~kwd2
#' @examples
#' 
#' library('opera')              # load the package
#' set.seed(1)                   
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
#' mod = mixture(y=Y, experts=X, aggregationRule = list(name="EWA", eta=1, loss.type='squareloss', loss.gradient=FALSE), awake=awake) 
#' # plot weights assigned to both experts (when an expert is not available its weight is 0)
#' matplot(mod$weights, type='l', main='EWA with fixed learning rate', col=2:3) 
#' cat('EWA mixture, rmse :', rmse(mod$prediction,Y), '\n')
#' 
#' # ewa algorithm with gradient loss function
#' mod = mixture(y=Y, experts=X, aggregationRule = list(name = "EWA", eta=1, loss.type='squareloss', loss.gradient=TRUE), awake=awake) 
#' matplot(mod$weights, type='l', main='EWA with gradient losses', col=2:3) 
#' cat('EWA mixture with gradient losses, rmse :', rmse(mod$prediction,Y), '\n')
#' 
#' # ewa algorithm with automatic calibration of the learning parameter
#' mod = mixture(y=Y, experts=X, aggregationRule = "EWA", awake = awake)
#' matplot(mod$weights, type='l', main = 'Automatic EWA', col=2:3) 
#' cat('EWA mixture with automatic tuning, rmse :', rmse(mod$prediction,Y), '\n')
#' 
#' # MLpol aggregation rule
#' mod = mixture(y=Y, experts=X, aggregationRule="MLpol", awake = awake)
#' mod$prediction = apply(mod$weights*X, 1, sum)
#' matplot(mod$weights, type='l', main = 'MLpol mixture', col=2:3, ylim = c(0,1))
#' cat('MLpol mixture, rmse :', rmse(mod$prediction,Y), '\n')
#' 
#' @export mixture
mixture <-
function(y, experts, 
  aggregationRule = "MLpol",  w0 = NULL,
  awake = NULL, href = 1, period = 1, delay = 0, y.ETR = NULL)
{
  if (is.character(aggregationRule)) {
    aggregationRule = list(name = aggregationRule)
  }
  if (aggregationRule$name == "Ridge") {
    if (is.null(aggregationRule$lambda)) {
      return(ridgeCalib(y = y, experts = experts, href = href, period = period, w0 = w0, delay = delay, y.ETR = y.ETR))
    } else {
      return(ridgeHour(y, experts, aggregationRule$lambda, href, period, w0, delay, y.ETR))
    }
  } else {

    if (is.null(aggregationRule$loss.type)) {aggregationRule$loss.type = "squareloss"}
    if (is.null(aggregationRule$loss.gradient)) {aggregationRule$loss.gradient = TRUE}
    if (is.null(aggregationRule$tau)) {aggregationRule$tau = 0.5}

    if (aggregationRule$name == "MLpol" || aggregationRule$name == "MLprod") {
      if (!is.null(w0)) { stop(paste(aggregationRule$name, "cannot handle non-uniform prior weight vector"))}
      algo <- eval(parse(text = aggregationRule$name))
      return(algo(y, experts, awake = awake, loss.type = aggregationRule$loss.type, loss.gradient = aggregationRule$loss.gradient, period = period, href = href, delay = delay, y.ETR = y.ETR, tau = aggregationRule$tau))
    }

    if (aggregationRule$name == "EWA") {
      if (is.null(aggregationRule$eta)) {
        return(ewaCalib(y = y, experts = experts, awake = awake, loss.type = aggregationRule$loss.type, loss.gradient = aggregationRule$loss.gradient, w0 = w0, href = href, period = period, delay = delay, y.ETR = y.ETR))
      } else {
        return(ewaHour(y = y, experts = experts, eta = aggregationRule$eta, awake = awake, loss.type = aggregationRule$loss.type, loss.gradient = aggregationRule$loss.gradient, w0 = w0, href = href, period = period, delay = delay, y.ETR = y.ETR))
      }
    }

    if ((aggregationRule$name == "BOA")|| aggregationRule$name == "MLewa") {
      if (!is.null(w0)) { stop(paste(aggregationRule$name, "cannot handle non-uniform prior weight vector"))}
      if ((delay > 0) || (!is.null(y.ETR))) {
        stop(paste(aggregationRule$name, "cannot handle delayed signal"))
      }
      algo <- eval(parse(text = aggregationRule$name))
      return(algo(y, experts, awake = awake, loss.type = aggregationRule$loss.type, loss.gradient = aggregationRule$loss.gradient, period = period, href = href))
    }

    if ((aggregationRule$name == "FS")) {
      if (!is.null(y.ETR) || (delay > 0)) {
        warning("Fixed-share aggregation rule cannot handle delayed signal well")
      }
      if (is.null(aggregationRule$eta) || is.null(aggregationRule$alpha)) {
        return(fixedshareCalib(y = y, experts = experts, awake = awake, loss.type = aggregationRule$loss.type, loss.gradient = aggregationRule$loss.gradient, w0 = w0, href = href, period = period))
      } else {
        return(fixedshareHour(y = y, experts = experts, eta = aggregationRule$eta, alpha = aggregationRule$alpha, awake = awake, loss.type = aggregationRule$loss.type, loss.gradient = aggregationRule$loss.gradient, w0 = w0, href = href, period = period))
      }
    }

    if ((aggregationRule$name == "gamMixture")) {
      warning("This aggregation rule is not stable in the current version")
      if (is.null(aggregationRule$lambda)) {
        stop("gamMixture cannot handle automatic calibration")
      }
      if (is.null(aggregationRule$z)){
        stop("A matrix of exogeneous variables z must be given")
      }
      if (is.null(aggregationRule$nknots)){ aggregationRule$nknots = 5}
      if (is.null(aggregationRule$degree)){ aggregationRule$degree = 3}
      if (is.null(aggregationRule$loss.type)){ aggregationRule$loss.type = "squareloss"}
      if (is.null(aggregationRule$uniform)){ aggregationRule$uniform = FALSE}
      if (is.null(aggregationRule$knots)){ aggregationRule$knots = NULL}
      if (is.null(aggregationRule$tau)){ aggregationRule$tau = 0.5}

      return(gamMixture(y = y, experts = experts, z = , aggregationRule$lambda, nknots = aggregationRule$nknots, degree = aggregationRule$degree, loss.type = aggregationRule$loss.type, href = href, period = period, uniform = aggregationRule$uniform, knots = aggregationRule$knots, tau = aggregationRule$tau))

    }

    if ((aggregationRule$name == "pinball")) {
      if (is.null(aggregationRule$lambda)) {
        stop("pinball cannot handle automatic calibration")
      }
      if (is.null(aggregationRule$tau)){ aggregationRule$tau = 0.5}
      pinballHour(y = y, experts = experts, lambda = aggregationRule$lambda, href = href, period = period, w0 = w0, tau = aggregationRule$tau)
    }


  }

}