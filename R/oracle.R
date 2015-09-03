#' Compute an aggregation rule
#' 
#' The function \code{oracle} performs a strategie that cannot be defined online 
#' (in contrast to \link{mixture}). It requires in advance the knowledge of the whole 
#' data set \code{y} and the expert advice to be well defined.
#' Example of oracles are the best fixed expert, the best fixed convex 
#' combination rule, the best linear combination rule, or the best expert
#' that can shift a few times. 
#'
#' @param y  A vector containing the observations
#' to be predicted.
#' @param experts A matrix containing the experts
#' forecasts. Each column corresponds to the predictions proposed by an expert
#' to predict \code{Y}.  It has as many columns as there are experts.
#' @param oracle Either a character string specifying the oracle to use or a list with a component \code{name} specifying the oracle and any additional parameter needed. 

#' Currently available oracles are:
#' \describe{
#'    \item{"expert"}{The best fixed expert oracle.}
#'    \item{""}{Fixed-share aggregation rule. As for \code{ewa}, a learning rate \strong{eta} can be chosen by the user or calibrated online. The main difference with \code{ewa} aggregation
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
#' @return  \item{weights }{ A matrix of dimension \code{c(T,N)}, with
#' \code{T} the number of instances to be predicted and \code{N} the number of
#' experts.  Each row contains the convex combination to form the predictions }
#' \item{prediction }{ A vector of length \code{T} that contains the
#' predictions outputted by the aggregation rule.  } \item{cumulativeLoss }{ The
#' cumulated loss suffered by the aggregation rule.  } \item{regret }{ An array
#' that contains the cumulated regret suffered by the aggregation rule against
#' each expert.  }
#' @author Pierre Gaillard <pierre-p.gaillard@@edf.fr>

oracle <-
function(y, experts, oracle = "convex", awake=NULL, ...)
{
  if (is.character(oracle)) {
    oracle = list(name = oracle)
  }
  if (is.null(oracle$loss.type)) {oracle$loss.type = "squareloss"}

  # if we are looking for the best convex combination of experts
  if (oracle$name == "convex") {
    if (is.null(oracle$niter)) {oracle$niter = 1}
    return(bestConvex(y, experts, awake = awake,
      loss.type = oracle$loss.type, niter = oracle$niter, ...))
  } 

  if (oracle$name == "linear") {
    if (!is.null(awake)) {stop('Sleeping not allowed for linear oracle!')}
    if (oracle$loss.type != "squareloss") {stop('Sorry, the linear oracle is only implemented for square loss yet.')}
    if (is.null(oracle$lambda)) {oracle$lambda = 0}
    return(bestLinear(y,experts,lambda=oracle$lambda))
  }

  if (oracle$name = "shifting") {
    return(bestShifts(y,experts,awake=oracle$awake,loss.type=oracle$loss.type))
  }
}