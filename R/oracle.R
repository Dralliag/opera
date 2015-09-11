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
#'    \item{"expert"}{The best fixed (constant over time) expert oracle.}
#'    \item{"convex"}{The best fixed convex combination (vector of non-negative weights that sum to 1)}
#'    \item{"linear"}{The best fixed linear combination of expert}
#'    \item{"shifting"}{It computes for all number $m$ of stwitches the
#' sequence of experts with at most $m$ shifts that would have performed the
#' best to predict the sequence of observations in \code{y}.}
#' }
#' Possible optional additional parameters are:
#' \describe{
#'    \item{loss.type}{(not possible for "linear" oracle which is currently restrained to square loss) a string specifying
#' the loss function considered to evaluate the performance.  It can be
#' "squareloss", "mae", "mape", or "pinballloss". See \code{\link{loss}} for
#' more details. If "pinballloss" is chosen, the quantile to be predicted can be set 
#' with parameter \code{tau} in \code{(0,1)} is possible (the default value is 0.5 to predict the median).
#' }
#'    \item{lambda}{For "linear" oracle. A possible $L_2$ regularization parameter for computing the linear oracle (if the design matrix is not identifiable)}
#'    \item{niter}{For "convex" oracle. Number of optimization steps to process in order to approximate the best convex combination rule if it is hard to compute in a straighforward way.}
#' }
#' @param awake A matrix specifying the
#' activation coefficients of the experts. Its entries lie in \code{[0,1]}.
#' Needed if some experts are specialists and do not always form and suggest
#' prediction.  If the expert number \code{k} at instance \code{t} does not
#' form any prediction of observation \code{Y_t}, we can put
#' \code{awake[t,k]=0} so that the mixture does not consider expert \code{k} in
#' the mixture to predict \code{Y_t}.
#' @param ... If oracle = "convex". Additional parameters
#' that are passed to \code{\link{optim}} function is order to perform convex optimization.
#'
#' @return
#' \item{loss}{ The average loss suffered by the oracle. For the "shifting" oracle,
#' it is a vector of length \code{T} where
#' \code{T} is the number of instance to be predicted (i.e., the length of the
#' sequence \code{y}). The value of $loss(m)$ is the loss
#' (determined by the parameter \code{loss.type}) suffered by the
#' best sequence of expert with at
#' most $m-1$ shifts.
#' }
#' \item{weights}{ Not for the "shifting" oracle. A vector containing the best weight vector corresponding to the oracle. }
#' \item{prediction}{ Not for the "shifting" oracle. A vecor containing the
#' predictions of the oracle.  }
#' \item{rmse}{If loss.type is the square loss (default) only.
#' The root mean square error (i.e., it is the square root of \code{loss}.}
#' @author Pierre Gaillard <pierre-p.gaillard@@edf.fr>
#' @export oracle
oracle <-
  function(y, experts, oracle = "convex", awake = NULL, ...)
  {
    if (is.character(oracle)) {
      oracle = list(name = oracle)
    }
    if (is.null(oracle$loss.type)) {
      oracle$loss.type = "squareloss"
    }
    if (is.null(oracle$tau)) {
      oracle$tau = 0.5
    }
    
    # if we are looking for the best convex combination of experts
    if (oracle$name == "convex") {
      if (is.null(oracle$niter)) {
        oracle$niter = 1
      }
      return(
        bestConvex(
          y, experts, awake = awake,
          loss.type = oracle$loss.type, niter = oracle$niter, 
          tau = oracle$tau,...
        )
      )
    }
    
    if (oracle$name == "linear") {
      if (!is.null(awake)) {
        stop('Sleeping not allowed for linear oracle!')
      }
      if (oracle$loss.type != "squareloss") {
        stop('Sorry, the linear oracle is only implemented for square loss yet.')
      }
      if (is.null(oracle$lambda)) {
        oracle$lambda = 0
      }
      return(bestLinear(y,experts,lambda = oracle$lambda))
    }
    
    if (oracle$name == "shifting") {
      return(bestShifts(
        y, experts, awake = oracle$awake,
        loss.type = oracle$loss.type
      ))
    }
    
    if (oracle$name == "expert") {
      if (!is.null(awake)) {stop("Sleeping not allowed for best expert oracle.")} 
      if (!is.null(oracle$lambda)) {warning("unused lambda parameter")} 
      if (!is.null(oracle$niter)) {warning("unused niter parameter")} 
      
      #browser()
      loss.experts = apply(apply(experts, 2, function (x) {
        loss(x,y,loss.type = oracle$loss.type,tau = oracle$tau)}), 2, mean)
      best.loss = min(loss.experts)
      weights =  (loss.experts == best.loss) / sum(loss.experts == best.loss)
      best.expert = which(weights > 0)[1]
      res = list(loss = best.loss,
                 weights = weights, 
                 prediction = experts[,best.expert])
      if (oracle$loss.type == "squareloss") {
        res$rmse = sqrt(res$loss)
      }
      return(res)
    }
    
    stop('oracle parameter is wrong')
  }