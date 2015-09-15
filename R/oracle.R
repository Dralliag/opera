#' Compute oracle predictions
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
#' @param model Either a character string specifying the oracle to use or a list with a component \code{name} specifying the oracle and any additional parameter needed.
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
#' "square", "absolute", "percentage", or "pinball". See \code{\link{loss}} for
#' more details. If "pinball" is chosen, the quantile to be predicted can be set 
#' with parameter \code{tau} in \code{(0,1)} is possible (the default value is 0.5 to predict the median).
#' }
#'    \item{lambda}{For "linear" oracle. A possible $L_2$ regularization parameter for computing the linear oracle (if the design matrix is not identifiable)}
#'    \item{niter}{For "convex" and "linear" oracle (if direct computation of the oracle is not possible). 
#'      Number of optimization steps to process in order to approximate the oracle. 
#'      (default value is 3).}
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
#' \item{coefficients}{ Not for the "shifting" oracle. A vector containing the best weight vector corresponding to the oracle. }
#' \item{prediction}{ Not for the "shifting" oracle. A vector containing the
#' predictions of the oracle.  }
#' \item{rmse}{If loss.type is the square loss (default) only.
#' The root mean square error (i.e., it is the square root of \code{loss}.}
#' @author Pierre Gaillard <pierre@@gaillard.me>
#' @export oracle

oracle <- function(y, experts, model = "convex", awake = NULL, ...) UseMethod("oracle")


#' @export 
oracle.default <-
  function(y, experts, model = "convex", awake = NULL, ...)
  {
    if (is.character(model)) {model = list(name = model)}
    if (is.null(model$loss.type)) {model$loss.type = "square"}
    if (!is.null(model$tau) && model$loss.type != "pinball") {
      warning("Unused parameter tau (loss.type != 'pinball')")
    }
    if (is.null(model$tau)) {model$tau = 0.5}
    if (!is.null(model$lambda) && model$name != "linear") {
      warning("Unused lambda parameter (model != linear)")}
    if (is.null(model$lambda)) model$lambda = 0
    if (!is.null(model$niter) && model$name!= "convex" && model$name != "linear") {
      warning("Unused niter parameter (model should be convex or linear)")} 
    if (is.null(model$niter)) model$niter = 3
    if ((!is.null(awake) || sum(is.na(experts)>0)) && model$name != "convex" && model$name != "shifting") {
      stop(paste("Sleeping or missing values not allowed for best", model$name, "oracle."))}  

    if (!(model$name %in% c("convex","linear","shifting","expert"))) {
      stop("Wrong model specification") 
    }
    if (min(y) <= 0 && model$loss.type == "percentage") {
      stop("Y should be non-negative for percentage loss function")
    }
    # if we are looking for the best convex combination of experts
    if (model$name == "convex") {
      res <- bestConvex(
          y, experts, awake = awake,
          loss.type = model$loss.type, niter = model$niter, 
          tau = model$tau,...)
    }
    
    if (model$name == "linear") {
      res <- bestLinear(y, experts, lambda = model$lambda, 
              loss.type = model$loss.type, tau = model$tau)
    }
    
    if (model$name == "shifting") {
      res <- bestShifts(y, experts, awake = model$awake,
              loss.type = model$loss.type)
    }
    
    if (model$name == "expert") {
      
      loss.experts = apply(apply(experts, 2, function (x) {
        loss(x,y,loss.type = model$loss.type,tau = model$tau)}), 2, mean)
      best.loss = min(loss.experts)
      coefficients =  (loss.experts == best.loss) / sum(loss.experts == best.loss)
      best.expert = which(coefficients > 0)[1]
      res = list(loss = best.loss,
                 coefficients = coefficients, 
                 prediction = experts[,best.expert])
      if (model$loss.type == "square") {
        res$rmse = sqrt(res$loss)
      }
    }
    
    res$model <- model$name
    res$residuals <- y - res$prediction
    res$call <- match.call()
    res$Y <- y
    res$experts <- experts
    res$loss.type <- model$loss.type
    
    class(res) <- "oracle"
    return(res)
  }
