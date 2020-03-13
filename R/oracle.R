#' Compute oracle predictions
#'
#' The function \code{oracle} performs a strategie that cannot be defined online
#' (in contrast to \link{mixture}). It requires in advance the knowledge of the whole
#' data set \code{Y} and the expert advice to be well defined.
#' Examples of oracles are the best fixed expert, the best fixed convex
#' combination rule, the best linear combination rule, or the best expert
#' that can shift a few times.
#'
#' @param Y  A vector containing the observations
#' to be predicted.
#' @param experts A matrix containing the experts
#' forecasts. Each column corresponds to the predictions proposed by an expert
#' to predict \code{Y}.  It has as many columns as there are experts.
#' @param model A character string specifying the oracle to use or a list with a component \code{name} specifying the oracle and any additional parameter needed.
#' Currently available oracles are:
#' \describe{
#'    \item{'expert'}{The best fixed (constant over time) expert oracle.}
#'    \item{'convex'}{The best fixed convex combination (vector of non-negative weights that sum to 1)}
#'    \item{'linear'}{The best fixed linear combination of expert}
#'    \item{'shifting'}{It computes for all number $m$ of stwitches the
#' sequence of experts with at most $m$ shifts that would have performed the
#' best to predict the sequence of observations in \code{Y}.}
#' }
#' @param loss.type A string or a list with a component 'name' specifying
#' the loss function considered to evaluate the performance. It can be
#' 'square', 'absolute', 'percentage', or 'pinball'. In the case of the pinball loss, the quantile 
#' can be provided by assigning to loss.type a list of two elements: 
#' \describe{
#'      \item{name}{A string defining the name of the loss function (i.e., 'pinball')}
#'      \item{tau}{ A number in \code{[0,1]} defining the quantile to be predicted. The default value is 0.5 to predict the median.}
#' } 
#' 
#' @param lambda A positive number used by the 'linear' oracle only. 
#' A possible $L_2$ regularization parameter for computing the linear oracle 
#' (if the design matrix is not identifiable)
#' @param niter A positive integer for 'convex' and 'linear' oracles 
#' if direct computation of the oracle is not implemented. 
#' It defines the number of optimization steps to perform in 
#' order to approximate the oracle (default value is 3).
#' 
#' @param awake A matrix specifying the
#' activation coefficients of the experts. Its entries lie in \code{[0,1]}.
#' Possible if some experts are specialists and do not always form and suggest
#' prediction. If the expert number \code{k} at instance \code{t} does not
#' form any prediction of observation \code{Y_t}, we can put
#' \code{awake[t,k]=0} so that the mixture does not consider expert \code{k} in
#' the mixture to predict \code{Y_t}. Remark that to compute the best expert oracle, 
#' the performance of unactive (or partially active) experts is computed by using 
#' the prediction of the uniform average of active experts.
#' 
#' @param ... Additional parameters
#' that are passed to \code{\link{optim}} function is order to perform convex optimization 
#' (see parameter \code{niter}).
#'
#' @return An object of class 'oracle' that contains:
#' \item{loss}{ The average loss suffered by the oracle. For the 'shifting' oracle,
#' it is a vector of length \code{T} where
#' \code{T} is the number of instance to be predicted (i.e., the length of the
#' sequence \code{Y}). The value of $loss(m)$ is the loss
#' (determined by the parameter \code{loss.type}) suffered by the
#' best sequence of expert with at
#' most $m-1$ shifts.
#' }
#' \item{coefficients}{ Not for the 'shifting' oracle. A vector containing the best weight vector corresponding to the oracle. }
#' \item{prediction}{ Not for the 'shifting' oracle. A vector containing the
#' predictions of the oracle.  }
#' \item{rmse}{If loss.type is the square loss (default) only.
#' The root mean square error (i.e., it is the square root of \code{loss}.}
#' 
#' @author Pierre Gaillard <pierre@@gaillard.me>
#' @export oracle

oracle <- function(Y, experts, model = "convex", loss.type = "square", awake = NULL, 
  lambda = NULL, niter = NULL, ...) UseMethod("oracle")


#' @export 
oracle.default <- function(Y, experts, model = "convex", loss.type = "square", awake = NULL, 
  lambda = NULL, niter = NULL, ...) {
  
  # Test that Y and experts have correct dimensions
  if (is.null(Y) || is.null(experts)) {
    stop("Y and experts should not be null")
  }
  if (length(Y) == 1) {
    experts <- as.matrix(experts)
    if (nrow(experts) == 1 || ncol(experts) == 1) {
      experts <- matrix(experts, nrow = 1)
    } else {
      stop("Bad dimensions: length(Y) should be equal to nrow(experts)")
    }
  }
  d = 1
  # We convert the data to 1-dimensional data if needed
  if (length(dim(Y)) == 2 && length(dim(experts)) == 3) {
    d = dim(Y)[2]
    T = dim(Y)[1]
    
    Y = blockToSeries(Y)
    experts = blockToSeries(experts)
    if (!is.null(awake)) {
      awake = blockToSeries(awake)
    }
  }
  
  if (!(length(Y) == nrow(experts))) {
    stop("Bad dimensions: length(Y) should be equal to nrow(experts)")
  }
  
  if (is.null(loss.type)) {
    loss.type <- list(name = "square")
  }
  if (!is.list(loss.type)) {
    loss.type <- list(name = loss.type)
  }
  if (!(loss.type$name %in% c("pinball", "square", "percentage", "absolute"))) {
    stop("loss.type should be one of these: 'absolute', 'percentage', 'square', 'pinball'")
  }
  if (!is.null(loss.type$tau) && loss.type$name != "pinball") {
    warning("Unused parameter tau (loss.type != 'pinball')")
  }
  if (!is.null(lambda) && model != "linear") {
    warning("Unused lambda parameter (model != 'linear')")
  }
  if (is.null(lambda) && model == "linear") 
    lambda <- 0 
  
  if (!is.null(niter) && model != "convex" && model != "linear") {
    warning("Unused niter parameter (model should be 'convex' or 'linear')")
  }
  if (is.null(niter)) 
    niter <- 3
  
  if ((!is.null(awake) || sum(is.na(experts) > 0)) && model != "convex" && model != 
    "shifting") {
    if (model != "expert") {
      stop(paste("Sleeping or missing values not allowed for best", model, "oracle."))
    }
    else {
      warning("When experts are unactive (or sleeping), their prediction are replaced with the uniform average of active experts")
    }
  }
  
  
  if (!(model %in% c("convex", "linear", "shifting", "expert"))) {
    stop("Wrong model specification")
  }
  if (min(Y) <= 0 && loss.type$name == "percentage") {
    stop("Y should be non-negative for percentage loss function")
  }
  names.experts <- colnames(experts)
  experts <- matrix(as.numeric(as.matrix(experts)), nrow = length(Y))
  colnames(experts) <- names.experts
  
  # if we are looking for the best convex combination of experts
  if (model == "convex") {
    res <- bestConvex(Y, experts, awake = awake, loss.type = loss.type, niter = niter, 
      ...)
  }
  
  if (model == "linear") {
    res <- tryCatch(
        bestLinear(Y, experts, lambda = lambda, loss.type = loss.type),
        error = function(err) {
          bestLinear(Y, experts, lambda = 1e-14, loss.type = loss.type) 
        })
  }
  
  if (model == "shifting") {
    res <- bestShifts(Y, experts, awake = awake, loss.type = loss.type)
  }
  
  if (!is.null(awake)) {
    pond <- apply(awake,1,mean)
    pred.unif <- apply(experts * awake, 1,mean) /pond
    experts.pred <- experts * awake + pred.unif * (1-awake)
  } else {
    experts.pred <- experts
  }
  
  loss.experts <- apply(apply(experts.pred, 2, function(x) {
    loss(x, Y, loss.type = loss.type)
  }), 2, mean)
  
  if (model == "expert") {
    best.loss <- min(loss.experts)
    coefficients <- (loss.experts == best.loss)/sum(loss.experts == best.loss)
    best.expert <- which(coefficients > 0)[1]
    res <- list(loss = best.loss, coefficients = coefficients, prediction = experts.pred[, 
      best.expert], loss.experts = loss.experts)
  }
  
  res$d <- d
  res$loss.experts <- loss.experts
  res$model <- model
  res$loss.type <- loss.type
  res$call <- match.call()
  if (model != "shifting") {
    res$residuals <- Y - res$prediction
    res$loss <- mean(loss(res$prediction, Y, loss.type))
    res$rmse <- sqrt(mean(loss(res$prediction, Y, "square")))
  }
  
  # we convert the data back to d-dimensional series if needed
  if (d>1){
    Y <- seriesToBlock(Y,d)
    experts <- seriesToBlock(experts,d)
    if (!is.null(res$residuals)) {
      res$residuals <- seriesToBlock(res$residuals,d)
    }
    res$prediction <- seriesToBlock(res$prediction,d)
    if (!is.null(awake)) {
      awake <- seriesToBlock(awake,d)
    }
  }
  res$Y <- Y
  res$experts <- experts
  res$awake <- awake
  
  
  class(res) <- "oracle"
  return(res)
} 
