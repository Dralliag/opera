 #' Follow the penalized leader aggregation rule for pinball loss
#' 
#'  The
#' function \code{pinballHour} is a mixing aggregation rule for quantile
#' regression.  At each instance, it forms the mixture by performing a convex
#' minimisation. It chooses the mixture that minimizes among the past a
#' penalized criterion based on cumulated pinball loss.
#' 
#' 
#' @param y  A vector containing the observations
#' to be predicted.
#' @param experts A matrix containing the
#' experts forecasts. Each column corresponds to the predictions proposed by an
#' expert to predict \code{Y}. It has as many columns as there are experts.
#' @param lambda A positive penalty
#' coefficient.
#' @param href A number in \code{[1,period]}
#' specifying the instant in the day when the aggregation rule can update its
#' weights.  It should lie in the interval \code{c(1,period)}.
#' @param period The number of instants in
#' each day.
#' @param w0 A vector containing the prior
#' weights of the experts.
#' @param tau A number in \code{[0,1]}
#' describing the quantile to be predicted. Used only if \code{loss.type =
#' "pinballloss"}.
#' @return  \item{weights}{ A matrix of dimension \code{c(T,N)}, with
#' \code{T} the number of instances to be predicted and \code{N} the number of
#' experts.  Each row contains the convex combination to form the predictions.
#' } \item{prediction}{ A vector of length \code{T} that contains the quantiles
#' predictions outputted by the aggregation rule.  }
#' @author Pierre Gaillard <pierre-p.gaillard@@edf.fr>
#' @seealso 
#' \code{\link{loss}}
#' @keywords ~kwd1 ~kwd2

pinballHour <-
function(y, experts, lambda, href = 1, period = 1, w0 = NULL, tau = 0.5) {
   experts = as.matrix(experts)
   N <- ncol(experts)
   T <- nrow(experts)
   if (is.null(w0)) {w0 <- rep(1/N,N)} # Uniform intial weight vector if unspecified
   if (sum(is.na(experts)) > 0) {warning("There are not allowed NA's in expert advice")}

   
   pinballloss = function(beta, Y, X) {
      beta = as.matrix(beta, ncol = 1)
      X = as.matrix(X)
      pred = X%*%beta
      return(mean((tau - (Y > pred)) * (pred - Y)) + lambda * sum(beta^2))
   }
   
   if (is.null(w0)) {w0 <- rep(1/N,N)} # Uniform intial weight vector if unspecified
   
   weights <- matrix(0, ncol = N, nrow = T)
   w = w0
   for (t in 1:T){  
      weights[t,] = w
      h <- (((t - 1) %% period) + 1)
      if (h == href) {
         cat('Iteration :',t,'\n')
         w = optim(w, fn=function(beta) {pinballloss(beta, y[1:t],experts[1:t,])})$par
         }
   }
   prediction <- apply(experts * weights, 1, sum)
   return(list(weights = weights, prediction = prediction))
}
