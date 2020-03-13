# Instantaneous loss suffered by a prediction The function \code{lossPred}
# computes the loss of a prediction \code{x} of an obstervation \code{y}
# knowing that we want to correct the error commited by the prediction
# \code{pred}.  It is used in the mixing aggregation rule (see
# \code{\link{mixture}}) to compute the loss of the experts at each instance.
# We can choose \code{pred} to be the prediction of \code{y} outputed by the
# aggregation rule.  @param x A vector of \code{N} prediction of the observation
# \code{y} to be evaluated.  @param y A number containing the observation.
# @param pred A reference prediction that the predictions in \code{x} aim to
# correct.  @param loss.type A string specifying the loss function considered to
# evaluate the performance. It can be 'square', 'absolute', 'percentage', or
# 'pinball'. See \code{\link{loss}} for more details.  @param loss.gradient A
# boolean. If TRUE (default) the aggregation rule will not be directly applied to
# the loss function at hand but to a gradient version of it. The aggregation rule
# is then similar to gradient descent aggregation rule.  @param tau A number in
# \code{[0,1]} describing the quantile to be predicted. Used only if
# \code{loss.type = 'pinball'}.  @return A vector containing the loss suffered
# by the \code{N} predictions in \code{x}.  @author Pierre Gaillard
# <pierre@@gaillard.me> @seealso \code{\link{loss}} @keywords ~kwd1 ~kwd2
lossPred <- function(x, y, pred = NULL, loss.type = "square", loss.gradient = FALSE) {
  
  if (!is.list(loss.type)) {
    loss.type <- list(name = loss.type)
  }
  if (is.null(loss.type$tau) && loss.type$name == "pinball") {
    loss.type$tau <- 0.5
  }
  npred <- length(pred)
  nx <- length(x)
  if (npred > 1 && nx > 1) {
    if (!loss.gradient) {
      if (loss.type$name == "square") 
        l <- matrix(rep((x - y)^2, npred), ncol = npred) 
      else if (loss.type$name == "absolute") 
        l <- matrix(rep(abs(x - y), npred), ncol = npred) 
      else if (loss.type$name == "percentage") 
        l <- matrix(rep(abs(x - y)/y, npred), ncol = npred) 
      else if (loss.type$name == "log") 
        l <- matrix(rep(-log(x), npred), ncol = npred) 
      else if (loss.type$name == "pinball") 
        l <- matrix(rep(((y < x) - loss.type$tau) * (x - y), npred), ncol = npred)
    } else {
      if (loss.type$name == "square") 
        l <- 2 * t(matrix(rep(pred - y, nx), ncol = nx)) * matrix(rep(x, npred), ncol = npred) 
      else if (loss.type$name == "absolute") 
        l <- t(matrix(rep(sign(pred - y), nx), ncol = nx)) * matrix(rep(x, npred), ncol = npred) 
      else if (loss.type$name == "percentage") 
        l <- matrix(rep(x, npred), ncol = npred)/y * t(matrix(rep(sign(pred - y), nx), ncol = nx))
      else if (loss.type$name == "log") 
        l <- 2 * t(matrix(rep(-1/pred, nx), ncol = nx)) * matrix(rep(x, npred), ncol = npred) 
      else if (loss.type$name == "pinball") 
        l <- t(matrix(rep((y < pred) - loss.type$tau, nx), ncol = nx)) * matrix(rep(x, npred), ncol = npred)
    }
  } else {
    if (!loss.gradient) {
      if (loss.type$name == "square") 
        l <- (c(x - y))^2 
      else if (loss.type$name == "absolute") 
        l <- abs(c(x - y))
      else if (loss.type$name == "percentage") 
        l <- c(abs(x - y)/y)
      if (loss.type$name == "log") 
        l <- -log(c(x))
      else if (loss.type$name == "pinball") 
        l <- c((y < x) - loss.type$tau) * c(x - y)
    } else {
      if (loss.type$name == "square") 
        l <- 2 * c(pred - y) * x 
      else if (loss.type$name == "absolute") 
        l <- sign(c(pred - y)) * x 
      else if (loss.type$name == "percentage") 
        l <- x/y * sign(c(pred - y))
      if (loss.type$name == "log") 
        l <- -c(x/pred)
      else if (loss.type$name == "pinball") 
        l <- c((y < pred) - loss.type$tau) * x
    }
  }
  return(l)
} 
