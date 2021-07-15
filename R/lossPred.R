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
lossPred <- function(x, y, pred = NULL, loss.type = list(name = "square"), loss.gradient = FALSE) {
  
  npred <- length(pred)
  nx <- length(x)
  
  args <- list("x" = if (! is.function(loss.gradient) && loss.gradient == FALSE) {x} else {pred},
               "y" = y)
  
  if (! class(loss.type) == "function") {
    args <- c(args, if (length(loss.type) > 1) {loss.type[setdiff(names(loss.type), "name")]}  else {NULL})
    
    if (loss.gradient) {
      loss.gradient <- get(paste0("gradient_", loss.type$name))
    }
    loss.type <- get(paste0("loss_", loss.type$name))
  }
  
  l <- tryCatch({
    if (npred > 1 && nx > 1) {
      if (! is.function(loss.gradient) && loss.gradient == FALSE) {
        matrix(rep(do.call(loss.type, args), npred), ncol = npred) 
      } else {
        t(matrix(rep(do.call(loss.gradient, args), nx), ncol = nx)) * matrix(rep(x, npred), ncol = npred)
      }
    } else {
      if (! is.function(loss.gradient) && loss.gradient == FALSE) {
        c(do.call(loss.type, args))
      } else {
        c(do.call(loss.gradient, args)) * x
      }
    }
  }, 
  error = function(e) {
    stop("Error when trying to apply custom loss function : \n",
         e$message)
  })
  
  return(l)
} 



### loss funs + gradients

# SQUARE LOSS
loss_square <- function(x, y) {
  (x - y)^2
}
gradient_square <- function(x, y) {
  2 * (x - y)
}

# ABSOLUTE LOSS
loss_absolute <- function(x, y) {
  abs(x - y)
}
gradient_absolute <- function(x, y) {
  sign(c(x - y))
}

# PERCENTAGE LOSS
loss_percentage <- function(x, y) {
  abs(x - y) / y
}
gradient_percentage <- function(x, y) {
  1 / y * sign(c(x - y))
}

# PINBALL LOSS
loss_pinball <- function(x, y, tau) {
  c((y < x) - tau) * c(x - y)
}
gradient_pinball <- function(x, y, tau) {
  (y < x) - tau
}
