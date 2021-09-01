#' Errors suffered by a sequence of predictions 
#' 
#' The
#' function \code{loss} computes the sequence of instantaneous losses suffered
#' by the predictions in \code{x} to predict the observation in \code{y}.
#' 
#' @param x \code{numeric}. A vector of length \code{T} containing the sequence of prediction to be evaluated.
#' @param y \code{numeric}. A vector of length \code{T} that contains the observations to be predicted.
#' @param pred \code{numeric}. A vector of length \code{T} containing the sequence of real values.
#' @param loss.type \code{character, list or function} ("square").
#' \itemize{
#'      \item{character}{ Name of the loss to be applied ('square', 'absolute', 'percentage', or 'pinball');}
#'      \item{list}{ List with field \code{name} equal to the loss name. If using pinball loss, field \code{tau} 
#'      equal to the required quantile in [0,1];}
#'      \item{function}{ A custom loss as a function of two parameters.}
#' }
#' @param loss.gradient \code{boolean, function} (TRUE). 
#' \itemize{
#'      \item{boolean}{ If TRUE, the aggregation rule will not be directly applied to the loss function at hand,
#'      but to a gradient version of it. The aggregation rule is then similar to gradient descent aggregation rule. }
#'      \item{function}{ If loss.type is a function, the derivative should be provided to be used (it is not automatically 
#'      computed).}
#' }
#' 
#' @return  A vector of length \code{T} containing the sequence of
#' instantaneous losses suffered by the expert previsions (x) or the gradient computed on the aggregated previsions (pred).
#' 
#' @author Pierre Gaillard <pierre@@gaillard.me>
#' @export
loss <- function(x, y, pred = NULL, loss.type = list(name = "square"), loss.gradient = FALSE) {
  
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
