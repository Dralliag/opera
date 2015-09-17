#' Errors suffered by a sequence of predictions
#' 
#'  The
#' function \code{loss} computes the sequence of instantaneous losses suffered
#' by the predictions in \code{x} to predict the observation in \code{y}.
#' 
#' 
#' @param x A vector of length \code{T}
#' containing the sequence of prediction to be evaluated.
#' @param y  A vector of length \code{T} that
#' contains the observations to be predicted.
#' @param loss.type A string specifying
#' the loss function considered to evaluate the performance. It can be
#' "square", "absolute", "percentage", or "pinball". See \code{\link{loss}} for
#' more details.
#' @param tau A number in \code{[0,1]}
#' describing the quantile to be predicted. Used only if \code{loss.type =
#' "pinball"}.
#' @return  A vector of length \code{T} containing the sequence of
#' instantaneous losses suffered by the prediction \code{x}.
#' @author Pierre Gaillard <pierre@@gaillard.me>
#' @keywords ~kwd1 ~kwd2
#' @export loss
loss <-
function(x,y,loss.type = 'square') {

   if (!is.list(loss.type)) {
      loss.type = list(name = loss.type)
   }
   if (is.null(loss.type$tau) && loss.type$name == 'pinball') {
      loss.type$tau = 0.5
   }

   if (loss.type$name == 'square')
      l <- (x-y)^2
   else if (loss.type$name == 'absolute')
      l <- abs(x-y)
   else if (loss.type$name == 'percentage')
      l <- abs(x-y)/y
   else if (loss.type$name == 'pinball')
      l <- (loss.type$tau - (y<x)) * (y-x)
   else
     stop("loss.type should be one of these: 'absolute', 'percentage', 'square', 'pinball'")
   return(l)
}
