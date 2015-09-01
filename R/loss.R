#' Error suffered by a sequence of prediction
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
#' "squareloss", "mae", "mape", or "pinballloss". See \code{\link{loss}} for
#' more details.
#' @param tau A number in \code{[0,1]}
#' describing the quantile to be predicted. Used only if \code{loss.type =
#' "pinballloss"}.
#' @return  A vector of length \code{T} containing the sequence of
#' instantaneous losses suffered by the prediction \code{x}.
#' @author Pierre Gaillard <pierre-p.gaillard@@edf.fr>
#' @seealso 
#' \code{\link{mape}}, \code{\link{rmse}}, \code{\link{lossConv}}
#' @keywords ~kwd1 ~kwd2
#' @export loss
loss <-
function(x,y,loss.type = 'squareloss', tau = 0.1) {
   if (loss.type == 'squareloss')
      l <- (x-y)^2
   else if (loss.type == 'mae')
      l <- abs(x-y)
   else if (loss.type == 'mape')
      l <- abs(x-y)/y
   else if (loss.type == 'pinballloss')
      l <- (tau - (y>x)) * (x - y)
   return(l)
}
