#' Instantaneous loss suffered by a prediction
#' 
#'  The
#' function \code{lossPred} computes the loss of a prediction \code{x} of an
#' obstervation \code{y} knowing that we want to correct the error commited by
#' the prediction \code{pred}.  It is used in the mixing aggregation rule (see
#' \code{\link{mixture}}) to compute the loss of the experts at each
#' instance.  We can choose \code{pred} to be the prediction of \code{y}
#' outputed by the aggregation rule.
#' 
#' 
#' @param x A vector of \code{N} prediction of
#' the observation \code{y} to be evaluated.
#' @param y  A number containing the observation.
#' @param pred A reference prediction that the
#' predictions in \code{x} aim to correct.
#' @param loss.type A string specifying
#' the loss function considered to evaluate the performance. It can be
#' "squareloss", "mae", "mape", or "pinballloss". See \code{\link{loss}} for
#' more details.
#' @param loss.gradient A boolean. If
#' TRUE (default) the aggregation rule will not be directly applied to the loss
#' function at hand but to a gradient version of it. The aggregation rule is
#' then similar to gradient descent aggregation rule.
#' @param tau A number in \code{[0,1]}
#' describing the quantile to be predicted. Used only if \code{loss.type =
#' "pinballloss"}.
#' @return  A vector containing the loss suffered by the \code{N}
#' predictions in \code{x}.
#' @author Pierre Gaillard <pierre-p.gaillard@@edf.fr>
#' @seealso \code{\link{loss}}
#' @keywords ~kwd1 ~kwd2
#' @export lossPred
lossPred <-
function(x, y, pred = NULL, loss.type = 'squareloss', loss.gradient = FALSE, tau = 0.1) {
   npred <- length(pred)
   nx <- length(x)
   if (npred > 1 && nx > 1) {
      if (!loss.gradient) {
         if (loss.type == 'squareloss')
            l = matrix(rep((x-y)^2, npred), ncol = npred)
         else if (loss.type == 'mae')
            l = matrix(rep(abs(x - y), npred), ncol = npred)
         else if (loss.type == 'mape')
            l = matrix(rep(abs(x - y) / y, npred), ncol = npred)
         else if (loss.type == 'pinballloss')
            l = matrix(rep(((y<x)-tau) * (x - y), npred), ncol = npred)
      } else {
         if (loss.type == 'squareloss') 
            l = 2 * t(matrix(rep(pred - y, nx), ncol = nx)) * matrix(rep(x, npred), ncol = npred)
         else if (loss.type == 'mae')
            l = t(matrix(rep(sign(pred - y), nx), ncol = nx)) * matrix(rep(x, npred), ncol = npred)
         else if (loss.type == 'mape')
            l =  matrix(rep(x, npred), ncol = npred) / y * t(matrix(rep(sign(pred - y), nx), ncol = nx))
         else if (loss.type == 'pinballloss')
            l = t(matrix(rep((y < pred)-tau, nx), ncol = nx)) * matrix(rep(x, npred), ncol = npred)
      }
   } else {
      if (!loss.gradient) {
         if (loss.type == 'squareloss')
            l = (x-y)^2
         else if (loss.type == 'mae')
            l = abs(x - y)
         else if (loss.type == 'mape')
            l = abs(x - y) / y
         else if (loss.type == 'pinballloss')
            l = ((y<x)-tau) * (x - y)
      } else {
         if (loss.type == 'squareloss') 
            l = 2 * (pred - y) * x
         else if (loss.type == 'mae')
            l = sign(pred - y) * x
         else if (loss.type == 'mape')
            l = x / y * sign(pred - y)  
         else if (loss.type == 'pinballloss') 
            l = ((y < pred)-tau) * x
      }
   }
   return(l)
}
