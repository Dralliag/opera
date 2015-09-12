#' Mean Absolute Percentage Error
#' 
#'  The
#' function \code{mape} computes the mean absolute percentage error of a
#' sequence of predictions. It is an equivalent way of running
#' \code{mean(loss(x,y,loss.type='percentage'))} (see \code{\link{loss}}).
#' 
#' 
#' @param x A vector containing the sequence of
#' predictions to be evaluated.
#' @param y  A vector containing the sequence of
#' observations to be predicted.
#' @param na.rm logical. Should missing values (including NaN) be removed?
#' @return  The mean absolute percentage error of the sequence of
#' prediction \code{x}.
#' @author Pierre Gaillard <pierre-p.gaillard@@edf.fr>
#' @seealso 
#' \code{\link{mae}}, \code{\link{rmse}}, \code{\link{loss}}
#' @keywords ~kwd1 ~kwd2
#' @export mape
mape <-
function(x,y, na.rm=TRUE) {mean(abs(x-y)/y, na.rm=na.rm)}
