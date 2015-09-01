#' Mean Absolute Error
#' 
#'  The
#' function \code{mae} computes the average absolute error of a sequence.
#' 
#' 
#' @param x A vector containing the sequence of
#' predictions to be evaluated.
#' @param y  A vector containing the sequence of
#' observations to be predicted.
#' @param awake A vector specifying the
#' activation coefficients of the predictions in \code{x}. Its entries lie in
#' \code{[0,1]}.
#' @param na.rm logical. Should missing values (including NaN) be removed?
#' @return  Average absolute error of the sequence
#' @author Pierre Gaillard <pierre-p.gaillard@@edf.fr>
#' @seealso 
#' \code{\link{rmse}}, \code{\link{mape}}, \code{\link{loss}}
#' @keywords ~kwd1 ~kwd2
#' @export mae
mae <-
function (x,y,awake=NULL, na.rm=TRUE) {
  if (is.null(awake)) {
    mean(abs(x-y), na.rm=na.rm)
  }
  else {
    sum(abs(x-y) * awake, na.rm=na.rm) / sum(awake * (!is.na(x-y)), na.rm=na.rm)
  }
}
