#' Root Mean Square Error
#' 
#'  The
#' function \code{rmse} computes the root mean square error of a sequence of
#' predictions. It is an equivalent way of running
#' \code{sqrt(mean(loss(x,y,loss.type='square')))} (see
#' \code{\link{loss}}).
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
#' @return  The root mean square error of the sequence of prediction
#' \code{x}.
#' @author Pierre Gaillard <pierre@@gaillard.me>
#' @keywords ~kwd1 ~kwd2
#' @export rmse
rmse <- function(x, y, awake = NULL, na.rm = TRUE) {
  if (is.null(awake)) {
    sqrt(mean((x - y)^2, na.rm = na.rm))
  } else {
    sqrt(sum((x - y)^2 * awake, na.rm = na.rm)/sum(awake * (!is.na(x - y))))
  }
} 
