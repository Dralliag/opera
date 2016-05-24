#' Convert a 1-dimensional series to blocks
#' 
#' The functions \code{seriesToBlock} and \code{blockToSeries} convert 1-dimensional series into series of higher dimension.
#' For instance, suppose you have a time-series that consists of \eqn{T = 100} days of \eqn{d = 24} hours. 
#' The function seriesToBlock converts the time-series X of \eqn{Td = 2400} observations into a matrix of size \code{c(T=100,d =24)}, 
#' where each line corresponds to a specific day. This function is usefull if you need to perform the prediction day by day, instead of hour by hour.
#' The function can also be used to convert a matrix of expert prediction of dimension \code{c(dT,K)} where K is the number of experts,
#' into an array of dimension \code{c(T,d,K)}. The new arrays of observations and of expert predictions can be
#' given to the aggregation rule procedure to perform \code{d}-dimensional predictions (i.e., day predictions).
#' 
#' The function blockToSeries performs the inverse operation.
#' @param X An array or a vector to be converted.
#' @param d A positive integer defining the block size.
#' 
#' 
#' @export seriesToBlock
seriesToBlock <- function(X, d) {
  f <- function(Y){
    matrix(Y, byrow = TRUE, ncol = d)
  }
  if (is.null(dim(X)) || dim(X)[2] == 1) {
    if ((length(X) %% d) != 0) {
      stop("length(X) must be a multiple of d")
    }
    return(f(X))
  } else {
    n <- dim(X)[1]
    K <- dim(X)[2]
    if ((n %% d) != 0) {
      stop("dim(X)[1] should be a multiple of d")
    }
    M <- array(apply(X,2,f),dim = c(n/d,d,K))
    return(M)
  }
}

#' @rdname seriesToBlock 
#' @export blockToSeries
blockToSeries <- function(X){
  if (is.null(dim(X)) || length(dim(X)) > 3){
    stop("X must be of dimension 2 or 3")
  } 
  # if X is a matrix, we convert it to a vector 
  if (length(dim(X)) == 2) {
    return(c(t(X)))
  } 
  if (length(dim(X)) == 3) {
    return(array(c(aperm(X,c(2,1,3))),dim = c(dim(X)[1]*dim(X)[2],dim(X)[3])))
  }
}
