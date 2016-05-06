## the opposite : take a series of d-dimensional blocks and concatenates them.
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