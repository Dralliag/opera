## A function to convert a R-series of concatenated blocks into an array of d-dimensional elements
seriesToBlock <- function(X, d) {
  f <- function(Y){
    matrix(Y, byrow = TRUE, ncol = d)
  }
  if (is.null(dim(X))) {
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

