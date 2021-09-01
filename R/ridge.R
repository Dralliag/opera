
ridge <- function(y, experts, lambda, w0 = NULL, training = NULL,
                  use_cpp = getOption("opera_use_cpp", default = TRUE), quiet = FALSE) {
  
  new = getOption("opera_use_new", default = FALSE)
  
  experts <- as.matrix(experts)
  N <- ncol(experts)
  T <- nrow(experts)
  
  # Uniform initial weight vector if unspecified
  if (is.null(w0)) {
    w0 <- matrix(1/N, ncol = N)
  }
  if (sum(is.na(experts)) > 0) {
    warning("There are not allowed NA's in expert advice")
  }
  
  w <- matrix(0, ncol = N, nrow = T)
  
  if (!is.null(training)) {
    At <- training$At
    bt <- training$bt
  } else {
    At <- 1/lambda * diag(1, N)
    bt <- matrix(lambda * w0, nrow = N)
  }
  
  if (use_cpp){
    error_code<-computeRidgeCPP(experts,w,At,bt,y, quiet = quiet)
    if (error_code != 0){
      stop("matrix is not invertible")
    }
  }
  else {
    if (!quiet) steps <- init_progress(T)
    
    for (t in 1:T) {
      if (!quiet) update_progress(t, steps)
      
      w[t, ] <- At %*% bt
      a <- At %*% experts[t, ]
      At <- At - a %*% t(a) / c(1 + experts[t,] %*% a)
      bt <- bt + y[t] * experts[t, ]
    }
    if (! quiet) end_progress()
  }
  
  object <- list(model = "Ridge", loss.type = list(name = "square"), coefficients = At%*%bt)
  
  object$parameters <- list(lambda = lambda)
  object$weights <- w
  object$prediction <- rowSums(experts * w)
  
  object$training <- list(At = At, bt = bt)
  
  return(object)
}
