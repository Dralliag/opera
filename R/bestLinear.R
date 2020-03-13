# best linear oracle
#' @importFrom stats rnorm optim 
bestLinear <- function(y, experts, lambda = 0, loss.type = list(name = "square"), 
  niter = 1, ...) {
  experts <- as.matrix(experts)
  N <- ncol(experts)
  
  coefficients <- NULL
  if (loss.type$name == "square") {
    coefficients <- 
      tryCatch(
        solve(lambda * diag(1, ncol(experts)) + t(experts) %*% experts, t(experts) %*% y),
        error = function(err){
          lambda = 1e-14
          solve(lambda * diag(1, ncol(experts)) + t(experts) %*% experts, t(experts) %*% y)
          warning("Ill conditioned problem. Regularized with lambda = 1e-14.")
        })
      
  } else if (loss.type$name == "pinball") {
    if (is.null(loss.type$tau)) {
      loss.type$tau <- 0.5
    }
    if (!requireNamespace("quantreg", quietly = TRUE)) {
        warning("The quantreg package must be installed to use this functionality")
        #Either exit or do something without quantreg
        return(NULL)
      } else {
        coefficients <- tryCatch({
          quantreg::rq(y ~ experts - 1, tau = loss.type$tau)$coefficients
        }, error = function(e) {
          NULL
        })
      }
  }
  if (is.null(coefficients)) {
    warning("The best linear oracle is only approximated (using optim).")
    lossu <- function(u) {
      return(mean(loss(x = experts %*% matrix(u, nrow = ncol(experts)), y = y, 
        loss.type = loss.type)))
    }
    
    best_u <- rep(0, N)
    bestLoss <- exp(700)
    
    for (i in 1:niter) {
      # Random initialization
      u <- rnorm(N, 0, 1)
      
      # Convex initialization
      w <- optim(u, lossu, gr = NULL, ...)$par
      l <- lossu(w)
      if (bestLoss > l) {
        bestLoss <- l
        best_u <- w
      }
    }
    coefficients <- matrix(best_u, nrow = N)
  }
  
  prev <- experts %*% coefficients
  return(list(coefficients = c(coefficients), prediction = prev))
} 
