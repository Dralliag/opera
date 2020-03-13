# best convex oracle
#' @importFrom stats runif optim 
bestConvex <- function(y, experts, awake = NULL, loss.type = list(name = "square"), 
  niter = 1, ...) {
  experts <- matrix(as.numeric(as.matrix(experts)), nrow = length(y))
  N <- ncol(experts)
  
  # if there are no NA and if awake is null we can perform an exact resolution for
  # the square loss
  idx.na <- which(is.na(experts))
  if (length(idx.na) == 0 && is.null(awake) && loss.type$name == "square") {
    y.na <- is.na(y)
    y <- y[!y.na]
    x <- experts[!y.na, ]
    eq <- paste("y ~ x-1")
    
    Q <- crossprod(x)
    c <- crossprod(x, y)
    A <- cbind(1, diag(nrow(Q)))
    b <- c(1, rep(0, nrow(Q)))
    m <- 1
    res <- tryCatch({
      if (!requireNamespace("quadprog", quietly = TRUE)) {
        warning("The quadprog package must be installed to use this functionality")
        #Either exit or do something without quadprog
        return(NULL)
      } else {
        quadprog::solve.QP(Dmat = Q, dvec = c, Amat = A, bvec = b, meq = m)
      }
    }, error = function(e) {
      NULL
    })
    if (!is.null(res)) {
      coefficients <- matrix(res$solution, ncol = N)
      prediction <- experts %*% t(coefficients)
      bestLoss <- mean(loss(x = prediction, y))
    }
  } else {
    res <- NULL
  }
  if (is.null(res)) {
    warning("The best convex oracle is only approximated (using optim).")
    if (is.null(awake)) {
      awake <- as.matrix(array(1, dim(experts)))
    }
    awake[idx.na] <- 0
    experts[idx.na] <- 0
    
    lossp <- function(p) {
      return(lossConv(p, y, experts, awake, loss.type))
    }
    
    best_p <- rep(0, N)
    bestLoss <- exp(700)
    
    for (i in 1:niter) {
      # Random initialization
      p <- runif(N, 0, 1)
      p <- p/sum(p)
      
      # Convex optimization
      w <- optim(p, lossp, gr = NULL, lower = 1e-20, method = "L-BFGS-B", ...)
      
      # Projection on the simplex
      w <- pmax(w$par, 0)
      l <- lossp(w)
      if (bestLoss > l) {
        bestLoss <- l
        best_p <- w
      }
    }
    coefficients <- matrix(best_p, ncol = N)
    coefficients <- coefficients/apply(coefficients, 1, sum)
    pond <- awake %*% t(coefficients)
    prediction <- ((as.numeric(experts) * awake) %*% t(coefficients))/pond
  }
  res <- list(coefficients = coefficients, prediction = prediction)
  return(res)
} 
