
BOA <- function(y, experts, awake = NULL, loss.type = "square", loss.gradient = TRUE, 
                w0 = NULL, training = NULL, use_cpp = getOption("opera_use_cpp", default = FALSE), quiet = FALSE) {
  
  experts <- as.matrix(experts)
  N <- ncol(experts)
  T <- nrow(experts)
  
  # Uniform initial weight vector if unspecified
  if (is.null(w0)) {
    w0 <- rep(1, N)
  }
  
  awake <- as.matrix(awake)
  idx.na <- which(is.na(experts))
  awake[idx.na] <- 0
  experts[idx.na] <- 0
  
  R <- rep(0, N)
  R.reg <- rep(0, N)
  # /!\ caution with *Copy on Write* before using RCPP
  w <- w0[]
  weights <- matrix(0, ncol = N, nrow = T)
  prediction <- rep(0, T)
  eta <- matrix(1, ncol = N, nrow = T + 1)
  V <- rep(0, N)
  B <- rep(2^{-20}, N)
  
  if (!is.null(training)) {
    w0 <- training$w0
    R <- training$R
    R.reg <- training$R.reg
    R.aux <- log(w0) + training$eta * R.reg
    w <- exp(R.aux - max(R.aux))
    eta[1, ] <- training$eta
    B <- training$B
    V <- training$V
  }
  
  if (use_cpp){
    loss_tau <- ifelse(! is.null(loss.type$tau), loss.type$tau, 0) 
    loss_name <- loss.type$name
    computeBOAEigen(awake,eta,experts,weights,y,prediction,
                    w,w0,R,R.reg,B,V,loss_name,loss_tau,loss.gradient, quiet = quiet);
  }
  else{
    if (! quiet) steps <- init_progress(T)
    
    for (t in 1:T) {
      if (! quiet) update_progress(t, steps)
      
      idx <- awake[t,] > 0
      R.aux <- log(eta[t,]) + log(w0) + eta[t, ] * R.reg
      R.max <- max(R.aux[idx])
      w <- numeric(N)
      w[idx] <- exp(R.aux[idx] - R.max)
      
      p <- awake[t, ]  * w /sum(awake[t, ] * w) 
      pred <- experts[t, ] %*% p
      
      weights[t, ] <- p
      prediction[t] <- pred
      
      lpred <- loss(pred, y[t], pred, loss.type = loss.type, loss.gradient = loss.gradient)
      lexp <- loss(experts[t, ], y[t], pred, loss.type = loss.type, loss.gradient = loss.gradient)
      
      # Instantaneous regret
      r <-  awake[t, ] * c(c(lpred) - lexp)
      
      # Update the learning rates
      B <- pmax(B, abs(r))
      B2 <- 2^ceiling(log(B,2))
      V <- V + r^2
      eta[t + 1, ] <- pmin(1/B2, sqrt(log(1/w0)/V))
      
      # Update the regret and the regularized regret used by BOA
      r.reg <- 1/2 * (r - eta[t+1, ] * r^2 + B2 * (eta[t+1,] * r > 1/2))
      R <- R + r
      R.reg <- R.reg + r.reg
    }
    if (! quiet) end_progress()
  }
  
  R.aux <- log(eta[T+1,]) + log(w0) + eta[T + 1, ] * R.reg
  R.max <- max(R.aux)
  w <-  exp(R.aux - R.max) 
  
  object <- list(model = "BOA", loss.type = loss.type, loss.gradient = loss.gradient, 
                 coefficients = w/sum(w))
  
  object$parameters <- list(eta = eta[1:T, ])
  object$weights <- weights
  object$prediction <- prediction
  
  object$training <- list(eta = eta[T + 1, ], R = R, w0 = w0, R.reg = R.reg, V= V, B=B)
  class(object) <- "mixture"
  
  return(object)
} 
