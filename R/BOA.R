
BOA <- function(y, experts, awake = NULL, loss.type = "square", loss.gradient = TRUE, 
                w0 = NULL, training = NULL, quiet = FALSE) {
  
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
  eta_inv2 <- matrix(0, ncol = N, nrow = T + 1)
  r.reg <- numeric(N)
  
  
  if (!is.null(training)) {
    w0 <- training$w0
    R <- training$R
    R.reg <- training$R.reg
    eta_inv2[1, ] <- training$eta_inv2
  }
  
  idx_nonzero <- eta_inv2[1, ] > 0
  empty <- sum(idx_nonzero) == 0
    
  if (! quiet) steps <- init_progress(T)

  for (t in 1:T) {
    if (! quiet) update_progress(t, steps)
    
    idx <- awake[t,] > 0

    w <- w0
    if (!empty) {
      R.aux <- -log(eta_inv2[t,idx_nonzero])/2 + log(w0[idx_nonzero]) + R.reg[idx_nonzero] / sqrt(eta_inv2[t, idx_nonzero])
      R.max <- max(R.aux)
      w[idx & idx_nonzero] <- sum(w0[idx & idx_nonzero]) *  exp(R.aux - R.max) / sum(exp(R.aux - R.max) )      
    }

    p <- awake[t, ]  * w /sum(awake[t, ] * w) 
    pred <- experts[t, ] %*% p
    
    weights[t, ] <- p
    prediction[t] <- pred
    
    lpred <- loss(pred, y[t], pred, loss.type = loss.type, loss.gradient = loss.gradient)
    lexp <- loss(experts[t, ], y[t], pred, loss.type = loss.type, loss.gradient = loss.gradient)
    
    # Instantaneous regret
    r <-  awake[t, ] * c(c(lpred) - lexp)
    

    # Update the learning rates
    eta_inv2[t + 1, ] <- eta_inv2[t, ] + 2.2 * r^2
    idx_nonzero <- eta_inv2[t+1, ] > 0
    if (empty) {
      empty <- sum(idx_nonzero) == 0
    }
    # Update the regret and the regularized regret used by BOA
    if (!empty) {
        r.reg[idx_nonzero] <- r[idx_nonzero] - r[idx_nonzero]^2 /sqrt(eta_inv2[t+1, idx_nonzero])
    }
    R <- R + r
    R.reg <- R.reg + r.reg
  }
  if (! quiet) end_progress()
  
  w <- w0
  if (!empty) {
    R.aux <- -log(eta_inv2[T+1,idx_nonzero])/2 + log(w0[idx_nonzero]) + R.reg[idx_nonzero] / sqrt(eta_inv2[T+1, idx_nonzero])
    R.max <- max(R.aux)
    w[idx_nonzero] <- sum(w0[idx_nonzero]) *  exp(R.aux - R.max) / sum(exp(R.aux - R.max) )      
  }
  
  object <- list(model = "BOA", loss.type = loss.type, loss.gradient = loss.gradient, 
                 coefficients = w/sum(w))
  
  object$parameters <- list(eta_inv2 = eta_inv2[1:T, ])
  object$weights <- weights
  object$prediction <- prediction
  
  object$training <- list(eta_inv2 = eta_inv2[T + 1, ], R = R, w0 = w0, R.reg = R.reg)
  class(object) <- "mixture"
  
  return(object)
} 
