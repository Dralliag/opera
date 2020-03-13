
BOA <- function(y, experts, awake = NULL, loss.type = "square", loss.gradient = TRUE, 
  w0 = NULL, training = NULL) {
  
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
  R.reg <- R
  weights <- matrix(0, ncol = N, nrow = T)
  prediction <- rep(0, T)
  w <- w0
  eta <- matrix(exp(350), ncol = N, nrow = T + 1)
  V <- 0
  B <- 0
  
  if (!is.null(training)) {
    w0 <- training$w0
    R <- training$R
    R.reg <- training$R.reg
    w <- truncate1(exp(log(w0) + training$eta * R.reg))
    eta[1, ] <- training$eta
    B <- training$B
    V <- training$V
  }
  
  for (t in 1:T) {
    p <- awake[t, ] * w/sum(awake[t, ] * w)
    pred <- experts[t, ] %*% p
    
    
    weights[t, ] <- p
    prediction[t] <- pred
    
    lpred <- lossPred(pred, y[t], pred, loss.type = loss.type, loss.gradient = loss.gradient)
    lexp <- lossPred(experts[t, ], y[t], pred, loss.type = loss.type, loss.gradient = loss.gradient)
    
    # Instantaneous regret
    r <-  awake[t, ] * c(c(lpred) - lexp)
    
    # Update the learning rates
    B <- pmax(B, abs(r))
    V <- V + r^2
    eta[t + 1, ] <- pmin(pmin(1/(2 * B), sqrt(log(N)/V)),exp(350))
    
    if (max(eta[t+1, ]) > exp(300)) {
      # if some losses still have not been observed
      r.reg <- r
    } else {
      r.reg <- r - eta[t+1, ] * r^2
    }
    
    # Update the regret and the regularized regret used by BOA
    R <- R + r
    R.reg <- R.reg + r.reg
    
    w <- truncate1(exp(log(w0) + eta[t + 1, ] * R.reg))
    
  }
  
  object <- list(model = "BOA", loss.type = loss.type, loss.gradient = loss.gradient, 
    coefficients = w/sum(w))
  
  object$parameters <- list(eta = eta[1:T, ])
  object$weights <- weights
  object$prediction <- prediction
  
  object$training <- list(eta = eta[T + 1, ], R = R, w0 = w0, R.reg = R.reg, V= V, B=B)
  class(object) <- "mixture"
  return(object)
} 
