
BOA <- function(y, experts, awake = NULL, loss.type = "square", loss.gradient = TRUE, w0 = NULL, training = NULL) {
  
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
  weights <- matrix(0, ncol = N, nrow = T)
  prediction <- rep(0, T)
  w <- w0
  eta <- matrix(exp(350), ncol = N, nrow = T + 1)
  
  if (!is.null(training)) {
    w0 <- training$w0
    R <- training$R
    w <- truncate1(exp(log(w0) + training$eta * R))
    eta[1, ] <- training$eta
  }
  
  for (t in 1:T) {
    p <- awake[t, ] * w/sum(awake[t, ] * w)
    pred <- experts[t, ] %*% p
    
    weights[t, ] <- p
    prediction[t] <- pred
    
    lpred <- lossPred(pred, y[t], pred, loss.type = loss.type, loss.gradient = loss.gradient)
    lexp <- lossPred(experts[t, ], y[t], pred, loss.type = loss.type, loss.gradient = loss.gradient)
    
    if (max(eta[t, ]) > exp(300)) {
      # if some losses still have not been observed
      r <- awake[t, ] * (lpred - lexp)
    } else {
      r <- awake[t, ] * (lpred - lexp) - eta[t, ] * (awake[t, ] * (lpred - lexp))^2
    }
    R <- R + r
    eta[t + 1, ] <- sqrt(log(N)/(log(N)/eta[t, ]^2 + (awake[t, ] * (lpred - lexp))^2))
    w <- truncate1(exp(log(w0) + eta[t + 1, ] * R))
  }
  
  object <- list(model = "BOA", loss.type = loss.type, loss.gradient = loss.gradient, coefficients = w/sum(w))
  
  object$parameters <- list(eta = eta[1:T, ])
  object$weights <- weights
  object$prediction <- prediction
  
  object$training <- list(eta = eta[T + 1, ], R = R, w0 = w0)
  class(object) <- "mixture"
  return(object)
} 
