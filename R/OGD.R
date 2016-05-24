OGD <- function(y, experts, loss.type = "square", training = NULL, alpha, simplex, w0 = NULL) {
  
  experts <- as.matrix(experts)
  N <- ncol(experts)
  T <- nrow(experts)
  
  # weight assigned to each expert
  weights <- matrix(0, ncol = N, nrow = T)
  prediction <- rep(0, T)
  
  # Initialization (number of previous instances already trained)
  t0 <- 0
  eta <- numeric(T+1)
  
  # Initial weights
  if (is.null(w0)) {
    w0 <- rep(0, N)
  }
  w <- w0
  if (simplex) {
    w <- simplexProj(w0)
  }
  B <- 0
  
  # Previous training ?
  if (!is.null(training)) {
    t0 <- training$t0
    w <- training$w
    B <- training$B
  } else {
    training <- list()
  }
  
  for (t in 1:T) {
    pred <- experts[t, ] %*% w
    
    # save the mixture and the prediction
    weights[t, ] <- w
    prediction[t] <- pred
    
    # Observe losses
    lexp <- lossPred(experts[t, ], y[t], pred, loss.type = loss.type, loss.gradient = TRUE)
    B <- max(B, sqrt(sum(lexp^2)))
    
    # Update the learning rate
    eta[t+1] <-  (t+t0)^(-alpha) / B
    
    # Gradient step
    w <- w - eta[t + 1] * lexp 
    if (simplex) {
      w <- simplexProj(w)
    }
  }


  object <- list(model = "OGD", loss.type = loss.type, loss.gradient = TRUE, 
                 coefficients = w)
  
  object$parameters <- list(alpha = alpha,
                            simplex = simplex,
                            w0 = w0)
  object$weights <- weights
  object$prediction <- prediction
  
  object$training <- list(w = w, t0 = t0 + T, B = B)
  class(object) <- "mixture"
  return(object)
} 
