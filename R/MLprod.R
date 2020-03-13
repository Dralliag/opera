MLprod <- function(y, experts, awake = NULL, loss.type = "square", loss.gradient = TRUE, 
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
  
  weights <- matrix(0, ncol = N, nrow = T)
  R <- log(w0)
  
  L <- rep(1, N)
  maxloss <- 0
  
  eta <- matrix(exp(700), ncol = N, nrow = T + 1)
  prediction <- rep(0, T)
  
  if (!is.null(training)) {
    eta[1, ] <- training$eta
    R <- training$R
    L <- training$L
    maxloss <- training$maxloss
  }
  
  for (t in 1:T) {
    
    # Update weights
    w <- truncate1(exp(R))
    w <- eta[t, ] * w/sum(eta[t, ] * w)
    p <- awake[t, ] * w/sum(awake[t, ] * w)
    pred <- experts[t, ] %*% p  # Predict
    
    weights[t, ] <- p
    prediction[t] <- pred
    
    # Observe losses
    lpred <- lossPred(pred, y[t], pred, loss.type = loss.type, loss.gradient = loss.gradient)
    lexp <- lossPred(experts[t, ], y[t], pred, loss.type = loss.type, loss.gradient = loss.gradient)
    
    r <- awake[t, ] * (c(c(lpred) - lexp))
    L <- L + r^2
    maxloss <- pmax(maxloss, abs(r))
    neweta <- pmin(1/(2 * maxloss), sqrt(log(N)/L))
    
    # Update regret and learning parameter
    R <- neweta/eta[t, ] * R + log(1 + awake[t, ] * neweta * (c(c(lpred) - lexp)))
    eta[t + 1, ] <- neweta
    
    if (is.na(sum(R))) {
      browser("Nan in R")
    }
  }
  
  w <- truncate1(exp(R))
  w <- eta[T + 1, ] * w/sum(eta[T + 1, ] * w)
  
  object <- list(model = "MLprod", loss.type = loss.type, loss.gradient = loss.gradient, 
    coefficients = w)
  
  object$parameters <- list(eta = eta[1:T, ])
  object$weights <- weights
  object$prediction <- prediction
  
  object$training <- list(eta = eta[T + 1, ], R = R, L = L, maxloss = maxloss)
  class(object) <- "mixture"
  
  return(object)
} 
