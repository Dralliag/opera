
MLewa <- function(y, experts, awake = NULL, loss.type = "square", loss.gradient = TRUE, 
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
  
  R <- rep(0, N)  # Regret vector
  w <- w0
  
  # weight assigned to each expert
  weights <- matrix(0, ncol = N, nrow = T)
  prediction <- rep(0, T)
  
  # Initialization or the learning parameter
  eta <- matrix(exp(350), ncol = N, nrow = T + 1)
  
  if (!is.null(training)) {
    w0 <- training$w0
    eta[1, ] <- training$eta
    R <- training$R
    # Update weights
    R.aux <- log(w0) + eta[1, ] * R
    R.max <- max(R.aux)
    w <- exp(R.aux - R.max)
  }
  
  if (! quiet) steps <- init_progress(T)
  
  for (t in 1:T) {
    if (! quiet) update_progress(t, steps)
    
    idx <- awake[t,] > 0
    R.aux <- log(w0) + eta[t, ] * R
    R.max <- max(R.aux[idx])
    w <- numeric(N)
    w[idx] <- exp(R.aux[idx] - R.max)
    
    # form the each-instant updated mixture and prediction
    p <- awake[t, ] * w/sum(awake[t, ] * w)
    pred <- experts[t, ] %*% p
    
    # form the operational mixture and the prediction of the aggregation rule
    weights[t, ] <- p
    prediction[t] <- pred
    
    # observe losses
    lpred <- loss(pred, y[t], pred, loss.type = loss.type, loss.gradient = loss.gradient)
    lexp <- loss(experts[t, ], y[t], pred, loss.type = loss.type, loss.gradient = loss.gradient)
    
    # update regret and weights
    r <- awake[t, ] * (c(c(lpred) - lexp))
    R <- R + r
    eta[t + 1, ] <- sqrt(log(N)/(log(N)/eta[t, ]^2 + r^2))
    idx <- awake[t,] > 0
    R.aux <- log(w0) + eta[t + 1, ] * R
    R.max <- max(R.aux[idx])
    w[idx] <- exp(R.aux[idx] - R.max)
  }
  if (! quiet) end_progress()
  
  w <- w/sum(w)
  
  object <- list(model = "MLewa", loss.type = loss.type, loss.gradient = loss.gradient, 
                 coefficients = w/sum(w))
  
  object$parameters <- list(eta = eta[1:T, ])
  object$weights <- weights
  object$prediction <- prediction
  
  object$training <- list(eta = eta[T + 1, ], R = R, w0 = w0)
  class(object) <- "mixture"
  return(object)
} 
