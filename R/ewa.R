ewa <- function(y, experts, eta, awake = NULL, loss.type = "square", loss.gradient = TRUE, 
                w0 = NULL, training = NULL, use_cpp = getOption("opera_use_cpp", default = TRUE), quiet = FALSE) {
  experts <- as.matrix(experts)
  
  N <- ncol(experts)  # Number of experts
  T <- nrow(experts)  # Number of instants
  
  # Uniform initial weight vector if unspecified
  if (is.null(w0)) {
    w0 <- rep(1, N)
  }
  
  awake <- as.matrix(awake)
  idx.na <- which(is.na(experts))
  awake[idx.na] <- 0
  experts[idx.na] <- 0
  
  R.w0 <- log(w0)/eta  # Regret vector
  pred <- rep(0, T)  # Prediction vector
  cumulativeLoss <- 0  # Cumulative losses of the mixture
  weights <- matrix(0, ncol = N, nrow = T)  # Matrix of weights formed by the mixture
  
  if (!is.null(training)) {
    w0 <- training$w0
    R.w0 <- training$R + log(w0)/eta
    cumulativeLoss <- training$cumulativeLoss
  }
  
  if (use_cpp){
    loss_tau <- ifelse(! is.null(loss.type$tau), loss.type$tau, 0)
    loss_name <- loss.type$name
    cumulativeLoss<-computeEWAEigen(
      awake, experts,weights,y,pred, 
      R.w0,eta,cumulativeLoss,
      loss_name,loss_tau,loss.gradient, quiet = quiet);
  } 
  else {
    if (! quiet) steps <- init_progress(T)
    
    for (t in 1:T) {
      if (! quiet) update_progress(t, steps)
      
      # Weight update
      idx = awake[t,] > 0 # index of active experts
      R.aux <- eta * R.w0
      R.max = max(R.aux[idx]) # max of active experts
      weights[t,idx] <- t(exp(R.aux[idx] - R.max)) * t(awake[t,idx]) 
      weights[t,idx] <- weights[t,idx]/sum(weights[t,idx])
      
      # Prediction and losses
      pred[t] <- experts[t, ] %*% weights[t, ]
      cumulativeLoss <- cumulativeLoss + loss(pred[t], y[t], loss.type)
      lpred <- lossPred(pred[t], y[t], pred[t], loss.type, loss.gradient)
      lexp <- lossPred(experts[t, ], y[t], pred[t], loss.type, loss.gradient)
      
      # Regret update
      R.w0 <- R.w0 + awake[t, ] * (c(c(lpred) - lexp))
    }
    if (! quiet) end_progress()
  }
  R.aux <- eta * R.w0
  R.max <- max(R.aux)
  w <- t(exp(R.aux - R.max)) / sum(exp(R.aux - R.max))
  
  object <- list(model = "EWA", loss.type = loss.type, loss.gradient = loss.gradient, 
                 coefficients = w)
  
  R <- R.w0 - log(w0) / eta
  object$parameters <- list(eta = eta)
  object$weights <- weights
  object$prediction <- pred
  
  object$training <- list(R = R, w0 = w0, cumulativeLoss = cumulativeLoss)
  
  return(object)
} 
