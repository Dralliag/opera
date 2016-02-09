
ewaCalib <- function(y, experts, grid.eta = 1, awake = NULL, loss.type = "square", 
  loss.gradient = TRUE, w0 = NULL, trace = F, gamma = 2, training = NULL) {
  experts <- as.matrix(experts)
  
  N <- ncol(experts)  # Number of experts
  T <- nrow(experts)  # Number of instants
  
  # Uniform initial weight vector if unspecified
  if (is.null(w0)) {
    w0 <- rep(1, N)
  }
  
  if (is.null(gamma)) {
    gamma <- 2
  }
  
  awake <- as.matrix(awake)
  idx.na <- which(is.na(experts))
  awake[idx.na] <- 0
  experts[idx.na] <- 0
  
  T0 <- 0  # number of observations in previous runs
  neta <- length(grid.eta)  # Initial number of learning parameter in the grid to be optimized
  besteta <- floor(neta)/2 + 1  # We start with the parameter eta in the middle of the grid
  eta <- rep(grid.eta[besteta], T)  # Vector of calibrated learning rates (will be filled online by the algorithm)
  cumulativeLoss <- rep(0, neta)  # Cumulative loss suffered by each learning rate
  
  R <- array(0, c(N, neta))  # Matrix of regret suffered by each learning rate (columns) against each expert (rows)
  weta <- array(w0, dim = c(N, neta))  # Weight matrix assigned by each algorithm EWA(eta) 
  weights <- matrix(0, ncol = N, nrow = T)  # Matrix of weights formed by the mixture
  prediction <- rep(0, T)  # Predictions formed by the mixture
  
  if (!is.null(training)) {
    weta <- training$weta
    w0 <- training$w0
    R <- training$R
    cumulativeLoss <- training$cumulativeLoss
    besteta <- training$besteta
    eta[1] <- grid.eta[besteta]
    T0 <- training$T
  }
  
  R.w0 <- R
  for (k in 1:neta) {
    R.w0[, k] <- R[, k] + log(w0)/grid.eta[k]
  }  # We initialize the regrets so that w0 is the initial weight vector
  
  
  for (t in 1:T) {
    # Display the state of progress of the algorithm
    if (!(t%%floor(T/10)) && trace) 
      cat(floor(10 * t/T) * 10, "% -- ")
    
    # Weights, prediction formed by EWA(eta[t]) where eta[t] is the learning rate
    # calibrated online
    weights[t, ] <- weta[, besteta] * awake[t, ]/sum(weta[, besteta] * awake[t, 
      ])
    prediction[t] <- experts[t, ] %*% weights[t, ]
    eta[t] <- grid.eta[besteta]
    
    # Weights, predictions formed by each EWA(eta) for eta in the grid 'grid.eta'
    pred <- experts[t, ] %*% t(t(weta * awake[t, ])/apply(weta * awake[t, ], 
      2, sum))
    cumulativeLoss <- cumulativeLoss + loss(pred, y[t], loss.type)  # cumulative loss without gradient trick
    lpred <- diag(lossPred(pred, y[t], pred, loss.type, loss.gradient))  # gradient loss suffered by each eta on the grid
    lexp <- lossPred(experts[t, ], y[t], pred, loss.type, loss.gradient)  # gradient loss suffered by each expert
    
    
    # Regret update
    R.w0 <- R.w0 + awake[t, ] * t(lpred - t(lexp))
    
    # Update of the best parameter
    besteta <- order(cumulativeLoss)[1]
    
    # We increase the size of the grid if the best parameter lies in an extremity
    if (besteta == neta) {
      if (trace) 
        cat(" + ")
      neweta <- grid.eta[besteta] * gamma^(1:3)
      grid.eta <- c(grid.eta, neweta)
      neta <- neta + length(neweta)
      R.w0 <- cbind(R.w0, array(0, dim = c(N, length(neweta))))
      for (k in 1:length(neweta)) {
        perfneweta <- ewa(c(training$oldY, y[1:t]), rbind(training$oldexperts, 
          matrix(experts[1:t, ], ncol = N)), neweta[k], awake = rbind(training$oldawake, 
          matrix(awake[1:t, ], ncol = N)), loss.type = loss.type, loss.gradient = loss.gradient, 
          w0 = w0)
        cumulativeLoss <- c(cumulativeLoss, perfneweta$training$cumulativeLoss)
        R.w0[, besteta + k] <- perfneweta$training$R + log(w0) / neweta[k]
      }
    }
    
    if (besteta == 1) {
      if (trace) 
        cat(" - ")
      neweta <- grid.eta[besteta]/gamma^(1:3)
      neta <- neta + length(neweta)
      besteta <- besteta + length(neweta)
      R.w0 <- cbind(array(0, dim = c(N, length(neweta))), R.w0)
      for (k in 1:length(neweta)) {
        grid.eta <- c(neweta[k], grid.eta)
        perfneweta <- ewa(y = c(training$oldY, y[1:t]), experts = rbind(training$oldexperts, 
          matrix(experts[1:t, ], ncol = N)), eta = neweta[k], awake = rbind(training$oldawake, 
          matrix(awake[1:t, ], ncol = N)), loss.type = loss.type, loss.gradient = loss.gradient, 
          w0 = w0)
        cumulativeLoss <- c(perfneweta$training$cumulativeLoss, cumulativeLoss)
        R.w0[, besteta - k] <- perfneweta$training$R + log(w0) / neweta[k]
      }
    }
    weta <- truncate1(exp(t(t(matrix(R.w0, ncol = neta)) * grid.eta)))
  }
  
  # Next weights
  w <- weta[, besteta]/sum(weta[, besteta])
  
  R <- R.w0
  for (k in 1:neta) {
    R[, k] <- R.w0[, k] - log(w0)/grid.eta[k]
  }  
  
  object <- list(model = "EWA", loss.type = loss.type, loss.gradient = loss.gradient, 
    coefficients = w)
  
  object$parameters <- list(eta = eta[1:T], grid.eta = grid.eta)
  object$weights <- weights
  object$prediction <- prediction
  
  object$training <- list(T = T0 + T, weta = weta, w0 = w0, R = R, cumulativeLoss = cumulativeLoss, 
    besteta = besteta, grid.loss = cumulativeLoss/(T0 + T), oldexperts = rbind(training$oldexperts, 
      experts), oldY = c(training$oldY, y), oldawake = rbind(training$oldawake, 
      awake))
  
  class(object) <- "mixture"
  if (trace) 
    cat("\n")
  return(object)
} 
