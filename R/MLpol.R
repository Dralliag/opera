
MLpol <-
  function(y, experts, awake = NULL,
           loss.type = 'squareloss', loss.gradient = TRUE,
           period = 1, href = 1,
           delay = 0, y.ETR = NULL,
           tau = 0.5) 
  {
    experts <- as.matrix(experts)
    N <- ncol(experts)
    T <- nrow(experts)
    
    if (is.null(awake)) {awake = matrix(1, nrow = T, ncol = N)} # Full activation if unspecified
    if (is.null(y.ETR) || (delay == 0)) {y.ETR = y}
    
    awake <- as.matrix(awake)
    idx.na <- which(is.na(experts))
    awake[idx.na] <- 0
    experts[idx.na] <- 0
    
    # weight assigned to each expert
    weights <- matrix(0, ncol = N, nrow = T)
    prediction <- rep(0, T)
    
    # Initialization or the learning parameter
    etaop <- matrix(exp(700),ncol=N,nrow=T+1)
    eta <- rep(exp(700),N)
    
    # regret suffered by each expert
    R <- rep(0,N)
    Rop <- rep(0,N)
    
    # Losses
    lpred <- numeric(T)     # Losses of the mixture
    lpred.ETR <- numeric(T) # Real time estimate of mixture losses
    lexp <- matrix(0, ncol = N, nrow = T)     # Expert losses
    lexp.ETR <- matrix(0, ncol = N, nrow = T) # Real time estimate of the expert losses
    
    for (t in 1:T) {
      # We check if there is at least one expert with positive weight
      if (max(awake[t,] * Rop) > 0) {
        wop <- etaop[t,] * pmax(Rop,0) / sum(etaop[t,] * pmax(Rop,0))
      } else {
        wop <- rep(1,N)
      }
      
      # operational weight vectors, predictions
      weights[t,] <- awake[t,] * wop / sum(awake[t,]*wop ) 
      prediction[t] <- experts[t,] %*% weights[t,]
      
      # We check if there is at least one expert with positive weight
      if (max(awake[t,] * R) > 0) {
        w <- eta * pmax(R,0) / sum(eta * pmax(R,0))
      } else {
        w <- rep(1,N)
      }
      
      # form the mixture and the prediction
      p <- awake[t,] * w / sum(awake[t,]*w ) 
      pred <- experts[t,] %*% p
      
      # Observe losses
      lpred.ETR[t] <- lossPred(pred, y.ETR[t], pred, tau = tau,
                               loss.type = loss.type, loss.gradient = loss.gradient)
      lexp.ETR[t,] <- lossPred(experts[t,], y.ETR[t], pred, tau = tau,
                               loss.type = loss.type, loss.gradient = loss.gradient)
      
      # Update regret 
      r <- awake[t,] * (lpred.ETR[t] - lexp.ETR[t,])
      R <- R + r
      
      # Update the learning rate
      eta <- 1 / (1/eta + r^2)   
      
      # We readjust by using the true observations
      if (delay > 0) {
        lpred[t] <- lossPred(pred, y[t], pred, loss.type, loss.gradient)
        lexp[t,] <- lossPred(experts[t,], y[t],  pred, loss.type, loss.gradient)
        
        if (t > delay) {
          r.ETR <- awake[t-delay,] *  (lpred.ETR[t-delay] - lexp.ETR[t-delay,])
          r <- awake[t-delay,] * (lpred[t-delay] - lexp[t-delay,]) 
          R <- R + r - r.ETR
          eta <- 1 / (1/eta + r^2 - r.ETR^2)                                     
        }
      }      
      
      # If it is time to update...
      h <- (((t - 1) %% period) + 1)
      if (h == href) {
        Rop <- R
        etaop[t+1,] <- eta
      } else {
        etaop[t+1,] <- etaop[t,]
      }
    }    
    w = eta * pmax(R,0) / sum(eta * pmax(R,0))
    return(list(weights = weights, prediction = prediction, eta = etaop, last.weights = w))
  }
