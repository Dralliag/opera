
fixedshareHour <-
  function(y, experts, eta, alpha, awake = NULL, 
           loss.type = 'squareloss', loss.gradient = TRUE, w0 = NULL, 
           href = 1, period = 1,
           delay = 0, y.ETR = NULL)
  {
    experts <- as.matrix(experts)
    
    N <- ncol(experts)  # Number of experts
    T <- nrow(experts)  # Number of instants
    
    if (is.null(w0)) {w0 <- rep(1,N)} # Uniform intial weight vector if unspecified
    if (is.null(awake)) {awake = matrix(1, nrow = T, ncol = N)} # Full activation if unspecified
    if (is.null(y.ETR) || (delay == 0)) {y.ETR = y}

    awake <- as.matrix(awake)
    idx.na <- which(is.na(experts))
    awake[idx.na] <- 0
    experts[idx.na] <- 0
    
    prediction <- rep(0, T)    # Prediction vector
    cumulativeLoss <- 0        # Cumulative losses of the mixture
    weights <- matrix(0, ncol = N, nrow = T)
    
    w <- w0             # Poids (hors sleeping) mis à jour chaque t
    wop <- w0           # Poids (hors sleeping) mis à jour seulement à href
    

    # Losses
    lpred <- numeric(T)     # Losses of the mixture
    lpred.ETR <- numeric(T) # Real time estimate of mixture losses
    lexp <- matrix(0, ncol = N, nrow = T)     # Expert losses
    lexp.ETR <- matrix(0, ncol = N, nrow = T) # Real time estimate of the expert losses
  
    for(t in 1:T){
      # operational weight vectors, predictions, and losses
      weights[t,] <- wop * awake[t,] / sum(wop * awake[t,])
      prediction[t] <- experts[t,] %*% weights[t,]
      cumulativeLoss <- cumulativeLoss + loss(prediction[t],y[t],loss.type)
      
      # Non operational predictions and losses (update every t)
      pred <- experts[t,] %*% ((w * awake[t,]) / sum(w * awake[t,]))
      lpred.ETR[t] <- lossPred(pred, y.ETR[t], pred, loss.type, loss.gradient)
      lexp.ETR[t,] <- lossPred(experts[t,], y.ETR[t], pred, loss.type, loss.gradient)
      
      # Regret update (translation invariant)
      R <- log(w)/eta + awake[t,] * (lpred.ETR[t] - lexp.ETR[t,])    

      # On ajuste les regrets en fonction du signal réalisé
      if (delay > 0) {
        lpred[t] <- lossPred(pred, y[t], pred, loss.type, loss.gradient)
        lexp[t,] <- lossPred(experts[t,], y[t],  pred, loss.type, loss.gradient)
        
        if (t > delay) {
          R <- R  + awake[t-delay,] * ((lpred[t-delay] - lexp[t-delay,]) - (lpred.ETR[t-delay] - lexp.ETR[t-delay,]))
        }
      }
      # Penser à ajouter un oubli !!! (C'est compliqué !!!)
      
      v <- truncate1(exp(eta*R))
      v <- v / sum(v)
      
      # Mise à jour des poids
      w <- alpha/N + (1-alpha) * v
      
      # Si c'est l'heure de mise à jour on met à jour wop
      h <- (((t - 1) %% period) + 1)
      if (h == href) {wop <- w} else {wop <- alpha/N * sum(wop) + (1-alpha) * wop}
    }
    # Renvoi de la matrice de poids 
    return(list(weights = weights, prediction = prediction, cumulativeLoss = cumulativeLoss, lastweight = w))
  }
