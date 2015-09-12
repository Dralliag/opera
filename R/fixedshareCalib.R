fixedshareCalib <-
  function(y, experts, grid.eta = 1, grid.alpha = 10^(-4:-1), awake = NULL,
           loss.type = 'square', loss.gradient = TRUE, w0 = NULL, trace = F, gamma = 2,
           tau = 0.5)
  {
    experts <- as.matrix(experts)
    
    N <- ncol(experts)  # Number of experts
    T <- nrow(experts)  # Number of instants
    
    if (is.null(w0)) {w0 <- rep(1,N)} # Uniform intial weight vector if unspecified
    if (is.null(awake)) {awake = matrix(1, nrow = T, ncol = N)} # Full activation if unspecified
    
    awake = as.matrix(awake) 
    idx.na <- which(is.na(experts))
    awake[idx.na] <- 0
    experts[idx.na] <- 0
    
    neta <- length(grid.eta)       # Number of learning rates in the grid
    nalpha <-length(grid.alpha)    # Number of mixing rates in the grid
    bestpar <- c(floor(neta)/2,floor(nalpha)/2)+1 # We start with the parameters in the middle of the grids
    par <- NULL
    
    cumulativeLoss <- array(0, dim = c(neta, nalpha))
    
    
    wpar <- array(w0, dim = c(N,neta,nalpha))   # Weight array (in 3 dimensions) assigned by each algorithm FixedShare(eta,alpha) 
    weights <- matrix(0, ncol = N, nrow = T)    # Matrix of weights formed by the mixture
    prediction <- rep(0, T)                     # Predictions formed by the mixture
    
    for(t in 1:T){
      # Display the state of progress of the algorithm
      if (!(t %% floor(T/10)) && trace) cat(floor(10 * t/T)*10, '% -- ')
      
      # Weights, prediction forme by FixedShare(eta[t],alpha[t]) where par[t,] = c(eta[t],alpha[t]) are the parameters calibrated online
      weights[t,] <- wpar[,bestpar[1],bestpar[2]] * awake[t,] / sum(wpar[,bestpar[1],bestpar[2]] * awake[t,])
      prediction[t] <- experts[t,] %*% weights[t,]
      par <- rbind(par, data.frame(eta = grid.eta[bestpar[1]],
                                   alpha = grid.alpha[bestpar[2]]))
      
      # Loop over the mixing rates alpha in the grid "grid.alpha"
      for (k in 1:nalpha) {
        
        # Weights, prediction, and losses formed by FixedShare(eta,alpha) for (eta,alpha) in the grid
        waux <-  t(t(wpar[,,k] * awake[t,]) / apply(as.matrix(wpar[,,k]*awake[t,]), 2, sum))
        pred <- experts[t,] %*% waux
        cumulativeLoss[,k] <- cumulativeLoss[,k] + loss(pred, y[t], loss.type,tau=tau) # Non gradient cumulative losses
        lpred <- diag(lossPred(pred, y[t], pred, loss.type, loss.gradient,tau=tau))
        lexp <- lossPred(experts[t,], y[t], pred, loss.type, loss.gradient,tau=tau)
        
        # Regret update
        R <- t(t(log(wpar[,,k]))/grid.eta) + awake[t,] * t(lpred - t(lexp))
        
        # Weight update
        v <- truncate1(exp(t(t(matrix(R, ncol=neta)) * grid.eta)))
        v <- t(t(v) / apply(v,2,sum)) # Renormalization of each column
        wpar[,,k] <- grid.alpha[k]/N + (1-grid.alpha[k]) * v
      }
      
      # Grid update
      # ***************
      # find the best parameter in the grid
      bestpar <- which(cumulativeLoss == min(cumulativeLoss), arr.ind = TRUE)[1,]  
      
      # Expand the grid if the best parameter lies on an extremity (only for the first component)
      if (bestpar[1] == neta) {
        if (trace) cat(' + ')
        neweta <- grid.eta[neta] * gamma^(1:3)
        grid.eta <- c(grid.eta, neweta)
        neta <- neta + length(neweta)
        wparaux <- wpar
        wpar <- array(dim = c(N,neta,nalpha))
        wpar[,1:bestpar[1],] <- wparaux
        for (j in 1:length(neweta)) {
          cumulativeLoss <- rbind(cumulativeLoss, 0)
          for (k in 1:nalpha) {
            perfnewpar <- fixedshare(y[1:t], matrix(experts[1:t,],ncol=N), neweta[j], grid.alpha[k], 
                                     awake =  matrix(awake[1:t,],ncol=N),
                                     loss.type = loss.type, loss.gradient = loss.gradient, 
                                     w0 = w0, tau = tau)
            cumulativeLoss[bestpar[1]+j,k] <- perfnewpar$loss * t
            wpar[,bestpar[1]+j,k] <- perfnewpar$weights.forecast
          }
        }
      }
      if (bestpar[1] == 1) {
        if (trace) cat(' - ')
        neweta <- grid.eta[1] / gamma^(1:3)
        neta <- neta + length(neweta)
        bestpar[1] <- bestpar[1] + length(neweta)
        wparaux <- wpar
        wpar <- array(dim = c(N, neta, nalpha))
        wpar[,bestpar[1]:neta,] <- wparaux
        for (j in 1:length(neweta)) {
          grid.eta <- c(neweta[j], grid.eta)
          cumulativeLoss <- rbind(0, cumulativeLoss)
          for (k in 1:nalpha) {
            perfnewpar <- fixedshare(y[1:t], matrix(experts[1:t,],ncol=N),
                                     neweta[j], grid.alpha[k],
                                     awake = matrix(awake[1:t,],ncol=N),
                                     loss.type = loss.type, loss.gradient = loss.gradient, 
                                     w0 = w0, tau = tau)
            cumulativeLoss[1, k] <- perfnewpar$loss * t
            wpar[, bestpar[1]-j, k] <- perfnewpar$weights.forecast
          }
        }
      }
    }
    
    # Next weights
    w <- wpar[,bestpar[1],bestpar[2]]  / sum(wpar[,bestpar[1],bestpar[2]])
    
    # Losses
    l <-  mean(loss(prediction, y ,loss.type = loss.type, tau = tau)) # Average loss suffered by the algorithm
    mloss <- cumulativeLoss / T         # Average loss of each learning rate on the grid
    

    rownames(mloss) = grid.eta
    colnames(mloss) = grid.alpha
    res = list(weights = weights, prediction = prediction,
                eta = par$eta, grid.eta = grid.eta, 
                alpha = par$alpha, grid.alpha = grid.alpha,
                loss = l, grid.loss = mloss, weights.forecast = w)
    if (loss.type == 'square') {
      res$grid.rmse <- sqrt(mloss)
      res$rmse <- sqrt(l)
    }
    if (trace) cat('\n')
    return(res)
  }
