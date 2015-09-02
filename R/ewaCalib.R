
ewaCalib <-
function(y, experts, grideta = 1, awake = NULL,
        loss.type = 'squareloss', loss.gradient = TRUE, 
        w0 = NULL, trace = F, gamma = 2)
{
  experts <- as.matrix(experts)
  
  N <- ncol(experts)  # Number of experts
  T <- nrow(experts)  # Number of instants
  
  if (is.null(w0)) {w0 <- rep(1,N)} # Uniform intial weight vector if unspecified
  if (is.null(awake)) {awake = matrix(1, nrow = T, ncol = N)} # Full activation if unspecified

  awake <- as.matrix(awake)
  idx.na <- which(is.na(experts))
  awake[idx.na] <- 0
  experts[idx.na] <- 0
  

  neta <- length(grideta)         # Initial number of learning parameter in the grid to be optimized
  besteta <- floor(neta)/2 + 1    # We start with the parameter eta in the middle of the grid
  eta <- rep(grideta[besteta],T)  # Vector of calibrated learning rates (will be filled online by the algorithm)
  cumulativeLoss <- rep(0,neta)   # Cumulative loss suffered by each learning rate
  
  R <- array(0,c(N,neta))         # Matrix of regret suffered by each learning rate (columns) against each expert (rows)
  for (k in 1:neta) {R[,k] <- log(w0)/grideta[k]} # We initialize the regrets so that w0 is the initial weight vector
  
  weta <- array(w0, dim = c(N,neta))          # Weight matrix assigned by each algorithm EWA(eta) 
  weights <- matrix(0, ncol = N, nrow = T)    # Matrix of weights formed by the mixture
  prediction <- rep(0, T)                     # Predictions formed by the mixture
  
  pred.eta.nonop <- matrix(0, ncol = neta, nrow = T) # Prediction of mixture algorithm with different learning rates eta

  for(t in 1:T){
    # Display the state of progress of the algorithm
    if (!(t %% floor(T/10)) && trace) cat(floor(10 * t/T)*10, '% -- ')

    # Weights, prediction forme by EWA(eta[t]) where eta[t] is the learning rate calibrated online
    weights[t,] <- weta[,besteta] * awake[t,] / sum(weta[,besteta] * awake[t,])
    prediction[t] <- experts[t,] %*% weights[t,]
    eta[t] <- grideta[besteta]

    # Weights, predictions formed by each EWA(eta) for eta in the grid "grideta"
    pred.eta[t,] <- experts[t,] %*% t(t(weta * awake[t,]) / apply(weta * awake[t,],2,sum))
    lpred <- loss(pred.eta[t,], y[t], loss.type)
    cumulativeLoss <- cumulativeLoss + lpred
    
    # Regret update
    R <- R + awake[t,] * t(lpred.ETR - t(lexp.ETR))

    # Update of the best parameter
    besteta <- order(cumulativeLoss)[1]

    # We increase the size of the grid if the best parameter lies in an extremity
    if (besteta == neta) { 
      if (trace) cat(' + ')
        neweta <- grideta[besteta] * gamma^(1:3)
        grideta <- c(grideta, neweta)
        neta <- neta + length(neweta)
        R <- cbind(R, array(0, dim = c(N,length(neweta))))
        for (k in 1:length(neweta)) {
          perfneweta <- ewa(y[1:t], matrix(experts[1:t,],ncol=N), neweta[k], matrix(awake[1:t,],ncol=N), 
                                  loss.type = loss.type, loss.gradient = loss.gradient, w0 = w0)
          pred.eta <- cbind(pred.eta, c(perfneweta$prediction, rep(0, (T-t))))          
          cumulativeLoss <- c(cumulativeLoss, perfneweta$cumulativeLoss)
          R[,besteta+k] <- perfneweta$regret
        }
    }

    if (besteta == 1) {
      if (trace) cat(' - ')
      neweta <- grideta[besteta] / gamma^(1:3)
      neta <- neta + length(neweta)
      besteta <- besteta + length(neweta)
      R <- cbind(array(0, dim = c(N,length(neweta))), R)        
      for (k in 1:length(neweta)) {
        grideta <- c(neweta[k],grideta)
        perfneweta <- ewa(y[1:t], matrix(experts[1:t,],ncol=N), neweta[k], matrix(awake[1:t,],ncol=N), 
          loss.type = loss.type, loss.gradient = loss.gradient, w0 = w0)
        pred.eta <- cbind(c(perfneweta$prediction, rep(0, (T-t))), pred.eta)
        cumulativeLoss <- c(perfneweta$cumulativeLoss,cumulativeLoss)
        R[,besteta-k] <- perfneweta$regret
      }
    }

    weta <- truncate1(exp(t(t(matrix(R, ncol=neta)) * grideta)))
  }
  l <-  mean(loss(prediction, y, loss.type=loss.type))
  mloss <- cumulativeLoss / T
  if (loss.type == 'squareloss') {
    mloss <- sqrt(mloss)
    l <- sqrt(l)
  }
  if (trace) cat('\n')
  return(list(weights = weights, prediction = prediction, 
              eta = eta, grid = grideta, 
              loss = l, gridloss = mloss, last.weights = weta[,besteta] / sum(weta[,besteta] )))
}
