
ewaCalib <-
function(y, experts, grideta = 1, awake = NULL,
        loss.type = 'squareloss', loss.gradient = TRUE, 
        w0 = NULL, href = 1, period = 1, 
        delay = 0, y.ETR = NULL, trace = F)
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
  
  # Grille du paramètre eta et perte cumulée
  neta <- length(grideta)
  besteta <- floor(neta)/2 + 1 # We start with the parameter eta in the middle of the grid
  eta <- rep(grideta[besteta],T)
  cumulativeLoss <- rep(0,neta)
  
  # Poids, regret et prévision
  R <- array(0,c(N,neta)) 
  for (k in 1:neta) {R[,k] <- log(w0)/grideta[k]}
  
  weta <- array(w0, dim = c(N,neta))    # Matrice des poids donnés par chaque ewa(eta) mis à jour chaque instant
  wetaop <- array(w0, dim = c(N,neta))  # Matrice des poids donnés par chaque ewa(eta) avec mise à jour opérationnelle 
  weights <- matrix(0, ncol = N, nrow = T)    # Matrix of weights formed by the mixture
  prediction <- rep(0, T)                # Prévisions du mélange
  
  pred.eta.op <- matrix(0, ncol = neta, nrow = T) # Prediction of mixture algorithm with different learning rates eta
  pred.eta.nonop <- matrix(0, ncol = neta, nrow = T)

  for(t in 1:T){
    # Affichage de l'avancement de l'aggregation rulee
    if (!(t %% floor(T/10)) && trace) cat(floor(10 * t/T)*10, '% -- ')

    # Poids, prévision et perte opérationnels
    weights[t,] <- wetaop[,besteta] * awake[t,] / sum(wetaop[,besteta] * awake[t,])
    prediction[t] <- experts[t,] %*% weights[t,]
    eta[t] <- grideta[besteta]

    # Prévisions et pertes opérationnelles pour chaque eta de la grille
    pred.eta.op[t,] <- experts[t,] %*% t(t(wetaop * awake[t,]) / apply(wetaop * awake[t,],2,sum))
    lpredop <- loss(pred.eta.op[t,], y.ETR[t], loss.type)
    cumulativeLoss <- cumulativeLoss + lpredop
    
    # Predictions and losses of each parameter of the grid with update at each instant with update at each instant
    pred.eta.nonop[t,] <- experts[t,] %*% t(t(weta * awake[t,]) / apply(weta * awake[t,], 2, sum))
    pred <- matrix(pred.eta.nonop[t,], nrow = 1)
    lpred.ETR <- diag(lossPred(pred, y.ETR[t], pred, loss.type, loss.gradient))
    lexp.ETR <- lossPred(experts[t,], y.ETR[t], pred, loss.type, loss.gradient)
    
    # Non operational regret update
    R <- R + awake[t,] * t(lpred.ETR - t(lexp.ETR))

    if ( (delay > 0) && (t > delay) ) {
      # Regret update according to the difference between y.ETR[t-delay] and y[t-delay]
      pred <- matrix(pred.eta.nonop[t-delay,], nrow = 1)
      lpred.ETR.delay <- diag(lossPred(pred, y.ETR[t-delay], pred, loss.type, loss.gradient))
      lexp.ETR.delay <- lossPred(experts[t-delay,], y.ETR[t-delay], pred, loss.type, loss.gradient)
      lpred.delay <- diag(lossPred(pred, y[t-delay], pred, loss.type, loss.gradient))
      lexp.delay <- lossPred(experts[t-delay,], y[t-delay], pred, loss.type, loss.gradient)

      R <- R  + awake[t-delay,] * ((lpred.delay - t(lexp.delay)) 
                                   - (lpred.ETR.delay - t(lexp.ETR.delay)))

      # Update of the accumulated loss of each EWA(eta) according to the difference between y.ETR[t-delay] and y[t-delay]
      cumulativeLoss <- cumulativeLoss + loss(pred.eta.op[t-delay,], y[t-delay], loss.type) - loss(pred.eta.op[t-delay,], y.ETR[t-delay], loss.type)
    }

    # If it is time, Update of wop, the learning rate and the grid
    h <- (((t - 1) %% period) + 1)
    if (h == href) {
      # Update of the best parameter
      besteta <- order(cumulativeLoss)[1]

      # We increase the size of the grid if the best parameter lies in an extremity
      if (besteta == neta) { 
        if (trace) cat(' + ')
        neweta <- grideta[besteta] * 2^(1:3)
        grideta <- c(grideta, neweta)
        neta <- neta + length(neweta)
        R <- cbind(R, array(0, dim = c(N,length(neweta))))
        for (k in 1:length(neweta)) {
          perfneweta <- ewaHour(y[1:t], matrix(experts[1:t,],ncol=N), neweta[k], matrix(awake[1:t,],ncol=N), href = href, period = period, loss.type = loss.type, loss.gradient = loss.gradient, w0 = w0, delay = delay, y.ETR = y.ETR[1:t])
          pred.eta.op <- cbind(pred.eta.op, c(perfneweta$prediction, rep(0, (T-t))))
          pred.eta.nonop <- cbind(pred.eta.nonop, c(perfneweta$prediction.non.operational, rep(0, (T-t))))
          cumulativeLoss <- c(cumulativeLoss, perfneweta$cumulativeLoss)
          R[,besteta+k] <- perfneweta$regret
        }
      }

      if (besteta == 1) {
        if (trace) cat(' - ')
        neweta <- grideta[besteta] / 2^(1:3)
        neta <- neta + length(neweta)
        besteta <- besteta + length(neweta)
        R <- cbind(array(0, dim = c(N,length(neweta))), R)        
        for (k in 1:length(neweta)) {
          grideta <- c(neweta[k],grideta)
          perfneweta <- ewaHour(y[1:t], matrix(experts[1:t,],ncol=N), neweta[k], matrix(awake[1:t,],ncol=N), href = href, period = period, loss.type = loss.type, loss.gradient = loss.gradient, w0 = w0, delay = delay, y.ETR = y.ETR[1:t])
          pred.eta.op <- cbind(c(perfneweta$prediction, rep(0, (T-t))), pred.eta.op)
          pred.eta.nonop <- cbind(c(perfneweta$prediction.non.operational, rep(0, (T-t))), pred.eta.nonop)

          cumulativeLoss <- c(perfneweta$cumulativeLoss,cumulativeLoss)
          R[,besteta-k] <- perfneweta$regret
        }
      }

      weta <- truncate1(exp(t(t(matrix(R, ncol=neta)) * grideta)))
      wetaop <- weta
    }
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
              loss = l, gridloss = mloss, last.weights = weta[,besteta] / sum(wetaop[,besteta] )))
}
