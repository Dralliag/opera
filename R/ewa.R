ewa <-
function(y, experts, eta, awake = NULL, loss.type = 'squareloss', 
                loss.gradient = TRUE, w0 = NULL)
{
  experts <- as.matrix(experts)
  
  N <- ncol(experts)  # Nombre d'experts
  T <- nrow(experts)  # Nombre d'instants
  
  if (is.null(w0)) {w0 <- rep(1,N)} # Poids initial uniforme si non spécifié
  if (is.null(awake)) {awake = matrix(1, nrow = T, ncol = N)} # Activation 1 si non spécifiée
  awake = as.matrix(awake)
  idx.na <- which(is.na(experts))
  awake[idx.na] <- 0
  experts[idx.na] <- 0
  
  R = log(w0)/eta     # Regret des experts
  pred <- rep(0,T)    # Vecteur des prévisions du mélange
  cumulatedloss <- 0  # Perte cumulée de l'algo
  weights <- matrix(0,ncol=N,nrow=T)    # Matrice des poids du mélange

  for(t in 1:T){
    # Mise à jour du vecteur de poids de l'algo
    weights[t,] <- t(truncate1(exp(eta*R)) * t(awake[t,]))
    weights[t,] <- weights[t,] / sum(weights[t,])
    
    # Prévision et perte non loss.gradient de l'algo
    pred[t] <- experts[t,] %*% weights[t,]
    cumulatedloss <- cumulatedloss + loss(pred[t],y[t],loss.type)
    
    # Perte de l'algo et des experts (peut être loss.gradient)
    lpred <- lossPred(pred[t], y[t], pred[t], loss.type, loss.gradient)
    lexp <- lossPred(experts[t,], y[t], pred[t], loss.type, loss.gradient)
    
    # Mise à jour du vecteur de regret
    R <- R + awake[t,] * (lpred - lexp)
  }
  return(list(weights = weights, prediction = pred, cumulatedloss = cumulatedloss, regret = R))
}
