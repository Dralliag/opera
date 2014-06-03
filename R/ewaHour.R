ewaHour <-
function(y, experts, eta, awake = NULL, 
                 loss.type = 'squareloss', loss.gradient = TRUE, w0 = NULL, href = 1, period = 1)
{
  experts <- as.matrix(experts)
  
  N <- ncol(experts)      # Nombre d'experts
  T <- nrow(experts)      # Nombre d'instants
  
  if (is.null(w0)) {w0 <- rep(1,N)} # Poids initial uniforme si non spécifié
  if (is.null(awake)) {awake = matrix(1, nrow = T, ncol = N)} # Activation 1 si non spécifiée
  awake = as.matrix(awake)
  idx.na <- which(is.na(experts))
  awake[idx.na] <- 0
  experts[idx.na] <- 0
  
  R = log(w0)/eta       # Regret des experts
  weights <- matrix(0,ncol=N,nrow=T)    # Matrice des poids du mélange
  prediction <- rep(0,T)  # Prévisions de l'algo
  cumulatedloss <- 0      # Perte cumulée de l'algo
  
  w <- w0                 # Vecteur de poids mis à jour tous les instants
  wop <- w0               # Vecteur de poids mis à jour seulement à href
  
  for(t in 1:T){
    # Poids, prévision et perte opérationels
    weights[t,] <- wop * awake[t,] / sum(wop * awake[t,])
    prediction[t] <- experts[t,] %*% weights[t,]
    cumulatedloss <- cumulatedloss + loss(prediction[t],y[t],loss.type)
      
    # Prévision et pertes avec mise à jours tous les instants
    pred <- experts[t,] %*% ((w * awake[t,]) / sum(w * awake[t,]))
    lpred <- lossPred(pred, y[t], pred, loss.type, loss.gradient)
    lexp <- lossPred(experts[t,], y[t], pred, loss.type, loss.gradient)
    
    # Mise à jour des regrets et des poids
    R <- R + awake[t,] * (lpred - lexp)
    w <- truncate1(exp(eta * R))
    
    # Si c'est l'heure de mise à jour on met à jour wop
    h <- (((t-1)%%period)+1)
    if (h == href) {wop <- w}
  }
  return(list(weights = weights, prediction = prediction, cumulatedloss = cumulatedloss, regret = R))
}
