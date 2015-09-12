ewa <-
function(y, experts, eta, awake = NULL, loss.type = 'square', 
                loss.gradient = TRUE, w0 = NULL, tau = 0.5)
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
  
  R = log(w0)/eta     # Regret vector
  pred <- rep(0, T)    # Prediction vector
  cumulativeLoss <- 0  # Cumulative losses of the mixture
  weights <- matrix(0, ncol = N, nrow = T)    # Matrix of weights formed by the mixture

  for(t in 1:T){
    # Mise à jour du vecteur de poids de l'algo
    weights[t,] <- t(truncate1(exp(eta*R)) * t(awake[t,]))
    weights[t,] <- weights[t,] / sum(weights[t,])
    
    # Prévision et perte non loss.gradient de l'algo
    pred[t] <- experts[t,] %*% weights[t,]
    cumulativeLoss <- cumulativeLoss + loss(pred[t],y[t],loss.type,tau = tau)
    
    # Perte de l'algo et des experts (peut être loss.gradient)
    lpred <- lossPred(pred[t], y[t], pred[t], loss.type, loss.gradient, tau = tau)
    lexp <- lossPred(experts[t,], y[t], pred[t], loss.type, loss.gradient, tau = tau)
    
    # Mise à jour du vecteur de regret
    R <- R + awake[t,] * (lpred - lexp)
  }
  w = t(truncate1(exp(eta*R))) / sum(t(truncate1(exp(eta*R))))
  
  return(list(weights = weights, prediction = pred, 
              loss = cumulativeLoss/T, regret = R,
              weights.forecast =  w))
}
