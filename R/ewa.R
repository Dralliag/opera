ewa <-
function(y, experts, eta, awake = NULL, loss.type = 'square', 
                loss.gradient = TRUE, w0 = NULL, training = NULL)
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

  if (!is.null(training)) {
    w0 <- training$w0
    R <- training$R
    cumulativeLoss <- training$cumulativeLoss
  }

  for(t in 1:T){
    # Mise à jour du vecteur de poids de l'algo
    weights[t,] <- t(truncate1(exp(eta*R)) * t(awake[t,]))
    weights[t,] <- weights[t,] / sum(weights[t,])
    
    # Prévision et perte non loss.gradient de l'algo
    pred[t] <- experts[t,] %*% weights[t,]
    cumulativeLoss <- cumulativeLoss + loss(pred[t],y[t],loss.type)
    
    # Perte de l'algo et des experts (peut être loss.gradient)
    lpred <- lossPred(pred[t], y[t], pred[t], loss.type, loss.gradient)
    lexp <- lossPred(experts[t,], y[t], pred[t], loss.type, loss.gradient)
    
    # Mise à jour du vecteur de regret
    R <- R + awake[t,] * (lpred - lexp)
  }
  w = t(truncate1(exp(eta*R))) / sum(t(truncate1(exp(eta*R))))
  
 object <- list(model = "EWA", loss.type = loss.type, loss.gradient = loss.gradient,
    coefficients = w)

  object$parameters <- list(eta = eta)
  object$weights <- weights
  object$prediction <- pred

  object$training = list(
    R = R,
    w0 = w0,
    cumulativeLoss = cumulativeLoss)

  return(object)
}
