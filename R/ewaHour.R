
ewaHour <-
function(y, experts, eta, awake = NULL, 
         loss.type = 'squareloss', loss.gradient = TRUE, 
         w0 = NULL, href = 1, period = 1,
         delay = 0, y.ETR = NULL)
{
  experts <- as.matrix(experts)
  
  N <- ncol(experts)      # Number of experts
  T <- nrow(experts)      # Number of instants
  
  if (is.null(w0)) {w0 <- rep(1,N)} # Poids initial uniforme si non spécifié
  if (is.null(awake)) {awake = matrix(1, nrow = T, ncol = N)} # Full activation if unspecified
  if (is.null(y.ETR) || (delay == 0)) {y.ETR = y}
  
  awake <- as.matrix(awake)
  idx.na <- which(is.na(experts))
  awake[idx.na] <- 0
  experts[idx.na] <- 0
  
  R = log(w0)/eta       # Regret vector
  weights <- matrix(0, ncol = N, nrow = T)    # Matrix of weights formed by the mixture
  prediction <- rep(0, T)  # Predictions of the mixture
  cumulativeLoss <- 0   # Cumulative loss of the mixture
  pred <- rep(0, T)     # Non operational prediction (update at each instant)

  w <- w0                 # Weight vector updated at each instant
  wop <- w0               # Weight vector updated only at href

  # Losses
  lpred <- numeric(T)     # Losses of the mixture
  lpred.ETR <- numeric(T) # Real time estimate of mixture losses
  lexp <- matrix(0, ncol = N, nrow = T)     # Expert losses
  lexp.ETR <- matrix(0, ncol = N, nrow = T) # Real time estimate of the expert losses
  
  for(t in 1:T){
    # operational weight vectors, predictions, and loss
    weights[t,] <- wop * awake[t,] / sum(wop * awake[t,])
    prediction[t] <- experts[t,] %*% weights[t,]
    cumulativeLoss <- cumulativeLoss + loss(prediction[t],y.ETR[t],loss.type)
    
    # Predictions and losses
    pred[t] <- experts[t,] %*% ((w * awake[t,]) / sum(w * awake[t,]))
    lpred.ETR[t] <- lossPred(pred[t], y.ETR[t], pred[t], loss.type, loss.gradient)
    lexp.ETR[t,] <- lossPred(experts[t,], y.ETR[t], pred[t], loss.type, loss.gradient)
    
    # Update regret
    R <- R + awake[t,] * (lpred.ETR[t] - lexp.ETR[t,])
    
    # On ajuste les regrets en fonction du signal réalisé
    if (delay > 0) {
      lpred[t] <- lossPred(pred[t], y[t], pred[t], loss.type, loss.gradient)
      lexp[t,] <- lossPred(experts[t,], y[t],  pred[t], loss.type, loss.gradient)
      
      if (t > delay) {
        R <- R  + awake[t-delay,] * ((lpred[t-delay] - lexp[t-delay,]) 
                                     - (lpred.ETR[t-delay] - lexp.ETR[t-delay,]))
        cumulativeLoss <- cumulativeLoss + loss(prediction[t-delay],y[t-delay],loss.type) - loss(prediction[t-delay],y.ETR[t-delay],loss.type)
      }
    }
    
    w <- truncate1(exp(eta * R))
    
    # If it is time to update
    h <- (((t - 1) %% period) + 1)
    if (h == href) {wop <- w}
  }
  return(list(weights = weights, prediction = prediction, cumulativeLoss = cumulativeLoss, regret = R, prediction.non.operational = pred))
}
