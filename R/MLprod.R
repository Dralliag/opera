MLprod <-
function(y, experts, awake = NULL, 
             loss.type = 'squareloss', loss.gradient = TRUE, tau = 0.5, w0 = NULL)
{
  
  experts <- as.matrix(experts)
  N <- ncol(experts)
  T <- nrow(experts)
  
  if (is.null(w0)) {w0 <- rep(1,N)} # Uniform intial weight vector if unspecified
  if (is.null(awake)) {awake = matrix(1, nrow = T, ncol = N)} # Activation 1 si non spécifiée
  awake = as.matrix(awake)
  idx.na <- which(is.na(experts))
  awake[idx.na] <- 0
  experts[idx.na] <- 0

  weights <- matrix(0,ncol=N,nrow=T)
  R <- log(w0)
  
  L <- rep(1, N)
  maxloss <- 0

  eta <- matrix(exp(700),ncol=N,nrow=T+1)
  prediction <- rep(0,T)
  
  for(t in 1:T){
    
    # Update weights
    w <- truncate1(exp(R))
    w <- eta[t,] * w / sum(eta[t,] * w)
    p <- awake[t,] * w  / sum(awake[t,] * w)
    pred <- experts[t,] %*% p       # Predict

    weights[t,] <- p
    prediction[t] <- pred
     
    # Observe losses
    lpred <- lossPred(pred,  y[t], pred, tau = tau,
                 loss.type = loss.type, loss.gradient = loss.gradient)
    lexp <- lossPred(experts[t,], y[t], pred, tau = tau,
                loss.type = loss.type, loss.gradient = loss.gradient)

    r = awake[t,] * (lpred - lexp)
    L = L + r^2
    maxloss = pmax(maxloss,abs(r))
    neweta <- pmin(1/(2* maxloss),sqrt(log(N)/L))
    
    # Update regret and learning parameter
    R <- neweta/eta[t,] * R + log (1 + awake[t,] * neweta * (lpred - lexp))
    eta[t+1,] <- neweta
    
    if (is.na(sum(R)) ) {
      browser("Nan in R")
    }
  }
  w <- truncate1(exp(R))
  w <- eta[T+1,] * w / sum(eta[T+1,] * w)
  return(list(weights = weights, prediction = prediction, eta = eta,
          weights.forecast = w))
}
