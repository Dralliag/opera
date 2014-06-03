MLprod <-
function(y, experts, awake = NULL, 
                   loss.type = 'squareloss', loss.gradient = TRUE,
                   href = 1, period = 1, tau = 0.5)
{
   
   experts <- as.matrix(experts)
   N <- ncol(experts)
   T <- nrow(experts)
   
   if (is.null(awake)) {awake = matrix(1, nrow = T, ncol = N)} # Activation 1 si non spécifiée
   awake = as.matrix(awake)
   idx.na <- which(is.na(experts))
   awake[idx.na] <- 0
   experts[idx.na] <- 0

   weights <- matrix(0,ncol=N,nrow=T)
   R <- rep(0,N)
   Rop <- rep(0,N)
   
   L <- rep(1, N)
   maxloss <- 0
   gamma <- 1/2
   eta <- rep(exp(700),N)
   etaop <- matrix(exp(700),ncol=N,nrow=T+1)
   
   prediction <- rep(0,T)
   
   for(t in 1:T){
      # Operational update
      wop <- truncate1(exp(Rop))
      wop <- etaop[t,] * wop / sum(etaop[t,] * wop)
      weights[t,] <- awake[t,] * wop  / sum(awake[t,] * wop)
      prediction[t] <- experts[t,] %*% weights[t,]
      
      # Update weights
      w <- truncate1(exp(R));
      w <- eta * w / sum(eta * w)
      p <- awake[t,] * w  / sum(awake[t,] * w) 
      
      # Predict
      pred <- experts[t,] %*% p  
      
      # Observe losses
      lpred <- lossPred(pred,  y[t], pred, tau = tau,
                         loss.type = loss.type, loss.gradient = loss.gradient)
      lexp <- lossPred(experts[t,], y[t], pred, tau = tau,
                        loss.type = loss.type, loss.gradient = loss.gradient)
      
      r = awake[t,]*(lpred - lexp)
      L = L + r^2
      maxloss = pmax(maxloss,abs(r))
      neweta <- pmin(1/(2* maxloss),sqrt(log(N)/L))
      
      # Update regret and learning parameter
      R <- neweta/eta * R + log (1 + awake[t,] * neweta * (lpred - lexp))
      eta <- neweta
      
      if (is.na(sum(R)) ) {
         browser("Nan in R")
      }
      
      # Si c'est l'heure de mise à jour on met à jour weights
      h <- (((t-1)%%period)+1)
      if (h == href) {
         Rop <- R
         etaop[t+1,] <- eta
      } else {
         etaop[t+1,] <- etaop[t,]
      }
   }
   return(list(weights = weights, prediction = prediction, eta = etaop))
}
