
MLprod <-
function(y, experts, awake = NULL,
                  loss.type = 'squareloss', loss.gradient = TRUE,
                  period = 1, href = 1,
                  delay = 0, y.ETR = NULL,
                  tau = 0.5) 
{
   
   experts <- as.matrix(experts)
   N <- ncol(experts)
   T <- nrow(experts)
   
   if (is.null(awake)) {awake = matrix(1, nrow = T, ncol = N)} # Full activation if unspecified
   if (is.null(y.ETR) || (delay == 0)) {y.ETR = y}

   awake <- as.matrix(awake)
   idx.na <- which(is.na(experts))
   awake[idx.na] <- 0
   experts[idx.na] <- 0

   weights <- matrix(0, ncol = N, nrow = T)
   prediction <- rep(0, T)

   R <- rep(0,N)
   Rop <- rep(0,N)
   
   L <- rep(1, N)
   maxloss <- 0
   gamma <- 1/2
   eta <- rep(exp(700),N)
   etaop <- matrix(exp(700), ncol = N, nrow = T + 1)
   
   # Losses
   lpred <- numeric(T)     # Losses of the mixture
   lpred.ETR <- numeric(T) # Real time estimate of mixture losses
   lexp <- matrix(0, ncol = N, nrow = T)     # Expert losses
   lexp.ETR <- matrix(0, ncol = N, nrow = T) # Real time estimate of the expert losses

   for(t in 1:T){
      # Operational update
      wop <- truncate1(exp(Rop))
      wop <- etaop[t,] * wop / sum(etaop[t,] * wop)

      weights[t,] <- awake[t,] * wop  / sum(awake[t,] * wop)
      prediction[t] <- experts[t,] %*% weights[t,]
      
      # Weight update
      w <- truncate1(exp(R));
      w <- eta * w / sum(eta * w)
      p <- awake[t,] * w  / sum(awake[t,] * w) 
      
      # Predict
      pred <- experts[t,] %*% p  
      
      # Observe losses
      lpred.ETR[t] <- lossPred(pred, y.ETR[t], pred, tau = tau,
                         loss.type = loss.type, loss.gradient = loss.gradient)
      lexp.ETR[t,] <- lossPred(experts[t,], y.ETR[t], pred, tau = tau,
                        loss.type = loss.type, loss.gradient = loss.gradient)
      
      
      r <- awake[t,] * (lpred.ETR[t] - lexp.ETR[t,])
      L <- L + r^2
      maxloss <- pmax(maxloss, abs(r))
      neweta <- pmin(1/(2* maxloss), sqrt(log(N)/L))
      
      # Update regret and learning parameter
      R <- neweta/eta * R + log (1 + neweta * r)
      eta <- neweta
      
      # We readjust by using the true observations
      if (delay > 0) {
        lpred[t] <- lossPred(pred, y[t], pred, loss.type, loss.gradient)
        lexp[t,] <- lossPred(experts[t,], y[t],  pred, loss.type, loss.gradient)
         
        if (t > delay) {
          r.ETR <- awake[t-delay,] *  (lpred.ETR[t-delay] - lexp.ETR[t-delay,])
          r <- awake[t-delay,] * (lpred[t-delay] - lexp[t-delay,]) 
          L <- L + r^2 - r.ETR^2
          maxloss <- pmax(maxloss, abs(r))
          neweta <- pmin(1 / (2 * maxloss), sqrt(log(N)/L))
          R <- neweta/eta * R + log (1 + neweta * r) - log(1 + neweta * r.ETR)
          eta <- neweta
        }
      }      

      if (is.na(sum(R)) ) {
         browser("Nan in R")
      }
      
      # Si c'est l'heure de mise à jour on met à jour weights
      h <- (((t - 1) %% period) + 1)
      if (h == href) {
         Rop <- R
         etaop[t+1,] <- eta
      } else {
         etaop[t+1,] <- etaop[t,]
      }
   }
   return(list(weights = weights, prediction = prediction, eta = etaop))
}
