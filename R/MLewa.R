
MLewa <-
function(y, experts, awake = NULL, 
                  loss.type='squareloss', loss.gradient = TRUE,
                  period = 1, href = 1)
{
   experts <- as.matrix(experts)
   N <- ncol(experts)
   T <- nrow(experts)
   
   # Experts are always active if awake is unspecified
   if (is.null(awake)) {awake = matrix(1, nrow = T, ncol = N)} 
   awake <- as.matrix(awake)
   idx.na <- which(is.na(experts))
   awake[idx.na] <- 0
   experts[idx.na] <- 0

   R <- rep(0,N)
   
   # weight assigned to each expert
   weights <- matrix(0, ncol = N, nrow = T)
   prediction <- rep(0, T)
   
   w <- rep(1,N)/N
   wop <- rep(1,N)/N
   
   # Initialization or the learning parameter
   etaop <- matrix(exp(350),ncol=N,nrow=T+1)
   eta <- rep(exp(350),N)
   
   for(t in 1:T){
      # form the operational mixture and the prediction of the aggregation rule
      weights[t,] <- awake[t,] * wop / sum(awake[t,] * wop) 
      prediction[t] <- experts[t,] %*% weights[t,]
      
      # form the each-instant updated mixture and prediction
      p <- awake[t,] * w / sum(awake[t,] * w) 
      pred <- experts[t,] %*% p 
      
      # observe losses
      lpred <- lossPred(pred, y[t], pred, 
                         loss.type = loss.type, loss.gradient = loss.gradient)
      lexp <- lossPred(experts[t,], y[t], pred, 
                        loss.type = loss.type, loss.gradient = loss.gradient)
      
      # update regret and weights
      r <- awake[t,] * (lpred - lexp)
      R <- R + r
      eta = sqrt(log(N)/(log(N)/eta^2+r^2))
      w <- truncate1(exp(eta*R)) 
      
      # If it is time to update...
      h <- (((t - 1) %% period) + 1)
      if (h == href) {
         wop <- w
         etaop[t+1,] <- eta
      } else {
         etaop[t+1,] <- etaop[t,]
      }
   }
   return(list(weights = weights, prediction = prediction, eta = etaop))
}
