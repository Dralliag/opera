
MLewa <-
function(y, experts, awake = NULL, 
                  loss.type='squareloss', loss.gradient = TRUE, w0 = NULL)
{
   experts <- as.matrix(experts)
   N <- ncol(experts)
   T <- nrow(experts)
   
   if (is.null(w0)) {w0 <- rep(1,N)} # Uniform intial weight vector if unspecified

   # Experts are always active if awake is unspecified
   if (is.null(awake)) {awake = matrix(1, nrow = T, ncol = N)} 
   awake <- as.matrix(awake)
   idx.na <- which(is.na(experts))
   awake[idx.na] <- 0
   experts[idx.na] <- 0

   R = rep(0,N)     # Regret vector
   w <- w0
   
   # weight assigned to each expert
   weights <- matrix(0, ncol = N, nrow = T)
   prediction <- rep(0, T)
   
   # Initialization or the learning parameter
   eta <- matrix(exp(350),ncol=N,nrow=T+1)
   for(t in 1:T){
      # form the each-instant updated mixture and prediction
      p <- awake[t,] * w / sum(awake[t,] * w) 
      pred <- experts[t,] %*% p     

      # form the operational mixture and the prediction of the aggregation rule
      weights[t,] <- p
      prediction[t] <- pred

      # observe losses
      lpred <- lossPred(pred, y[t], pred, 
                         loss.type = loss.type, loss.gradient = loss.gradient)
      lexp <- lossPred(experts[t,], y[t], pred, 
                        loss.type = loss.type, loss.gradient = loss.gradient)
      
      # update regret and weights
      r <- awake[t,] * (lpred - lexp)
      R <- R + r
      eta[t+1,] = sqrt(log(N)/(log(N)/eta[t,]^2+r^2))
      w <- truncate1(exp(log(w0) + eta[t+1,]*R)) 
   }
   w = w / sum(w)

   return(list(weights = weights, prediction = prediction, eta = eta,
              weights.forecast =  w))
}
