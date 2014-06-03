pinballHour <-
function(y, experts, lambda, href=1, period = 1, w0 = NULL, tau = 0.5) {
   experts = as.matrix(experts)
   N <- ncol(experts)
   T <- nrow(experts)
   if (is.null(w0)) {w0 <- rep(1/N,N)} # Poids initial uniforme si non spécifié
   if (sum(is.na(experts)) > 0) {warning("There are not allowed NA's in expert advice")}

   
   pinballloss = function(beta, Y, X) {
      beta = as.matrix(beta, ncol = 1)
      X = as.matrix(X)
      pred = X%*%beta
      return(mean((tau - (Y > pred)) * (pred - Y)) + lambda * sum(beta^2))
   }
   
   if (is.null(w0)) {w0 <- rep(1/N,N)} # Poids initial uniforme si non spécifié
   
   weights <- matrix(0,ncol=N,nrow=T)
   w = w0
   for (t in 1:T){  
      weights[t,] = w
      h <- (((t-1)%%period)+1)
      if (h == href) {
         cat('Iteration :',t,'\n')
         w = optim(w, fn=function(beta) {pinballloss(beta, y[1:t],experts[1:t,])})$par
         }
   }
   prediction <- apply(experts * weights, 1, sum)
   return(list(weights = weights, prediction = prediction))
}
