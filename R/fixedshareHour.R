fixedshareHour <-
function(y, experts, eta, alpha, awake = NULL, 
                         loss.type = 'squareloss', loss.gradient = TRUE, w0 = NULL, 
                         href = 1, period = 1)
{
   experts <- as.matrix(experts)
   
   N <- ncol(experts)  # Nombre d'experts
   T <- nrow(experts)  # Nombre d'instants
   
   if (is.null(w0)) {w0 <- rep(1,N)} # Poids initial uniforme si non spécifié
   if (is.null(awake)) {awake = matrix(1, nrow = T, ncol = N)} # Activation 1 si non spécifiée
   awake = as.matrix(awake)
   idx.na <- which(is.na(experts))
   awake[idx.na] <- 0
   experts[idx.na] <- 0

   prediction <- rep(0,T)    # Vecteur des prévisions du mélange
   cumulatedloss <- 0        # Perte cumulée de l'algo
   weights <- matrix(0,ncol=N,nrow=T)
   
   w <- w0             # Poids (hors sleeping) mis à jour chaque t
   wop <- w0           # Poids (hors sleeping) mis à jour seulement à href
   
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
      R <- log(w)/eta + awake[t,] * (lpred - lexp)    
      v <- truncate1(exp(eta*R))
      v <- v / sum(v)
      
      # Mise à jour des poids
      w <- alpha/N + (1-alpha) * v
      
      # Si c'est l'heure de mise à jour on met à jour wop
      h <- (((t-1)%%period)+1)
      if (h == href) {wop <- w} else {wop <- alpha/N * sum(wop) + (1-alpha) * wop}
   }
   # Renvoi de la matrice de poids 
   return(list(weights = weights, prediction = prediction, cumulatedloss = cumulatedloss, lastweight = w))
}
