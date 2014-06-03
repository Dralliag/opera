fixedshareCalib <-
function(y, experts, grideta = 1, gridalpha = 10^(-4:-1), awake = NULL,
                            loss.type = 'squareloss', loss.gradient = TRUE, w0 = NULL, href = 1, period = 1, trace = F)
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
   
   # Grille du paramètre eta et alpha, on commence par le milieu de la grille
   neta <- length(grideta)
   nalpha <-length(gridalpha)
   bestpar <- c(floor(neta)/2,floor(nalpha)/2)+1
   par <- NULL
   
   cumulatedloss <- array(0, dim = c(neta, nalpha))
   
   # Poids et regret
   wpar <- array(w0, dim = c(N,neta,nalpha))
   wparop <- array(w0, dim = c(N,neta,nalpha))
   
   weights <- matrix(0,ncol=N,nrow=T)    # Matrice des poids du mélange
   prediction <- rep(0,T)                # Prévisions du mélange
   
   test = rep(0,T)
   
   for(t in 1:T){
      # Affichage de l'avancement de l'aggregation rulee
      if (!(t %% floor(T/10)) && trace) cat(floor(10 * t/T)*10, '% -- ')
      
      # Poids prévisions et pertes opérationnels
      weights[t,] <- wparop[,bestpar[1],bestpar[2]] * awake[t,] / sum(wparop[,bestpar[1],bestpar[2]] * awake[t,])
      prediction[t] <- experts[t,] %*% weights[t,]
      par <- rbind(par, data.frame(eta = grideta[bestpar[1]],
                                   alpha = gridalpha[bestpar[2]]))
      
      # boucle sur les paramètres alpha
      for (k in 1:nalpha) {
         
         # Prévision et pertes opérationnelles pour chaque (eta,alpha) de la grille
         waux <-  t(t(wparop[,,k] * awake[t,]) / apply(as.matrix(wparop[,,k]*awake[t,]), 2, sum))
         predop <- experts[t,] %*% waux
         lpredop <- loss(predop, y[t], loss.type)
         cumulatedloss[,k] <- cumulatedloss[,k] + lpredop
         
         # Prévision et pertes non opérationnelles des algo et des experts
         waux <- t(t(awake[t,] * wpar[,,k]) / apply(as.matrix(awake[t,] * wpar[,,k]), 2, sum))
         pred <- experts[t,] %*% waux
         lpred <- diag(lossPred(pred, y[t], pred, loss.type, loss.gradient))
         lexp <- lossPred(experts[t,], y[t], pred, loss.type, loss.gradient)
         
         # Mise à jour des regrets et des poids
         R <- t(t(log(wpar[,,k]))/grideta) + awake[t,] * t(lpred - t(lexp))
         v <- truncate1(exp(t(t(matrix(R, ncol=neta)) * grideta)))
         v <- t(t(v) / apply(v,2,sum)) # on renormalise chaque colonne
         
         # Mise à jour des poids
         wpar[,,k] <- gridalpha[k]/N + (1-gridalpha[k]) * v
      }
      
      
      # Si c'est l'heure de mise à jour...
      h <- (((t-1)%%period)+1)
      if (h == href) {
         # Mise à jour de la grille
         bestpar <- which(cumulatedloss == min(cumulatedloss), arr.ind = TRUE)[1,]
         # On augmente la grille si on est à une extrémité
         if (bestpar[1] == neta) {
            if (trace) cat(' + ')
            neweta <- grideta[neta] * 2^(1:3)
            grideta <- c(grideta, neweta)
            neta <- neta + length(neweta)
            wparaux <- wpar
            wpar <- array(dim = c(N,neta,nalpha))
            wpar[,1:bestpar[1],] <- wparaux
            for (j in 1:length(neweta)) {
               cumulatedloss <- rbind(cumulatedloss, 0)
               for (k in 1:nalpha) {
                  perfnewpar <- fixedshareHour(y[1:t], matrix(experts[1:t,],ncol=N), neweta[j], gridalpha[k], 
                                             awake =  matrix(awake[1:t,],ncol=N),
                                             loss.type = loss.type, loss.gradient = loss.gradient, 
                                             href = href, period = period, w0 = w0)
                  cumulatedloss[bestpar[1]+j,k] <- perfnewpar$cumulatedloss
                  wpar[,bestpar[1]+j,k] <- perfnewpar$lastweight
               }
            }
         }
         if (bestpar[1] == 1) {
            if (trace) cat(' - ')
            neweta <- grideta[1] / 2^(1:3)
            neta <- neta + length(neweta)
            bestpar[1] <- bestpar[1] + length(neweta)
            wparaux <- wpar
            wpar <- array(dim = c(N,neta,nalpha))
            wpar[,bestpar[1]:neta,] <- wparaux
            for (j in 1:length(neweta)) {
               grideta <- c(neweta[j], grideta)
               cumulatedloss <- rbind(0, cumulatedloss)
               for (k in 1:nalpha) {
                  perfnewpar <- fixedshareHour(y[1:t], matrix(experts[1:t,],ncol=N),
                                             neweta[j], gridalpha[k],
                                             awake = matrix(awake[1:t,],ncol=N),
                                             loss.type = loss.type, loss.gradient = loss.gradient,
                                             href = href, period = period, w0 = w0)
                  cumulatedloss[1,k] <- perfnewpar$cumulatedloss
                  wpar[,bestpar[1]-j,k] <- perfnewpar$lastweight
               }
            }
         }
         wparop <- wpar
      } else {
         for (k in 1:nalpha) {
            wparop[,,k] <- gridalpha[k]/N + (1-gridalpha[k]) * wparop[,,k]
         }
      }
   }
  l <-  mean(loss(prediction,y,loss.type=loss.type))
   mloss <- cumulatedloss / T
   if (loss.type == 'squareloss') {
      mloss <- sqrt(mloss)
      l <- sqrt(l)
   }
   rownames(mloss) = grideta
   colnames(mloss) = gridalpha
   
   if (trace) cat('\n')
   return(list(weights = weights, prediction = prediction, 
               par = par, grideta = grideta, gridalpha = gridalpha,
               loss = l, gridloss = mloss))
}
