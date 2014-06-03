ewaCalib <-
function(y, experts, grideta = 1, awake = NULL,
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
  
  # Grille du paramètre eta et perte cumulée
  neta <- length(grideta)
  besteta <- floor(neta)/2 + 1 # Pour commencer on choisi le eta au milieu de la grille
  eta <- rep(grideta[besteta],T)
  cumulatedloss <- rep(0,neta)
  
  # Poids, regret et prévision
  R <- array(0,c(N,neta)) 
  for (k in 1:neta) {R[,k] <- log(w0)/grideta[k]}
  
  weta <- array(w0, dim = c(N,neta))    # Matrice des poids donnés par chaque ewa(eta) mis à jour chaque instant
  wetaop <- array(w0, dim = c(N,neta))  # Matrice des poids donnés par chaque ewa(eta) avec mis à jour opérationnelle 
  weights <- matrix(0,ncol=N,nrow=T)    # Matrice des poids du mélange
  prediction <- rep(0,T)                # Prévisions du mélange
    
  for(t in 1:T){
    # Affichage de l'avancement de l'aggregation rulee
    if (!(t %% floor(T/10)) && trace) cat(floor(10 * t/T)*10, '% -- ')
    # Poids, prévision et perte opérationnels
    weights[t,] <- wetaop[,besteta] * awake[t,] / sum(wetaop[,besteta] * awake[t,])
    prediction[t] <- experts[t,] %*% weights[t,]
    eta[t] <- grideta[besteta]
    # Prévisions et pertes opérationnelles pour chaque eta de la grille
    predop <- experts[t,] %*% t(t(wetaop * awake[t,]) / apply(wetaop * awake[t,],2,sum))
    lpredop <- loss(predop,y[t],loss.type)
    cumulatedloss <- cumulatedloss + lpredop
    
    # Prévisions et pertes pour chaque eta de la grille avec mise à jour à chaque instant
    pred <- experts[t,] %*% t(t(weta * awake[t,]) / apply(weta * awake[t,], 2, sum))
    lpred <- diag(lossPred(pred, y[t], pred, loss.type, loss.gradient))
    lexp <- lossPred(experts[t,], y[t], pred, loss.type, loss.gradient)
    
    # Mise à jour non opérationnelle
    R <- R + awake[t,] * t(lpred - t(lexp))
    weta <- truncate1(exp(t(t(matrix(R, ncol=neta)) * grideta)))
    
    # Si c'est l'heure de mise à jour on met à jour wop et le paramètre eta et la grille
    h <- (((t-1)%%period)+1)
    if (h == href) {
      # Mise à jour de la grille
      besteta <- order(cumulatedloss)[1]

      # On aumente la grille si on est sur une extrémité
      if (besteta == neta) { 
        if (trace) cat(' + ')
        neweta <- grideta[besteta] * 2^(1:3)
        grideta <- c(grideta, neweta)
        neta <- neta + length(neweta)
        R <- cbind(R, array(0, dim = c(N,length(neweta))))
        for (k in 1:length(neweta)) {
          perfneweta <- ewaHour(y[1:t], matrix(experts[1:t,],ncol=N), neweta[k], matrix(awake[1:t,],ncol=N),
                              href = href, period = period, 
                              loss.type = loss.type, loss.gradient = loss.gradient,
                              w0 = w0)
          cumulatedloss <- c(cumulatedloss,perfneweta$cumulatedloss)
          R[,besteta+k] <- perfneweta$regret
        }
      }
      if (besteta == 1) {
        if (trace) cat(' - ')
        neweta <- grideta[besteta] / 2^(1:3)
        neta <- neta + length(neweta)
        besteta <- besteta + length(neweta)
        R <- cbind(array(0, dim = c(N,length(neweta))), R)        
        for (k in 1:length(neweta)) {
          grideta <- c(neweta[k],grideta)
          perfneweta <- ewaHour(y[1:t], matrix(experts[1:t,],ncol=N), neweta[k], matrix(awake[1:t,],ncol=N),
                              href = href, period = period, 
                              loss.type = loss.type, loss.gradient = loss.gradient, w0 = w0)
          cumulatedloss <- c(perfneweta$cumulatedloss,cumulatedloss)
          R[,besteta-k] <- perfneweta$regret
        }
      }
      weta <- truncate1(exp(t(t(matrix(R, ncol=neta)) * grideta)))
      wetaop <- weta
    }
  }
  l <-  mean(loss(prediction,y,loss.type=loss.type))
  mloss <- cumulatedloss / T
  if (loss.type == 'squareloss') {
    mloss <- sqrt(mloss)
    l <- sqrt(l)
  }
  if (trace) cat('\n')
  return(list(weights = weights, prediction = prediction, 
              eta = eta, grid = grideta, 
              loss = l, gridloss = mloss))
}
