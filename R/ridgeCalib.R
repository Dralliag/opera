ridgeCalib <-
function(y, experts, gridlambda = 1, 
                       href = 1, period = 1, w0 = NULL, trace = F)
{
  experts <- as.matrix(experts)
  
  N <- ncol(experts)  # Nombre d'experts
  T <- nrow(experts)  # Nombre d'instants
  
  if (is.null(w0)) {w0 <- matrix(1/N,,ncol = N)} # Poids initial uniforme si non spécifié
  if (sum(is.na(experts)) > 0) {warning("There are not allowed NA's in expert advice")}

  # Grille du paramètre lambda et perte cumulée
  nlambda <- length(gridlambda)
  bestlambda <- floor(nlambda)/2 + 1 # Pour commencer on choisi le lambda au milieu de la grille
  gridlambda = matrix(gridlambda, nrow=nlambda)
  
  lambda <- rep(gridlambda[bestlambda],T)
  cumulatedloss <- rep(0,nlambda)
  
  wlambda <- array(w0, dim = c(N,nlambda))    # Matrice des poids donnés par chaque ewa(lambda) mis à jour chaque instant
  wlambdaop <- array(w0, dim = c(N,nlambda))  # Matrice des poids donnés par chaque ewa(lambda) avec mis à jour opérationnelle 

  weights <- matrix(0,ncol=N,nrow=T)    # Matrice des poids du mélange
  prediction <- rep(0,T)                # Prévisions du mélange

  At = diag(0,N)
  bt = matrix(t(w0)%*%gridlambda, nrow=N, ncol = nlambda)

  for(t in 1:T){
    # Affichage de l'avancement de l'aggregation rulee
    if (!(t %% floor(T/10)) && trace) cat(floor(10 * t/T)*10, '% -- ')
    
    # Poids, prévision et perte opérationnels
    weights[t,] <- wlambdaop[,bestlambda] 
    prediction[t] <- experts[t,] %*% weights[t,]
    lambda[t] <- gridlambda[bestlambda]
    
    # Prévisions et pertes opérationnelles pour chaque lambda de la grille
    predop <- experts[t,] %*% wlambdaop 
    lpredop <- (predop-y[t])^2
    cumulatedloss <- cumulatedloss + lpredop
    
    # Mise à jour
    At = At + experts[t,] %*% t(experts[t,])
    bt = bt + y[t] * experts[t,]
    
    # Si c'est l'heure de mise à jour on met à jour wop et le paramètre lambda et la grille
    h <- (((t-1)%%period)+1)
    if (h == href) {
      # Mise à jour de la grille
      bestlambda <- order(cumulatedloss)[1]

      # On aumente la grille si on est sur une extrémité
      if (bestlambda == nlambda) { 
        if (trace) cat(' + ')
        newlambda <- gridlambda[bestlambda] * 10^(1:3)
        gridlambda <- c(gridlambda, newlambda)
        nlambda <- nlambda + length(newlambda)
        for (k in 1:length(newlambda)) {
          perfnewlambda <- tryCatch(ridgeHour(y[1:t], matrix(experts[1:t,],ncol=N), newlambda[k],
                                   href = href, period = period, w0 = w0),
                                    error = function(e) {list(prediction = rep(0,t))})
           newcumulatedloss <- sum((perfnewlambda$prediction - y[1:t])^2)
          cumulatedloss <- c(cumulatedloss, newcumulatedloss)
        }
      }
      if (bestlambda == 1) {
        if (trace) cat(' - ')
        newlambda <- gridlambda[bestlambda] / 10^1
        nlambda <- nlambda + length(newlambda)
        bestlambda <- bestlambda + length(newlambda)
        for (k in 1:length(newlambda)) {
          gridlambda <- c(newlambda[k],gridlambda)
          perfnewlambda <- tryCatch(ridgeHour(y[1:t], matrix(experts[1:t,],ncol=N), newlambda[k],
                                   href = href, period = period, w0 = NULL),
                                    error = function(e) {list(prediction = rep(0,t))})
          newcumulatedloss <- sum((perfnewlambda$prediction - y[1:t])^2)
          cumulatedloss <- c(newcumulatedloss, cumulatedloss)
        }
      }
      wlambdaop <- matrix(0,nrow = N, ncol = nlambda)
      for (k in 1:nlambda) {
        wlambdaop[,k] = tryCatch(solve(gridlambda[k]*diag(1,N) + At,bt),
                                 error = function(e) {0})
      }
    }
  }
  l <- rmse(prediction,y)
  rmse <- sqrt(cumulatedloss / T)
  if (trace) cat('\n')
  return(list(weights = weights, prediction = prediction, 
              lambda = lambda, grid = gridlambda, 
              loss = l, gridloss = rmse))
}
