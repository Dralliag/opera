simplexProj <- function(u, U = 1) {
  N <- length(u)
  v = sort(u,decreasing = TRUE)
  rho = max(which((v + (U-cumsum(v)) / (1:N)) >0))
  if (rho >= 1) {
    lambda = 1/rho * (U - sum(v[1:rho]))
    return(pmax(u+ lambda,0))
  } else {
    return(rep(1/N,N))
  }
}

norme1 = function(x){sum(abs(x))}

projectionL1 = function(x,U=1) {
  if (norme1(x)<=U){
    return(x)
  } else {
    return(sign(x)*simplexProj(abs(x),U))
  }
}


# algorithme 1 du papier pour la perte carrée
SarseAccEG = function(Y,X,U,d0,alpha,B,delta) {
  n <- length(Y)
  d <- ncol(X)
  
  epsilon <-  16 * d0 * B * (2 * sqrt(2*log(d)) + sqrt(2*log(n/delta))) / alpha
  
  i.start <- ceiling(2*log(2*epsilon/U,base = 2))
  i.end <- floor(log(n,base = 2))
  if (i.start < i.end) {
    i.seq <- c(0,i.start:i.end)
  } else {
    i.seq <- 0
  }
  ti <- c(2^i.seq,n+1)
  
  # les sorties de l'algo
  theta.bar <- matrix(0, nrow = d, ncol = length(ti)) # theta.bar.i
  theta.bar.d0 <- theta.bar # theta.bar tronqué
  y.hat <- numeric(n) # prevision
  theta.hat <- matrix(0, ncol = d, nrow = n)
  
  for (i in 1:(length(ti)-1)) {
    t.next <- ti[i]:(ti[i+1]-1) # prochains instants de prevision
    id = order(abs(theta.bar[,i]),decreasing = TRUE)[1:d0] # composantes de theta.bar à tronquer
    theta.bar.d0[id,i] <- theta.bar[id,i] # troncature
    
    # les coins de la petite boule centrée en theta.bar.d0 dans laquelle theta.star est garanti
    corners <- theta.bar.d0[,i] + epsilon * 2^(-i/2) * cbind(diag(d),-diag(d))
    
    # On projette ces coins sur la boule centrée en 0 de rayon U
    for (j in 1:ncol(corners)){
      corners[,j] <- projectionL1(x = corners[,j], U = U)
    }
    
    experts.forecasts <- X[t.next,] %*% corners
    # On fait BOA
    m <- mixture(Y = Y[t.next], experts = experts.forecasts[t.next,], model = "BOA")
    # les previsions de BOA
    theta.hat[t.next,] <- m$weights %*% t(corners)
    y.hat[t.next] <- m$prediction
    theta.bar[,i+1] <- apply(theta.hat[t.next,],1,mean)
  }
  return(list(y.hat,theta.hat,theta.bar,theta.bar.d0))
}
