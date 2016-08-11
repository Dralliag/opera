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


