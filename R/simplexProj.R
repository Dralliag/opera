simplexProj <- function(u) {
  N <- length(u)
  v = sort(u,decreasing = TRUE)
  rho = max(which((v + (1-cumsum(v)) / (1:N)) >0))
  if (rho >= 1) {
    lambda = 1/rho * (1 - sum(v[1:rho]))
    return(pmax(u+ lambda,0))
  } else {
    return(rep(1/N,N))
  }
  
}