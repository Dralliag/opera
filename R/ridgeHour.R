ridgeHour <-
function(y, experts, lambda, href=1, period = 1, w0 = NULL) {
  experts = as.matrix(experts)
  N <- ncol(experts)
  T <- nrow(experts)
  
  if (is.null(w0)) {w0 <- rep(1/N,N)} # Poids initial uniforme si non spécifié
  if (sum(is.na(experts)) > 0) {warning("There are not allowed NA's in expert advice")}
  
  weights <- matrix(0,ncol=N,nrow=T)
  At = lambda * diag(1,N)
  bt = matrix(lambda * w0, nrow=N)
  
  w = w0
  for (t in 1:T){  
    weights[t,] = w
    At = At + experts[t,] %*% t(experts[t,])
    bt = bt + y[t] * experts[t,]
    h <- (((t-1)%%period)+1)
    if (h == href) {w = solve(At,bt)}
  }
  prediction <- apply(experts * weights, 1, sum)
  return(list(weights = weights, prediction = prediction))
}
