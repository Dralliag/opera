
ridge <-
function(y, experts, lambda, w0 = NULL) {
  experts = as.matrix(experts)
  N <- ncol(experts)
  T <- nrow(experts)
  
  if (is.null(w0)) {w0 <- matrix(1/N,ncol=N)} # Uniform intial weight vector if unspecified
  if (sum(is.na(experts)) > 0) {warning("There are not allowed NA's in expert advice")}

  R <- rep(0,N)
  w <- matrix(0, ncol = N, nrow = T)
  
  At = lambda * diag(1,N)
  bt = matrix(lambda * w0, nrow=N)
  for (t in 1:T){
    w[t,] = solve(At,bt)
    At = At + experts[t,] %*% t(experts[t,])
    bt = bt + y[t] * experts[t,]
  }
  #w[1,] = w0
  
  weights.forecast = solve(At,bt)
  prediction <- apply(experts * w, 1, sum)
  l <- mean(loss(prediction,y))
  return(list(weights = w, prediction = prediction, 
    weights.forecast = weights.forecast,
    loss =  l))
}
