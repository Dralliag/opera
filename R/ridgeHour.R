ridgeHour <-
function(y, experts, lambda, href = 1, period = 1, w0 = NULL,
                  delay = 0, y.ETR = NULL) 
{
  experts <- as.matrix(experts)
  N <- ncol(experts)
  T <- nrow(experts)
  
  if (is.null(w0)) {w0 <- rep(1/N,N)} # Uniform intial weight vector if unspecified
  if (sum(is.na(experts)) > 0) {warning("There are not allowed NA's in expert advice")}
  if (is.null(y.ETR) || (delay == 0)) {y.ETR = y}

  weights <- matrix(0, ncol = N, nrow = T)
  At <- lambda * diag(1,N)
  bt <- matrix(lambda * w0, nrow = N)
  
  w <- w0
  for (t in 1:T){  
    weights[t,] = w
    At <- At + experts[t,] %*% t(experts[t,])
    bt <- bt + y.ETR[t] * experts[t,]
    
    if ((delay > 0) && (t > delay)) {
      bt <- bt +  (y[t-delay] - y.ETR[t-delay]) * experts[t,]
    }

    h <- (((t - 1) %% period) + 1)
    if (h == href) {w <- solve(At, bt)}
  }
  prediction <- apply(experts * weights, 1, sum)
  return(list(weights = weights, prediction = prediction))
}
