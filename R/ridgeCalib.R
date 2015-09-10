# Ridge aggregation rule with automatic calibration of smoothing parameters
ridgeCalib <-
function(y, experts, gridlambda = 1, w0 = NULL, trace = F, gamma = 2)
{
  experts <- as.matrix(experts)
  
  N <- ncol(experts)  # Number of experts
  T <- nrow(experts)  # Number of instants
  
  if (is.null(w0)) {w0 <- matrix(1/N, ncol = N)} # Uniform intial weight vector if unspecified
  if (sum(is.na(experts)) > 0) {warning("There are not allowed NA's in expert advice")}
  
  # Smoothing parameter grid 
  nlambda <- length(gridlambda)
  bestlambda <- floor(nlambda)/2 + 1 # We start with the parameter in the middle of the grid
  gridlambda <- matrix(gridlambda, nrow = nlambda)
  
  lambda <- rep(gridlambda[bestlambda],T)
  cumulativeLoss <- rep(0,nlambda)
  
  wlambda <- array(w0, dim = c(N, nlambda))   # Weight matrix proposed by each Ridge(lambda) where lambda is a parameter of the grid
  weights <- matrix(0, ncol = N, nrow = T)    # Matrix of weights formed by the mixture
  prediction <- rep(0, T)                     # Vector of predictions formed by the mixing algorithm
  pred.lambda <- matrix(0, ncol = nlambda, nrow = T) # Prediction of mixture algorithm with different learning rates eta

  At <- diag(0, N)
  bt <- matrix(t(w0) %*% gridlambda, nrow = N, ncol = nlambda)

  for(t in 1:T){
    # Display the state of progress of the algorithm
    if (!(t %% floor(T/10)) && trace) cat(floor(10 * t/T)*10, '% -- ')
    
    # Weights, prediction forme by Ridge(lambda[t]) where lambda[t] is the learning rate calibrated online
    weights[t,] <- wlambda[,bestlambda] 
    prediction[t] <- experts[t,] %*% weights[t,]
    lambda[t] <- gridlambda[bestlambda]
    
    # Weights, predictions formed by each Ridge(lambda) for lambda in the grid "gridlambda"
    pred.lambda[t,] <- experts[t,] %*% wlambda 
    cumulativeLoss <- cumulativeLoss + (pred.lambda[t,] - y[t])^2
    
    # Mise à jour
    At <- At + experts[t,] %*% t(experts[t,])
    bt <- bt + y[t] * experts[t,]

    # Grid update
    # **************
    bestlambda <- order(cumulativeLoss)[1] # find the best smoothing rate lambda in the grid

    # Expand the grid if the best parameter lies on an extremity
    if (bestlambda == nlambda) { 
      if (trace) cat(' + ')
      newlambda <- gridlambda[bestlambda] * gamma^(1:3)
      gridlambda <- c(gridlambda, newlambda)
      nlambda <- nlambda + length(newlambda)
      for (k in 1:length(newlambda)) {
        perfnewlambda <- tryCatch(
          ridge(y[1:t], matrix(experts[1:t,],ncol=N), newlambda[k], w0 = w0), 
          error = function(e) {list(prediction = rep(0, T))})
        newcumulativeLoss <- sum((perfnewlambda$prediction - y[1:t])^2)
        cumulativeLoss <- c(cumulativeLoss, newcumulativeLoss)
        pred.lambda <- cbind(pred.lambda, c(perfnewlambda$prediction, rep(0, (T-t))))
      }
    }
    if (bestlambda == 1) {
      if (trace) cat(' - ')
      newlambda <- gridlambda[bestlambda] / gamma^(1:3)
      nlambda <- nlambda + length(newlambda)
      bestlambda <- bestlambda + length(newlambda)
      for (k in 1:length(newlambda)) {
        gridlambda <- c(newlambda[k],gridlambda)
        perfnewlambda <- tryCatch(
          ridge(y[1:t], matrix(experts[1:t,],ncol = N), newlambda[k], w0 = NULL), 
          error = function(e) {list(prediction = rep(0, T))})
        newcumulativeLoss <- sum((perfnewlambda$prediction - y[1:t])^2)
        cumulativeLoss <- c(newcumulativeLoss, cumulativeLoss)
        pred.lambda <- cbind(c(perfnewlambda$prediction, rep(0, (T-t))), pred.lambda)
      }
    }
    wlambda <- matrix(0,nrow = N, ncol = nlambda)
    for (k in 1:nlambda) {
      wlambda[,k] = tryCatch(solve(gridlambda[k]*diag(1,N) + At,bt),
                               error = function(e) {0})
    }
  }
  l <- rmse(prediction,y)
  rmse <- sqrt(cumulativeLoss / T)
  if (trace) cat('\n')
  return(list(weights = weights, prediction = prediction, 
              lambda = lambda, grid = gridlambda, 
              loss = l, gridloss = rmse, weights.forecast = wlambda[,bestlambda]))
}
