# Ridge aggregation rule with automatic calibration of smoothing parameters
ridgeCalib <-
function(y, experts, gridlambda = 1, w0 = NULL, trace = F, gamma = 2)
{

  experts <- as.matrix(experts)
  
  N <- ncol(experts)  # Number of experts
  T <- nrow(experts)  # Number of instants
  
  if (is.null(w0)) {w0 <- matrix(1/N, ncol = N)} # Uniform intial weight vector if unspecified
  
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
  bt <- rep(0, N)

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
          error = function(e) {list(prediction = rep(0, t))})
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
          ridge(y[1:t], matrix(experts[1:t,],ncol = N), newlambda[k], w0 = w0), 
          error = function(e) {list(prediction = rep(NA, t))})
        newcumulativeLoss <- sum((perfnewlambda$prediction - y[1:t])^2)
        cumulativeLoss <- c(newcumulativeLoss, cumulativeLoss)
        pred.lambda <- cbind(c(perfnewlambda$prediction, rep(NA, (T-t))), pred.lambda)
      }
    }
    wlambda <- matrix(0,nrow = N, ncol = nlambda)
    for (k in 1:nlambda) {
      wlambda[,k] = tryCatch(solve(gridlambda[k]*diag(1,N) + At, matrix(gridlambda[k]*w0, nrow=N) + bt),
                               error = function(e) {NA})
    }
    # the smoothing parameter has to big large enough
    # in order to have invertible design matrix
    lambda.min = which(!(is.na(wlambda[1,])))[1]
    bestlambda = max(lambda.min, bestlambda)
  }
  l <- mean(loss(prediction,y))
  grid.loss <- cumulativeLoss / T
  l.rmse <- sqrt(l)
  grid.rmse <- sqrt(grid.loss)
  if (trace) cat('\n')
  return(list(weights = weights, prediction = prediction, 
              lambda = lambda, grid.lambda = gridlambda, 
              loss = l, grid.loss = grid.loss, weights.forecast = wlambda[,bestlambda],
              rmse = l.rmse, grid.rmse = grid.rmse))
}
