# Ridge aggregation rule with automatic calibration of smoothing parameters
ridgeCalib <- function(y, experts, grid.lambda = 1, w0 = NULL, trace = FALSE, gamma = 2, 
  training = NULL) {
  experts <- as.matrix(experts)
  
  N <- ncol(experts)  # Number of experts
  T <- nrow(experts)  # Number of instants
  
  # Uniform intial weight vector if unspecified
  if (is.null(w0)) {
    w0 <- matrix(1/N, ncol = N)
  }
  if (is.null(grid.lambda)) {
    grid.lambda <- 1
  }
  if (is.null(gamma)) {
    gamma <- 2
  }
  
  # Smoothing parameter grid
  nlambda <- length(grid.lambda)
  grid.lambda <- matrix(grid.lambda, nrow = nlambda)
  cumulativeLoss <- rep(0, nlambda)
  
  weights <- matrix(0, ncol = N, nrow = T)  # Matrix of weights formed by the mixture
  prediction <- rep(0, T)  # Vector of predictions formed by the mixing algorithm
  pred.lambda <- matrix(0, ncol = nlambda, nrow = T)  # Prediction of mixture algorithm with different learning rates eta
  
  if (!is.null(training)) {
    At <- training$At
    bt <- training$bt
    bestlambda <- training$bestlambda
    wlambda <- training$wlambda
    w0 <- training$w0
    cumulativeLoss <- training$cumulativeLoss
    T0 <- training$T
  } else {
    At <- diag(0, N)
    bt <- rep(0, N)
    bestlambda <- floor(nlambda)/2 + 1  # We start with the parameter in the middle of the grid
    wlambda <- array(w0, dim = c(N, nlambda))  # Weight matrix proposed by each Ridge(lambda) where lambda is a parameter of the grid
    T0 <- 0
  }
  
  lambda <- rep(grid.lambda[bestlambda], T)
  for (t in 1:T) {
    # Display the state of progress of the algorithm
    if (!(t%%floor(T/10)) && trace) 
      cat(floor(10 * t/T) * 10, "% -- ")
    
    # Weights, prediction forme by Ridge(lambda[t]) where lambda[t] is the learning
    # rate calibrated online
    weights[t, ] <- wlambda[, bestlambda]
    prediction[t] <- experts[t, ] %*% weights[t, ]
    lambda[t] <- grid.lambda[bestlambda]
    
    # Weights, predictions formed by each Ridge(lambda) for lambda in the grid
    # 'grid.lambda'
    pred.lambda[t, ] <- experts[t, ] %*% wlambda
    cumulativeLoss <- cumulativeLoss + (pred.lambda[t, ] - y[t])^2
    
    # Parameter update
    At <- At + experts[t, ] %*% t(experts[t, ])
    bt <- bt + y[t] * experts[t, ]
    
    # Grid update **************
    bestlambda <- order(cumulativeLoss)[1]  # find the best smoothing rate lambda in the grid
    
    # Expand the grid if the best parameter lies on an extremity
    if (bestlambda == nlambda) {
      if (trace) 
        cat(" + ")
      newlambda <- grid.lambda[bestlambda] * gamma^(1:3)
      grid.lambda <- c(grid.lambda, newlambda)
      nlambda <- nlambda + length(newlambda)
      for (k in 1:length(newlambda)) {
        perfnewlambda <- tryCatch(ridge(y = c(training$oldY, y[1:t]), experts = rbind(training$oldexperts, 
          matrix(experts[1:t, ], ncol = N)), lambda = newlambda[k], w0 = w0), 
          error = function(e) {
          list(prediction = rep(0, t))
          })
        newcumulativeLoss <- sum((perfnewlambda$prediction - y[1:t])^2)
        cumulativeLoss <- c(cumulativeLoss, newcumulativeLoss)
        pred.lambda <- cbind(pred.lambda, c(perfnewlambda$prediction, rep(0, 
          (T - t))))
      }
    }
    if (bestlambda == 1) {
      if (trace) 
        cat(" - ")
      newlambda <- grid.lambda[bestlambda]/gamma^(1:3)
      nlambda <- nlambda + length(newlambda)
      bestlambda <- bestlambda + length(newlambda)
      for (k in 1:length(newlambda)) {
        grid.lambda <- c(newlambda[k], grid.lambda)
        perfnewlambda <- tryCatch(y = ridge(c(training$oldY, y[1:t]), experts = rbind(training$oldexperts, 
          matrix(experts[1:t, ], ncol = N)), lambda = newlambda[k], w0 = w0), 
          error = function(e) {
          list(prediction = rep(NA, t))
          })
        newcumulativeLoss <- sum((perfnewlambda$prediction - y[1:t])^2)
        cumulativeLoss <- c(newcumulativeLoss, cumulativeLoss)
        pred.lambda <- cbind(c(perfnewlambda$prediction, rep(NA, (T - t))), 
          pred.lambda)
      }
    }
    wlambda <- matrix(0, nrow = N, ncol = nlambda)
    for (k in 1:nlambda) {
      wlambda[, k] <- tryCatch(solve(grid.lambda[k] * diag(1, N) + At, matrix(grid.lambda[k] * 
        w0, nrow = N) + bt), error = function(e) {
        NA
      })
    }
    # the smoothing parameter has to big large enough in order to have invertible
    # design matrix
    lambda.min <- which(!(is.na(wlambda[1, ])))[1]
    bestlambda <- max(lambda.min, bestlambda)
  }
  
  object <- list(model = "Ridge", loss.type = list(name = "square"), coefficients = wlambda[, 
    bestlambda])
  
  object$parameters <- list(lambda = c(lambda[1:T]), grid.lambda = c(grid.lambda))
  object$weights <- weights
  object$prediction <- prediction
  
  object$training <- list(T = T0 + T, wlambda = wlambda, w0 = w0, At = At, bt = bt, 
    bestlambda = bestlambda, cumulativeLoss = cumulativeLoss, grid.loss = cumulativeLoss/(T0 + 
      T), oldexperts = rbind(training$oldexperts, experts), oldY = c(training$oldY, 
      y))
  
  if (trace) 
    cat("\n")
  return(object)
} 
