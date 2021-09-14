# Ridge aggregation rule with automatic calibration of smoothing parameters
ridgeCalib <- function(y, experts, grid.lambda = 1, w0 = NULL, gamma = 2, 
                       training = NULL, use_cpp = getOption("opera_use_cpp", default = FALSE), quiet = FALSE) {
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
  
  if (!is.null(training)) {
    At <- training$At
    bt <- training$bt
    bestlambda <- training$bestlambda
    wlambda <- training$wlambda
    w0 <- training$w0
    cumulativeLoss <- training$cumulativeLoss
    T0 <- training$T
  } else {
    At <- list()
    for (k in 1:nlambda){
      At[[k]] <- diag(0,N) / grid.lambda[k]
    }
    bt <- rep(0, N)
    bestlambda <- floor(nlambda)/2 + 1  # We start with the parameter in the middle of the grid
    wlambda <- array(w0, dim = c(N, nlambda))  # Weight matrix proposed by each Ridge(lambda) where lambda is a parameter of the grid
    T0 <- 0
  }
  
  
  lambda <- rep(grid.lambda[bestlambda], T)
  
  if (! quiet) steps <- init_progress(T)
  
  for (t in 1:T) {
    if (! quiet) update_progress(t, steps)
    
    # Weights, prediction formed by Ridge(lambda[t]) where lambda[t] is the learning
    # rate calibrated online
    if (use_cpp){
      bestlambda<-RidgeCalibStep1(t,bestlambda,
                                  experts, weights,
                                  wlambda, w0,
                                  bt, 
                                  grid.lambda,  
                                  y, lambda,
                                  cumulativeLoss, prediction)
    }
    else{
      weights[t, ] <- wlambda[, bestlambda]
      prediction[t] <- experts[t, ] %*% weights[t, ]
      lambda[t] <- grid.lambda[bestlambda]
      
      # Weights, predictions formed by each Ridge(lambda) for lambda in the grid
      # 'grid.lambda'
      cumulativeLoss <- cumulativeLoss + c(experts[t, ] %*% wlambda - y[t])^2
      # Grid update **************
      bestlambda <- order(cumulativeLoss)[1]  # find the best smoothing rate lambda in the grid
      bt <- bt + y[t] * experts[t, ]
    }
    
    # Parameter update
    for (k in 1:nlambda){
      a <- At[[k]] %*% experts[t, ]
      At[[k]] <- At[[k]] - a %*% t(a) / c(1 + experts[t,] %*% a)
    }
    
    # Expand the grid if the best parameter lies on an extremity
    if (bestlambda == nlambda) {
      newlambda <- grid.lambda[bestlambda] * gamma^(1:3)
      grid.lambda <- c(grid.lambda, newlambda)
      for (k in 1:length(newlambda)) {
        perfnewlambda <- tryCatch(ridge(y = c(training$oldY, y[1:t]), experts = rbind(training$oldexperts, 
                                                                                      matrix(experts[1:t, ], ncol = N)), lambda = newlambda[k], w0 = w0, use_cpp = use_cpp, quiet = TRUE), 
                                  error = function(e) {
                                    list(prediction = rep(0, t))
                                  })
        newcumulativeLoss <- sum((perfnewlambda$prediction - c(training$oldY, y[1:t]))^2)
        At[[nlambda + k]] <- perfnewlambda$training$At
        cumulativeLoss <- c(cumulativeLoss, newcumulativeLoss)
      }
      nlambda <- nlambda + length(newlambda)
    }
    if (bestlambda == 1) {
      newlambda <- grid.lambda[bestlambda]/gamma^(1:3)
      nlambda <- nlambda + length(newlambda)
      bestlambda <- bestlambda + length(newlambda)
      for (k in 1:length(newlambda)) {
        grid.lambda <- c(newlambda[k], grid.lambda)
        perfnewlambda <- tryCatch(ridge(c(training$oldY, y[1:t]), experts = rbind(training$oldexperts, 
                                                                                  matrix(experts[1:t, ], ncol = N)), lambda = newlambda[k], w0 = w0, use_cpp = use_cpp, quiet = TRUE), 
                                  error = function(e) {
                                    list(prediction = rep(NA, t))
                                  })
        newcumulativeLoss <- sum((perfnewlambda$prediction - c(training$oldY, y[1:t]))^2)
        cumulativeLoss <- c(newcumulativeLoss, cumulativeLoss)
        At = c(list(perfnewlambda$training$At),At)
      }
    }
    if (nlambda!=ncol(wlambda)){
      wlambda <- matrix(0, nrow = N, ncol = nlambda)
    }
    for (k in 1:nlambda) {
      wlambda[, k] <-  At[[k]] %*% (grid.lambda[k] * w0 + bt)
    }
    # the smoothing parameter has to big large enough in order to have invertible
    # design matrix
    lambda.min <- which(!(is.na(wlambda[1, ])))[1]
    bestlambda <- max(lambda.min, bestlambda)
  }
  if (! quiet) end_progress()
  
  object <- list(model = "Ridge", loss.type = list(name = "square"), coefficients = wlambda[, 
                                                                                            bestlambda])
  
  object$parameters <- list(lambda = c(lambda[1:T]), grid.lambda = c(grid.lambda))
  object$weights <- weights
  object$prediction <- prediction
  
  object$training <- list(T = T0 + T, wlambda = wlambda, w0 = w0, At = At, bt = bt, 
                          bestlambda = bestlambda, cumulativeLoss = cumulativeLoss, grid.loss = cumulativeLoss/(T0 + 
                                                                                                                  T), oldexperts = rbind(training$oldexperts, experts), oldY = c(training$oldY, 
                                                                                                                                                                                 y))
  
  return(object)
} 
