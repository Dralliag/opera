

#' @export 
predict.mixture <- function(object, newexperts = NULL, newY = NULL, awake = NULL, 
  online = TRUE, type = c("model", "response", "weights", "all"), z = NULL, ...) {
  
  type = match.arg(type)
  
  # Number of instant and number of experts
  if (!is.null(newY)) {
    T <- length(newY)
    newexperts <- matrix(newexperts, nrow = T)
    N <- ncol(newexperts)
  } else if (object$coefficients[1] != "Uniform") {
    N <- length(object$coefficients)
    newexperts <- matrix(newexperts, ncol = N)
    T <- nrow(newexperts)
  } else {
    warning("You should provide observations to train non trivial model")
    N <- ncol(newexperts)
    T <- nrow(newexperts)

    if (is.null(newexperts)) {
      result <- switch(type,
                model = object,
                response = NULL,
                weights = NULL,
                all = list(model = object, response = NULL, weights = NULL))
      return(result)
    }
  }

  if (!is.null(awake)) {
    awake <- matrix(awake, nrow = T)
    if (!identical(dim(awake),dim(newexperts))) {
      stop("Bad dimensions: awake and newexperts should have same dimensions")
    }
  } else {
    awake <- matrix(1, nrow = T, ncol = N)
  }
  idx.na <- which(is.na(newexperts))
  awake[idx.na] <- 0
  newexperts[idx.na] <- 0

  # test possible warnings and errors
  if (is.null(object$training) && (object$coefficients != "Uniform") && (object$model == "MLpol")) {
    stop(paste(object$model, "cannot handle non-uniform prior weight vector"))
  }

  if (object$coefficients[1] == "Uniform") {
    object$coefficients <- rep(1/N, N)
  }

  if (length(object$coefficients) != N) {
    stop("Bad number of experts: (length(object$coefficients) != nrow(newexperts))")
  }

  if (!is.null(object$loss.type$tau) && object$loss.type$name != "pinball") {
    warning("Unused parameter tau (loss.type$name != 'pinball)")
  }
  
  if (object$loss.type$name != "square" && (object$model == "Ridge" || object$model == "gamMixture")) {
    stop(paste("Square loss is require for", object$model, "model."))
  }

  # if no expert advice is provided, it returns the fitted object
  if (is.null(newexperts)) {
    if (!is.null(newY)) {
      stop("Expert advice should be provided if newY is non null")
    }
    if (object$loss.type$name == "gamMixture") {
      stop("Method predict is currently not implemented for 'gamMixture' aggregation rule")
    }
    result <- switch(type,
            model = object,
            response = object$prediction,
            weights = object$weights,
            all = list(model = object, response = object$prediction, weights = object$weights))
    return(result)
  }

  if (!is.null(newexperts) && is.null(newY) && online) {
    stop("newY cannot be null to perform online prediction. Provide newY or set online = FALSE")
  }
  
  # Batch prediction and weights
  if (!online) {
    w <- matrix(object$coefficients, ncol = 1)
    pond <- c(awake %*% w)
    newpred <- c(((newexperts * awake) %*% w)/pond)
    newweights <- t(matrix(rep(w,T),ncol=T)) / pond
  }

  # Online prediction and weights if newY is provided
  if (!is.null(newY)) {

    if (min(newY) <= 0 && object$loss.type$name == "percentage") {
      stop("Y should be non-negative for percentage loss function")
    }

    if (object$model == "Ridge") {
      if (is.null(object$parameters$lambda)) {
        newobject <- ridgeCalib(y = newY, experts = newexperts, w0 = object$coefficients, gamma = object$parameters$gamma, 
          grid.lambda = object$parameters$grid.lambda, training = object$training)
        newobject$parameters$lambda <- c(object$parameters$lambda, newobject$parameters$lambda)
      } else {
        newobject <- ridge(y = newY, experts = newexperts, lambda = object$parameters$lambda, w0 = object$coefficients, 
          training = object$training)
      }
    }
    
    if (object$model == "MLpol") {
      newobject <- MLpol(y = newY, experts = newexperts, awake = awake, loss.type = object$loss.type, loss.gradient = object$loss.gradient, 
        training = object$training)
      newobject$parameters <- list(eta = rbind(object$parameters$eta, newobject$parameters$eta))
    }
    
    if ((object$model == "BOA") || (object$model == "MLewa") || (object$model == "MLprod")) {
      algo <- eval(parse(text = object$model))
      newobject <- algo(y = newY, experts = newexperts, awake = awake, loss.type = object$loss.type, loss.gradient = object$loss.gradient, 
        w0 = object$coefficients, training = object$training)
      newobject$parameters <- list(eta = rbind(object$parameters$eta, newobject$parameters$eta))
    }
    
    if (object$model == "EWA") {
      if (is.null(object$parameters$eta) || !is.null(object$parameters$grid.eta)) {
        if (is.null(object$parameters$grid.eta)) {
          object$parameters$grid.eta <- 1
        }
        
        newobject <- ewaCalib(y = newY, experts = newexperts, awake = awake, loss.type = object$loss.type, 
          loss.gradient = object$loss.gradient, w0 = object$coefficients, gamma = object$parameters$gamma, 
          grid.eta = sort(object$parameters$grid.eta), training = object$training)
        newobject$parameters$eta <- c(object$parameters$eta, newobject$parameters$eta)
      } else {
        newobject <- ewa(y = newY, experts = newexperts, eta = object$parameters$eta, awake = awake, loss.type = object$loss.type, 
          loss.gradient = object$loss.gradient, w0 = object$coefficients, training = object$training)
      }
    }
    
    if ((object$model == "FS")) {
      if (is.null(object$parameters$eta) || is.null(object$parameters$alpha) || !is.null(object$parameters$grid.eta) || 
        !is.null(object$parameters$grid.alpha)) {
        if (is.null(object$parameters$grid.eta)) {
          object$parameters$grid.eta <- 1
        }
        if (is.null(object$parameters$grid.alpha)) {
          object$parameters$grid.alpha <- 10^(-4:-1)
        }
        newobject <- fixedshareCalib(y = newY, experts = newexperts, awake = awake, loss.type = object$loss.type, 
          loss.gradient = object$loss.gradient, w0 = object$coefficients, gamma = object$parameters$gamma, 
          grid.eta = object$parameters$grid.eta, grid.alpha = object$parameters$grid.alpha, training = object$training)
      } else {
        newobject <- fixedshare(y = newY, experts = newexperts, eta = object$parameters$eta, alpha = object$parameters$alpha, 
          awake = awake, loss.type = object$loss.type, loss.gradient = object$loss.gradient, w0 = object$coefficients, 
          training = object$training)
      }
    }
    
    if ((object$model == "gamMixture")) {
      warning("This aggregation rule is not stable in the current version")
      if (is.null(parameters$lambda)) {
        stop("gamMixture cannot handle automatic calibration")
      }
      if (is.null(z)) {
        stop("A matrix (or vector) of exogeneous variables z must be provided")
      }
      if (is.null(parameters$nknots)) {
        parameters$nknots <- 5
      }
      if (is.null(parameters$degree)) {
        parameters$degree <- 3
      }
      if (is.null(parameters$uniform)) {
        parameters$uniform <- FALSE
      }
      if (is.null(parameters$knots)) {
        parameters$knots <- NULL
      }
      
      newobject <- gamMixture(y = newY, experts = newexperts, z = z, parameters$lambda, nknots = parameters$nknots, 
        degree = parameters$degree, loss.type = object$loss.type$name, uniform = parameters$uniform, 
        knots = parameters$knots)
    }

    newobject$Y <- c(object$Y, newY)
    newobject$experts <- rbind(object$experts, newexperts)
    newobject$awake <- rbind(object$awake, awake)
  
    # If online is true, we use online prediction
    if (online) {
      newpred <- newobject$prediction
      newweights <- newobject$weights
    } 

    newobject$prediction <- c(object$prediction, newpred)
    newobject$weights <- rbind(object$weights, newweights)
    newobject$loss <- mean(loss(newobject$prediction, newobject$Y, loss.type = newobject$loss.type))

  } else {
    newobject <- object
  }
  class(newobject) <- "mixture"

  result <- switch(type, 
              model = newobject,
              response = newpred,
              weights = newweights,
              all = list(model = newobject, response = newpred, weights = newweights)
            )

  return(result)
} 
