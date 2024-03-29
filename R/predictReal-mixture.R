predictReal <- function(object, newexperts = NULL, newY = NULL, awake = NULL, 
                        online = TRUE, type = c("model", "response", "weights", "all"), quiet = FALSE, ...) {
  
  type <- match.arg(type)
  
  
  # Number of instant and number of experts
  if (!is.null(newY)) {
    T <- length(newY)
    if (is.null(object$names.experts)) {
      object$names.experts <- colnames(newexperts)
    }
    newexperts <- matrix(newexperts, nrow = T)
    N <- ncol(newexperts)
  } else if (object$coefficients[1] != "Uniform") {
    N <- length(object$coefficients)
    newexperts <- matrix(newexperts, ncol = N)
    T <- nrow(newexperts)
  } else {
    # warning("You should provide observations to train non trivial model")
    N <- ncol(newexperts)
    T <- nrow(newexperts)
    
    if (is.null(newexperts)) {
      result <- switch(type, model = object, response = NULL, weights = NULL, 
                       all = list(model = object, response = NULL, weights = NULL))
      return(result)
    }
  }
  
  if (!is.null(awake)) {
    awake <- matrix(awake, nrow = T)
    if (!identical(dim(awake), dim(newexperts))) {
      stop("Bad dimensions: awake and newexperts should have same dimensions")
    }
  } else {
    awake <- matrix(1, nrow = T, ncol = N)
  }
  idx.na <- which(is.na(newexperts))
  awake[idx.na] <- 0
  newexperts[idx.na] <- 0
  
  init = FALSE
  if (object$coefficients[1] == "Uniform") {
    object$coefficients <- rep(1/N, N)
    init = TRUE
  }
  
  if (length(object$coefficients) != N) {
    stop("Bad number of experts: (length(object$coefficients) != nrow(newexperts))")
  }
  
  if (!is.null(awake) && !identical(awake, matrix(1, nrow = T, ncol = N)) && (object$model == 
                                                                              "Ridge" || object$model == "OGD")) {
    stop(paste("Sleeping or missing values not allowed for", object$model, "model."))
  }
  
  
  # if no expert advice is provided, it returns the fitted object
  if (is.null(newexperts)) {
    if (!is.null(newY)) {
      stop("Expert advice should be provided if newY is non null")
    }
    result <- switch(type, model = object, response = object$prediction, weights = object$weights, 
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
    newweights <- (t(t(awake) * c(w)) / pond)[seq(1,T,by=object$d),]
  }  
  
  # Online prediction and weights if newY is provided
  if (!is.null(newY)) {
    ## If averaged is true, the models do not use coefficient as the next weight
    if (!is.null(object$parameters$averaged) && object$parameters$averaged && !is.null(object$training)) {
      object$coefficients <- object$training$next.weights
    }
    
    if (object$model == "Ridge") {
      if (is.null(object$parameters$lambda) || !is.null(object$parameters$grid.lambda)) {
        newobject <- ridgeCalib(y = newY, experts = newexperts, w0 = object$coefficients, 
                                gamma = object$parameters$gamma, grid.lambda = object$parameters$grid.lambda, 
                                training = object$training, quiet = quiet)
        newobject$parameters$lambda <- c(object$parameters$lambda, newobject$parameters$lambda)
      } else {
        newobject <- ridge(y = newY, experts = newexperts, lambda = object$parameters$lambda, 
                           w0 = object$coefficients, training = object$training, quiet = quiet)
      }
      newobject$loss.gradient = FALSE
    }
    
    if (object$model == "MLpol") {
      newobject <- MLpol(y = newY, experts = newexperts, awake = awake, loss.type = object$loss.type, 
                         loss.gradient = object$loss.gradient, training = object$training, quiet = quiet)
      newobject$parameters <- list(eta = rbind(object$parameters$eta, newobject$parameters$eta))
    }
    
    if (object$model == "OGD") {
      if (is.null(object$parameters$alpha)) {object$parameters$alpha = 0.5}
      if (is.null(object$parameters$simplex)) {object$parameters$simplex = TRUE}
      newobject <- OGD(y = newY, experts = newexperts, loss.type = object$loss.type, 
                       training = object$training, alpha = object$parameters$alpha, simplex = object$parameters$simplex,
                       w0 = object$coefficients, quiet = quiet)
    }
    
    if ((object$model == "BOA") || (object$model == "MLewa") || (object$model == 
                                                                 "MLprod")) {
      algo <- eval(parse(text = object$model))
      if (object$model != "MLewa")  {
        newobject <- algo(y = newY, experts = newexperts, awake = awake, loss.type = object$loss.type, 
                          loss.gradient = object$loss.gradient, w0 = object$coefficients, 
                          training = object$training, quiet = quiet) 
      } else {
        newobject <- algo(y = newY, experts = newexperts, awake = awake, loss.type = object$loss.type, 
                          loss.gradient = object$loss.gradient, w0 = object$coefficients, 
                          training = object$training, quiet = quiet)
      }
      newobject$parameters <- list(eta = rbind(object$parameters$eta, newobject$parameters$eta))
    }
    
    if (object$model == "EWA") {
      if (is.null(object$parameters$eta) || !is.null(object$parameters$grid.eta)) {
        newobject <- ewaCalib(y = newY, experts = newexperts, awake = awake, 
                              loss.type = object$loss.type, loss.gradient = object$loss.gradient, 
                              w0 = object$coefficients, gamma = object$parameters$gamma, grid.eta = sort(object$parameters$grid.eta), 
                              training = object$training, quiet = quiet)
        newobject$parameters$eta <- c(object$parameters$eta, newobject$parameters$eta)
      } else {
        newobject <- ewa(y = newY, experts = newexperts, eta = object$parameters$eta, 
                         awake = awake, loss.type = object$loss.type, loss.gradient = object$loss.gradient, 
                         w0 = object$coefficients, training = object$training, quiet = quiet)
      }
    }
    
    if (object$model == "FS") {
      if (is.null(object$parameters$eta) || is.null(object$parameters$alpha) || 
          !is.null(object$parameters$grid.eta) || !is.null(object$parameters$grid.alpha)) {
        if (is.null(object$parameters$grid.alpha)) {
          object$parameters$grid.alpha <- 10^(-4:-1)
        }
        newobject <- fixedshareCalib(y = newY, experts = newexperts, awake = awake, 
                                     loss.type = object$loss.type, loss.gradient = object$loss.gradient, 
                                     w0 = object$coefficients, gamma = object$parameters$gamma, grid.eta = object$parameters$grid.eta, 
                                     grid.alpha = object$parameters$grid.alpha, training = object$training, quiet = quiet)
        newobject$parameters$eta <- c(object$parameters$eta, newobject$parameters$eta)
        newobject$parameters$alpha <- c(object$parameters$alpha, newobject$parameters$alpha)
      } else {
        newobject <- fixedshare(y = newY, experts = newexperts, eta = object$parameters$eta, 
                                alpha = object$parameters$alpha, awake = awake, loss.type = object$loss.type, 
                                loss.gradient = object$loss.gradient, w0 = object$coefficients, 
                                training = object$training, quiet = quiet)
      }
    }
    
    if (object$model == "FTRL") {
      if (is.null(object$training) && ! any(c("fun_reg", "constr_ineq", "constr_eq") %in% names(object$parameters))) {
        default <- TRUE
      } else {
        default <- FALSE
      }
      if (init) {
        object$coefficients = NULL
      }
      newobject <- FTRL("y" = newY, "experts" = newexperts, 
                        "eta" = object$parameters$eta,
                        "fun_reg" = object$parameters$fun_reg, "fun_reg_grad" = object$parameters$fun_reg_grad, 
                        "constr_eq" = object$parameters$constr_eq, "constr_eq_jac" = object$parameters$constr_eq_jac, 
                        "constr_ineq" = object$parameters$constr_ineq, "constr_ineq_jac" = object$parameters$constr_ineq_jac, 
                        "max_iter" = object$parameters$max_iter,
                        "obj_tol" = object$parameters$obj_tol,
                        "loss.type" = object$loss.type, "loss.gradient" = object$loss.gradient, 
                        "w0" = object$coefficients,
                        "training" = object$training,
                        "default" =  default, "quiet" = quiet)
    }
    
    
    newobject$Y <- rbind(object$Y, matrix(newY, ncol = object$d))
    newobject$experts <- rbind(object$experts, newexperts)
    newobject$names.experts <- object$names.experts
    if (is.null(newobject$names.experts)) {
      if (!is.null(colnames(newexperts))) {
        newobject$names.experts <- colnames(newexperts)
      } else {
        if (!is.null(names(newexperts))) {
          newobject$names.experts <- colnames(newexperts)
        } else {
          newobject$names.experts <- paste("X",1:N,sep="")
        }
      }
    }
    newobject$awake <- rbind(object$awake, awake)
    
    colnames(newobject$experts) <- newobject$names.experts
    colnames(newobject$weights) <- newobject$names.experts
    colnames(newobject$awake) <- newobject$names.experts
    
    # Averaging of the weights if asked by the averaged parameter
    if (is.null(object$parameters$averaged)) {
      newobject$parameters$averaged = FALSE
    } else {
      newobject$parameters$averaged = object$parameters$averaged
    }
    
    newobject$weights = newobject$weights[seq(1,T,by=object$d),]
    
    if (newobject$parameters$averaged) {
      if (object$T == 0) {
        newweights.avg <- apply(newobject$weights, 2, cumsum) / (1:(T/object$d))
      } else {
        newweights.avg <- (object$training$sumweights + apply(newobject$weights, 2, cumsum)) / (object$T + 1:(T/object$d))
      }
      
      newobject$training$sumweights <- (object$T + T/object$d) * newweights.avg[T/object$d,] + newobject$coefficients
      newobject$training$next.weights <- newobject$coefficients
      newobject$coefficients <- newobject$training$sumweights / (object$T + T/object$d + 1)
    }
    
    # If online is true, we use online prediction
    if (online) {
      if (newobject$parameters$averaged) {
        newweights <- newweights.avg
        newpred <- rowSums(newweights.avg * newexperts)
      } else {
        newweights <- newobject$weights
        newpred <- newobject$prediction
      }
    }
    
    newobject$prediction <- rbind(object$prediction, matrix(newpred, ncol = object$d))
    newobject$weights <- rbind(object$weights, newweights)
    rownames(newobject$weights) <- NULL
    newobject$loss <- mean(loss(x = c(newobject$prediction), y = c(newobject$Y), loss.type = newobject$loss.type)) 
    newobject$T <- object$T + T/object$d
    newobject$d <- object$d
  } else {
    newobject <- object
  }
  class(newobject) <- "mixture"
  
  
  
  result <- switch(type, model = newobject, response = matrix(newpred, ncol = object$d), weights = newweights, 
                   all = list(model = newobject, response = newpred, weights = newweights))
  
  return(result)
} 
