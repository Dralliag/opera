

#' @export 
predict.mixture <- function(object, newexperts = NULL, newY = NULL, awake = NULL,
                            online = FALSE, type=c("object$model","response"),...)
  {
    # test possible warnings and errors
    if (is.null(newexperts)) {
      cat("Renvoyer les prévisions du modèle fitté")
    }
    K = ncol(newexperts)
    if (is.null(object$training) && (object$coefficients != "Uniform") && (object$model == "MLpol")) {
      stop(paste(object$model, "cannot handle non-uniform prior weight vector"))
    }
    if (object$coefficients[1] == "Uniform") {object$coefficients <- rep(1/K,K)}
    if (!is.null(object$loss.type$tau) && object$loss.type$name != "pinball") {
      warning("Unused parameter tau (loss.type != 'pinball')")
    }

    if (object$loss.type$name != "square" && (object$model == "Ridge" || object$model == "gamMixture")) {
      stop(paste("Square loss is require for", object$model, "aggregationRule."))
    }
    if (min(newY) <= 0 && object$loss.type$name == "percentage") {
      stop("Y should be non-negative for percentage loss function")
    }
    
    if (object$model == "Ridge") {
      if (is.null(object$parameters$lambda)) {
        res <- ridgeCalib(y = newY, experts = newexperts, w0 = object$coefficients, gamma = object$parameters$gamma,
                          grid.lambda = object$parameters$grid.lambda, training = object$training)
        res$parameters$lambda <- c(object$parameters$lambda,res$parameters$lambda)
      } else {
        res <- ridge(y = newY, experts = newexperts, lambda = object$parameters$lambda, w0 = object$coefficients, 
          training = object$training)
      }
    }
    
    if (object$model == "MLpol") {
      res <- MLpol(y = newY, experts = newexperts, awake = awake, loss.type = object$loss.type, 
                  loss.gradient = object$loss.gradient, training = object$training)
      res$parameters <- list(eta = rbind(object$parameters$eta,res$parameters$eta))
    }

    if ((object$model == "BOA") || (object$model == "MLewa") || (object$model == "MLprod")) {
      algo <- eval(parse(text = object$model))
      res <- algo(y = newY, experts = newexperts, awake = awake, loss.type = object$loss.type, 
                  loss.gradient = object$loss.gradient, w0 = object$coefficients, training = object$training)
      res$parameters <- list(eta = rbind(object$parameters$eta,res$parameters$eta))
    }

    if (object$model == "EWA") {
      if (is.null(object$parameters$eta)  || !is.null(object$parameters$grid.eta)) {
        if (is.null(object$parameters$grid.eta)) {object$parameters$grid.eta = 1}
        
        res <- ewaCalib(y = newY, experts = newexperts, 
                        awake = awake, loss.type = object$loss.type, 
                        loss.gradient = object$loss.gradient, w0 = object$coefficients,  
                        gamma = object$parameters$gamma,
                        grid.eta = sort(object$parameters$grid.eta), training = object$training)
        res$parameters$eta <- c(object$parameters$eta,res$parameters$eta)
      } else {
        res <- ewa(y = newY, experts = newexperts, eta = object$parameters$eta, 
                   awake = awake, loss.type = object$loss.type, 
                   loss.gradient = object$loss.gradient, w0 = object$coefficients, 
                   training = object$training)
      }
    }
    
    if ((object$model == "FS")) {
      if (is.null(object$parameters$eta) || is.null(object$parameters$alpha) || !is.null(object$parameters$grid.eta) || !is.null(object$parameters$grid.alpha)) {
        if (is.null(object$parameters$grid.eta)) {object$parameters$grid.eta = 1}
        if (is.null(object$parameters$grid.alpha)) {object$parameters$grid.alpha = 10^(-4:-1)}
        res <- fixedshareCalib(y = newY, experts = newexperts, awake = awake, 
                               loss.type = object$loss.type, 
                               loss.gradient = object$loss.gradient, w0 = object$coefficients,  gamma = object$parameters$gamma,
                               grid.eta = object$parameters$grid.eta, grid.alpha = object$parameters$grid.alpha, training = object$training)
      } else {
        res <- fixedshare(y = newY, experts = newexperts, eta = object$parameters$eta, alpha = object$parameters$alpha, awake = awake, loss.type = object$loss.type, 
                          loss.gradient = object$loss.gradient, w0 = object$coefficients, training = object$training)
      }
    }
    
    if ((object$model == "gamMixture")) {
      warning("This aggregation rule is not stable in the current version")
      if (is.null(aggregationRule$lambda)) {
        stop("gamMixture cannot handle automatic calibration")
      }
      if (is.null(aggregationRule$z)){
        stop("A matrix (or vector) of exogeneous variables z must be given in aggregationRule")
      }
      if (is.null(aggregationRule$nknots)){ aggregationRule$nknots = 5}
      if (is.null(aggregationRule$degree)){ aggregationRule$degree = 3}
      if (is.null(aggregationRule$uniform)){ aggregationRule$uniform = FALSE}
      if (is.null(aggregationRule$knots)){ aggregationRule$knots = NULL}
      
      res <- gamMixture(y = y, experts = experts, z = aggregationRule$z, aggregationRule$lambda, 
        nknots = aggregationRule$nknots, degree = aggregationRule$degree, 
        loss.type = object$loss.type$name, uniform = aggregationRule$uniform, 
        knots = aggregationRule$knots)
    }
    #res$call <- match.call()
    res$prediction <- c(object$prediction,res$prediction)
    res$weights <- rbind(object$weights,res$weights)
    res$Y <- c(object$Y, newY)
    res$experts <- rbind(object$experts, newexperts)
    res$awake <- rbind(object$awake, awake)
    res$loss <- mean(loss(res$prediction,res$Y,loss.type=res$loss.type))

    class(res) <- "mixture"
    return(res)
  }
