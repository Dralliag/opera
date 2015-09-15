

#' @export 
predict.mixture <- function(object, newexpert = NULL, newY = NULL, 
                            online = FALSE, type=c("object$model","response"),...)
  {

    # test possible warnings and errors
    if (is.null(object$parameters$gamma)) {
      object$parameters$gamma = 2
    }
    if (is.null(newexpert)) {
      cat("Renvoyer les prévisions du modèle fitté")
    }
    K = ncol(newexpert)
    if (is.null(object$training) && (object$coefficients != "Uniform") && (object$model == "MLpol")) {
      stop(paste(object$model, "cannot handle non-uniform prior weight vector"))
    }
    if (object$coefficients == "Uniform") {object$coefficients <- rep(1/K,K)}
    if (!is.null(object$loss.type$tau) && object$loss.type$name != "pinball") {
      warning("Unused parameter tau (loss.type != 'pinball')")
    }
    if (is.null(object$loss.type$tau)) {tau = 0.5} else {tau=object$loss.type$tau} # voir comment le rajouter dans ...

    if (object$loss.type$name != "square" && (object$model == "Ridge" || object$model == "gamMixture")) {
      stop(paste("Square loss is require for", object$model, "aggregationRule."))
    }
    if (min(newY) <= 0 && object$loss.type$name == "percentage") {
      stop("Y should be non-negative for percentage loss function")
    }
    
    if (object$model == "Ridge") {
      if (is.null(object$parameters$lambda)) {
        res <- ridgeCalib(y = newY, experts = newexpert, w0 = object$coefficients, gamma = object$parameters$gamma,
                          grid.lambda = object$parameters$grid.lambda, training = object$training)
      } else {
        res <- ridge(y = newY, experts = newexpert, lambda = object$parameters$lambda, w0 = object$coefficients, training = object$training)
      }
    }
    
    if (object$model == "MLpol") {
      res <- MLpol(y = newY, experts = newexpert, awake = object$awake, loss.type = object$loss.type$name, 
                  loss.gradient = object$loss.gradient, tau = tau, 
                  training = object$training)
    }

    if ((object$model == "BOA") || (object$model == "MLewa") || (object$model == "MLprod")) {
      algo <- eval(parse(text = object$model))
      res <- algo(y, experts, awake = awake, loss.type = object$loss.type$name, 
                  loss.gradient = object$loss.gradient,
                  tau = tau, w0 = object$coefficients, training = object$training)
    }

    if (object$model == "EWA") {
      if (is.null(aggregationRule$eta)) {
        if (is.null(aggregationRule$grid.eta)) {aggregationRule$grid.eta = 1}
        res <- ewaCalib(y = y, experts = experts, 
                        awake = awake, loss.type = object$loss.type$name, 
                        loss.gradient = object$loss.gradient, w0 = w0, 
                        tau = tau, gamma = aggregationRule$gamma,
                        grid.eta = sort(aggregationRule$grid.eta))
      } else {
        res <- ewa(y = y, experts = experts, eta = aggregationRule$eta, 
                   awake = awake, loss.type = object$loss.type$name, 
                   loss.gradient = object$loss.gradient, w0 = w0,
                   tau = tau)
      }
    }
    
    if ((object$model == "FS")) {
      if (is.null(aggregationRule$eta) || is.null(aggregationRule$alpha)) {
        if (is.null(aggregationRule$grid.eta)) {aggregationRule$grid.eta = 1}
        if (is.null(aggregationRule$grid.alpha)) {aggregationRule$grid.alpha = 10^(-4:-1)}
        res <- fixedshareCalib(y = y, experts = experts, awake = awake, 
                               loss.type = object$loss.type$name, 
                               loss.gradient = object$loss.gradient, w0 = w0,
                               tau = tau, gamma = aggregationRule$gamma,
                               grid.eta = aggregationRule$grid.eta, grid.alpha = aggregationRule$grid.alpha)
      } else {
        res <- fixedshare(y = y, experts = experts, eta = aggregationRule$eta, alpha = aggregationRule$alpha, awake = awake, loss.type = object$loss.type$name, 
                          loss.gradient = object$loss.gradient, w0 = w0,
                          tau = tau)
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
      if (is.null(object$loss.type$name)){ object$loss.type$name = "square"}
      if (is.null(aggregationRule$uniform)){ aggregationRule$uniform = FALSE}
      if (is.null(aggregationRule$knots)){ aggregationRule$knots = NULL}
      
      res <- gamMixture(y = y, experts = experts, z = aggregationRule$z, aggregationRule$lambda, 
        nknots = aggregationRule$nknots, degree = aggregationRule$degree, 
        loss.type = object$loss.type$name, uniform = aggregationRule$uniform, 
        knots = aggregationRule$knots)
    }
    res$object$model <- object$model
    res$residuals <- newY - res$prediction
    res$call <- match.call()
    res$Y <- newY
    res$experts <- newexpert

    class(res) <- "mixture"
    return(res)
  }
