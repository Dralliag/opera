#' Compute an aggregation rule (method)
#' @inheritParams mixture
#' @export mixture.default

mixture.default <-
  function(y, experts, 
           aggregationRule = "MLpol",  w0 = NULL, awake = NULL)
  {
    # test possible warning and errors
    if (is.character(aggregationRule)) {
      aggregationRule = list(name = aggregationRule)
    } else {
      if (!is.list(aggregationRule)) {
        stop("Bad type for aggregationRule: it must be either a list or a character.")
      }
    }
    if (is.null(aggregationRule$name)){
      aggregationRule$name  = "MLpol"
      warning("Aggregation rule MLpol is chosen")
    }
    if (is.null(aggregationRule$gamma)) {
      aggregationRule$gamma = 2
    }
    if (!is.null(w0) && (aggregationRule$name == "MLpol")) {
      stop(paste(aggregationRule$name, "cannot handle non-uniform prior weight vector"))
    }
    if (is.null(aggregationRule$loss.type)) {aggregationRule$loss.type = "square"}
    if (is.null(aggregationRule$loss.gradient)) {aggregationRule$loss.gradient = TRUE} 
    if (!is.null(aggregationRule$tau) && aggregationRule$loss.type != "pinball") {
      warning("Unused parameter tau (loss.type != 'pinball')")
    }
    if (is.null(aggregationRule$tau)) {aggregationRule$tau = 0.5}

    if (aggregationRule$loss.type != "square" && (aggregationRule$name == "Ridge" || aggregationRule$name == "gamMixture")) {
      stop(paste("Square loss is require for", aggregationRule$name, "aggregationRule."))
    }
    if (min(y) <= 0 && model$loss.type == "percentage") {
      stop("Y should be non-negative for percentage loss function")
    }
    
    if (aggregationRule$name == "Ridge") {
      if (is.null(aggregationRule$lambda)) {
        if (is.null(aggregationRule$grid.lambda)) {aggregationRule$grid.lambda = 1}
        res <- ridgeCalib(y = y, experts = experts, w0 = w0, gamma = aggregationRule$gamma,
                          grid.lambda = aggregationRule$grid.lambda)
      } else {
        res <- ridge(y, experts, aggregationRule$lambda, w0)
      }
    } 
    
    if (aggregationRule$name == "MLpol") {
      res <- MLpol(y, experts, awake = awake, loss.type = aggregationRule$loss.type, 
                  loss.gradient = aggregationRule$loss.gradient, tau = aggregationRule$tau)
    }

    if ((aggregationRule$name == "BOA") || (aggregationRule$name == "MLewa") || (aggregationRule$name == "MLprod")) {
      algo <- eval(parse(text = aggregationRule$name))
      res <- algo(y, experts, awake = awake, loss.type = aggregationRule$loss.type, 
                  loss.gradient = aggregationRule$loss.gradient,
                  tau = aggregationRule$tau, w0 = w0)
    }

    if (aggregationRule$name == "EWA") {
      if (is.null(aggregationRule$eta)) {
        if (is.null(aggregationRule$grid.eta)) {aggregationRule$grid.eta = 1}
        res <- ewaCalib(y = y, experts = experts, 
                        awake = awake, loss.type = aggregationRule$loss.type, 
                        loss.gradient = aggregationRule$loss.gradient, w0 = w0, 
                        tau = aggregationRule$tau, gamma = aggregationRule$gamma,
                        grid.eta = sort(aggregationRule$grid.eta))
      } else {
        res <- ewa(y = y, experts = experts, eta = aggregationRule$eta, 
                   awake = awake, loss.type = aggregationRule$loss.type, 
                   loss.gradient = aggregationRule$loss.gradient, w0 = w0,
                   tau = aggregationRule$tau)
      }
    }
    
    if ((aggregationRule$name == "FS")) {
      if (is.null(aggregationRule$eta) || is.null(aggregationRule$alpha)) {
        if (is.null(aggregationRule$grid.eta)) {aggregationRule$grid.eta = 1}
        if (is.null(aggregationRule$grid.alpha)) {aggregationRule$grid.alpha = 10^(-4:-1)}
        res <- fixedshareCalib(y = y, experts = experts, awake = awake, 
                               loss.type = aggregationRule$loss.type, 
                               loss.gradient = aggregationRule$loss.gradient, w0 = w0,
                               tau = aggregationRule$tau, gamma = aggregationRule$gamma,
                               grid.eta = aggregationRule$grid.eta, grid.alpha = aggregationRule$grid.alpha)
      } else {
        res <- fixedshare(y = y, experts = experts, eta = aggregationRule$eta, alpha = aggregationRule$alpha, awake = awake, loss.type = aggregationRule$loss.type, 
                          loss.gradient = aggregationRule$loss.gradient, w0 = w0,
                          tau = aggregationRule$tau)
      }
    }
    
    if ((aggregationRule$name == "gamMixture")) {
      warning("This aggregation rule is not stable in the current version")
      if (is.null(aggregationRule$lambda)) {
        stop("gamMixture cannot handle automatic calibration")
      }
      if (is.null(aggregationRule$z)){
        stop("A matrix (or vector) of exogeneous variables z must be given in aggregationRule")
      }
      if (is.null(aggregationRule$nknots)){ aggregationRule$nknots = 5}
      if (is.null(aggregationRule$degree)){ aggregationRule$degree = 3}
      if (is.null(aggregationRule$loss.type)){ aggregationRule$loss.type = "square"}
      if (is.null(aggregationRule$uniform)){ aggregationRule$uniform = FALSE}
      if (is.null(aggregationRule$knots)){ aggregationRule$knots = NULL}
      
      res <- gamMixture(y = y, experts = experts, z = aggregationRule$z, aggregationRule$lambda, 
        nknots = aggregationRule$nknots, degree = aggregationRule$degree, 
        loss.type = aggregationRule$loss.type, uniform = aggregationRule$uniform, 
        knots = aggregationRule$knots)
    }
    
    res$model <- aggregationRule$name
    res$residuals <- y - res$prediction
    res$call <- match.call()
    res$Y <- y
    res$experts <- experts

    class(res) <- "mixture"
    return(res)
  }