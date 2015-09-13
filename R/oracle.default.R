#' Compute an oracle (method)
#' @inheritParams oracle
#' @export oracle.default

oracle.default <-
  function(y, experts, model = "convex", awake = NULL, ...)
  {
    if (is.character(model)) {model = list(name = model)}
    if (is.null(model$loss.type)) {model$loss.type = "square"}
    if (!is.null(model$tau) && model$loss.type != "pinball") {
      warning("Unused parameter tau (loss.type != 'pinball')")
    }
    if (is.null(model$tau)) {model$tau = 0.5}
    if (!is.null(model$lambda) && model$name != "linear") {
      warning("Unused lambda parameter (model != linear)")}
    if (is.null(model$lambda)) model$lambda = 0
    if (!is.null(model$niter) && model$name!= "convex" && model$name != "linear") {
      warning("Unused niter parameter (model should be convex or linear)")} 
    if (is.null(model$niter)) model$niter = 3
    if ((!is.null(awake) || sum(is.na(experts)>0)) && model$name != "convex" && model$name != "shifting") {
      stop(paste("Sleeping or missing values not allowed for best", model$name, "oracle."))}  

    if (!(model$name %in% c("convex","linear","shifting","expert"))) {
      stop("Wrong model specification") 
    }
    if (min(y) <= 0 && model$loss.type == "percentage") {
      stop("Y should be non-negative for percentage loss function")
    }
    # if we are looking for the best convex combination of experts
    if (model$name == "convex") {
      res <- bestConvex(
          y, experts, awake = awake,
          loss.type = model$loss.type, niter = model$niter, 
          tau = model$tau,...)
    }
    
    if (model$name == "linear") {
      res <- bestLinear(y, experts, lambda = model$lambda, 
              loss.type = model$loss.type, tau = model$tau)
    }
    
    if (model$name == "shifting") {
      res <- bestShifts(y, experts, awake = model$awake,
              loss.type = model$loss.type)
    }
    
    if (model$name == "expert") {
      
      loss.experts = apply(apply(experts, 2, function (x) {
        loss(x,y,loss.type = model$loss.type,tau = model$tau)}), 2, mean)
      best.loss = min(loss.experts)
      coefficients =  (loss.experts == best.loss) / sum(loss.experts == best.loss)
      best.expert = which(coefficients > 0)[1]
      res = list(loss = best.loss,
                 coefficients = coefficients, 
                 prediction = experts[,best.expert])
      if (model$loss.type == "square") {
        res$rmse = sqrt(res$loss)
      }
    }
    
    res$model <- model$name
    res$residuals <- y - res$prediction
    res$call <- match.call()
    res$Y <- y
    res$experts <- experts
    res$loss.type <- model$loss.type
    
    class(res) <- "oracle"
    return(res)
  }