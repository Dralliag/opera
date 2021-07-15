#' Summary of an aggregation procedure
#' 
#' @param object An object of class mixture
#' @param ... Additional parameters
#' 
#' @export 
#' 
#' @rdname mixture-opera
#' 
summary.mixture <- function(object, ...) {
  if (is.null(object$Y)) {
    K <- "Unknown"
    T <- 0
    d <- "Unknown"
    TAB <- c("No losses yet")
  } else {
    T <- object$T
    K <- length(object$coefficients)
    d <- object$d
    
    rmse.algo <- sqrt(mean(loss(x = c(object$prediction), y = c(object$Y), loss.type = "square")))
    mape.algo <- mean(loss(x = c(object$prediction), y = c(object$Y), loss.type = "percentage"))
    rmse.unif <- sqrt(lossConv(rep(1/K, K), c(t(object$Y)), object$experts, awake = object$awake))
    mape.unif <- lossConv(rep(1/K, K), c(t(object$Y)), object$experts, awake = object$awake, 
                          loss.type = list(name = "percentage"))
    
    TAB <- data.frame(rmse = c(rmse.algo, rmse.unif), mape = c(mape.algo, mape.unif))
    rownames(TAB) <- c(object$model, "Uniform")
  }
  
  res <- list(object = object, coefficients = object$coefficients, losses = TAB, 
              n.experts = K, n.observations = T, n.dimension = d)
  class(res) <- "summary.mixture"
  res
}

#' @export 
print.summary.mixture <- function(x, ...) {
  print(x$object)
  cat("\nNumber of experts: ", x$n.experts)
  cat("\nNumber of observations: ", x$n.observations)
  cat("\nDimension of the data: ", x$n.dimension, "\n\n")
  
  if (!is.null(dim(x$losses))) {
    print(signif(x$losses, digits = 3))
  }
}


