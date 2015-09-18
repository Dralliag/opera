#' @export 
print.mixture <- function(x, ...) {
  cat("Aggregation rule: ")
  cat(x$model, "\n")
  cat("Loss function: ", x$loss.type$name, "loss", "\n")
  cat("Gradient trick: ", x$loss.gradient, "\n")
  cat("Coefficients: ")
  print(x$coefficients)
}

#' @export 
summary.mixture <- function(object, ...) {
  if (is.null(object$Y)) {
    K <- "Unknown"
    T <- 0
    TAB <- c("No losses yet")
  } else {
    T <- length(Y)
    K <- length(object$parameters)

    rmse.algo <- rmse(object$prediction, object$Y)
    mape.algo <- mean(loss(object$prediction, object$Y, loss.type = "percentage"))
    rmse.unif <- sqrt(lossConv(rep(1/K,K),object$Y,object$experts, awake = object$awake))
    mape.unif <- lossConv(rep(1/K,K),object$Y,object$experts, awake = object$awake, loss.type = "percentage")

    TAB <- data.frame(rmse = c(rmse.algo, rmse.unif), mape = c(mape.algo, mape.unif))
    rownames(TAB) <- c(object$model,"Uniform")
  }

  res <- list(object = object, 
    coefficients = object$coefficients, 
    losses = TAB,
    n.experts = K,
    n.observations = T)
  class(res) <- "summary.mixture"
  res
}

#' @export 
print.summary.mixture <- function(x, ...) {
  print(x$object)
  cat("\nNumber of experts: ", x$n.experts)
  cat("\nNumber of observations: ", x$n.observations,"\n\n")

  if (!is.null(dim(x$losses))) {
    print(x$losses)
  }
}

#' @export 
plot.mixture <- function(x, ...) {
  x$experts <- data.frame(x$experts)
  
  K <- ncol(x$experts)
  if (is.null(names(x$experts))) {
    names(x$experts) <- colnames(x$experts)
  }
  if (is.null(names(x$experts))) {
    names(x$experts) <- paste("Expert", 1:K)
  }
  if (x$model == "gamMixture") {
    stop("Plot method is not yet available for gamMixture")
  }
  if (x$model == "Ridge") {
    # Linear aggregation rule
    matplot(x$weights, type = "l", xlab = "Time steps", ylab = "Weights", lty = 1:5, col = 1:8, ...)
    legend("topright", names(x$experts), lty = 1:5, col = 1:8, ...)
  } else {
    # Convex aggregation rule Mettre le plot en polygon des poids

    
  }
} 
