#' @export 
print.oracle <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  
  if (x$d > 1) {
    x$experts <- blockToSeries(x$experts)
  }
  if (x$model != "shifting") {
    cat("\nCoefficients:\n")
    x$coefficients <- data.frame(matrix(as.numeric(as.matrix(x$coefficients)), 
                                        nrow = 1))
    names(x$coefficients) <- colnames(x$experts)
    rownames(x$coefficients) <- ""
    print(signif(x$coefficients, digits = 3))
  }
  
  cat("\n")
  print(signif(summary(x)$losses, digits = 3))
}

#' @importFrom stats quantile
#' @export
summary.oracle <- function(object, ...) {
  
  if (object$d > 1) {
    object$experts <- blockToSeries(object$experts)
    object$Y <- blockToSeries(object$Y)
  }
  
  if (object$model == "expert") {
    
    T <- length(object$Y)
    K <- ncol(object$experts)
    
    rmse.algo <- sqrt(mean(loss(object$prediction, object$Y)))
    mape.algo <- mean(loss(object$prediction, object$Y, loss.type = "percentage"))
    rmse.unif <- sqrt(lossConv(rep(1/K, K), object$Y, object$experts, awake = object$awake))
    mape.unif <- lossConv(rep(1/K, K), object$Y, object$experts, awake = object$awake, 
                          loss.type = "percentage")
    
    
    TAB <- data.frame(rmse = c(rmse.algo, rmse.unif), mape = c(mape.algo, mape.unif))
    rownames(TAB) <- c(paste("Best", object$model, "oracle: "), "Uniform combination: ")
  }
  
  if (object$model == "linear" || object$model == "convex") {
    
    x <- summary(oracle(object$Y, object$experts, model = "expert", loss.type = object$loss.type, 
                 awake = object$awake))
    
    rmse.algo <- sqrt(mean(loss(object$prediction, object$Y)))
    mape.algo <- mean(loss(object$prediction, object$Y, loss.type = "percentage"))
    
    TAB.lin <- data.frame(rmse = rmse.algo, mape = mape.algo)
    rownames(TAB.lin) <- paste("Best", object$model, "oracle: ")
    TAB <- rbind(x$losses, TAB.lin)
    
  }
  
  if (object$model == "shifting") {
    
    K <- nrow(object$experts)
    n.shifts <- round(quantile(1:K))
    TAB <- matrix(object$loss[n.shifts], nrow = 1)
    colnames(TAB) <- paste(n.shifts - 1, "shifts")
    rownames(TAB) <- paste("Average", object$loss.type, "loss:")
    
    if (object$loss.type == "square") {
      TAB <- sqrt(TAB)
      rownames(TAB) <- "rmse:"
    } else if (object$loss.type == "absolute") {
      rownames(TAB) <- "mean absolute error:"
    } else if (object$loss.type == "percentage") {
      rownames(TAB) <- "mape:"
    }
  }
  
  res <- list(call = object$call, coefficients = object$coefficients, losses = TAB, 
              model = object$model)
  class(res) <- "summary.oracle"
  res
}

#' @export 
print.summary.oracle <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  
  if (x$model != "shifting") {
    cat("\nCoefficients:\n")
    print(signif(x$coefficients, digits = 3))
  }
  
  cat("\n")
  print(signif(x$losses, digits = 3))
}

#' @export
predict.oracle <- function(object, newexpert = NULL, ...) {
  if (missing(newexpert) || is.null(newexpert)) {
    stop("You should enter new expert predictions")
  }
  
  K <- length(object$coefficients)
  if (object$d == 1) {
    newexpert <- as.matrix(newexpert)
    if (ncol(newexpert) == 1 && nrow(newexpert) > 1) {
      newexpert <- t(newexpert)
    }
  }
  w <- matrix(object$coefficients, ncol = 1)
  pred <- newexpert %*% w
  
  return(pred)
} 
