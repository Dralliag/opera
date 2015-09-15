#' @export 
print.mixture <- function(x, ...)
{
    cat("Call:\n")
    print(x$call)
    cat("\nCoefficients:\n")
    print(x$coefficients)
}

#' @export 
summary.mixture <- function(object, ...)
{
    rmse <- rmse(object$prediction, object$Y)
    mape <- mean(loss(object$prediction, object$Y,loss.type="percentage"))
    
    # mettre les rmse des oracles pour comparer
    TAB <- cbind(rmse = rmse, mape = mape)
    rownames(TAB) = "Average losses:"

    res <- list(call=object$call,
                coefficients = object$coefficients,
                losses=TAB)
    class(res) <- "summary.mixture"
  res 
}

#' @export 
print.summary.mixture <- function(x,...) 
{
    cat("Call:\n")
    print(x$call)
    cat("\n")

    cat("\nCoefficients:\n")
    print(x$coefficients)
    print(x$losses)
}

#' @export 
plot.mixture <- function(x, ...) 
{
  x$experts = data.frame(x$experts)

  K = ncol(x$experts)
  if (is.null(names(x$experts))) {names(x$experts) = colnames(x$experts)}
  if (is.null(names(x$experts))) {names(x$experts) = paste("Expert",1:K)}
  if (x$model == "gamMixture") {
  	stop("Plot method is not yet available for gamMixture")
  }
  if (x$model == "Ridge") { # Linear aggregation rule
    matplot(x$weights, type = 'l', xlab = "Time steps", ylab = "Weights", lty = 1:5, col =1:8, ...)
    legend("topright", names(x$experts), lty = 1:5, col = 1:8, ...)
  }
    else { # Convex aggregation rule
        # Mettre le plot en polygon des poids
    }
}

# #' @export
# predict.mixture <- function(object, newexpert = NULL, newY = NULL, 
# 														online = FALSE, type=c("model","response"),...)
# {
# 	if (object$model == "gamMixture") {
# 	  stop("Predict method is not yet available for gamMixture")
# 	}

# 	if (is.null(newexpert)) {
# 		prediction = NULL
# 		new.model = object
# 	}
# 	else {

# 		newexpert <- as.matrix(newexpert)
# 		if (ncol(newexpert) == 1 && nrow(newexpert) > 1 ) {
# 			newexpert = t(newexpert)
# 		}

# 		if (is.null(newY) || (!online)) {
# 			if (online) {
# 				warning("Non null observations newY is needed to perform online prediction.")
# 			}
# 			w <- matrix(object$coefficients, ncol = 1)	
# 	  	if (!online) {
# 	  		prediction = newexpert %*% w
# 	  	}
# 		}
# 	}
# 	return(prediction)
# }