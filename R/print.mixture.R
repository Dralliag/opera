#' Print mixture method
#' @method print mixture
#' @export print.mixture
print.mixture <- function(x, ...)
{
    cat("Call:\n")
    print(x$call)
    cat("\nCoefficients:\n")
    print(x$coefficients)
}

#' @method summary mixture
#' @export summary.mixture
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


#' @method print sumary.mixture
#' @export print.summary.mixture
print.summary.mixture <- function(x,...) 
{
    cat("Call:\n")
    print(x$call)
    cat("\n")

    cat("\nCoefficients:\n")
    print(x$coefficients)
    print(x$losses)
}



#' @method plot mixture
#' @export plot.mixture
plot.mixture <- function(object, ...) 
{
    object$experts = data.frame(object$experts)

    K = ncol(object$experts)
    if (is.null(names(object$experts))) {names(object$experts) = colnames(object$experts)}
    if (is.null(names(object$experts))) {names(object$experts) = paste("Expert",1:K)}
	if (object$model == "gamMixture") {

	}
	if (object$model == "Ridge") { # Linear aggregation rule
		matplot(object$weights, type = 'l', xlab = "Time steps", ylab = "Weights", lty = 1:5, col =1:8, ...)
		legend("topright", names(object$experts), lty = 1:5, col = 1:8, ...)
	}
    else { # Convex aggregation rule
        # Mettre le plot en polygon des poids
    }
}


