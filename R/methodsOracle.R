#' @export 
print.oracle <- function(x, ...)
{
    cat("Call:\n")
    print(x$call)

    if (x$model != "shifting") {
        cat("\nCoefficients:\n")
        print(x$coefficients)   
    }

    cat('\n')
    print(summary(x)$losses)
}
  
#' @export
summary.oracle <- function(object, ...)
{

    if (object$model != "shifting") {   
        rmse <- rmse(object$prediction, object$Y)
        mape <- mean(loss(object$prediction, object$Y,loss.type="percentage"))
    
        # mettre les rmse des oracles pour comparer
        TAB <- cbind(rmse = rmse, mape = mape)
        rownames(TAB) = paste("Best",object$model,"oracle: ")

    } else {

        K = nrow(object$experts)
        n.shifts = round(quantile(1:K))
        TAB = matrix(object$loss[n.shifts], nrow = 1)
        colnames(TAB) = paste(n.shifts-1,"shifts")
        rownames(TAB) = paste("Average", object$loss.type, "loss:")

        if (object$loss.type == "square") {
            TAB = sqrt(TAB)
            rownames(TAB) = "rmse:"
        } else if (object$loss.type == "absolute") {
            rownames(TAB) = "mean absolute error:"
        } else if (object$loss.type == "percentage") {
            rownames(TAB) = "mape:"
        }
    }

    res <- list(call=object$call,
                coefficients = object$coefficients,
                losses=TAB,
                model = object$model)
    class(res) <- "summary.oracle"
	res 
}

#' @export 
print.summary.oracle <- function(x,...) 
{
    cat("Call:\n")
    print(x$call)

    if (x$model != "shifting") {
        cat("\nCoefficients:\n")
        print(x$coefficients)   
    }

    cat('\n')
    print(x$losses)
}

#' @export 
plot.oracle <- function(x, ...) 
{
    if (x$model != "shifting") {
        stop("No plot method for", x$model, "oracle.")
    } else {
        L = x$loss
        if (x$loss.type == "square") {
            L = sqrt(x$loss)
            y.lab = "rmse"
        } else if (x$loss.type == "absolute") {
            y.lab = "mean absolute error"
        } else if (x$loss.type == "percentage") {
            y.lab = "mape"
        }
        plot(0:(length(L)-1), L, xlab = "Number of shifts", ylab = y.lab, type = "o", pch = 20, cex = .6, lwd = 2,
            main = "Error suffered by the shifting oracle")
    }
}

