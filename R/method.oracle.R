#' Print oracle method
#' @method print oracle
#' @export print.oracle
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
  

#' @method summary oracle
#' @export summary.oracle
summary.oracle <- function(x, ...)
{

    if (x$model != "shifting") {   
        rmse <- rmse(x$prediction, x$Y)
        mape <- mean(loss(x$prediction, x$Y,loss.type="percentage"))
    
        # mettre les rmse des oracles pour comparer
        TAB <- cbind(rmse = rmse, mape = mape)
        rownames(TAB) = paste("Best",x$model,"oracle: ")

    } else {

        K = nrow(x$experts)
        n.shifts = round(quantile(1:K))
        TAB = matrix(m$loss[n.shifts], nrow = 1)
        colnames(TAB) = paste(n.shifts-1,"shifts")
        rownames(TAB) = paste("Average", x$loss.type, "loss:")

        if (x$loss.type == "square") {
            TAB = sqrt(TAB)
            rownames(TAB) = "rmse:"
        } else if (x$loss.type == "absolute") {
            rownames(TAB) = "mean absolute error:"
        } else if (x$loss.type == "percentage") {
            rownames(TAB) = "mape:"
        }
    }

    res <- list(call=x$call,
                coefficients = x$coefficients,
                losses=TAB,
                model = x$model)
    class(res) <- "summary.oracle"
	res 
}

#' @method print sumary.oracle
#' @export print.summary.oracle
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

#' @method plot oracle
#' @export plot.oracle
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

