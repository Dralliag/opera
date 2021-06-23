#' Print an aggregation procedure
#' 
#' @param x An object of class mixture
#' 
#' @export 
#' 
#' @rdname mixture-opera
print.mixture <- function(x, ...) {
  cat("Aggregation rule: ")
  cat(x$model, "\n")
  cat("Loss function: ", x$loss.type$name, "loss", "\n")
  cat("Gradient trick: ", x$loss.gradient, "\n")
  cat("Coefficients: ")
  if (x$coefficients[1] != "Uniform") {
    cat("\n")
    x$coefficients <- data.frame(signif(matrix(as.matrix(x$coefficients), nrow = 1)))
    names(x$coefficients) <- colnames(x$experts)
    rownames(x$coefficients) <- ""
    print(signif(x$coefficients, digits = 3))
  } else {
    print("Uniform mixture")
  }
}