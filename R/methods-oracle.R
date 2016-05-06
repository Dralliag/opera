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
    
    x <- summary(oracle(object$Y, object$experts, model = "expert", loss.type = object$loss.type), 
                 awake = object$awake)
    
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
plot.oracle <- function(x, sort = TRUE, col = NULL, ...) {
  def.par <- par(no.readonly = TRUE)
  
  if (x$d > 1) {
    x$experts <- blockToSeries(x$experts)
    x$Y <- blockToSeries(x$Y)
  }
  x$experts <- data.frame(x$experts)
  
  T <- nrow(x$experts)
  K <- ncol(x$experts)
  if (K <= 9) {
    col.palette <- RColorBrewer::brewer.pal(n = K,name = "Set1")  
  } else {
    col.palette <- RColorBrewer::brewer.pal(n = K,name = "Paired")  
  }
  if (is.null(col)) {
    col <- col.palette
  }
  
  
  if (is.null(names(x$experts))) {
    names(x$experts) <- colnames(x$experts)
  }
  
  if (x$model == "expert") {
    err.unif <- lossConv(rep(1/K, K), x$Y, x$experts, awake = x$awake, loss.type = x$loss.type)
    if (sort) {
      idx.sorted <- order(c(x$loss.experts, err.unif))
      i.min <- 1
    } else {
      idx.sorted = c(K+1,1:K)
      i.min <- order(x$loss.experts)[1]
    }
    my.col <- rep(1, K + 1)
    my.col[which(idx.sorted == K + 1)] <- col[1]
    my.col[which(idx.sorted != K + 1)[i.min]] <- col[2]
    
    par(mar = c(4.5, 4, 2, 2))
    plot(c(x$loss.experts, err.unif)[idx.sorted], xlab = "", ylab = paste(x$loss.type$name, 
                                                                          "loss"), main = "Average loss suffered by the experts", axes = F, pch = 3, 
         col = my.col, lwd = 2,...)
    axis(1, at = 1:(K + 1), labels = FALSE)
    mtext(at = 1:(K + 1), text = c(names(x$experts), "Uniform")[idx.sorted], 
          side = 1, las = 2, col = my.col, line = 0.8)
    axis(2)
    box()
    
  }
  
  if (x$model == "convex" || x$model == "linear") {
    err.unif <- lossConv(rep(1/K, K), x$Y, x$experts, awake = x$awake, loss.type = x$loss.type)
    idx.sorted <- order(c(x$loss.experts, err.unif, x$loss))
    my.col <- rep(1, K + 2)
    my.col[which(idx.sorted == K + 1)] <- col[1]
    my.col[which(idx.sorted == K + 2)] <- col[2]
    my.col[which(!(idx.sorted %in% c(K + 1, K + 2)))[1]] <- col[3]
    y.max <- c(x$loss.experts, err.unif, x$loss)[idx.sorted]
    
    par(mar = c(4.5, 4, 2, 2))
    plot(c(x$loss.experts, err.unif, x$loss)[idx.sorted], xlab = "", ylab = paste(x$loss.type$name, 
                                                                                  "loss"), main = "Average loss suffered by the experts", axes = F, pch = 3, 
         col = my.col, lwd = 2)
    axis(1, at = 1:(K + 2), labels = FALSE)
    mtext(at = 1:(K + 2), text = c(names(x$experts), "Uniform", 
                                   x$model)[idx.sorted], 
          side = 1, las = 2, col = my.col, line = 0.8)
    axis(2)
    box()
  }
  
  if (x$model == "shifting") {
    L <- x$loss
    if (x$loss.type == "square") {
      L <- sqrt(x$loss)
      y.lab <- "rmse"
    } else if (x$loss.type == "absolute") {
      y.lab <- "mean absolute error"
    } else if (x$loss.type == "percentage") {
      y.lab <- "mape"
    }
    
    plot(0:(length(L) - 1), L, xlab = "Number of shifts", ylab = y.lab, type = "o", 
         pch = 20, cex = 0.6, lwd = 2, main = "Error suffered by the shifting oracle")
  }
  par(def.par)
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
