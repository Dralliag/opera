#' Plot an aggregation procedure
#' @describeIn oracle \code{plot}. It has one optional arguments. 
#' @param x An object of class \code{oracle}. 
#' @param sort if set to TRUE (default), it sorts the experts by performance before the plots.
#' @param col colors
#' @importFrom graphics axis box mtext par plot
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

  if (is.null(col)) {
    if(!requireNamespace("RColorBrewer", quietly = TRUE)) {
      print("The RColorBrewer package must be installed to get better colors\n")
      col.palette <- 2:min((K+1),7)
    } else{
      col.palette <- RColorBrewer::brewer.pal(n = 3,name = "Set1")    
    }
    col <- col.palette
  }
  
  if (!is.null(x$names.experts)) {
     names(x$experts) <- x$names.experts
  } else {
    if (is.null(names(x$experts))) {
      names(x$experts) <- x$names.experts <- paste("X", 1:K,sep="")
    }
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