#' @export 
print.mixture <- function(x, ...) {
  cat("Aggregation rule: ")
  cat(x$model, "\n")
  cat("Loss function: ", x$loss.type$name, "loss", "\n")
  cat("Gradient trick: ", x$loss.gradient, "\n")
  cat("Coefficients:")
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

#' @export 
summary.mixture <- function(object, ...) {
  if (is.null(object$Y)) {
    K <- "Unknown"
    T <- 0
    TAB <- c("No losses yet")
  } else {
    T <- length(object$Y)
    K <- length(object$coefficients)
    
    rmse.algo <- sqrt(mean(loss(object$prediction, object$Y, loss.type = "square")))
    mape.algo <- mean(loss(object$prediction, object$Y, loss.type = "percentage"))
    rmse.unif <- sqrt(lossConv(rep(1/K, K), object$Y, object$experts, awake = object$awake))
    mape.unif <- lossConv(rep(1/K, K), object$Y, object$experts, awake = object$awake, 
                          loss.type = "percentage")
    
    TAB <- data.frame(rmse = c(rmse.algo, rmse.unif), mape = c(mape.algo, mape.unif))
    rownames(TAB) <- c(object$model, "Uniform")
  }
  
  res <- list(object = object, coefficients = object$coefficients, losses = TAB, 
              n.experts = K, n.observations = T)
  class(res) <- "summary.mixture"
  res
}

#' @export 
print.summary.mixture <- function(x, ...) {
  print(x$object)
  cat("\nNumber of experts: ", x$n.experts)
  cat("\nNumber of observations: ", x$n.observations, "\n\n")
  
  if (!is.null(dim(x$losses))) {
    print(signif(x$losses, digits = 3))
  }
}

#' @export 
plot.mixture <- function(x, pause = FALSE, losses = FALSE, ...) {
  
  x$experts <- data.frame(x$experts)
  if (is.null(names(x$experts))) {
    names(x$experts) <- colnames(x$experts)
  }
  
  op <- par(no.readonly = TRUE) # save default, for resetting...
  if (!losses) {
    if (!pause) {
      layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE))  
    }
    par(mar = c(3, 3, 1.6, 0.1), mgp = c(2, 0.5, 0))
    x$experts <- data.frame(x$experts)
    x$weights <- data.frame(x$weights)
    T <- nrow(x$experts)
    K <- ncol(x$experts)
    
    if (is.null(names(x$experts))) {
      names(x$experts) <- colnames(x$experts)
    }
    if (is.null(names(x$experts))) {
      names(x$experts) <- paste("Expert", 1:K)
    }
    names(x$weights) <- names(x$experts)
    l.names <- max(nchar(names(x$experts))) / 3 + 1.7
    
    if (x$model == "Ridge") {
      # Linear aggregation rule
      par(mar = c(3, 3, 1.6, l.names/2), mgp = c(1, 0.5, 0))
      matplot(x$weights, type = "l", xlab = "", ylab = "", lty = 1:5, 
              col = 2:(K+1), main = "Weights associated with the experts",...)
      mtext(side = 2, text = "Weights", line = 1.8, cex = 1)
      mtext(side = 1, text = "Time steps", line = 1.8, cex = 1)
      mtext(side = 4, text = names(x$experts), at = x$weights[T,], las = 2, col = 2:(K+1), cex= 0.5, line = 0.3)
    } else {
      # Convex aggregation rule
      par(mar = c(3, 3, 1.6, l.names/2), mgp = c(1, 0.5, 0))
      plot(c(1), type = "l", col = 1:8, lwd = 2, axes = F, xlim = c(1, T), ylim = c(0, 
                                                                                    1), ylab = "", xlab = "", main = "Weights associated with the experts")
      mtext(side = 2, text = "Weights", line = 1.8, cex = 1)
      mtext(side = 1, text = "Time steps", line = 1.8, cex = 1)
      x.idx <- c(1, 1:T, T:1)
      for (i in 1:K) {
        w.summed <- apply(matrix(as.matrix(x$weights[, i:K]), nrow = T), 1, sum)
        y.idx <- c(0, w.summed, rep(0, T))
        polygon(x = x.idx, y = y.idx, col = i + 1)
      }
      axis(1)
      axis(2)
      box()
      mtext(side = 4, text = names(x$experts), 
            at = (1-cumsum(c(x$weights[T,])))  + x$weights[T,]/2, las = 2, col = 2:(K+1), cex= 0.5, line = 0.3)
    }
    
    if (pause) { # break
      cat("Hit <Return> to see next plot:")
      plop <- readLines(n = 1)
    }
    
    # Box plot
    par(mar = c(l.names, 3, 1.6, 0.1))
    boxplot(x$weights, main = "Weights associated with the experts", col = 2:(K+1), axes = FALSE)
    mtext(side = 2, text = "Weights", line = 1.8, cex = 1)
    axis(1, at = 1:(K + 1), labels = FALSE)
    mtext(at = 1:K, text = names(x$weights), side = 1, las = 2, col = 2:(K+1), line = 0.8)
    axis(2)
    box()
    
    if (pause) { # break
      cat("Hit <Return> to see next plot:")
      plop <- readLines(n = 1)
    }
    
    #note: always pass alpha on the 0-255 scale
    makeTransparent<-function(someColor, alpha=100)
    {
      newColor<-col2rgb(someColor)
      apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
                                                  blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
    }
    
    # Cumulative loss
    par(mar = c(3, 3, 1.6, l.names/2), mgp = c(1, 0.5, 0))
    cumul.losses <- apply(loss(x$experts, x$Y, x$loss.type), 2, cumsum)
    matplot(cumul.losses, type = "l", lty = 1, xlab = "", ylab = "", 
            main = paste("Cumulative", x$loss.type$name, "loss"), col = makeTransparent(2:(K+1)))
    lines(cumsum(loss(x$prediction, x$Y, x$loss.type)), col = 1, lwd = 2)
    mtext(side = 2, text = "Cumulative loss", line = 1.8, cex = 1)
    mtext(side = 1, text = "Time steps", line = 1.8, cex = 1)
    mtext(side = 4, text = names(x$experts), at = cumul.losses[T,], las = 2, col = 2:(K+1), cex= 0.5, line = 0.3)
    legend("topleft", c("Experts", x$model), bty = "n", lty = 1, col = c("gray", 1), lwd = c(1,2))
    
    if (pause) { # break
      cat("Hit <Return> to see next plot:")
      plop <- readLines(n = 1)
    }
    
    # Cumulative residuals
    par(mar = c(3, 3, 1.6,l.names/2), mgp = c(1, 0.5, 0))
    cumul.residuals <- apply(x$Y - x$experts, 2, cumsum)
    matplot(cumul.residuals, type = "l", lty = 1, xlab = "", ylab = "", 
            main = paste("Cumulative residuals"), col = makeTransparent(2:(K+1)))
    lines(cumsum(x$Y - x$prediction), col = 1, lwd = 2)
    mtext(side = 2, text = "Cumulative residuals", line = 1.8, cex = 1)
    mtext(side = 1, text = "Time steps", line = 1.8, cex = 1)
    if (max(cumul.residuals) > abs(min(cumul.residuals))) {
      place = "topleft"
    } else {
      place = "bottomleft"
    }
    mtext(side = 4, text = names(x$experts), at = cumul.residuals[T,], las = 2, col = 2:(K+1), cex= 0.5, line = 0.3)
    legend(place, c("Experts", x$model), bty = "n", lty = 1, col = c("gray", 1), lwd = c(1,2))
  } else {
    x$loss.experts <- oracle(x$Y, x$experts, model = "expert", loss.type = x$loss.type)$loss.experts
    err.unif <- lossConv(rep(1/K, K), x$Y, x$experts, awake = x$awake, loss.type = x$loss.type)
    err.mixt <- x$loss
    idx.sorted <- order(c(x$loss.experts, err.unif, err.mixt))
    i.min <- 1
    my.col <- rep(1, K + 2)
    my.col[which(idx.sorted == K + 1)] <- 2
    my.col[which(idx.sorted == K + 2)] <- "forestgreen"
    my.col[which(idx.sorted != K + 1)[i.min]] <- 4
    
    par(mar = c(4.5, 4, 2, 2), mgp = c(3, 1, 0))
    plot(c(x$loss.experts, err.unif, err.mixt)[idx.sorted], xlab = "Experts", ylab = paste(x$loss.type$name, 
                                                                                 "loss"), main = "Average loss suffered by the experts", axes = F, pch = 3, 
         col = my.col, lwd = 2)
    axis(1, at = 1:(K + 2), labels = FALSE)
    mtext(at = 1:(K + 2), text = c(names(x$experts), "Uniform", x$model)[idx.sorted], 
          side = 1, las = 2, col = my.col, line = 0.8)
    axis(2)
    box()
  }
  par(op)
} 
