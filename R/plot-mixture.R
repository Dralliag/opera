#' Plot an aggregation procedure
#' @describeIn mixture \code{plot}. It has two optional arguments. 
#' The argument pause = TRUE displays the plots separately.
#' The argument losses = TRUE prints only the performance achieved by the aggregation procedures and the experts.
#' @export 
plot.mixture <- function(x, pause = FALSE, losses = FALSE, col = NULL, ...) {
  def.par <- par(no.readonly = TRUE) # save default, for resetting...
  x$experts <- data.frame(x$experts)
  K <- length(x$experts)
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
  
  if (!losses) {
    if (!pause) {
      layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE))  
    }
    par(mar = c(3, 3, 1.6, 0.1), mgp = c(2, 0.5, 0))
    x$experts <- data.frame(x$experts)
    x$Y <- c(t(x$Y))
    x$prediction <- c(t(x$prediction))
    x$weights <- data.frame(x$weights)
    T <- x$T
    d <- x$d
    
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
      
      matplot(x$weights, type = "l", xlab = "", ylab = "", lty = 1:5, main = "Weights associated with the experts", col = col,...)
      mtext(side = 2, text = "Weights", line = 1.8, cex = 1)
      mtext(side = 1, text = "Time steps", line = 1.8, cex = 1)
      mtext(side = 4, text = names(x$experts), at = x$weights[T,], las = 2, col = col, cex= 0.5, line = 0.3)
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
        polygon(x = x.idx, y = y.idx, col = col[i])
      }
      axis(1)
      axis(2)
      box()
      mtext(side = 4, text = names(x$experts), 
            at = (1-cumsum(c(x$weights[T,])))  + x$weights[T,]/2, las = 2, col = col, cex= 0.5, line = 0.3)
    }
    
    if (pause) { # break
      cat("Hit <Return> to see next plot:")
      plop <- readLines(n = 1)
    }
    
    # Box plot
    par(mar = c(l.names, 3, 1.6, 0.1))
    boxplot(x$weights, main = "Weights associated with the experts", col = col, axes = FALSE)
    mtext(side = 2, text = "Weights", line = 1.8, cex = 1)
    axis(1, at = 1:(K + 1), labels = FALSE)
    mtext(at = 1:K, text = names(x$weights), side = 1, las = 2, col = col, line = 0.8)
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
    cumul.losses <- apply(loss(x$experts, x$Y, x$loss.type), 2, cumsum)[seq(d,T*d,by=d),]
    cumul.exploss <- cumsum(loss(x$prediction, x$Y, x$loss.type))[seq(d,T*d,by=d)]
    matplot(cumul.losses, type = "l", lty = 1, xlab = "", ylab = "", 
            main = paste("Cumulative", x$loss.type$name, "loss"), col = makeTransparent(col), ylim = range(c(cumul.losses,cumul.exploss)))
    lines(cumul.exploss, col = 1, lwd = 2)
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
    cumul.residuals <- apply(x$Y - x$experts, 2, cumsum)[seq(d,T*d,by=d),]
    cumul.expres <- cumsum(x$Y - x$prediction)[seq(d,T*d,by=d)]
    matplot(cumul.residuals, type = "l", lty = 1, xlab = "", ylab = "", 
            main = paste("Cumulative residuals"), col = makeTransparent(col), ylim = range(c(cumul.residuals,cumul.expres)))
    lines(cumul.expres, col = 1, lwd = 2)
    mtext(side = 2, text = "Cumulative residuals", line = 1.8, cex = 1)
    mtext(side = 1, text = "Time steps", line = 1.8, cex = 1)
    if (max(cumul.residuals) > abs(min(cumul.residuals))) {
      place = "topleft"
    } else {
      place = "bottomleft"
    }
    mtext(side = 4, text = names(x$experts), at = cumul.residuals[T,], las = 2, col = col, cex= 0.5, line = 0.3)
    legend(place, c("Experts", x$model), bty = "n", lty = 1, col = c("gray", 1), lwd = c(1,2))
  } else {
    x$loss.experts <- oracle(x$Y, x$experts, model = "expert", loss.type = x$loss.type)$loss.experts
    err.unif <- lossConv(rep(1/K, K), x$Y, x$experts, awake = x$awake, loss.type = x$loss.type)
    err.mixt <- x$loss
    idx.sorted <- order(c(x$loss.experts, err.unif, err.mixt))
    i.min <- 1
    my.col <- rep(1, K + 2)
    my.col[which(idx.sorted == K + 1)] <- col[1]
    my.col[which(idx.sorted == K + 2)] <- col[2]
    my.col[which(idx.sorted != K + 1)[i.min]] <- col[3]
    
    par(mar = c(4.5, 4, 2, 2))
    plot(c(x$loss.experts, err.unif, err.mixt)[idx.sorted], xlab = "", ylab = paste(x$loss.type$name, 
                                                                                    "loss"), main = "Average loss suffered by the experts", axes = F, pch = 3, 
         col = my.col, lwd = 2)
    axis(1, at = 1:(K + 2), labels = FALSE)
    mtext(at = 1:(K + 2), text = c(names(x$experts), "Uniform", x$model)[idx.sorted], 
          side = 1, las = 2, col = my.col, line = 0.8)
    axis(2)
    box()
  }
  par(def.par)
} 