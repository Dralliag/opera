#' Plot an aggregation procedure
#' @describeIn mixture \code{plot}. It has two optional arguments. 
#' The argument pause = TRUE displays the plots separately.
#' The argument losses = TRUE prints only the performance achieved by the aggregation procedures and the experts.
#' @export 
plot.mixture <- function(x, pause = FALSE, losses = FALSE, col = NULL, ...) {
  def.par <- par(no.readonly = TRUE) # save default, for resetting...
  x$experts <- data.frame(x$experts)
  K <- length(x$experts)
  w.order <- order(apply(x$weights,2,mean),decreasing = TRUE)
  
  if (is.null(col)) {
    if (K <= 9) {
      col <- RColorBrewer::brewer.pal(n = K,name = "Set1")  
    } else {
      c0 <- RColorBrewer::brewer.pal(n = 9,name = "Set1") 
      c1 <- RColorBrewer::brewer.pal(n = 12,name = "Set3")  
      c2 <- RColorBrewer::brewer.pal(n = 8,name = "Dark2") 
      c3 <- RColorBrewer::brewer.pal(n = 12,name = "Paired")
      
      my.colors <- c(c0,c1,c2,c3)
      col <- numeric(K)
      col[w.order] <- rep(my.colors,ceiling(K/length(my.colors)))[1:K]
    }
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
      names(x$experts) <- paste("X", 1:K,sep="")
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
      w.summed <- rep(1,T)
      i.remaining = rep(TRUE,K)
      i.order <- rep(0,K)
      for (i in 1:K) {
        if (i <K){
          j <- which(i.remaining)[which.min(apply(x$weights[,i.remaining],2,function(p){sqrt(var(w.summed-p))}))]
        } else {
          j <- which(i.remaining)
        }
        i.order[i] <- j
        y.idx <- c(0, w.summed, rep(0, T))
        polygon(x = x.idx, y = y.idx, col = col[j])
        w.summed.old <- w.summed
        w.summed <- w.summed - x$weights[,j]
        i.remaining[j] <- FALSE
        writeLegend(f = w.summed.old,w.summed,tau = 0.05,name = names(x$experts)[j])
      }
      axis(1)
      axis(2)
      box()
      mtext(side = 4, text = names(x$experts)[i.order], 
            at = (1-cumsum(c(x$weights[T,i.order])))  + x$weights[T,i.order]/2, las = 2, col = col[i.order], cex= 0.5, line = 0.3)
    }
    
    if (pause) { # break
      cat("Hit <Return> to see next plot:")
      plop <- readLines(n = 1)
    }
    
    # Box plot
    i.order <- w.order[1:min(K,20)]
    par(mar = c(l.names, 3, 1.6, 0.1))
    boxplot(x$weights[,i.order], main = "Weights associated with the experts", col = col[i.order], axes = FALSE)
    mtext(side = 2, text = "Weights", line = 1.8, cex = 1)
    axis(1, at = 1:(min(K,20)), labels = FALSE)
    mtext(at = 1:min(K,20), text = names(x$weights)[i.order], side = 1, las = 2, col = col[i.order], line = 0.8)
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

writeLegend <- function(f,g, tau, name) {
  Tab = matrix(0,ncol = 2, nrow = 100)
  for (i in 1:100) {
    x = 0.01 * i
    sel = which(g < x & f > x + tau)
    temp <- cumsum(c(1, diff(sel) - 1))
    temp2 <- rle(temp)
    Tab[i,1] <- max(temp2$lengths)
    Tab[i,2] <- sel[which(temp == with(temp2, values[which.max(lengths)]))][1]
  }
  id = which.max(Tab[,1])
  x <- id * 0.01
  l <- Tab[id,1]
  v <- Tab[id,2]
  if (l > length(f)/20){
    j = floor(60 *l/length(f))
    text(v+l/2,x+tau/2,substr(name,1,j),cex = 0.8,)
  }
}