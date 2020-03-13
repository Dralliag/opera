#' Plot an object of class mixture
#' 
#' provides different diagnostic plots for an aggregation procedure.
#' @param x an object of class mixture. If awake is provided (i.e., some experts are unactive), 
#' their residuals and cumulative losses are computed by using the predictions of the mixture.
#' @param pause if set to TRUE (default) displays the plots separately, otherwise on a single page
#' @param col the color to use to represent each experts, if set to NULL (default) use R\code{RColorBrewer::brewer.pal(...,"Spectral"}
#' @param ... additional plotting parameters
#' 
#' 
#' @return plots representing: plot of weights of each expert in function of time, boxplots of these weights,
#' cumulative loss \eqn{L_T=\sum_{t=1}^T l_{i,t}} of each expert in function of time, cumulative residuals \eqn{\sum_{t=1}^T (y_t-f_{i,t})} of each 
#' expert's forecast in function of time, average loss suffered by the experts and the contribution of each expert to the aggregation 
#' \eqn{p_{i,t}f_{i,t}} in function of time.
#' 
#' @author Pierre Gaillard <pierre@@gaillard.me>
#' @author Yannig  Goude <yannig.goude@edf.fr>
#' @seealso See \code{\link{opera-package}} and opera-vignette for a brief example about how to use the package.
#' @importFrom grDevices col2rgb rgb
#' @importFrom graphics axis box boxplot layout legend lines matplot mtext par plot polygon text
#' @importFrom stats lowess var
#' @export 
#' 
#'
plot.mixture <- function(x, pause = FALSE, col = NULL, ...) {
  def.par <- par(no.readonly = TRUE) # save default, for resetting...
  if (pause) par(ask=TRUE)
  x$experts <- data.frame(x$experts)
  K <- length(x$experts)
  w.order <- order(apply(x$weights,2,mean),decreasing = TRUE)
  
  if (is.null(col)) {
    if(!requireNamespace("RColorBrewer", quietly = TRUE)) {
      print("The RColorBrewer package must be installed to get better colors\n")
      col <- 2:min((K+1),7)
    } else{
      col <- rev(RColorBrewer::brewer.pal(n = max(min(K,11),4),name = "Spectral"))[1:min(K,11)]
    }
  }

    
  my.colors <- col
  
  col <- numeric(K)
  if (K <= length(my.colors)) {
    col[w.order] <- my.colors[1:K]
  } else {
    col[w.order] <- c(my.colors, rep(my.colors[length(my.colors)],K-length(my.colors)))
  }
  
  if (!pause) {
    layout(matrix(c(1,2,3,4,5,6),nrow = 3,ncol =  2, byrow = TRUE))  
  }
  par(mar = c(3, 3, 1.6, 0.1), mgp = c(2, 0.5, 0))
  x$experts <- data.frame(x$experts)
  x$Y <- c(t(x$Y))
  x$prediction <- c(t(x$prediction))
  x$weights <- data.frame(x$weights)
  T <- x$T
  d <- x$d
  
  if (!is.null(x$names.experts)) {
    names(x$weights) <- names(x$experts) <- x$names.experts
  } else {
    if (is.null(names(x$experts))) {
      names(x$weights) <- names(x$experts) <- x$names.experts <- paste("X", 1:K,sep="")
    }
  }
  l.names <- max(nchar(names(x$experts))) / 3 + 1.7
  
  if (x$model == "Ridge") {
    # Linear aggregation rule
    par(mar = c(3, 3, 2, l.names/2), mgp = c(1, 0.5, 0))
    
    matplot(x$weights, type = "l", xlab = "", ylab = "", lty = 1:5, main = "Weights associated with the experts", col = col,...)
    mtext(side = 2, text = "Weights", line = 1.8, cex = 1)
    # mtext(side = 1, text = "Time steps", line = 1.8, cex = 1)
    mtext(side = 4, text = names(x$experts), at = x$weights[T,], las = 2, col = col, cex= 0.5, line = 0.3)
  } else {
    # Convex aggregation rule
    par(mar = c(3, 3, 2, l.names/2), mgp = c(1, 0.5, 0))
    plot(c(1), type = "l", col = 1:8, lwd = 2, axes = F, xlim = c(1, T), ylim = c(0, 
                                                                                  1), ylab = "", xlab = "", main = "Weights associated with the experts")
    mtext(side = 2, text = "Weights", line = 1.8, cex = 1)
    # mtext(side = 1, text = "Time steps", line = 1.8, cex = 1)
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
      polygon(x = x.idx, y = y.idx, col = col[j], border=NA)
      w.summed.old <- w.summed
      w.summed <- w.summed - x$weights[,j]
      i.remaining[j] <- FALSE
      writeLegend(f = w.summed.old,w.summed,name = names(x$experts)[j])
    }
    axis(1)
    axis(2)
    box()
    names.toWrite <- names(x$experts)
    names.toWrite[w.order[-(1:min(K,15))]] <- ""
    mtext(side = 4, text = names.toWrite[i.order], 
          at = (1-cumsum(c(x$weights[T,i.order])))  + x$weights[T,i.order]/2, las = 2, col = col[i.order], cex= 0.5, line = 0.3)
  }
  
  
  # Box plot
  if (!is.null(x$awake)) {
    pond <- apply(x$awake[d*(1:T),],1,sum)
    normalized.weights <- x$weights * pond / (K*x$awake[d*(1:T),])
    normalized.weights[x$awake[d*(1:T),] == pond] <- NaN
  } else {
    normalized.weights <- x$weights 
  }
  
  i.order <- w.order[1:min(K,20)]
  par(mar = c(l.names, 3, 1.6, 0.1))
  boxplot(normalized.weights[,i.order], main = "Weights associated with the experts", col = col[i.order], axes = FALSE, pch='.')
  mtext(side = 2, text = "Weights", line = 1.8, cex = 1)
  axis(1, at = 1:(min(K,20)), labels = FALSE)
  mtext(at = 1:min(K,20), text = names(x$weights)[i.order], side = 1, las = 2, col = col[i.order], line = 0.8)
  axis(2)
  box()
  
  
  #note: always pass alpha on the 0-255 scale
  makeTransparent<-function(someColor, alpha=220)
  {
    newColor<-col2rgb(someColor)
    apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
                                                blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
  }
  
  # Cumulative loss
  par(mar = c(1.5, 3, 2.5, l.names/2), mgp = c(1, 0.5, 0))
  pred.experts <- (x$experts * x$awake + x$prediction * (1-x$awake))
  
  cumul.losses <- apply(loss(pred.experts, x$Y, x$loss.type), 2, cumsum)[seq(d,T*d,by=d),]
  cumul.exploss <- cumsum(loss(x$prediction, x$Y, x$loss.type))[seq(d,T*d,by=d)]
  
  matplot(cumul.losses, type = "l", lty = 1, xlab = "", ylab = "", 
          main = paste("Cumulative", x$loss.type$name, "loss"), col = makeTransparent(col), ylim = range(c(cumul.losses,cumul.exploss)))
  lines(cumul.exploss, col = 1, lwd = 2)
  mtext(side = 2, text = "Cumulative loss", line = 1.8, cex = 1)
  # mtext(side = 1, text = "Time steps", line = 1.8, cex = 1)
  mtext(side = 4, text = x$names.experts, at = cumul.losses[T,], las = 2, col = makeTransparent(col), cex= 0.5, line = 0.3)
  legend("topleft", c("Experts", x$model), bty = "n", lty = 1, col = c("gray", 1), lwd = c(1,2))
  
  
  # Cumulative residuals
  par(mar = c(1.5, 3, 2.5,l.names/2), mgp = c(1, 0.5, 0))
  cumul.residuals <- apply(x$Y - pred.experts, 2, cumsum)[seq(d,T*d,by=d),]
  cumul.expres <- cumsum(x$Y - x$prediction)[seq(d,T*d,by=d)]
  matplot(cumul.residuals, type = "l", lty = 1, xlab = "", ylab = "", 
          main = paste("Cumulative residuals"), col = makeTransparent(col), ylim = range(c(cumul.residuals,cumul.expres)))
  lines(cumul.expres, col = 1, lwd = 2)
  mtext(side = 2, text = "Cumulative residuals", line = 1.8, cex = 1)
  # mtext(side = 1, text = "Time steps", line = 1.8, cex = 1)
  if (max(cumul.residuals) > abs(min(cumul.residuals))) {
    place = "topleft"
  } else {
    place = "bottomleft"
  }
  mtext(side = 4, text = x$names.experts, at = cumul.residuals[T,], las = 2, col = col, cex= 0.5, line = 0.3)
  legend(place, c("Experts", x$model), bty = "n", lty = 1, col = c("gray", 1), lwd = c(1,2))
  
  #losses
  l.names <- max(max(nchar(names(x$experts))) / 3 + 1.7,4)
  x$loss.experts <- apply(loss(x = pred.experts,y = x$Y,loss.type = x$loss.type),2,mean)
  err.unif <- lossConv(rep(1/K, K), x$Y, x$experts, awake = x$awake, loss.type = x$loss.type)
  err.mixt <- x$loss
  idx.sorted <- order(c(x$loss.experts, err.unif, err.mixt))
  my.col <- c(col,1,1)[idx.sorted]
  my.pch <- c(rep(20, K),8,8)[idx.sorted]
  
  par(mar = c(l.names, 3, 2.5,l.names/2), mgp = c(1, 0.5, 0))
  plot(c(x$loss.experts, err.unif, err.mixt)[idx.sorted], xlab = "", ylab = "", main = "Average loss suffered by the experts", axes = F, pch = my.pch, 
       col = my.col, lwd = 2,type='b')
  mtext(side = 2, text = paste(x$loss.type$name,"loss"), line = 1.8, cex = 1)
  axis(1, at = 1:(K + 2), labels = FALSE)
  mtext(at = 1:(K + 2), text = c(names(x$experts), "Uniform", x$model)[idx.sorted], 
        side = 1, las = 2, col = my.col, line = 0.8,cex = .7)
  axis(2)
  box()
  
  # cumulative plot of the series
  par(mar = c(2, 3, 2.5,l.names/2), mgp = c(1, 0.5, 0))
  if (x$d ==1) {
    cumulativePlot(W = x$weights,X = x$experts, Y = x$Y,smooth = TRUE,alpha = 0.01,plot.Y = TRUE, col.pal = rev(my.colors))
  } else {
    X <- apply(seriesToBlock(X = x$experts,d = x$d),c(1,3),mean)
    Y <- apply(seriesToBlock(x$Y,d = x$d),1,mean)
    colnames(X) <- x$names.experts
    cumulativePlot(W = x$weights,X = X, Y = Y,smooth = TRUE,alpha = 0.01,plot.Y = TRUE, col.pal = col)    
  }
  par(def.par)
} 



writeLegend <- function(f,g,name,Y.lim=c(0,1), ...) {
  tau = Y.lim[2]/20
  Tab = matrix(0,ncol = 2, nrow = 100)
  y.seq <- seq(Y.lim[1],Y.lim[2],length.out = 100)
  for (i in 1:100) {
    x = y.seq[i]
    sel = which(g < x & f > x + tau)
    temp <- cumsum(c(1, diff(sel) - 1))
    temp2 <- rle(temp)
    Tab[i,1] <- max(temp2$lengths)
    Tab[i,2] <- sel[which(temp == with(temp2, values[which.max(lengths)]))][1]
  }
  id = which.max(Tab[,1])
  x <- y.seq[id]
  l <- Tab[id,1]
  v <- Tab[id,2]
  if (l > length(f)/20){
    j = floor(60 *l/length(f))
    text(v+l/2,x+tau/2,substr(name,1,j),cex = 0.8,...)
  }
}

cumulativePlot<-function(W,X,Y,col.pal=NULL, smooth = FALSE, plot.Y = FALSE, alpha = 0.1)
{
  time<-c(1:nrow(X))
  active.experts<-which(colMeans(W)>0)
  W<-W[,active.experts]  
  X<-X[,active.experts]
  
  K <- ncol(X)
  
  
  if(is.null(col.pal)) col.pal <- RColorBrewer::brewer.pal(n = min(K,9),name = "Spectral")
  if (length(col.pal) < K) col.pal <- c(rep(col.pal[1],K-length(col.pal)),col.pal)
  
  
  o<-order(colSums(W),decreasing = F)
  mat<-W[,o]*X[,o]
  Agg<-apply(mat,1,sum)
  colnames(mat)<-colnames(X)[o]
  
  if (!smooth)Y.lim = range(Agg,Y,mat)
  if (smooth) 
  {
    y.lo<-lowess(x = time,y = Y,f = alpha)$y
    Agg.lo<-lowess(x = time,y = Agg,f = alpha)$y
    
    mat.lo<-apply(mat,2,function(z){lowess(x = time,y = z,f = alpha)$y})
    Y.lim = range(Agg.lo,mat.lo)
  }
  
  
  plot(x = NULL,y = NULL,col=col.pal[1], type='l', xaxt='n',ylim=Y.lim,lty='dotted',
       yaxt='n',xlab="",ylab="",lwd=3,xlim = range(time),
       main = paste("Contribution of each expert to prediction"))
  y.summed <- Agg
  for(i in rev(c(1:ncol(mat))))
  {
    if (!smooth) addPoly(time,y.summed,col=col.pal[i])
    if (smooth) addPoly(time,lowess(y.summed,f = alpha)$y,col=col.pal[i])
    y.summed.old <- y.summed
    y.summed <- y.summed - mat[,i]
    if (!smooth) writeLegend(f=y.summed.old,g= y.summed, name = colnames(mat)[i],Y.lim,col='black')
    if (smooth) writeLegend(f=lowess(y.summed.old,f=alpha/10)$y,g=lowess(y.summed,f=alpha/10)$y, name = colnames(mat)[i],Y.lim,col='black')
  }
  if (plot.Y && !smooth) lines(time,Y,col=1,lwd=2,lty='dotted')
  if (plot.Y && smooth) lines(lowess(x = time,y = Y,f = alpha)$y,col=1,lwd=2,lty='dotted')
  axis(1)
  axis(2)
}


addPoly<-function(x,y,col)
{
  xx <- c(x, rev(x))
  yy <- c(rep(0, length(x)), rev(y))
  polygon(xx, yy, col=col, border=NA)
}
