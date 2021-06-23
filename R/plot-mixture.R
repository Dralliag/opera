#' Plot an object of class mixture
#' 
#' provides different diagnostic plots for an aggregation procedure.
#' 
#' @param x an object of class mixture. If awake is provided (i.e., some experts are unactive), 
#' their residuals and cumulative losses are computed by using the predictions of the mixture.
#' @param pause if set to TRUE (default) displays the plots separately, otherwise on a single page
#' @param col the color to use to represent each experts, if set to NULL (default) use R\code{RColorBrewer::brewer.pal(...,"Spectral"}
#' @param alpha \code{numeric}. Smoothing parameter for contribution plot (parameter 'f' of function \code{\link[stats]{lowess}}).
#' @param dynamic \code{boolean}. If TRUE, graphs are generated with \code{rAmCharts}, else with base R.
#' @param type \code{char}.
#' \itemize{
#'      \item{'all'}{ Display all the graphs ;}
#'      \item{'plot_weight', 'boxplot_weight', 'dyn_avg_loss', 'cumul_res', 'avg_loss', 'contrib'}{ Display the selected graph alone.}
#' }
#' @param max_experts \code{integer}. Maximum number of experts to be displayed (only the more influencial).
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
#' 
#' @seealso See \code{\link{opera-package}} and opera-vignette for a brief example about how to use the package.
#' 
#' @importFrom grDevices col2rgb rgb
#' @importFrom graphics axis box boxplot layout legend lines matplot mtext par plot polygon text
#' @importFrom stats lowess var
#' @importFrom htmltools browsable tagList
#' 
#' 
#' @export 
#' 
#'
plot.mixture <- function(x, 
                         pause = FALSE, 
                         col = NULL, 
                         alpha = 0.01,
                         dynamic = T, 
                         type = c('all', 'plot_weight', 'boxplot_weight', 
                                  'dyn_avg_loss', 'cumul_res', 
                                  'avg_loss', 'contrib'
                         ), 
                         max_experts = 50,
                         ...) {
  
  type <- tryCatch({
    match.arg(type)
  }, error = function(e){
    warning("Invalid 'type' argument. Set to 'all'")
    'all'
  })
  ############# add checks on x$loss

  
  def.par <- par(no.readonly = TRUE) # save default, for resetting...
  if (pause) par(ask=TRUE)
  K <- ncol(x$experts)
  w.order <- order(apply(x$weights,2,mean),decreasing = FALSE)
  
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
  
  if (!pause && type == "all") {
    layout(matrix(c(1,2,3,4,5,6),nrow = 3,ncol =  2, byrow = TRUE))  
  }
  
  x$Y <- c(t(x$Y))
  x$prediction <- c(t(x$prediction))
  x$weights <- data.frame(x$weights)
  T <- x$T
  d <- x$d
  
  if (!is.null(x$names.experts)) {
    names(x$weights) <- colnames(x$experts) <- x$names.experts
  } else {
    if (is.null(colnames(x$experts))) {
      names(x$weights) <- colnames(x$experts) <- x$names.experts <- paste("X", 1:K,sep="")
    } 
  }
  
  if (ncol(x$weights) > max_experts + 2) {
    l.names <- max(nchar(c(colnames(x$experts), "worst_others"))) / 3 + 1.7
    col[1:ncol(x$weight) <= ncol(x$weights) - max_experts] <- "grey"
    
  } else {
    l.names <- max(nchar(colnames(x$experts))) / 3 + 1.7
  }
  
  x$weights <- x$weights[, w.order]
  x$experts <- x$experts[, w.order]
  
  if (dynamic) {
    list_plt <- list()
  } else {
    par(mar = c(3, 3, 1.6, 0.1), mgp = c(2, 0.5, 0))
  }
  
  if (x$model == "Ridge" && (type == "all" || type == "plot_weight")) {
    # Linear aggregation rule
    if (! dynamic) {
      if (type == "all") {
        par(mar = c(3, 3, 2, l.names/2), mgp = c(1, 0.5, 0)) 
      }
      
      if (ncol(x$weights) > max_experts + 2) {
        tmp_weights <- x$weights[, c(1, max_experts, (ncol(x$weights) - max_experts + 1):ncol(x$weights))]
        names(tmp_weights)[1:2] <- c("worst_others", "best_others")
      } else {
        tmp_weights <- x$weights[, max(1, ncol(x$weights) - max_experts):ncol(x$weights)]
      }
      
      matplot(tmp_weights, type = "l", xlab = "", ylab = "", lty = 1:5, main = "Weights associated with the experts", col = col,...)
      mtext(side = 2, text = "Weights", line = 1.8, cex = 1)
      # mtext(side = 1, text = "Time steps", line = 1.8, cex = 1)
      mtext(side = 4, text = colnames(tmp_weights), at = tmp_weights[T,], las = 2, col = col, cex= 0.5, line = 0.3)
      
    } else {
      list_plt[[length(list_plt) + 1]] <- 
        {
          html_p <- rAmCharts::controlShinyPlot(
            plot_ridge_weights(data = x, colors = col, 
                               max_experts = max_experts, 
                               round = 3)
          )
          html_p$height <- 280 + 10 * min(K, max_experts)
          html_p
        }
    }
    
  } else if (type == "all" || type == "plot_weight") {
    # Convex aggregation rule
    if (! dynamic) {
      if (type == "all") {
        par(mar = c(3, 3, 2, l.names/2), mgp = c(1, 0.5, 0)) 
      } 
      if (ncol(x$weights) > max_experts) {
        tmp_weights <- x$weights[]
        tmp_weights <- cbind(apply(tmp_weights[1:(ncol(tmp_weights) - max_experts)], 1, sum), 
                             tmp_weights[, (ncol(tmp_weights) - max_experts + 1):ncol(tmp_weights)])
        names(tmp_weights)[1] <- "others"
        tmp_K <- min(K, max_experts + 1)
        tmp_cols <- c(rev(col)[1:(tmp_K-1)], "grey")
        
      } else {
        tmp_weights <- x$weights
        tmp_cols <- rev(col)
        tmp_K <- K
      }
      
      plot(c(1), type = "l", col = 1:8, lwd = 2, axes = F, xlim = c(1, T), 
           ylim = c(0,1), ylab = "", xlab = "", main = "Weights associated with the experts")
      mtext(side = 2, text = "Weights", line = 1.8, cex = 1)
      # mtext(side = 1, text = "Time steps", line = 1.8, cex = 1)
      x.idx <- c(1, 1:T, T:1)
      w.summed <- rep(1,T)
      i.remaining = rep(TRUE, tmp_K)
      for (i in 1:tmp_K) {
        y.idx <- c(0, w.summed, rep(0, T))
        polygon(x = x.idx, y = y.idx, col = tmp_cols[i], border=NA)
        w.summed.old <- w.summed
        w.summed <- w.summed - tmp_weights[, rev(names(tmp_weights))][, i]
        i.remaining[i] <- FALSE
        writeLegend(f = w.summed.old,w.summed,name = rev(colnames(tmp_weights))[i])
      }
      axis(1)
      axis(2)
      box()
      names.toWrite <- rev(colnames(tmp_weights))
      names.toWrite[-(1:min(tmp_K,15))] <- ""
      mtext(side = 4, text = names.toWrite,
            at = (1-cumsum(c(tmp_weights[, rev(names(tmp_weights))][T,])))  + 
              tmp_weights[, rev(names(tmp_weights))][T,]/2, las = 2, col = tmp_cols, cex= 0.5, line = 0.3)
      
    } else {
      list_plt[[length(list_plt) + 1]] <- 
        {
          html_p <- rAmCharts::controlShinyPlot(
            plot_weights(data = x, 
                         colors = col, 
                         max_experts = max_experts, 
                         round = 3
            )
          )
          html_p$height <- 325 + 25 * (min(K, max_experts) - 3)
          html_p
        }
    }
  }
  
  
  # boxplot of weights
  if (!is.null(x$awake)) {
    pond <- apply(x$awake[d*(1:T),],1,sum)
    normalized.weights <- x$weights * pond / (K*x$awake[d*(1:T),])
    normalized.weights[x$awake[d*(1:T),] == pond] <- NaN
  } else {
    normalized.weights <- x$weights 
  }
  
  if (type == "all" || type == "boxplot_weight") {
    if (! dynamic) {
      i.order <- 1:min(c(K, 15, max_experts))
      if (type == "all") {
        par(mar = c(l.names+2, 3, 1.6, 0.1))
      }
      
      if (ncol(x$weights) > max_experts + 2) {
        if (ncol(x$weights) <= 15) {
          i.order <- c(i.order, max_experts + 1, ncol(x$weights))
          names(normalized.weights)[c(1, ncol(x$weights) - min(14, max_experts + 1) + 1)] <- c("worst_others", "best_others")
          tmp_col <- rev(col) ; tmp_col[1:ncol(x$weight) > 13] <- "grey"
          
        } else {
          i.order <- c(i.order[1:min(13, max_experts)], min(14, max_experts + 1), ncol(x$weights))
          names(normalized.weights)[c(1, ncol(x$weights) - min(14, max_experts + 1) + 1)] <- c("worst_others", "best_others")
          tmp_col <- rev(col) ; tmp_col[1:ncol(x$weight) > 13] <- "grey"
        }
      } else {
        tmp_col <- rev(col)
      }
      
      boxplot(normalized.weights[, rev(names(normalized.weights))][,i.order], main = "Weights associated with the experts", 
              col = tmp_col[i.order], axes = FALSE, pch='.')
      mtext(side = 2, text = "Weights", line = 1.8, cex = 1)
      axis(1, at = 1:(min(c(K, 15, max_experts))), labels = FALSE)
      mtext(at = 1:min(c(K, 15, max_experts + 2)), text = rev(names(normalized.weights))[i.order], 
            side = 1, las = 2, col = tmp_col[i.order], line = 0.8)
      axis(2)
      box()
      
    } else {
      list_plt[[length(list_plt) + 1]] <- 
        {
          html_p <- rAmCharts::controlShinyPlot(
            boxplot_weights(data = x, colors = col, 
                            max_experts = max_experts
            )
          )
          html_p$height <- 300 + 10 * max(c(nchar(names(max_experts)), 17*(ncol(x$weights) > max_experts)))
          html_p
        }
    } 
  }
  
  # note: always pass alpha on the 0-255 scale
  if (! dynamic) {
    makeTransparent<-function(someColor, alpha=220)
    {
      newColor<-col2rgb(someColor)
      apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
                                                  blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
    }
  }
  
  # dynamic average loss
  if (type == "all" || type == "dyn_avg_loss") {
    if (! dynamic) {
      pred.experts <- data.frame(x$experts * x$awake + x$prediction * (1-x$awake))
      cumul.losses <- apply(loss(pred.experts, x$Y, x$loss.type), 2, cumsum)[seq(d,T*d,by=d),] / 1:T
      cumul.exploss <- cumsum(loss(x$prediction, x$Y, x$loss.type))[seq(d,T*d,by=d)] / 1:T
      
      if (ncol(x$weights) > max_experts + 2) {
        cumul.losses <- cumul.losses[, -c(2:(ncol(x$weights) - max_experts - 1))]
        colnames(cumul.losses)[1:2] <- c("worst_others", "best_others")
        tmp_col <- col[-c(2:(ncol(x$weights) - max_experts - 1))]
      } else {
        cumul.losses <- cumul.losses[, max(1, ncol(cumul.losses) - max_experts + 1):ncol(cumul.losses)]
        tmp_col <- col
      }
      
      if (type == "all") {
        par(mar = c(1.5, 3, 2.5, l.names/2), mgp = c(1, 0.5, 0))
      }
      
      matplot(cumul.losses, type = "l", lty = 1, xlab = "", ylab = "",main = "Dynamic average loss", 
              col = makeTransparent(tmp_col), ylim = range(c(cumul.losses,cumul.exploss)))
      lines(cumul.exploss, col = 1, lwd = 2)
      mtext(side = 2, text = "Cumulative loss", line = 1.8, cex = 1)
      # mtext(side = 1, text = "Time steps", line = 1.8, cex = 1)
      mtext(side = 4, text = colnames(cumul.losses), 
            at = cumul.losses[T,], las = 2, col = makeTransparent(tmp_col), cex= 0.5, line = 0.3)
      legend("topleft", c("Experts", x$model), bty = "n", lty = 1, col = c("gray", 1), lwd = c(1,2))
      
    } else {
      list_plt[[length(list_plt) + 1]] <- 
        {
          html_p <- rAmCharts::controlShinyPlot(
            plot_dyn_avg_loss(data = x, colors = col, 
                              max_experts = max_experts, round = 3
            )
          )
          html_p$height <- 322 + 22 * (min(K, max_experts) - 3)
          html_p
        }
    } 
  }
  
  
  # cumulative residuals
  if (type == "all" || type == "cumul_res") {
    if (! dynamic) {
      pred.experts <- data.frame(x$experts * x$awake + x$prediction * (1-x$awake))
      cumul.residuals <- apply(x$Y - pred.experts, 2, cumsum)[seq(d,T*d,by=d),]
      cumul.expres <- cumsum(x$Y - x$prediction)[seq(d,T*d,by=d)]
      
      if (ncol(x$weights) > max_experts + 2) {
        cumul.residuals <- cumul.residuals[, -c(2:(ncol(x$weights) - max_experts - 1))]
        colnames(cumul.residuals)[1:2] <- c("worst_others", "best_others")
        tmp_col <- col[-c(2:(ncol(x$weights) - max_experts - 1))]
      } else {
        cumul.residuals <- cumul.residuals[, max(1, ncol(cumul.residuals) - max_experts + 1):ncol(cumul.residuals)]
        tmp_col <- col
      }
      
      if (type == "all") {
        par(mar = c(1.5, 3, 2.5,l.names/2), mgp = c(1, 0.5, 0))
      }
      
      matplot(cumul.residuals, type = "l", lty = 1, xlab = "", ylab = "",
              main = paste("Cumulative residuals"), col = makeTransparent(tmp_col), ylim = range(c(cumul.residuals,cumul.expres)))
      lines(cumul.expres, col = 1, lwd = 2)
      mtext(side = 2, text = "Cumulative residuals", line = 1.8, cex = 1)
      # mtext(side = 1, text = "Time steps", line = 1.8, cex = 1)
      if (max(cumul.residuals) > abs(min(cumul.residuals))) {
        place = "topleft"
      } else {
        place = "bottomleft"
      }
      mtext(side = 4, text = colnames(cumul.residuals), 
            at = cumul.residuals[T,], las = 2, col = tmp_col, cex= 0.5, line = 0.3)
      legend(place, c("Experts", x$model), bty = "n", lty = 1, col = c("gray", 1), lwd = c(1,2))
      
    } else {
      list_plt[[length(list_plt) + 1]] <- {
        html_p <- rAmCharts::controlShinyPlot(
          plot_cumul_res(
            data = x, colors = col, 
            max_experts = max_experts, round = 3
          )
        )
        html_p$height <- 322 + 22 * (min(K, max_experts) - 3)
        html_p
      }
    } 
  }
  
  
  # losses
  if (type == "all" || type == "avg_loss") {
    if (! dynamic) {
      pred.experts <- data.frame(x$experts * x$awake + x$prediction * (1-x$awake))
      x$loss.experts <- apply(loss(x = pred.experts,y = x$Y,loss.type = x$loss.type),2,mean)
      err.unif <- lossConv(rep(1/K, K), x$Y, x$experts, awake = x$awake, loss.type = x$loss.type)
      err.mixt <- x$loss
      
      if (ncol(x$weights) > max_experts + 2) {
        x$loss.experts <- c(x$loss.experts[-c(2:(ncol(x$weight) - max_experts - 1))], "uniform" = err.unif, "mixt" = err.mixt)
        names(x$loss.experts)[1:2] <- c("worst_others", "best_others")
        idx.sorted <- order(x$loss.experts)
        tmp_cols <- c("grey", "grey", col[-c(1:(ncol(x$weight) - max_experts))], "black", "black")[idx.sorted]
        my.pch <- c(rep(20, length(x$loss.experts)-2),8,8)[idx.sorted]
        
      } else {
        x$loss.experts <- c(x$loss.experts[max(1, length(x$loss.experts) - max_experts + 1):length(x$loss.experts)], "uniform" = err.unif, "mixt" = err.mixt)
        idx.sorted <- order(x$loss.experts)
        tmp_cols <- c(col[max(1, col(x$weights) - max_experts + 1):ncol(x$weights)], "black", "black")[idx.sorted]
        my.pch <- c(rep(20, length(x$loss.experts)-2),8,8)[idx.sorted]
      }
      
      l.names <- max(max(nchar(names(x$loss.experts))) / 3 + 1.7,4)
      
      if (type == "all") {
        par(mar = c(l.names, 3, 2.5,l.names/2), mgp = c(1, 0.5, 0))
      }
      
      plot(x$loss.experts[idx.sorted], xlab = "", ylab = "", main = "Average loss suffered by the experts", 
           axes = F, pch = my.pch, col = tmp_cols, lwd = 2, type='b')
      mtext(side = 2, text = paste(x$loss.type$name,"loss"), line = 1.8, cex = 1)
      axis(1, at = 1:(K + 2), labels = FALSE)
      mtext(at = 1:length(x$loss.experts), text = c(names(x$loss.experts), "Uniform", x$model)[idx.sorted],
            side = 1, las = 2, col = tmp_cols, line = 0.8,cex = .7)
      
      axis(2)
      box()
      
    } else {
      list_plt[[length(list_plt) + 1]] <- 
        {
          html_p <- rAmCharts::controlShinyPlot(
            plot_avg_loss(
              data = x, colors = col, 
              max_experts = max_experts, round = 3
            )
          )
          html_p$height <- 300
          html_p
        }
    }
  }
  
  
  # cumulative plot of the series
  if (type == "all" || type == "contrib") {
    if (! dynamic) {
      if (x$d ==1) {
        if (type == "all") {
          par(mar = c(2, 3, 2.5,l.names/2), mgp = c(1, 0.5, 0))
        }
        
        cumulativePlot(W = x$weights,X = x$experts, Y = x$Y,smooth = TRUE, alpha = alpha, 
                       plot.Y = TRUE, col.pal = col, max_experts = max_experts)
        
      } else {
        X <- apply(seriesToBlock(X = x$experts,d = x$d),c(1,3),mean)
        Y <- apply(seriesToBlock(x$Y,d = x$d),1,mean)
        colnames(X) <- names(x$weights)
        
        if (type == "all") {
          par(mar = c(2, 3, 2.5,l.names/2), mgp = c(1, 0.5, 0))
        }
        
        cumulativePlot(W = x$weights,X = X, Y = Y,smooth = TRUE,
                       alpha = alpha,plot.Y = TRUE, col.pal = rev(col),
                       max_experts = max_experts)
      }
    } else {
      list_plt[[length(list_plt) + 1]] <- 
        {
          html_p <- rAmCharts::controlShinyPlot(
            plot_contrib(
              data = x, colors = col, alpha = alpha, 
              max_experts = max_experts, round = 3
            )
          )
          html_p$height <- 325 + 25 * (min(K, max_experts) - 3)
          html_p
        }
    }
  }
  
  if (! dynamic) {
    par(def.par) 
    return(invisible(NULL))
  } else {
    res <- htmltools::browsable(
      htmltools::tagList(
        list_plt
      )
    )
    return(res)
  }
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
    j = max(4, floor(60 *l/length(f)))
    text(v+l/2,x+tau/2,substr(name,1,j),cex = 0.8,...)
  }
}

cumulativePlot<-function(W,X,Y,col.pal=NULL, smooth = FALSE, plot.Y = FALSE, alpha = 0.1, max_experts = 50)
{
  time<-c(1:nrow(X))
  active.experts<-which(colMeans(W)>0)
  W<-W[,active.experts]  
  X<-X[, names(W)][,active.experts]
  
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
  
  if (ncol(mat) > max_experts) {
    colnames(mat)[1] <- "others"
  }
  
  plot(x = NULL,y = NULL,col=col.pal[1], type='l', xaxt='n',ylim=Y.lim,lty='dotted',
       yaxt='n',xlab="",ylab="",lwd=3,xlim = range(time),
       main = paste("Contribution of each expert to prediction"))
  y.summed <- Agg
  for(i in rev(c(1:ncol(mat))))
  {
    if (!smooth) addPoly(time,y.summed,col=ifelse(i <= (ncol(mat) - max_experts), "grey", col.pal[i]))
    if (smooth) addPoly(time,lowess(y.summed,f = alpha)$y,col=ifelse(i <= (ncol(mat) - max_experts), "grey", col.pal[i]))
    y.summed.old <- y.summed
    y.summed <- y.summed - mat[,i]
    if (!smooth) writeLegend(f=y.summed.old,g= y.summed, 
                             name = colnames(mat)[i],Y.lim,col=ifelse(i <= (ncol(mat) - max_experts), "black", "black"))
    if (smooth) writeLegend(f=lowess(y.summed.old,f=alpha/10)$y,g=lowess(y.summed,f=alpha/10)$y, 
                            name = colnames(mat)[i],Y.lim,col=ifelse(i <= (ncol(mat) - max_experts), "black", "black"))
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



##########
### dynamic version of the plots

#' Functions to render dynamic mixture graphs using rAmCharts
#'
#' @param data \code{mixture object}. Displays graphs.
#' @param colors \code{character}. Colors of the lines and bullets.
#' @param max_experts \code{integer}. Maximum number of experts to be displayed (only the more influencial).
#' @param round \code{integer}. Precision of the displayed values.
#' @param alpha \code{numeric}. Smoothing parameter for contribution plot (parameter 'f' of function \code{\link[stats]{lowess}}).
#'
#' @return a \code{rAmCharts} plot
#' 
#' @import pipeR
#' @importFrom rAmCharts amSerialChart addValueAxis addGraph addTitle setExport setChartCursor setChartScrollbar setLegend 
#' amBoxplot setCategoryAxis controlShinyPlot
#' 
#' @rdname plot-opera-rAmCharts
#' 
plot_ridge_weights <- function(data,
                               colors = NULL,
                               max_experts = 50,
                               round = 3) {
  
  if (is.null(colors)) {
    colors <- RColorBrewer::brewer.pal(n = min(max(3, ncol(data)), 9), name = "Spectral")
  }
  
  data <- data$weights
  
  if (ncol(data) > max_experts + 2) {
    colors <- colors[-c(1:(ncol(data) - max_experts - 2))]
    data <- data[, c(1, (ncol(data) - max_experts):ncol(data))]
    names(data)[1:2] <- c("worst_others", "best_others")
  } else {
    colors <- colors[-c(1:(ncol(data) - max_experts))]
    data <- data[, max(1, ncol(data) - max_experts + 1):ncol(data)]
  }
  
  names_experts <- names(data)
  data$timestamp <- 1:nrow(data)
  
  plt <- amSerialChart(dataProvider = data,
                       categoryField = c("timestamp"), 
                       creditsPosition = "bottom-right",
                       thousandsSeparator = " ",
                       precision = round) %>>%
    rAmCharts::addValueAxis(title = "Weights")
  
  for (index in 1:length(names_experts)) {
    plt <- plt %>>%
      rAmCharts::addGraph(title = names_experts[index], id = names_experts[index],
                          valueField = names_experts[index], valueAxis = "timestamp", 
                          type = "line", lineColor = colors[index])
  }
  
  plt <- plt %>>%
    rAmCharts::addTitle(text = "Weights associated with the experts") %>>%
    rAmCharts::setExport(position = "bottom-right") %>>% 
    rAmCharts::setChartCursor() %>>% 
    # rAmCharts::setChartScrollbar(scrollbarHeight = 10, dragIconHeight = 26, offset = 8) %>>%
    rAmCharts::setLegend(useGraphSettings = F, valueText = "", position = "right", reversedOrder = T)
  
  plt@otherProperties$zoomOutButtonImageSize <- 0
  
  plt
}


#' @rdname plot-opera-rAmCharts
plot_weights <- function(data,
                         colors = NULL,
                         max_experts = 50,
                         round = 3) {
  
  if (is.null(colors)) {
    colors <- RColorBrewer::brewer.pal(n = min(max(3, ncol(data$experts)), 9), name = "Spectral")
  }
  
  data_weight <- data$weights
  
  if (ncol(data_weight) > max_experts + 2) {
    data_weight <- cbind(apply(data_weight[1:(ncol(data_weight) - max_experts)], 1, sum), data_weight[, (ncol(data_weight) - max_experts + 1):ncol(data_weight)])
    names(data_weight)[1] <- "others"
    colors <- colors[-c(2:(ncol(data$weights) - max_experts))]
  }
  
  names_weights <- colnames(data_weight)
  data_weight <- data.frame("timestamp" = 1:data$`T`, t(apply(data_weight, 1, cumsum)), round(data_weight, round))
  
  
  
  plt <- amSerialChart(dataProvider = data_weight,
                       categoryField = c("timestamp"), 
                       creditsPosition = "bottom-right",
                       thousandsSeparator = " ",
                       precision = round) %>>%
    rAmCharts::addValueAxis(title = "Weights", maximum = 1)
  
  for (index in 1:length(names_weights)) {
    if (index == 1) {
      plt <- plt %>>%
        rAmCharts::addGraph(title = names_weights[index], id = names_weights[index],
                            valueField = names_weights[index], valueAxis = "timestamp", 
                            type = "line", lineColor = colors[index],
                            fillToAxis = T, fillColors = colors[index], fillAlphas = .75,
                            balloonText = paste0("<b>", names_weights[index], " : </b>", "[[", names_weights[index], ".1]]")) 
      
    } else {
      plt <- plt %>>%
        rAmCharts::addGraph(title = names_weights[index], id = names_weights[index],
                            valueField = names_weights[index], valueAxis = "timestamp", 
                            type = "line", lineColor = colors[index],
                            fillToGraph = names_weights[index-1], fillColors = colors[index], fillAlphas = .75,
                            balloonText = paste0("<b>", names_weights[index], " : </b>", "[[", names_weights[index], ".1]]"))
    }
  }
  
  plt <- plt %>>%
    rAmCharts::addTitle(text = "Weights associated with the experts") %>>%
    rAmCharts::setExport(position = "bottom-right") %>>% 
    rAmCharts::setChartCursor() %>>% 
    # rAmCharts::setChartScrollbar(scrollbarHeight = 10, dragIconHeight = 26, offset = 8) %>>%
    rAmCharts::setLegend(useGraphSettings = F, valueText = "", position = "right", reversedOrder = T)
  
  plt@otherProperties$zoomOutButtonImageSize <- 0
  
  plt
}


#' @rdname plot-opera-rAmCharts
boxplot_weights <- function(data,
                            colors = NULL,
                            max_experts = 50) {
  
  if (is.null(colors)) {
    colors <- RColorBrewer::brewer.pal(n = min(ncol(data$experts), 9), name = "Spectral")
  }
  
  data_weight <- data$weight
  
  if (ncol(data_weight) > max_experts + 2) {
    data_weight <- data_weight[, -c(2:(ncol(data$weights) - max_experts - 1))]
    names(data_weight)[1:2] <- c("worst_others", "best_others")
    colors <- colors[-c(2:(ncol(data$weights) - max_experts - 1))]
  } else {
    data_weight <- data_weight[, max(1, ncol(data_weight) - max_experts):ncol(data_weight)]
  }
  
  plt <- rAmCharts::amBoxplot(data_weight[, rev(names(data_weight))], col = rev(colors),
                              ylab = "weights", creditsPosition = "bottom-right") %>>%
    rAmCharts::addTitle(text = "Weights associated with the experts") %>>%
    rAmCharts::setCategoryAxis(autoGridCount = FALSE, gridCount = ncol(data_weight), labelRotation = 90, labelOffset = 5) %>>%
    rAmCharts::setExport(position = "bottom-right") # %>>% 
  # rAmCharts::setLegend(useGraphSettings = TRUE, valueText = "", position = "right")
  
  plt
}


#' @rdname plot-opera-rAmCharts
plot_dyn_avg_loss <- function(data,
                              colors = NULL,
                              max_experts = 50,
                              round = 3) {
  
  if (is.null(colors)) {
    colors <- RColorBrewer::brewer.pal(n = min(ncol(data$experts), 9), name = "Spectral")
  }
  
  pred.experts <- data.frame(data$experts * data$awake + data$prediction * (1-data$awake))
  cumul.losses <- apply(loss(pred.experts, data$Y, data$loss.type), 2, cumsum)[seq(data$d, data$T*data$d, by = data$d), ] / 1:nrow(data$experts)
  cumul.exploss <- cumsum(loss(data$prediction, data$Y, data$loss.type))[seq(data$d, data$T*data$d, by = data$d)] / 1:nrow(data$experts)
  
  data_loss <- data.frame(cbind(cumul.losses, cumul.exploss))
  data_loss$timestamp <- 1:nrow(data_loss)
  data_loss[, c(names(data$weights), "cumul.exploss", "timestamp")]
  
  if (ncol(data$weight) > max_experts + 2) {
    data_loss <- data_loss[, -c(2:(ncol(data$weights) - max_experts - 1))]
    names(data_loss)[1:2] <- c("worst_others", "best_others")
    colors <- colors[-c(2:(ncol(data$weights) - max_experts - 1))]
  } else {
    data_loss <- data_loss[, max(1, ncol(data$weights) - max_experts + 1):ncol(data_loss)]
  }
  
  names_experts  <- setdiff(names(data_loss), c("cumul.exploss", "timestamp"))
  
  data_loss <- round(data_loss, round)
  
  plt <- amSerialChart(dataProvider = data_loss,
                       categoryField = "timestamp", 
                       creditsPosition = "bottom-right",
                       thousandsSeparator = " ",
                       precision = round) %>>%
    rAmCharts::addValueAxis(title = "Cumulative loss")
  
  for (index in 1:length(names_experts)) {
    plt <- plt %>>%
      rAmCharts::addGraph(title = names_experts[index], id = names_experts[index],
                          valueField = names_experts[index], valueAxis = "timestamp", 
                          type = "line", lineColor = colors[index],
                          balloonText = paste0("<b>", names_experts[index], " : </b>", "[[", names_experts[index], "]]")) 
  }
  
  plt <- plt %>>%
    rAmCharts::addGraph(title = data$model, id = "cumul.exploss",
                        valueField = "cumul.exploss", valueAxis = "timestamp", 
                        type = "line", lineThickness = 2, lineColor = "black",
                        balloonText = paste0("<b> cumul.exploss : </b>", "[[cumul.exploss]]"))
  
  plt <- plt %>>%
    rAmCharts::addTitle(text = "Dynamic average loss") %>>%
    rAmCharts::setExport(position = "bottom-right") %>>% 
    rAmCharts::setChartCursor() %>>% 
    # rAmCharts::setChartScrollbar(scrollbarHeight = 10, dragIconHeight = 26, offset = 8) %>>%
    rAmCharts::setLegend(useGraphSettings = F, valueText = "", position = "right", reversedOrder = T)
  
  plt@otherProperties$zoomOutButtonImageSize <- 0
  
  plt
}


#' @rdname plot-opera-rAmCharts
plot_cumul_res <- function(data,
                           colors = NULL,
                           max_experts = 50,
                           round = 3) {
  
  if (is.null(colors)) {
    colors <- RColorBrewer::brewer.pal(n = min(ncol(data$experts), 9), name = "Spectral")
  }
  
  pred.experts <- data.frame(data$experts * data$awake + data$prediction * (1-data$awake))
  cumul.residuals <- apply(data$Y - pred.experts, 2, cumsum)[seq(data$d, data$T*data$d, by = data$d),]
  cumul.expres <- cumsum(data$Y - data$prediction)[seq(data$d, data$T*data$d, by = data$d)]
  
  data_res <- data.frame(cbind(cumul.residuals, cumul.expres))
  data_res$timestamp <- 1:nrow(data_res)
  data_res[, c(names(data$weights), "cumul.expres", "timestamp")]
  
  if (ncol(data$weight) > max_experts + 2) {
    data_res <- data_res[, -c(2:(ncol(data$weights) - max_experts - 1))]
    names(data_res)[1:2] <- c("worst_others", "best_others")
    colors <- colors[-c(2:(ncol(data$weights) - max_experts - 1))]
  } else {
    data_res <- data_res[, max(1, ncol(data$weights) - max_experts + 1):ncol(data_res)]
  }
  
  names_experts  <- setdiff(names(data_res), c("cumul.expres", "timestamp"))
  
  data_res <- round(data_res, round)
  
  plt <- amSerialChart(dataProvider = data_res,
                       categoryField = "timestamp", 
                       creditsPosition = "bottom-right",
                       thousandsSeparator = " ",
                       precision = round) %>>%
    rAmCharts::addValueAxis(title = "Cumulative residuals")
  
  for (index in 1:length(names_experts)) {
    plt <- plt %>>%
      rAmCharts::addGraph(title = names_experts[index], id = names_experts[index],
                          valueField = names_experts[index], valueAxis = "timestamp", 
                          type = "line", lineColor = colors[index],
                          balloonText = paste0("<b>", names_experts[index], " : </b>", "[[", names_experts[index], "]]")) 
  }
  
  plt <- plt %>>%
    rAmCharts::addGraph(title = data$model, id = "cumul.expres",
                        valueField = "cumul.expres", valueAxis = "timestamp", 
                        type = "line", lineThickness = 2, lineColor = "black",
                        balloonText = paste0("<b> cumul.expres : </b>", "[[cumul.expres]]"))
  
  plt <- plt %>>%
    rAmCharts::addTitle(text = "Cumulative residuals") %>>%
    rAmCharts::setExport(position = "bottom-right") %>>% 
    rAmCharts::setChartCursor() %>>% 
    # rAmCharts::setChartScrollbar(scrollbarHeight = 10, dragIconHeight = 26, offset = 8) %>>%
    rAmCharts::setLegend(useGraphSettings = F, valueText = "", position = "right", reversedOrder = T)
  
  plt@otherProperties$zoomOutButtonImageSize <- 0
  
  plt
}


#' @rdname plot-opera-rAmCharts
plot_avg_loss <- function(data,
                          colors = NULL,
                          max_experts = 50,
                          round = 3) {
  
  K <- ncol(data$experts)
  pred.experts <- data.frame(data$experts * data$awake + data$prediction * (1-data$awake))
  data$loss.experts <- apply(loss(x = pred.experts, y = data$Y, loss.type = data$loss.type), 2, mean)
  err.unif <- lossConv(rep(1/K, K), data$Y, data$experts, awake = data$awake, loss.type = data$loss.type)
  err.mixt <- data$loss
  
  data_loss <- c(data$loss.experts, "uniform" = err.unif, "mixt" = data$loss)
  
  if (is.null(colors)) {
    colors <- RColorBrewer::brewer.pal(n = min(length(data_loss)-2, 9), name = "Spectral")
  }
  
  data_plot <- data.frame("values" = unname(data_loss), 
                          "bullet" = "diamond", "size" = 10, "cols" = "black", "names" = names(data_loss))
  data_plot$bullet = ifelse(data_plot$names %in% c("mixt", "uniform"), "diamond", "round")
  data_plot$size = ifelse(data_plot$names %in% c("mixt", "uniform"), 10, 8)
  colors_bis <- rep("black", length(data_loss)) ; colors_bis[! data_plot$names %in% c("mixt", "uniform")] <- colors
  data_plot$cols = colors_bis
  data_plot$names <- ifelse(data_plot$names == "mixt", data$model, data_plot$names)
  
  if (ncol(data$weights) > max_experts + 2) {
    data_plot <- data_plot[-c(2:(ncol(data$weight) - max_experts - 1)), ]
    data_plot[1:2, ]$names <- c("worst_others", "best_others")
    data_plot <- data_plot[order(data_plot$values), ]
  } else {
    data_plot <- data_plot[max(1, ncol(data$weight) - max_experts + 1):nrow(data_plot), ]
    data_plot <- data_plot[order(data_loss[max(1, ncol(data$weight) - max_experts + 1):nrow(data_plot)]), ]
  }
  
  plt <- amSerialChart(dataProvider = data_plot,
                       categoryField = "names", 
                       creditsPosition = "bottom-right",
                       thousandsSeparator = "",
                       precision = round) %>>%
    rAmCharts::addValueAxis(title = "Square loss") %>>%
    rAmCharts::addGraph(title = "lines", id = "lines",
                        valueField = "values", valueAxis = "names", 
                        type = "line", lineColor = "black",
                        showBalloon = F) %>>%
    rAmCharts::addGraph(title = "bullets", id = "bullets",
                        valueField = "values", valueAxis = "names", 
                        type = "line", lineAlpha = 0, 
                        bulletField = "bullet", bulletSizeField = "size", colorField = "cols") %>>%
    rAmCharts::addTitle(text = "Average loss suffered by the experts") %>>%
    rAmCharts::setExport(position = "bottom-right") %>>% 
    rAmCharts::setChartCursor() %>>%
    rAmCharts::setCategoryAxis(autoGridCount = FALSE, gridCount = nrow(data_plot), labelRotation = 90, labelColorField = "cols", labelOffset = 5)
  
  return(plt)
}


#' @rdname plot-opera-rAmCharts
plot_contrib <- function(data, 
                         colors = NULL, 
                         alpha = 0.1,
                         max_experts = 50,
                         round = 3) {
  
  W = data$weights
  K = ncol(data$experts)
  if (data$d == 1) {
    X = data$experts
    Y = data$Y
  } else {
    X <- apply(seriesToBlock(X = data$experts, d = data$d), c(1, 3), mean)
    Y <- apply(seriesToBlock(data$Y, d = data$d), 1, mean)
    colnames(X) <- names(data$weights)
  }
  
  if (is.null(colors)) {
    colors <- RColorBrewer::brewer.pal(n = min(K, 9), name = "Spectral")
  }
  if (length(colors) < K) colors <- c(rep(colors[1],K-length(colors)),colors)
  
  time<-c(1:nrow(X))
  active.experts<-which(colMeans(W)>0)
  W<-W[,active.experts]  
  X<-X[, names(W)][,active.experts]
  
  K <- ncol(X)
  
  o<-order(colSums(W),decreasing = F)
  mat<-W[,o]*X[,o]
  mat <- sapply(mat, function(x) lowess(x = 1:nrow(mat), y = x, f = alpha)$y)
  colnames(mat)<-colnames(X)[o]
  
  data_weight <- as.data.frame(mat)
  
  if (ncol(data_weight) > max_experts + 2) {
    data_weight <- cbind(apply(data_weight[1:(ncol(data_weight) - max_experts)], 1, sum), data_weight[, (ncol(data_weight) - max_experts + 1):ncol(data_weight)])
    names(data_weight)[1] <- "others"
    colors <- colors[-c(2:(ncol(mat) - max_experts))]
  }
  
  names_weights <- colnames(data_weight)
  data_weight <- data.frame("timestamp" = 1:nrow(data_weight), t(apply(data_weight, 1, cumsum)), round(data_weight, round))
  
  data_weight$pred <- round(lowess(x = time,y = Y,f = alpha)$y, round)
  
  plt <- amSerialChart(dataProvider = data_weight,
                       categoryField = c("timestamp"), 
                       creditsPosition = "bottom-right",
                       thousandsSeparator = " ",
                       precision = round) %>>%
    rAmCharts::addValueAxis(maximum = max(data_weight$pred), useScientificNotation = T)
  
  for (index in 1:length(names_weights)) {
    if (index == 1) {
      plt <- plt %>>%
        rAmCharts::addGraph(title = names_weights[index], id = names_weights[index],
                            valueField = names_weights[index], valueAxis = "timestamp", 
                            type = "line", lineColor = colors[index],
                            fillToAxis = T, fillColors = colors[index], fillAlphas = .75,
                            balloonText = paste0("<b>", names_weights[index], " : </b>", "[[", names_weights[index], ".1]]")) 
      
    } else {
      plt <- plt %>>%
        rAmCharts::addGraph(title = names_weights[index], id = names_weights[index],
                            valueField = names_weights[index], valueAxis = "timestamp", 
                            type = "line", lineColor = colors[index],
                            fillToGraph = names_weights[index-1], fillColors = colors[index], fillAlphas = .75,
                            balloonText = paste0("<b>", names_weights[index], " : </b>", "[[", names_weights[index], ".1]]"))
    }
  }
  
  plt <- plt %>>%
    rAmCharts::addGraph(title = "Prediction", id = "pred",
                        valueField = "pred", valueAxis = "timestamp", 
                        type = "line", dashLength = 5, lineThickness = 2, lineColor = "black",
                        balloonText = paste0("<b> Prediction : </b>", "[[pred]]"))
  
  plt <- plt %>>%
    rAmCharts::addTitle(text = "Contribution of each expert to the prediction") %>>%
    rAmCharts::setExport(position = "bottom-right") %>>% 
    rAmCharts::setChartCursor() %>>% 
    # rAmCharts::setChartScrollbar(scrollbarHeight = 10, dragIconHeight = 26, offset = 8) %>>%
    rAmCharts::setLegend(useGraphSettings = F, valueText = "", position = "right", reversedOrder = T)
  
  plt@otherProperties$zoomOutButtonImageSize <- 0
  
  plt
}
