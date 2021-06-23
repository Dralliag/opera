#' Plot an aggregation procedure
#' @description oracle \code{plot}. It has one optional arguments. 
#' @param x An object of class \code{oracle}. 
#' @param sort if set to TRUE (default), it sorts the experts by performance before the plots.
#' @param col colors
#' @importFrom graphics axis box mtext par plot
#' @export 
plot.oracle <- function(x, sort = TRUE, col = NULL, dynamic = T,  ...) {
  def.par <- par(no.readonly = TRUE)
  
  ############# add checks on x$loss
  
  if (x$d > 1) {
    x$experts <- blockToSeries(x$experts)
    x$Y <- blockToSeries(x$Y)
  }
  
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
     colnames(x$experts) <- x$names.experts
  } else {
    if (is.null(colnames(x$experts))) {
      colnames(x$experts) <- x$names.experts <- paste("X", 1:K,sep="")
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
    
    if (! dynamic) {
      par(mar = c(4.5, 4, 2, 2))
      plot(c(x$loss.experts, err.unif)[idx.sorted], xlab = "", ylab = paste(x$loss.type$name, 
                                                                            "loss"), main = "Average loss suffered by the experts", axes = F, pch = 3, 
           col = my.col, lwd = 2,...)
      axis(1, at = 1:(K + 1), labels = FALSE)
      mtext(at = 1:(K + 1), text = c(colnames(x$experts), "Uniform")[idx.sorted], 
            side = 1, las = 2, col = my.col, line = 0.8)
      axis(2)
      box() 
      
    } else {
      data <- c(x$loss.experts, err.unif) ; names(data)[length(data)] <- "Uniform"
      plt <- plt_oracle_convex(data = data, 
                               colors = my.col,
                               round = 3)
      
    }
  }
  
  if (x$model == "convex" || x$model == "linear") {
    err.unif <- lossConv(rep(1/K, K), x$Y, x$experts, awake = x$awake, loss.type = x$loss.type)

    idx.sorted <- order(c(x$loss.experts, err.unif, x$loss))
    my.col <- rep(1, K + 2)
    my.col[which(idx.sorted == K + 1)] <- col[1]
    my.col[which(idx.sorted == K + 2)] <- col[2]
    my.col[which(!(idx.sorted %in% c(K + 1, K + 2)))[1]] <- col[3]
    
    if (! dynamic) {
      y.max <- c(x$loss.experts, err.unif, x$loss)[idx.sorted]
      
      par(mar = c(4.5, 4, 2, 2))
      plot(c(x$loss.experts, err.unif, x$loss)[idx.sorted], xlab = "", ylab = paste(x$loss.type$name, 
                                                                                    "loss"), main = "Average loss suffered by the experts", axes = F, pch = 3, 
           col = my.col, lwd = 2)
      axis(1, at = 1:(K + 2), labels = FALSE)
      mtext(at = 1:(K + 2), text = c(colnames(x$experts), "Uniform", 
                                     x$model)[idx.sorted], 
            side = 1, las = 2, col = my.col, line = 0.8)
      axis(2)
      box()
      
    } else {
      data <- c(x$loss.experts, err.unif, x$loss) ; names(data)[(length(data)-1):length(data)] <- c("Uniform", "convex")
      plt <- plt_oracle_convex(data = data, 
                               colors = my.col,
                               round = 3)
    }
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
    
    if (! dynamic) {
      plot(0:(length(L) - 1), L, xlab = "Number of shifts", ylab = y.lab, type = "o", 
           pch = 20, cex = 0.6, lwd = 2, main = "Error suffered by the shifting oracle")
      
    } else {
      data <- data.frame(x = 0:(length(L) - 1),
                         y = L)
      plt <- plt_oracle_shift(data = data,
                       ylab = y.lab, 
                       round = 5)
    }
    
  }
  
  if (! dynamic) {
    res <- NULL
    par(def.par)  
  } else {
    res <- rAmCharts::plot(plt, height = 300)
  }
  
  return(res)
}



##########
### dynamic version of the plots 

#' Functions to render dynamic oracle graphs using rAmCharts
#'
#' @return a rAmCharts plot
#' 
#' @import rAmCharts pipeR
#' 
plt_oracle_convex <- function(data, 
                              colors,
                              round = 3) {
  
  if (is.null(colors)) {
    colors <- RColorBrewer::brewer.pal(n = min(length(data), 9), name = "Spectral")
  }
  
  data <- data[order(data)]
  
  data_plot <- data.frame("values" = unname(data), 
                          "bullet" = "diamond", 
                          "size" = 10, 
                          "cols" = colors[1:length(data)], 
                          "names" = names(data))
  
  plt <- amSerialChart(dataProvider = data_plot,
                       categoryField = "names", 
                       creditsPosition = "top-right") %>>%
    rAmCharts::addValueAxis(title = "Square loss") %>>%
    rAmCharts::addGraph(title = "lines", id = "lines",
                        valueField = "values", valueAxis = "names", 
                        type = "line", lineColor = "black",
                        showBalloon = F) %>>%
    rAmCharts::addGraph(title = "bullets", id = "bullets",
                        valueField = "values", valueAxis = "names", 
                        type = "line", lineAlpha = 0, 
                        bulletField = "bullet", bulletSizeField = "size", colorField = "cols",
                        precision = round) %>>%
    rAmCharts::addTitle(text = "Average loss suffered by the experts") %>>%
    rAmCharts::setExport(position = "top-right") %>>% 
    rAmCharts::setChartCursor() %>>%
    rAmCharts::setCategoryAxis(labelRotation = 90, labelColorField = "cols", labelOffset = 5)
  
  return(plt)
}


plt_oracle_shift <- function(data, 
                             ylab,
                             round = 5) {
  
  plt <- amSerialChart(dataProvider = data,
                       categoryField = "x", 
                       creditsPosition = "top-right") %>>%
    rAmCharts::addValueAxis(title = ylab) %>>%
    rAmCharts::addGraph(title = "shift", id = "shift",
                        valueField = "y", valueAxis = "x", 
                        type = "line", lineColor = "black",
                        precision = round) %>>%
    rAmCharts::addTitle(text = "Error suffered by the shifting oracle") %>>%
    rAmCharts::setExport(position = "top-right") %>>% 
    rAmCharts::setChartCursor() %>>%
    rAmCharts::setCategoryAxis(title = "Number of shifts", labelRotation = 90, labelColorField = "cols", labelOffset = 5)
  
  return(plt)
}