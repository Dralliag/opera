# best convex oracle
bestConvex <-
function(y, experts, awake=NULL, loss.type='square', niter = 1, tau = 0.5, ...)
{
   experts <- as.matrix(experts)
   N <- ncol(experts)
   
   # if there are no NA and if awake is null 
   # we can perform an exact resolution for the square loss
   idx.na <- which(is.na(experts))
   if (length(idx.na) == 0 && is.null(awake) && loss.type == "square") {
      y.na = is.na(y)
      y = y[!y.na]
      x = experts[!y.na,]
      eq = paste("y ~ x-1")

      Q <- crossprod(x)
      c <- crossprod(x, y)
      A <- cbind(1,diag(nrow(Q)))
      b <- c(1, rep(0,nrow(Q)))
      m <- 1
      res <- tryCatch({
        quadprog::solve.QP(Dmat = Q, dvec = c, Amat = A, bvec = b, meq = m)},
        error = function(e) {
          NULL
        })
      if (!is.null(res)) {
        weights = matrix(res$solution, ncol = N)
        prediction = experts %*% t(weights)
        bestLoss = mean(loss(x = prediction, y))
      }
   } 
   else {
     res = NULL
   }
   if (is.null(res)) {
      warning("The best convex oracle is only approximated (using optim).")
      if (is.null(awake)) {awake = as.matrix(array(1,dim(experts)))}
      awake[idx.na] <- 0
      experts[idx.na] <- 0
      
      lossp <- function(p)   {
         return(lossConv(p, y, experts, awake, loss.type, tau)) 
      }
      
      best_p <- rep(0,N)
      bestLoss <- exp(700)
      
      for (i in 1:niter)
      {
         # Optimisation convexe avec choix aléatoire de la condition initiale   
         p <- runif(N,0,1)
         p <- p/sum(p)
         w <- optim(p,lossp, gr = NULL, lower = 1e-20, method = "L-BFGS-B", ...)
         # Projection sur le simplex
         w <- pmax(w$par,0)
         l <- lossp(w)
         if (bestLoss > l) {
            bestLoss = l
            best_p = w  
         }
         #print(c(i, l))  
      }
      weights = matrix(best_p, ncol = N)
      weights = weights / apply(weights,1,sum)
      pond <- awake %*% t(weights)
      prediction <- ((experts* awake) %*% t(weights)) / pond
   }
   res = list(loss = bestLoss, weights = weights, prediction = prediction)
   if (loss.type == "square") {
      res$rmse = sqrt(res$loss)
   }
   return(res)
}
