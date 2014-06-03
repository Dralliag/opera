bestConvex <-
function(y,experts, awake=NULL, loss.type='squareloss', niter = 1, method = "L-BFGS-B", control = list(maxit = 100, trace=T))
{
   experts <- as.matrix(experts)
   N <- ncol(experts)
   eps <- 1e-20
   
   if (is.null(awake)) {
      awake = as.matrix(array(1,dim(experts)))
   }
   idx.na <- which(is.na(experts))
   awake[idx.na] <- 0
   experts[idx.na] <- 0
   
   lossp <- function(p)   {
      return(lossConv(p, y, experts, awake, loss.type)) 
   }
   
   best_p <- rep(0,N)
   bestLoss <- exp(700)
   
   for (i in 1:niter)
   {
      # Optimisation convexe avec choix aléatoire de la condition initiale   
      p <- runif(N,0,1)
      p <- p/sum(p)
      w <- optim(p,lossp, gr = NULL, method = method, lower = eps)
      # Projection sur le simplex
      w <- pmax(w$par,0)
      l <- lossp(w)
      if (bestLoss > l) {
         bestLoss = l
         best_p = w  
      }
      print(c(i, l))  
   }
   weights = matrix(best_p, ncol = N)
   weights = weights / apply(weights,1,sum)
   return(list(loss = bestLoss, weights = weights, prediction = experts %*% t(weights)))
}
