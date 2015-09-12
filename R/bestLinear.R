# best linear oracle
bestLinear <-
function(y, experts, lambda = 0, loss.type = "square", 
	tau = 0.5, niter = 1, ...)
{  
  experts <- as.matrix(experts)
  N = ncol(experts)

  weights <- NULL
  if (loss.type == "square") {
  	weights <- solve(lambda * diag(1,ncol(experts)) + t(experts) %*% experts,t(experts)%*%y)
  	
  } else if (loss.type == "pinball") {
  		weights <- tryCatch({
  						quantreg::rq(y~experts-1,tau = tau)$coefficients},
  						error = function(e) { NULL 
  					})
  } 
  if (is.null(weights)) {
  		warning("The best linear oracle is only approximated (using optim).")
  	 lossu <- function(u)   {
         return(mean(loss(
         		x = experts %*% matrix(u,nrow=ncol(experts)), 
         		y = y, 
         		loss.type = loss.type, 
         		tau = tau))) 
      }
      
      best_u <- rep(0,N)
      bestLoss <- exp(700)
      
      for (i in 1:niter)
      {
         # Optimisation convexe avec choix aléatoire de la condition initiale   
         u <- rnorm(N,0,1)
         w <- optim(u,lossu, gr = NULL, ...)$par
         l <- lossu(w)
         if (bestLoss > l) {
            bestLoss = l
            best_u = w
         }
      }
      weights = matrix(best_u, nrow = N)
  }
 
	prev <- experts %*% weights	
  loss <- mean(loss(prev, y, loss.type = loss.type, tau = tau))
  return(list(loss=loss,weights=c(weights),prediction=prev,
  	rmse=sqrt(mean(loss(prev,y)))))
}
