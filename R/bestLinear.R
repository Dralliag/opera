# best linear oracle
bestLinear <-
function(y, experts, lambda = 0, loss.type = "square", 
	tau = 0.5, niter = 1, ...)
{  
  experts <- as.matrix(experts)
  N = ncol(experts)

  coefficients <- NULL
  if (loss.type == "square") {
  	coefficients <- solve(lambda * diag(1,ncol(experts)) + t(experts) %*% experts,t(experts)%*%y)
  	
  } else if (loss.type == "pinball") {
  		coefficients <- tryCatch({
  						quantreg::rq(y~experts-1,tau = tau)$coefficients},
  						error = function(e) { NULL 
  					})
  } 
  if (is.null(coefficients)) {
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
      coefficients = matrix(best_u, nrow = N)
  }
 
	prev <- experts %*% coefficients	
  loss <- mean(loss(prev, y, loss.type = loss.type, tau = tau))
  return(list(loss=loss,coefficients=c(coefficients),prediction=prev,
  	rmse=sqrt(mean(loss(prev,y)))))
}
