#' best convex oracle
#' 
#'  The
#' function \code{bestConvex} computes the best convex combination oracle, or
#' in other words the fixed convex combination of experts in \code{experts}
#' that would have perfomed the best to predict the observations in \code{y}.
#' 
#' 
#' @param y  vector that contains the observations
#' to be predicted.
#' @param experts A matrix containing the
#' experts forecasts. Each column corresponds to the predictions proposed by an
#' expert to predict \code{Y}. It has as many columns as there are experts.
#' @param awake A matrix specifying the
#' activation coefficients of the experts. Its entries lie in \code{[0,1]}.
#' Needed if some experts are specialists and do not always form and suggest
#' prediction.  If the expert number \code{k} at instance \code{t} does not
#' form any prediction of observation \code{Y_t}, we can put
#' \code{awake[t,k]=0} so that the mixture does not consider expert \code{k} in
#' the mixture to predict \code{Y_t}.
#' @param loss.type A string specifying
#' the loss function considered to evaluate the performance. It can be
#' "squareloss", "mae", "mape", or "pinballloss". See \code{\link{loss}} for
#' more details.
#' @param niter Number of times the convex
#' optimisation process is repeated with different initial values. Note that if
#' the experts are always activated, it does not need to be greater than 1.
#' @param method The optimization method to be used. See \code{\link{optim}}.
#' @param control Additional parameters
#' that are passed to \code{optim} function.
#' @return \item{loss}{ The loss suffered by the best fixed convex combination.
#' } \item{weights}{ A vector containing the best fixed convex combination
#' chosen in hindsight.  } \item{prediction}{ A vecor containing the
#' predictions of the best fixed convex combination.  }
#' @author Pierre Gaillard <pierre-p.gaillard@@edf.fr>
#' @seealso 
#' \code{\link{bestLinear}}, \code{\link{bestShifts}}, \code{\link{loss}},
#' \code{\link{lossConv}}
#' @keywords ~kwd1 ~kwd2
#' @export bestConvex
bestConvex <-
function(y, experts, awake=NULL, loss.type='squareloss', niter = 1, method = "L-BFGS-B", control = list(maxit = 100, trace=T))
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
