#' Fixed-share aggregation rule with automatic tunning of the parameters
#' 
#'  The
#' function \code{fixedshareCalib} performs \code{fixedshareHour}
#' aggregation rule with automatic calibration of the learning parameter
#' \code{eta} and \code{alpha} by performing an optimization on a finite grid.
#' See \code{fixedshare} and \code{fixedshareHour} for more
#' details.
#' 
#' 
#' @param y  A vector containing the observations
#' to be predicted.
#' @param experts A matrix containing the experts
#' forecasts. Each column corresponds to the predictions proposed by an expert
#' to predict \code{Y}. It has as many columns as there are experts.
#' @param awake A matrix specifying the
#' activation coefficients of the experts. Its entries lie in \code{[0,1]}.
#' Needed if some experts are specialists and do not always form and suggest
#' prediction.  If the expert number \code{k} at instance \code{t} does not
#' form any prediction of observation \code{Y_t}, we can put
#' \code{awake[t,k]=0} so that the mixture does not consider expert \code{k} in
#' the mixture to predict \code{Y_t}.  A matrix containing the activation
#' coefficient of the experts. Its entries lie in \code{[0,1]}.  Needed if some
#' experts are specialists and do not always form and suggest prediction.  If
#' the expert number \code{k} at instance \code{t} does not form any prediction
#' of observation \code{Y_t}, we can put \code{awake[t,k]=0} so that the
#' mixture does not consider expert \code{k} in the mixture to predict
#' \code{Y_t}.
#' @param grideta A vector containing the
#' initial grid of allowed learning parameters to consider. It will be extended
#' if its borders perform well.
#' @param gridalpha A vector specifying the
#' allowed mixing parameters considered by the aggregation rule. It will not be
#' changed.
#' @param loss.type Loss function
#' considered to evaluate the performance. It can be "squareloss", "mae",
#' "mape", or "pinballloss". See \code{\link{loss}} for more details.
#' @param loss.gradient A boolean. If
#' TRUE (default) the aggregation rule will not be directly applied to the loss
#' function at hand but to a gradient version of it. The aggregation rule is
#' then similar to gradient descent aggregation rule.
#' @param w0 prior weights vector of the
#' experts.
#' @param href A number in \code{[1,period]}
#' specifying the instant in the day when the aggregation rule can update its
#' weights.  It should lie in the intervall \code{c(1,period)}.
#' @param period The number of instants in
#' each day.
#' @param trace A boolean. If TRUE, the
#' evolution of the aggregation rule is displayed. Usefull if the code is too
#' long.
#' @return  \item{weights }{a matrix of dimension \code{c(T,N)}, with
#' \code{T} the number of instances to be predicted and \code{N} the number of
#' experts. Each row contains the convex combination to form the predictions}
#' \item{prediction }{ A vector of length \code{T} that contains the
#' predictions outputted by the aggregation rule.  } \item{par}{the sequence of
#' parameters \code{eta} and \code{alpha} chosen by the aggregation rule}
#' \item{grideta}{the final grid of potential learning parameter \code{eta}
#' considered by the aggregation rule} \item{gridalpha}{the final grid of
#' potential mixing parameter \code{alpha} considered by the aggregation rule}
#' \item{loss}{the error suffered by the aggregation rule determined by
#' \code{loss.type}.  If \code{loss.type = 'squareloss'}, the \link{rmse} is
#' computed.} \item{gridloss}{errors suffered by the
#' \code{fixedshareHour} aggregation rule if it had picked the fixed
#' learning rates in \code{gridalpha} and \code{grideta}}
#' @author Pierre Gaillard <pierre-p.gaillard@@edf.fr>
#' @keywords ~kwd1 ~kwd2
fixedshareCalib <-
function(y, experts, grideta = 1, gridalpha = 10^(-4:-1), awake = NULL,
                            loss.type = 'squareloss', loss.gradient = TRUE, w0 = NULL, href = 1, period = 1, trace = F)
{
   experts <- as.matrix(experts)
   
   N <- ncol(experts)  # Number of experts
   T <- nrow(experts)  # Number of instants
   
   if (is.null(w0)) {w0 <- rep(1,N)} # Uniform intial weight vector if unspecified
   if (is.null(awake)) {awake = matrix(1, nrow = T, ncol = N)} # Full activation if unspecified

   awake = as.matrix(awake) 
   idx.na <- which(is.na(experts))
   awake[idx.na] <- 0
   experts[idx.na] <- 0
   
   # Grille du paramètre eta et alpha, on commence par le milieu de la grille
   neta <- length(grideta)
   nalpha <-length(gridalpha)
   bestpar <- c(floor(neta)/2,floor(nalpha)/2)+1
   par <- NULL
   
   cumulativeLoss <- array(0, dim = c(neta, nalpha))
   
   # Poids et regret
   wpar <- array(w0, dim = c(N,neta,nalpha))
   wparop <- array(w0, dim = c(N,neta,nalpha))
   
   weights <- matrix(0, ncol = N, nrow = T)    # Matrix of weights formed by the mixture
   prediction <- rep(0, T)                # Prévisions du mélange
   
   test = rep(0, T)
   
   for(t in 1:T){
      # Affichage de l'avancement de l'aggregation rulee
      if (!(t %% floor(T/10)) && trace) cat(floor(10 * t/T)*10, '% -- ')
      
      # Poids prévisions et pertes opérationnels
      weights[t,] <- wparop[,bestpar[1],bestpar[2]] * awake[t,] / sum(wparop[,bestpar[1],bestpar[2]] * awake[t,])
      prediction[t] <- experts[t,] %*% weights[t,]
      par <- rbind(par, data.frame(eta = grideta[bestpar[1]],
                                   alpha = gridalpha[bestpar[2]]))
      
      # boucle sur les paramètres alpha
      for (k in 1:nalpha) {
         
         # Prévision et pertes opérationnelles pour chaque (eta,alpha) de la grille
         waux <-  t(t(wparop[,,k] * awake[t,]) / apply(as.matrix(wparop[,,k]*awake[t,]), 2, sum))
         predop <- experts[t,] %*% waux
         lpredop <- loss(predop, y[t], loss.type)
         cumulativeLoss[,k] <- cumulativeLoss[,k] + lpredop
         
         # Prévision et pertes non opérationnelles des algo et des experts
         waux <- t(t(awake[t,] * wpar[,,k]) / apply(as.matrix(awake[t,] * wpar[,,k]), 2, sum))
         pred <- experts[t,] %*% waux
         lpred <- diag(lossPred(pred, y[t], pred, loss.type, loss.gradient))
         lexp <- lossPred(experts[t,], y[t], pred, loss.type, loss.gradient)
         
         # Mise à jour des regrets et des poids
         R <- t(t(log(wpar[,,k]))/grideta) + awake[t,] * t(lpred - t(lexp))
         v <- truncate1(exp(t(t(matrix(R, ncol=neta)) * grideta)))
         v <- t(t(v) / apply(v,2,sum)) # on renormalise chaque colonne
         
         # Mise à jour des poids
         wpar[,,k] <- gridalpha[k]/N + (1-gridalpha[k]) * v
      }
      
      
      # Si c'est l'heure de mise à jour...
      h <- (((t - 1) %% period) + 1)
      if (h == href) {
         # Mise à jour de la grille
         bestpar <- which(cumulativeLoss == min(cumulativeLoss), arr.ind = TRUE)[1,]
         # On augmente la grille si on est à une extrémité
         if (bestpar[1] == neta) {
            if (trace) cat(' + ')
            neweta <- grideta[neta] * 2^(1:3)
            grideta <- c(grideta, neweta)
            neta <- neta + length(neweta)
            wparaux <- wpar
            wpar <- array(dim = c(N,neta,nalpha))
            wpar[,1:bestpar[1],] <- wparaux
            for (j in 1:length(neweta)) {
               cumulativeLoss <- rbind(cumulativeLoss, 0)
               for (k in 1:nalpha) {
                  perfnewpar <- fixedshareHour(y[1:t], matrix(experts[1:t,],ncol=N), neweta[j], gridalpha[k], 
                                             awake =  matrix(awake[1:t,],ncol=N),
                                             loss.type = loss.type, loss.gradient = loss.gradient, 
                                             href = href, period = period, w0 = w0)
                  cumulativeLoss[bestpar[1]+j,k] <- perfnewpar$cumulativeLoss
                  wpar[,bestpar[1]+j,k] <- perfnewpar$lastweight
               }
            }
         }
         if (bestpar[1] == 1) {
            if (trace) cat(' - ')
            neweta <- grideta[1] / 2^(1:3)
            neta <- neta + length(neweta)
            bestpar[1] <- bestpar[1] + length(neweta)
            wparaux <- wpar
            wpar <- array(dim = c(N,neta,nalpha))
            wpar[,bestpar[1]:neta,] <- wparaux
            for (j in 1:length(neweta)) {
               grideta <- c(neweta[j], grideta)
               cumulativeLoss <- rbind(0, cumulativeLoss)
               for (k in 1:nalpha) {
                  perfnewpar <- fixedshareHour(y[1:t], matrix(experts[1:t,],ncol=N),
                                             neweta[j], gridalpha[k],
                                             awake = matrix(awake[1:t,],ncol=N),
                                             loss.type = loss.type, loss.gradient = loss.gradient,
                                             href = href, period = period, w0 = w0)
                  cumulativeLoss[1,k] <- perfnewpar$cumulativeLoss
                  wpar[,bestpar[1]-j,k] <- perfnewpar$lastweight
               }
            }
         }
         wparop <- wpar
      } else {
         for (k in 1:nalpha) {
            wparop[,,k] <- gridalpha[k]/N + (1-gridalpha[k]) * wparop[,,k]
         }
      }
   }
  l <-  mean(loss(prediction,y,loss.type=loss.type))
   mloss <- cumulativeLoss / T
   if (loss.type == 'squareloss') {
      mloss <- sqrt(mloss)
      l <- sqrt(l)
   }
   rownames(mloss) = grideta
   colnames(mloss) = gridalpha
   
   if (trace) cat('\n')
   return(list(weights = weights, prediction = prediction, 
               par = par, grideta = grideta, gridalpha = gridalpha,
               loss = l, gridloss = mloss))
}
