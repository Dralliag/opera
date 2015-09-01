#' Ridge aggregation rule with automatic calibration of smoothing parameters
#' 
#' The function \code{ridgeCalib} performs ridge aggregation rule with
#' automatic calibration of the smoothing parameter \code{lambda} by performing
#' an optimization on a finite grid. 
#' 
#' 
#' @param y  A vector containing the observations
#' to be predicted.
#' @param experts A matrix containing the
#' experts forecasts. Each column corresponds to the predictions proposed by an
#' expert to predict \code{Y}. It has as many columns as there are experts.
#' @param gridlambda A vector containing the initial possible values of the
#' smoothing parameter \code{lambda}.
#' @param href A number in \code{[1,period]}
#' specifying the instant in the day when the aggregation rule can update its
#' weights.  It should lie in the interval \code{c(1,period)}.
#' @param period The number of instants in
#' each day.
#' @param w0 A vector containing the prior
#' weights of the experts.
#' @param delay A positive number that indicates the number of instants before
#' the mixture has access to the true observations in \code{y}.  If \code{delay
#' > 0}, \code{y.ETR} can not be \code{NULL}.
#' @param y.ETR A vector containing real time estimations of \code{y} to be
#' used at each instant \code{t} instead of \code{y_t} to predict instants
#' \code{t+1,...,t+delay}.
#' @param trace  A boolean. If TRUE, the
#' evolution of the aggregation rule is displayed. Usefull if the code is too
#' long.
#' @return  \item{weights}{ A matrix of dimension \code{c(T,N)}, with
#' \code{T} the number of instances to be predicted and \code{N} the number of
#' experts.  Each row contains the convex combination to form the predictions.
#' } \item{prediction}{ A vector of length \code{T} that contains the quantiles
#' predictions outputted by the aggregation rule.  } \item{lambda}{ A vector of
#' length \code{T} containing the sequence of penalty coefficient chosen by the
#' aggregation rule.  } \item{grid}{ A vector containing the final grid of
#' potential penalty coefficient considered by the aggregation rule.  }
#' \item{loss}{ The error suffered by the aggregation rule determined by
#' \code{loss.type}.  If \code{loss.type = 'squareloss'}, the \link{rmse} is
#' computed.  } \item{gridloss}{ A vector of the same length as grid containing
#' the errors suffered by the \code{ridgeHour} aggregation rule if it
#' had picked the fixed learning rates in \code{grid}.  }
#' @author Pierre Gaillard <pierre-p.gaillard@@edf.fr>
#' @keywords ~kwd1 ~kwd2
ridgeCalib <-
function(y, experts, gridlambda = 1, 
                      href = 1, period = 1, w0 = NULL, 
                      delay = 0, y.ETR = NULL,
                      trace = F)
{
  experts <- as.matrix(experts)
  
  N <- ncol(experts)  # Number of experts
  T <- nrow(experts)  # Number of instants
  
  if (is.null(w0)) {w0 <- matrix(1/N, ncol = N)} # Uniform intial weight vector if unspecified
  if (sum(is.na(experts)) > 0) {warning("There are not allowed NA's in expert advice")}
  if (is.null(y.ETR) || (delay == 0)) {y.ETR = y}

  # Grille du paramètre lambda et perte cumulée
  nlambda <- length(gridlambda)
  bestlambda <- floor(nlambda)/2 + 1 # Pour commencer on choisi le lambda au milieu de la grille
  gridlambda <- matrix(gridlambda, nrow = nlambda)
  
  lambda <- rep(gridlambda[bestlambda],T)
  cumulativeLoss <- rep(0,nlambda)
  
  wlambda <- array(w0, dim = c(N, nlambda))    # Matrice des poids donnés par chaque ewa(lambda) mis à jour chaque instant
  wlambdaop <- array(w0, dim = c(N, nlambda))  # Matrice des poids donnés par chaque ewa(lambda) avec mis à jour opérationnelle 

  weights <- matrix(0, ncol = N, nrow = T)    # Matrix of weights formed by the mixture
  prediction <- rep(0, T)                # Prévisions du mélange

  pred.lambda.op <- matrix(0, ncol = nlambda, nrow = T) # Prediction of mixture algorithm with different learning rates eta

  At <- diag(0, N)
  bt <- matrix(t(w0) %*% gridlambda, nrow = N, ncol = nlambda)

  for(t in 1:T){
    # Affichage de l'avancement de l'aggregation rulee
    if (!(t %% floor(T/10)) && trace) cat(floor(10 * t/T)*10, '% -- ')
    
    # Poids, prévision et perte opérationnels
    weights[t,] <- wlambdaop[,bestlambda] 
    prediction[t] <- experts[t,] %*% weights[t,]
    lambda[t] <- gridlambda[bestlambda]
    
    # if (t == 2) { browser()}
    # Prévisions et pertes opérationnelles pour chaque lambda de la grille
    pred.lambda.op[t,] <- experts[t,] %*% wlambdaop 
    cumulativeLoss <- cumulativeLoss + (pred.lambda.op[t,] - y.ETR[t])^2
    
    # Mise à jour
    At <- At + experts[t,] %*% t(experts[t,])
    bt <- bt + y.ETR[t] * experts[t,]

    if ((delay > 0) && (t > delay)) {
      bt <- bt +  (y[t-delay] - y.ETR[t-delay]) * experts[t,]
      cumulativeLoss <- cumulativeLoss + (pred.lambda.op[t-delay,] - y[t-delay])^2 - (pred.lambda.op[t-delay,] - y.ETR[t-delay])^2
    }

    # Si c'est l'heure de mise à jour on met à jour wop et le paramètre lambda et la grille
    h <- (((t - 1) %% period) + 1)
    if (h == href) {
      # Mise à jour de la grille
      bestlambda <- order(cumulativeLoss)[1]

      # On aumente la grille si on est sur une extrémité
      if (bestlambda == nlambda) { 
        if (trace) cat(' + ')
        newlambda <- gridlambda[bestlambda] * 10^(1:3)
        gridlambda <- c(gridlambda, newlambda)
        nlambda <- nlambda + length(newlambda)
        for (k in 1:length(newlambda)) {
          perfnewlambda <- tryCatch(ridgeHour(y[1:t], matrix(experts[1:t,],ncol=N), newlambda[k], href = href, period = period, w0 = w0, delay = delay, y.ETR = y.ETR), error = function(e) {list(prediction = rep(0, T))})
          signal <-  y[1:t]
          if (delay > 0) { signal[(t-delay+1):t] = y.ETR[(t-delay+1):t] }
          newcumulativeLoss <- sum((perfnewlambda$prediction - signal)^2)
          cumulativeLoss <- c(cumulativeLoss, newcumulativeLoss)
          pred.lambda.op <- cbind(pred.lambda.op, c(perfnewlambda$prediction, rep(0, (T-t))))
        }
      }
      if (bestlambda == 1) {
        if (trace) cat(' - ')
        newlambda <- gridlambda[bestlambda] / 10^1
        nlambda <- nlambda + length(newlambda)
        bestlambda <- bestlambda + length(newlambda)
        for (k in 1:length(newlambda)) {
          gridlambda <- c(newlambda[k],gridlambda)
          perfnewlambda <- tryCatch(ridgeHour(y[1:t], matrix(experts[1:t,],ncol = N), newlambda[k], href = href, period = period, w0 = NULL, delay = delay, y.ETR = y.ETR), error = function(e) {list(prediction = rep(0, T))})
          signal <-  y[1:t]
          if (delay > 0) { signal[(t-delay+1):t] = y.ETR[(t-delay+1):t] }
          newcumulativeLoss <- sum((perfnewlambda$prediction - signal)^2)
          cumulativeLoss <- c(newcumulativeLoss, cumulativeLoss)
          pred.lambda.op <- cbind(c(perfnewlambda$prediction, rep(0, (T-t))), pred.lambda.op)
        }
      }
      wlambdaop <- matrix(0,nrow = N, ncol = nlambda)
      for (k in 1:nlambda) {
        wlambdaop[,k] = tryCatch(solve(gridlambda[k]*diag(1,N) + At,bt),
                                 error = function(e) {0})
      }
    }
  }
  l <- rmse(prediction,y)
  rmse <- sqrt(cumulativeLoss / T)
  if (trace) cat('\n')
  return(list(weights = weights, prediction = prediction, 
              lambda = lambda, grid = gridlambda, 
              loss = l, gridloss = rmse, last.weights = wlambdaop[,bestlambda] ))
}
