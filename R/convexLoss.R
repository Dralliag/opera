#' agerage error suffered by convex combination
#' 
#'  The
#' function \code{convexLoss} computes the loss incured by a fixed convex
#' combination of experts.
#' 
#' 
#' @param weights vector of weights
#' assigned to each experts. It must be the same length as the number of
#' columns of \code{experts}. It entries must be non-negative and sum to one.
#' @param y  A vector containing the observations
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
#' @return A dataframe containing several measures of losses suffered
#' by the mixture
#' @author Pierre Gaillard <pierre-p.gaillard@@edf.fr>
#' @seealso 
#' \code{\link{loss}}
#' @keywords ~kwd1 ~kwd2
#' @export convexLoss
convexLoss <-
function(weights,y,experts,awake=NULL) {
  N <- ncol(experts)
  T <- nrow(experts)

  # Experts are always active if awake is unspecified
  if (is.null(awake)) {awake = matrix(1, nrow = T, ncol = N)} 
  awake <- as.matrix(awake)
  idx.na <- which(is.na(experts))
  awake[idx.na] <- 0
  experts[idx.na] <- 0

  pond <- awake %*% weights
  pred <- (experts * awake) %*% weights / pond

  err1 <- loss(pred,y,'square') * pond
  err2 <- loss(pred,y,'absolute') * pond 
  err3 <- loss(pred,y,'percentage') * pond
  
  mpred <- mean(pred)
  mY <- mean(y)
  merr1 <- sum(err1)  / sum(pond)
  merr2 <- sum(err2)  / sum(pond)
  merr3 <- sum(err3)  / sum(pond)
  
  corr <- sum( (pred - mpred) * (y - mY)) / sqrt(sum((pred - mpred)^2) * sum((y - mY)^2))
  disp1 <- 1.96 / sqrt(T) * sqrt(mean((err1 - merr1)^2) / (4 * merr1))
  disp2 <- 1.96 / sqrt(T) * sd(err2)
  disp3 <- 1.96 / sqrt(T) * sd(err3)
  return(data.frame(rmse = sqrt(merr1), mae = merr2, mape = merr3, disp_rmse = disp1, disp_mae = disp2, disp_mape = disp3, corr = corr))
}
