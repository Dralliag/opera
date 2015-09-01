#' best linear oracle
#' 
#'  The
#' function \code{bestLinear} computes the best linear combination oracle or in
#' other words the fixed linear combination of experts in \code{experts} that
#' would have perfomed the best to predict the observations in \code{y}.
#' 
#' 
#' @param y  vector that contains the observations
#' to be predicted.
#' @param experts A matrix containing the
#' experts forecasts. Each column corresponds to the predictions proposed by an
#' expert to predict \code{Y}. It has as many columns as there are experts.
#' @param lambda Smoothing parameter. No
#' smoothing corresponds to \code{lambda = 0}, however the linearn system might
#' be singular.
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
#' @return \item{loss}{ The loss suffered by the best fixed linear combination.
#' } \item{u}{ A vector containing the best fixed linear combination chosen in
#' hindsight.  } \item{prediction}{ A vecor containing the predictions of the
#' best fixed linear combination.  }
#' @author Pierre Gaillard <pierre-p.gaillard@@edf.fr>
#' @seealso 
#' \code{\link{bestConvex}}, \code{\link{bestShifts}}, \code{\link{loss}}
#' @keywords ~kwd1 ~kwd2
#' @export bestLinear
bestLinear <-
function(y,experts, lambda = 1, awake=NULL, loss.type='squareloss')
{
  if (!is.null(awake)) {stop('Sleeping not allowed here!')}
  if (sum(is.na(experts)>0)) {warning("There are NA's in expert advice")}
  
  experts <- as.matrix(experts)
  u <- solve(lambda * diag(1,ncol(experts)) + t(experts) %*% experts,t(experts)%*%y)
  prev <- experts %*% u
  loss <- sqrt(mean((prev-y)^2))
  return(list(loss=loss,prediction=prev,u=u))
}
