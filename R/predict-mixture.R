#' Predict method for Mixture models
#' 
#' Performs sequential predictions and updates
#' of a mixture object based on new observations 
#' and expert advice.
#' 
#' @param object Object of class inheriting from 'mixture'
#' 
#' @param newexperts An optional matrix in which to look for expert advice with which
#' predict. If omitted, the past predictions of the object are returned and the
#' object is not updated.
#' 
#' @param newY An optional matrix with d columns (or vector if \eqn{d=1}) of observations to be predicted. If provided, it 
#' should have the same number of rows as the number of rows of \code{newexperts}.
#' If omitted, the object (i.e, the aggregation rule) is not updated.
#' 
#' @param awake An optional array specifying the
#' activation coefficients of the experts. It must have the same dimension as experts. Its entries lie in \code{[0,1]}.
#' Possible if some experts are specialists and do not always form and suggest
#' prediction. If the expert number \code{k} at instance \code{t} does not
#' form any prediction of observation \code{Y_t}, we can put
#' \code{awake[t,k]=0} so that the mixture does not consider expert \code{k} in
#' the mixture to predict \code{Y_t}.
#' 
#' @param online A boolean determining if the observations in newY are predicted
#' sequentially (by updating the object step by step) or not. If FALSE, 
#' the observations are predicting using the object (without using any past 
#' information in newY). If TRUE, newY and newexperts should not be null.
#' 
#' @param type Type of prediction. It can be 
#' \describe{
#'    \item{model}{return the updated version of object (using newY and newexperts).}
#'    \item{response}{return the forecasts. If type is 'model', forecasts can also 
#'    be obtained from the last values of object$prediction.}
#'    \item{weights}{return the weights assigned to the expert advice to 
#'    produce the forecasts. If type is 'model', forecasts can also 
#'    be obtained from the last rows of object$weights.}
#'    \item{all}{return a list containing 'model', 'response', and 'weights'.}
#'    }
#' 
#' @param ...  further arguments are ignored
#' 
#' @return \code{predict.mixture} produces a matrix of predictions 
#' (type = 'response'), an updated object (type = 'model'), or a matrix of
#' weights (type = 'weights').
#' 
#' @export 
predict.mixture <- function(object, newexperts = NULL, newY = NULL, awake = NULL, 
                            online = TRUE, type = c("model", "response", "weights", "all"), ...) {
  
  result <- object
  d <- object$d
  if ((d == 1) || (d == "unknown" && is.null(dim(newY)))) {
    object$d <- 1
    return(predictReal(object, newexperts, newY, awake, 
                online, type, ...))
  } else {
    if (d == "unknown") {
      d = dim(newY)[2]
      T = dim(newY)[1]
      # Bad dimension for experts
      if (T > 1 && length(dim(newexperts)) < 3) {
        stop("Bad dimensions: nrow(experts) should be equal to dim(experts)[3]")
      } 
      if (length(dim(newexperts)) == 3) {
        if ((dim(newexperts)[1] != T) || (dim(newexperts)[2] != d)){
          stop("Bad dimensions between Y and experts")
        }
      }
      if (T == 1) {
        if (length(dim(newexperts)) == 2) {
          if (dim(newexperts)[1] != d) {
            stop("Bad dimensions between Y and experts")
          } else {
            newexperts = array(newexperts, dim = c(1,dim(newexperts)))
          }
        }
      }
    }
    result$d <- d
    awakei <- NULL
    for (i in 1:nrow(newY)) {
      if (!online){
        stop("Batch prediction are currently not supported for dimension > 1")
      }
      if (!is.null(awake)){
        awakei <- as.matrix(awake[i,,])
      }
      result <- predictReal(result, newexperts = as.matrix(newexperts[i,,]), newY = c(newY[i,]), awake = awakei, 
                            online = FALSE, type, ...)
    }
  }
  result$weights <- matrix(result$weights, nrow = result$T)
  
  return(result)
}
