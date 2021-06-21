#' Errors suffered by a sequence of predictions
#' 
#'  The
#' function \code{loss} computes the sequence of instantaneous losses suffered
#' by the predictions in \code{x} to predict the observation in \code{y}.
#' 
#' 
#' @param x A vector of length \code{T}
#' containing the sequence of prediction to be evaluated.
#' @param y  A vector of length \code{T} that
#' contains the observations to be predicted.
#' @param loss.type \code{character, list or function}. 
#' \describe{
#'      \item{character}{ Name of the loss to be applied ('square', 'absolute', 'percentage', or 'pinball');}
#'      \item{list}{ When using pinball loss: list with field name equal to 'pinball' and field tau equal to the required quantile in [0,1];}
#'      \item{function}{ A custom loss as a function of two parameters.}
#' }
#' @return  A vector of length \code{T} containing the sequence of
#' instantaneous losses suffered by the prediction \code{x}.
#' @author Pierre Gaillard <pierre@@gaillard.me>
#' @export loss
loss <- function(x, y, loss.type = "square") {
  
  if (class(loss.type) == "function") {
    args <- formalArgs(loss.type)
    
    if (length(args) > 2) {
      stop("The provided loss function should contain exactly 2 arguments.")
    }
    
    l <- tryCatch({loss.type(x, y)}, 
                  error = function(e) {
                    stop("Error when trying to apply custom loss function : \n",
                         e$message)
                  })
    
  } else {
    if (!is.list(loss.type)) {
      loss.type <- list(name = loss.type) 
    }
    if (is.null(loss.type$tau) && loss.type$name == "pinball") {
      loss.type$tau <- 0.5
    }
    
    if (loss.type$name == "square") 
      l <- (x - y)^2 else if (loss.type$name == "absolute") 
        l <- abs(x - y) else if (loss.type$name == "percentage") 
          l <- abs(x - y)/y else if (loss.type$name == "pinball") 
            l <- (loss.type$tau - (y < x)) * (y - x) else if (loss.type$name == "log") 
              l <- -log(x) else stop("loss.type should be one of these: 'absolute', 'percentage', 'square', 'pinball'")
  }
  
  return(l)
} 
