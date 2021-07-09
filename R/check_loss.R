#' Function to check validy of provided loss function
#'
#' @param loss.type \code{character, list or function}. 
#' \describe{
#'      \item{character}{ Name of the loss to be applied ('square', 'absolute', 'percentage', or 'pinball');}
#'      \item{list}{ When using pinball loss: list with field name equal to 'pinball' and field tau equal to the required quantile in [0,1];}
#'      \item{function}{ A custom loss as a function of two parameters.}
#' }
#' @param loss.gradient \code{boolean, function}. 
#' \describe{
#'      \item{boolean}{ If TRUE, the aggregation rule will not be directly applied to the loss function at hand,
#'      but to a gradient version of it. The aggregation rule is then similar to gradient descent aggregation rule. }
#'      \item{function}{ If loss.type is a function, the derivative should be provided to be used (it is not automatically 
#'      computed).}
#' }
#' @param Y \code{numeric} (NULL). (Optional) Target values (to perform some checks).
#' @param model \code{character} (NULL). (Optional) Model used (to perform some checks).
#' @param use_cpp \code{boolean}. Whether or not to use cpp function to increase perf.
#' 
#' @importFrom methods formalArgs
#' 
#' @return loss.type
#'
check_loss <- function(loss.type, 
                       loss.gradient,
                       Y = NULL,
                       model = NULL,
                       use_cpp = getOption("opera_use_cpp", default = TRUE)) {
  
  if (! is.function(loss.type)) {
    if (! is.logical(loss.gradient)) {
      stop("loss.gradient should be a boolean when loss.type is not a function.")
    }
    if (is.vector(loss.type) && ! is.list(loss.type)) {
      loss.type <- list("name" = loss.type) 
    }
    if (! (loss.type$name %in% c("pinball", "square", "percentage", "absolute"))) {
      stop("Predifined loss.type should be one of these: ['absolute', 'percentage', 'square', 'pinball'].")
    }
    if (! is.null(loss.type$tau) && loss.type$name != "pinball") {
      warning("Unused parameter tau (loss.type != 'pinball').")
    }
    if (loss.type$name == "pinball" && is.null(loss.type$tau)) {
      loss.type$tau <- 0.5
    }
    if (! is.null(Y) && min(Y) <= 0 && loss.type$name == "percentage") {
      stop("Y should be non-negative for percentage loss function.")
    }
    if (loss.type$name != "square" && ! is.null(model) && model == "Ridge") {
      stop(paste("Square loss is require for Ridge model."))
    }
  } 
  else {
    if (use_cpp == TRUE && class(loss.type) == "function") {
      stop("Custom loss functions are not yet available when use_cpp == TRUE.") 
    }
    
    args_loss <- formalArgs(loss.type)
    if (! length(args_loss) == 2) {
      stop("loss.type should be a function of 2 arguments.")
    }
    
    if (! is.function(loss.gradient) && loss.gradient == TRUE) {
      stop("When loss.type is a function and you want to use the gradient version, you must provide it.")
      
    } else if (is.function(loss.gradient)) {
      args_loss_grad <- formalArgs(loss.gradient)
      if (! length(args_loss_grad) == 2) {
        stop("loss.gradient should be a function of 2 arguments.")
      } 
    } else if (! is.logical(loss.gradient)) {
      stop("when loss.type is a function, loss.gradient should either be a function or set to FALSE.")
    }
    # check if provided loss is vectorial
    is_vect <- tryCatch({
      ctrl_1 <- loss.type(1, 2)
      ctrl_2 <- loss.type(3.14, 2)
      ctrl_3 <- loss.type(c(1, 3.14), 2)
      all(ctrl_3 == c(ctrl_1, ctrl_2))
    }, error = function(e) {
      stop("Error when calling provided loss : \n", 
           e$message)
    })
    if (! all(is_vect)) {
      stop("The provided loss must be vectorial.")
    }
  }
  
  return(loss.type)
}