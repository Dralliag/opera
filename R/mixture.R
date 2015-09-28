#' Compute an aggregation rule
#' 
#' The function \code{mixture} builds an
#' aggregation rule chosen by the user. 
#' It can then be used to predict new observations Y sequentially.
#' If observations \code{Y} and expert advice \code{experts} are provided, 
#' \code{mixture} is trained by predicting the observations in \code{Y}
#' sequentially with the help of the expert advice in \code{experts}.  
#' At each time instance \code{t}, the mixture forms a prediction by assigning 
#' a weight to each expert and by combining the expert advice.
#' 
#' 
#' @param Y  A vector of length T (a non negative integer) containing the observations to be predicted sequentially
#' in order to train the aggregation rule.
#'  
#' @param experts A matrix containing the expert
#' forecasts. Each column corresponds to the predictions proposed by an expert
#' to predict \code{Y}.  It has as many columns as there are experts.
#' Its number of row should be \code{T}.
#' 
#' @param model A character string specifying the aggregation rule to use. 
#' Currently available aggregation rules are:
#' \describe{
#'    \item{'EWA'}{Exponentially weighted average aggregation rule. A positive learning rate \strong{eta} 
#' can be chosen by the user. The
#' bigger it is the faster the aggregation rule will learn from observations
#' and experts performances. However, too hight values lead to unstable weight
#' vectors and thus unstable predictions. If it is not specified, the learning rate is calibrated online. 
#' A finite grid of potential learning rates to be optimized online can be specified with \strong{grid.eta}.}
#'    \item{'FS'}{Fixed-share aggregation rule. As for \code{ewa}, a learning rate \strong{eta} 
#' can be chosen by the user or calibrated online. The main difference with \code{ewa} aggregation
#' rule rely in the mixing rate \strong{alpha}\eqn{\in [0,1]} wich considers at
#' each instance a small probability \code{alpha} to have a rupture in the
#' sequence and that the best expert may change. Fixed-share aggregation rule
#' can thus compete with the best sequence of experts that can change a few
#' times (see \code{\link{oracle}}), while \code{ewa} can only
#' compete with the best fixed expert. The mixing rate \strong{alpha} is either chosen by the user either calibrated online.
#' Finite grids of learning rates and mixing rates to be optimized can be specified with 
#' parameters \strong{grid.eta} and \strong{grid.alpha}.}
#'    \item{'Ridge'}{Ridge regression. It minimizes at
#' each instance a penalized criterion.  It forms at each instance linear
#' combination of the experts' forecasts and can assign negative weights that
#' not necessarily sum to one.  It is useful if the experts are biased or
#' correlated. It cannot be used with specialized experts. A positive regularization coefficient \strong{lambda} 
#' can either be chosen by the user or calibrated online. 
#' A finite grid of coefficient to be optimized can be specified with a parameter \strong{grid.lambda}.}
#'    \item{'MLpol'}{Polynomial Potential aggregation rule
#' with different learning rates for each expert.  The learning rates are
#' calibrated using theoretical values. There are similar aggregation rules 
#' like 'BOA' (Bernstein online Aggregation see [Wintenberger, 2014] 'MLewa', and 'MLprod' 
#' (see [Gaillard, Erven, and Stoltz, 2014])}
#' }
#' 
#' @param loss.type A string or a list with a component 'name' specifying
#' the loss function considered to evaluate the performance. It can be
#' 'square', 'absolute', 'percentage', or 'pinball'. In the case of the pinball loss, the quantile 
#' can be provided by assigning to loss.type a list of two elements: 
#' \describe{
#'      \item{name}{A string defining the name of the loss function (i.e., 'pinball')}
#'      \item{tau}{ A number in \code{[0,1]} defining the quantile to be predicted. 
#' The default value is 0.5 to predict the median.}
#' }. 'Ridge' aggregation rule is restricted to square loss.
#' 
#' @param loss.gradient A boolean. If
#' TRUE (default) the aggregation rule will not be directly applied to the loss
#' function at hand but to a gradient version of it.  The aggregation rule is
#' then similar to gradient descent aggregation rule. 
#' 
#' @param coefficients A vector containing the prior weights of the experts
#' (not possible for 'MLpol').
#' 
#' @param awake A matrix specifying the
#' activation coefficients of the experts. Its entries lie in \code{[0,1]}.
#' Possible if some experts are specialists and do not always form and suggest
#' prediction. If the expert number \code{k} at instance \code{t} does not
#' form any prediction of observation \code{Y_t}, we can put
#' \code{awake[t,k]=0} so that the mixture does not consider expert \code{k} in
#' the mixture to predict \code{Y_t}.
#' 
#' @param parameters A list that contains optional parameters for the aggregation rule. 
#' If no parameters are provided, the aggregation rule is fully calibrated
#' online. Possible parameters are:
#' \describe{
#'    \item{eta}{A positive number defining the learning rate. 
#'    Possible if model is either 'EWA' or 'FS'}
#'    \item{grid.eta}{A vector of positive numbers defining potential learning rates 
#'    for 'EWA' of 'FS'.
#'    The learning rate is then calibrated by sequentially optimizing the parameter in 
#'    the grid. The grid may be extended online if needed by the aggregation rule.}
#'    \item{gamma}{A positive number defining the exponential step of extension 
#'    of grid.eta when it is needed. The default value is 2.}
#'    \item{alpha}{A number in [0,1] defining the mixing rate for 'FS'.}
#'    \item{grid.alpha}{A vector of numbers in [0,1] defining potential mixing rates for 'FS'
#'    to be optimized online. The grid is fixed over time. The default value is \code{[0.0001,0.001,0.01,0.1]}.}
#'    \item{lambda}{A positive number defining the smoothing parameter of 'Ridge' aggregation rule.}
#'    \item{grid.lambda}{Similar to \code{grid.eta} for the parameter \code{lambda}.}
#' }
#' 
#' 
#' @return An object of class mixture that can be used to perform new predictions. 
#' It contains the parameters \code{model}, \code{loss.type}, \code{loss.gradient},
#' \code{experts}, \code{Y}, \code{awake}, and the fields
#' \item{coefficients}{A vector of coefficients 
#' assigned to each expert to perform the next prediction.}
#' 
#' \item{weights }{ A matrix of dimension \code{c(T,N)}, with
#' \code{T} the number of instances to be predicted and \code{N} the number of
#' experts.  Each row contains the convex combination to form the predictions }
#' \item{prediction }{ A vector of length \code{T} that contains the
#' predictions outputted by the aggregation rule.  } 
#' \item{loss}{ The average loss (as stated by parameter \code{loss.type}) suffered
#' by the aggregation rule.}
#' \item{parameters}{The learning parameters chosen by the aggregation rule.}
#' \item{training}{A list that contains useful temporary information of the 
#' aggregation rule to be updated and to perform predictions.}
#' @author Pierre Gaillard <pierre@@gaillard.me>
#' @keywords ~kwd1 ~kwd2
#' @seealso See \code{\link{opera-package}} and opera-vignette for a brief example about how to use the package.
#'  
#' @export mixture

mixture <- function(Y = NULL, experts = NULL, model = "MLpol", loss.type = "square", 
  loss.gradient = TRUE, coefficients = "Uniform", awake = NULL, parameters = list()) UseMethod("mixture")


#' @export 
mixture.default <- function(Y = NULL, experts = NULL, model = "MLpol", loss.type = "square", 
  loss.gradient = TRUE, coefficients = "Uniform", awake = NULL, parameters = list()) {
  
  if (!is.list(loss.type)) {
    loss.type <- list(name = loss.type)
  }
  if (!(loss.type$name %in% c("pinball", "square", "percentage", "absolute"))) {
    stop("loss.type should be one of these: 'absolute', 'percentage', 'square', 'pinball'")
  }
  
  object <- list(model = model, loss.type = loss.type, loss.gradient = loss.gradient, 
    coefficients = coefficients, parameters = parameters, Y = NULL, experts = NULL, 
    awake = NULL, training = NULL, names.experts = colnames(experts))
  class(object) <- "mixture"
  
  # Test that Y and experts have correct dimensions
  if ((is.null(Y) && !is.null(experts)) || (!is.null(Y) && is.null(experts))) {
    stop("Bad dimensions: length(Y) should be equal to nrow(experts)")
  }
  if (!is.null(Y)) {
    if (length(Y) == 1) {
      experts <- as.matrix(experts)
      if (nrow(experts) == 1 || ncol(experts) == 1) {
        experts <- matrix(as.numeric(experts), nrow = 1)
      } else {
        stop("Bad dimensions: length(Y) should be equal to nrow(experts)")
      }
    }
    if (!(length(Y) == nrow(experts))) {
      stop("Bad dimensions: length(Y) should be equal to nrow(experts)")
    }
    
    object <- predict(object, newY = Y, newexperts = experts, awake = awake, 
      type = "model")
    
  }
  return(object)
} 
