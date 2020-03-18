#' Compute an aggregation rule
#' 
#' The function \code{mixture} builds an
#' aggregation rule chosen by the user. 
#' It can then be used to predict new observations Y sequentially.
#' If observations \code{Y} and expert advice \code{experts} are provided, 
#' \code{mixture} is trained by predicting the observations in \code{Y}
#' sequentially with the help of the expert advice in \code{experts}.  
#' At each time instance \eqn{t=1,2,\dots,T}, the mixture forms a prediction of \code{Y[t,]} by assigning 
#' a weight to each expert and by combining the expert advice.
#' 
#' 
#' @param Y  A matrix with T rows and d columns. Each row \code{Y[t,]} contains a d-dimensional 
#' observation to be predicted sequentially.
#'  
#' @param experts An array of dimension \code{c(T,d,K)}, where \code{T} is the length of the data-set, 
#' \code{d} the dimension of the observations, and \code{K} is the number of experts. It contains the expert
#' forecasts. Each vector \code{experts[t,,k]} corresponds to the d-dimensional prediction of \code{Y[t,]} 
#' proposed by expert k at time \eqn{t=1,\dots,T}.
#' In the case of real prediction (i.e., \eqn{d = 1}), \code{experts} is a matrix with \code{T} rows and \code{K} columns.
#' 
#' @param model A character string specifying the aggregation rule to use. 
#' Currently available aggregation rules are:
#' \describe{
#'    \item{'EWA'}{Exponentially weighted average aggregation rule. A positive learning rate \strong{eta} 
#' can be chosen by the user. The
#' bigger it is the faster the aggregation rule will learn from observations
#' and experts performances. However, too high values lead to unstable weight
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
#' like 'BOA' (Bernstein online Aggregation see [Wintenberger, 2014] <doi:10.1007/s10994-016-5592-6>, 'MLewa', and 'MLprod' 
#' (see [Gaillard, Erven, and Stoltz, 2014])} 
#'  \item{'OGD'}{Online Gradient descent (see Zinkevich, 2003). The optimization is performed with a time-varying learning rate. 
#'  At time step \eqn{t \geq 1}, the learning rate is chosen to be \eqn{t^{-\alpha}}, where \eqn{\alpha} is provided by alpha in the parameters argument.
#'  The algorithm may or not perform a projection step into the simplex space (non-negative weights that sum to one) according to
#'  the value of the parameter 'simplex' provided by the user.}
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
#' } 'Ridge' aggregation rule is restricted to square loss.
#' 
#' @param loss.gradient A boolean. If
#' TRUE (default) the aggregation rule will not be directly applied to the loss
#' function at hand but to a gradient version of it.  The aggregation rule is
#' then similar to gradient descent aggregation rule. 
#' 
#' @param coefficients A probability vector of length K containing the prior weights of the experts
#' (not possible for 'MLpol'). The weights must be non-negative and sum to 1.
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
#'    \item{alpha}{A number in [0,1]. If the model is 'FS', it defines the mixing rate. 
#'    If the model is 'OGD', it defines the order of the learning rate: \eqn{\eta_t = t^{-\alpha}}.}
#'    \item{grid.alpha}{A vector of numbers in [0,1] defining potential mixing rates for 'FS'
#'    to be optimized online. The grid is fixed over time. The default value is \code{[0.0001,0.001,0.01,0.1]}.}
#'    \item{lambda}{A positive number defining the smoothing parameter of 'Ridge' aggregation rule.}
#'    \item{grid.lambda}{Similar to \code{grid.eta} for the parameter \code{lambda}.}
#'    \item{simplex}{A boolean that specifies if 'OGD' does a project on the simplex. In other words,
#'    if TRUE (default) the online gradient descent will be under the constraint that the weights sum to 1
#'    and are non-negative. If FALSE, 'OGD' performs an online gradient descent on K dimensional real space.
#'    without any projection step.}
#'    \item{averaged}{A boolean (default is FALSE). If TRUE the coefficients and the weights 
#'    returned (and used to form the predictions) are averaged over the past. It leads to more stability on the time evolution of the weights but needs
#'    more regularity assumption on the underlying process genearting the data (i.i.d. for instance). }
#' }
#' 
#' 
#' @return An object of class mixture that can be used to perform new predictions. 
#' It contains the parameters \code{model}, \code{loss.type}, \code{loss.gradient},
#' \code{experts}, \code{Y}, \code{awake}, and the fields
#' \item{coefficients}{A vector of coefficients 
#' assigned to each expert to perform the next prediction.}
#' 
#' \item{weights }{ A matrix of dimension \code{c(T,K)}, with
#' \code{T} the number of instances to be predicted and \code{K} the number of
#' experts.  Each row contains the convex combination to form the predictions }
#' \item{prediction }{ A matrix with \code{T} rows and \code{d} columns that contains the
#' predictions outputted by the aggregation rule.  } 
#' \item{loss}{ The average loss (as stated by parameter \code{loss.type}) suffered
#' by the aggregation rule.}
#' \item{parameters}{The learning parameters chosen by the aggregation rule or by the user.}
#' \item{training}{A list that contains useful temporary information of the 
#' aggregation rule to be updated and to perform predictions.}
#' @author Pierre Gaillard <pierre@@gaillard.me>
#' @keywords ~models ~ts
#' @seealso See \code{\link{opera-package}} and opera-vignette for a brief example about how to use the package.
#'  
#' @template example
#' @importFrom stats predict
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
    awake = NULL, training = NULL, names.experts = colnames(experts), T = 0, d = "unknown")
  
  class(object) <- "mixture"
  
  # Test that Y and experts have correct dimensions
  if ((is.null(Y) && !is.null(experts)) || (!is.null(Y) && is.null(experts))) {
    stop("Bad dimensions: length(Y) should be equal to nrow(experts)")
  }
  
  if (!is.null(Y)) {
    
    # Test the dimension of Y: if Y is a matrix, the number of columns is the space of prediction
    if (is.null(dim(Y))) {
      d = 1
      T = length(Y)
    } else {
      d = ncol(Y)
      T = nrow(Y)
      if (d > 1 && T > 1 && length(dim(experts)) < 3) {
        stop("Bad dimensions: nrow(experts) should be equal to dim(experts)[3]")
      } 
      if (length(dim(experts)) == 3) {
        if ((dim(experts)[1] != T) || (dim(experts)[2] != d)){
          stop("Bad dimensions between Y and experts")
        }
      }
      if (T == 1 && d>1) {
        if (length(dim(experts)) == 2) {
          if (dim(experts)[1] != d) {
            stop("Bad dimensions between Y and experts")
          }
        }
      }
    }
    
    if (T == 1 && d == 1) {
      experts <- as.matrix(experts)
      if (nrow(experts) == 1 || ncol(experts) == 1) {
        experts <- matrix(as.numeric(experts), nrow = 1)
      } else {
        stop("Bad dimensions: length(Y) should be equal to nrow(experts)")
      }
    }
    
    if (dim(experts)[1] != T) {
      stop("Bad dimensions: length(Y) should be equal to nrow(experts)")
    }
    object$d <- d
    object <- predict(object, newY = Y, newexperts = experts, awake = awake, 
      type = "model")
    
  }
  return(object)
} 
