#' Implementation of RFTL (Regularized Follow The Leader)
#'
#' @param y \code{vector}. Real observations.
#' @param experts \code{matrix}. Matrix of experts previsions.
#' @param eta \code{numeric}. Regularization parameter.
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
#' @param parameters \code{list of 3 items}. Parameters for the optimization part of the algorithm (using CVXR package).
#' \describe{
#'      \item{p}{ Object of class variable, as defined in package CVXR.}
#'      \item{reg}{ Regularization function to be applied on p.}
#'      \item{constr (NULL)}{ List of constraints to be applied on p, stored in char expr.}
#' }
#'
#' @return object of class mixture.
#'
#' @import alabama
#'
#' @examples
RFTL <- function(y, 
                 experts, 
                 eta, 
                 reg = function(x) sqrt(sum(x**2)), 
                 reg_grad = function(x) x * (x**2)**(-1/2),
                 heq, heq_jac = NULL,
                 hin, hin_jac = NULL,
                 loss.type = "square",
                 loss.gradient = TRUE, 
                 w0,
                 itmax = 50) {
  
  # checks
  if (is.null(eta)) {
    stop("eta must be provided as a numeric.")
  }
  if (is.null(heq) || ! is.function(heq)) {
    stop("heq must be provided as a function (see ?auglag).")
  }
  if (is.null(hin) || ! is.function(hin)) {
    stop("hin must be provided as a function (see ?auglag).")
  }
  if (! is.null(reg_grad) && ! is.function(reg_grad)) {
    stop("reg_grad must be a function (the gradient of the regulation function).") 
  }
  if (! is.null(heq_jac) && ! is.function(heq_jac)) {
   stop("heq_jac must be a function that returns a matrix (see ?auglag).") 
  }
  if (! is.null(hin_jac) && ! is.function(hin_jac)) {
    stop("hin_jac must be a function that returns a matrix (see ?auglag).") 
  }
  
  
  N <- ncol(experts)  # Number of experts
  T <- nrow(experts)  # Number of instants
  
  # init
  if (! loss.gradient == FALSE) {
    G <- numeric(N) 
  }
  weights <- matrix(0, nrow = T, ncol = N)
  prediction <- rep(0, T)
  if (identical(as.character(attributes(reg)$srcref), "function(x) sqrt(sum(x**2))")) {
    reg_grad <- function(x) x / sqrt(sum(x^2))
  }
  if (identical(as.character(attributes(heq)$srcref), "function(x) sum(x) - 1")) {
    heq_jac <- function(x) {matrix(1, ncol = N)}
  } else {
    heq_jac <- function(par, ...) jacobian(func = heq, x = par, method = "simple", ...)
  }
  if (identical(as.character(attributes(hin)$srcref), "function(x) x")) {
    hin_jac <- function(x) {diag(N)}
  } else {
    hin_jac <- function(par, ...) jacobian(func = hin, x = par, method = "simple", ...)
  }
  
  p0 <- rep(1/N, N)
  if (is.null(w0)) {
    result <- alabama::auglag(par = p0, fn = reg, heq = heq, hin = hin, control.outer = list(trace = F))
    weights[1, ] <- result$par
  } else {
    weights[1, ] <- w0
  }
  
  obj <- reg
  
  steps <- init_progress(T)
  for (t in 1:T) {
    update_progress(t, steps)
    
    if (! loss.gradient == FALSE) {
      
      pred <- prediction[t] <- weights[t,] %*% X[t,]
      
      # update gradient
      G <- G + lossPred(experts[t, ], Y[t], pred, loss.type = loss.type, loss.gradient = loss.gradient)
      
      # update obj function
      obj <- function(x) (reg(x) + eta * sum(G * x))
      
      # compute grad of obj function
      obj_grad <- function(x) {reg_grad(x) + eta * G}
      
      res_optim <- alabama::auglag(par = if (t == 1) {w0} else {weights[t-1, ]}, 
                      fn = obj, 
                      gr = obj_grad, 
                      heq = heq, 
                      heq.jac = heq_jac,
                      hin = hin, 
                      hin.jac = hin_jac,
                      control.outer = list(trace = F, itmax = itmax, kkt2.check = FALSE))
      
      # Q : check convergence / contraintes / .. et stop sinon
      if(res_optim$convergence == 0){
        weights[t + 1, ] <- res_optim$par
      } else {
        stop("No convergence...!")
      }

    } 
    # pas d'intéret car dès lors qu'on a le gradient, on peut faire la méthode la plus rapide
    # else {
    #   obj <- function(x) {
    #     res <- sapply(1:t, function(it) {
    #       eta * lossPred(sum(experts[it, ] * x), Y[it], loss.type = loss.type, loss.gradient = FALSE)
    #     })
    # 
    #     reg(x) + sum(res)
    #   }
    # 
    #   if (is.list(loss.type)) {
    #     loss_grad <- get_custom_loss_grad(loss.type)
    #     obj_grad <- function(x) {
    #       res <- sapply(1:t, function(it) {
    #         eta * loss_grad(sum(experts[it, ] * x), Y[it])
    #       })
    # 
    #       reg_grad(x) + sum(res)
    #     }
    #   }
    #   else {
    #     obj_grad <- function(par, ...) grad(func = obj, x = par, method = "simple", ...)
    #   }
    # 
    #   weights[t,] <- alabama::auglag(par = if (t == 1) {w0} else {weights[t-1, ]},
    #                                  fn = obj, gr = obj_grad,
    #                                  heq = heq, heq.jac = heq_jac,
    #                                  hin = hin, hin.jac = hin_jac,
    #                                  control.outer = list(trace = F, itmax = itmax))$par
    # }
  }
  end_progress()
  
  object <- list(model = "RFTL", loss.type = loss.type, loss.gradient = loss.gradient, 
                 coefficients = weights[T, ])
  class(object) <- "mixture"
  
  object$parameters <- list("eta" = eta, 
                            "reg" = reg,
                            "constr" = list("heq" = heq, "hin" = hin))
  object$weights <- weights
  object$prediction <- prediction
  
  return(object)
}



# get_custom_loss_grad <- function(loss.type) {
#   loss_grad <- NULL
#   
#   if (loss.type$name == "square") {
#    loss_grad <- function(x, y) 2*(x-y)
#   }
#   else if (loss.type$name == "percentage") {
#     loss_grad <- function(x, y) sign(x-y)/y
#   }
#   else if (loss.type$name == "absolute") {
#     loss_grad <- function(x, y) sign(x-y)
#   } 
#   else if (loss.type$name == "pinball") {
#     loss_grad <- function(x, y) -loss.type$tau
#   }
#   
#   return(loss_grad)
# }