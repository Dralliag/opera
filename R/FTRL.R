#' Implementation of FTRL (Follow The Regulaized Leader)
#'
#' @param y \code{vector}. Real observations.
#' @param experts \code{matrix}. Matrix of experts previsions.
#' @param eta \code{numeric}. Regularization parameter.
#' @param fun_reg \code{function} (NULL). Regularization function to be applied during the optimization. 
#' @param fun_reg_grad \code{function} (NULL). Gradient of the regularization function (to speed up the computations).
#' @param constr_eq \code{function} (NULL). Constraints (equalities) to be applied during the optimization.
#' @param constr_eq_jac \code{function} (NULL). Jacobian of the equality constraints (to speed up the computations). 
#' @param constr_ineq \code{function} (NULL). Constraints (inequalities) to be applied during the optimization (... > 0).
#' @param constr_ineq_jac \code{function} (NULL). Jacobian of the inequality constraints (to speed up the computations).
#' @param loss.type \code{character, list or function} ("square").  
#' \describe{
#'      \item{character}{ Name of the loss to be applied ('square', 'absolute', 'percentage', or 'pinball');}
#'      \item{list}{ When using pinball loss: list with field name equal to 'pinball' and field tau equal to the required quantile in [0,1];}
#'      \item{function}{ A custom loss as a function of two parameters.}
#' }
#' @param loss.gradient \code{boolean, function} (TRUE). 
#' \describe{
#'      \item{boolean}{ If TRUE, the aggregation rule will not be directly applied to the loss function at hand,
#'      but to a gradient version of it. The aggregation rule is then similar to gradient descent aggregation rule. }
#'      \item{function}{ If loss.type is a function, the derivative should be provided to be used (it is not automatically 
#'      computed).}
#' }
#' @param w0 \code{numeric} (NULL). Vector of initialization for the weights.
#' @param max_iter \code{integer} (50). Maximum number of iterations of the optimization algorithm per round.
#' @param default \code{boolean} (FALSE). Whether or not to use default parameters for fun_reg, constr_eq, constr_ineq and their grad/jac, 
#' which values are ALL ignored when TRUE.
#'
#' @return object of class mixture.
#'
#' @import alabama
#'
#' @examples
#' 
#' 
#' 
#' 
FTRL <- function(y, 
                 experts, 
                 eta, 
                 fun_reg = NULL, fun_reg_grad = NULL,
                 constr_eq = NULL, constr_eq_jac = NULL,
                 constr_ineq = NULL, constr_ineq_jac = NULL,
                 loss.type = "square",
                 loss.gradient = TRUE, 
                 w0 = NULL,
                 max_iter = 50,
                 obj_tol = 1e-2,
                 default = FALSE) {
  
  # checks
  if (is.null(eta)) {
    # stop("eta must be provided as a numeric in argument 'parameters'.")
    eta = 0.1 # to be initialized and dynamically updated when NULL
  }
  if (default == FALSE && (is.null(fun_reg) || ! is.function(fun_reg))) {
    stop("fun_reg must cannot be missing when other optimization parameters are provided (see ?auglag... fn).")
  }
  if (! is.null(fun_reg_grad) && ! is.function(fun_reg_grad)) {
    stop("fun_reg_grad must be a function (the gradient of the fun_reg function).")
  }
  if (! is.null(constr_eq) && ! is.function(constr_eq)) {
    stop("constr_eq must be provided as a function (see ?auglag... heq).")
  }
  if (! is.null(constr_ineq) && ! is.function(constr_ineq)) {
    stop("constr_ineq must be provided as a function (see ?auglag... hin).")
  }
  if (! is.null(constr_eq_jac) && ! is.function(constr_eq_jac)) {
   stop("constr_eq_jac must be a function that returns a matrix (see ?auglag... heq.jac).")
  }
  if (! is.null(constr_ineq_jac) && ! is.function(constr_ineq_jac)) {
    stop("constr_ineq_jac must be a function that returns a matrix (see ?auglag... hin.jac).")
  }
  if (! is.null(constr_eq_jac) && is.null(constr_eq)) {
    stop("constr_eq_jac is not null but contr_eq is missing.")
  }
  if (! is.null(constr_ineq_jac) && is.null(constr_ineq)) {
    stop("constr_ineq_jac is not null but contr_ineq is missing.")
  }
  if (is.null(max_iter)) {
    max_iter <- 50
  }
  if (is.null(obj_tol)) {
    obj_tol <- 1e-2 
  }
  
  # 
  N <- ncol(experts)  # Number of experts
  T <- nrow(experts)  # Number of instants
  
  if (default) {
    fun_reg <- function(x) sqrt(sum(x**2))
    fun_reg_grad <- function(x) x / sqrt(sum(x**2))
    constr_eq <- function(x) sum(x) - 1
    constr_eq_jac <- function(x) matrix(1, ncol = N)
    constr_ineq <- function(x) x
    constr_ineq_jac <- function(x) diag(N)
  }
  
  # inits
  if (! loss.gradient == FALSE) {
    G <- numeric(N) 
  }
  weights <- matrix(0, nrow = T, ncol = N)
  prediction <- rep(0, T)
  if (is.null(w0)) {
    result <- alabama::auglag(par = rep(1/N, N), 
                              fn = fun_reg, 
                              heq = constr_eq, 
                              hin = constr_ineq, 
                              control.outer = list(trace = FALSE, kkt2.check = FALSE, itmax = max_iter, eps = obj_tol))
    weights[1, ] <- result$par
  } else {
    weights[1, ] <- w0
  }
  
  steps <- init_progress(T)
  for (t in 1:T) {
    update_progress(t, steps)
    
    if (! loss.gradient == FALSE) {
      # compute current prevision
      prediction[t] <- weights[t,] %*% X[t,]
      
      # update gradient
      G <- G + lossPred(experts[t, ], Y[t], prediction[t], loss.type = loss.type, loss.gradient = loss.gradient)
      
      # update obj function
      obj <- function(x) (fun_reg(x) + eta * sum(G * x))
      
      # compute grad of obj function
      obj_grad <- if(is.null(fun_reg_grad)) {NULL} else {function(x) fun_reg_grad(x) + eta * G}
      
      # run optimization
      parms <- list("par" = if (t == 1) {w0} else {weights[t-1, ]},
                    "fn" = obj,
                    "gr" = obj_grad,
                    "heq" = constr_eq,
                    "heq.jac" = constr_eq_jac,
                    "hin" = constr_ineq,
                    "hin.jac" = constr_ineq_jac,
                    "control.outer" = list(trace = FALSE, itmax = max_iter, kkt2.check = FALSE))
      
      parms <- parms[! sapply(parms, is.null)]
      # 
      if (is.null(parms$heq) && is.null(parms$hin)) {
        parms <- c(parms[intersect(c("par", "fn", "gr"), names(parms))], "control" = list(list("trace" = 0, maxit = max_iter, abstol = obj_tol)))
        res_optim <- do.call(stats::optim, parms)
      }
      else {
        res_optim <- do.call(alabama::auglag, parms) 
      }
      
      # update weights
      if (t < T) {
        weights[t + 1, ] <- res_optim$par
      } else {
        coeffs <- res_optim$par
      }
      # check convergency
      if (! res_optim$convergence == 0){
        warning(paste0("Optimization didn't converge at step ", t, ". Your decision space might not be compact, 
                       or your regularisation function not convex."))
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
  
  # round to 0 when coefficients are slightly negative
  if (all(weights > -1e-5)) {
    coeffs <- pmax(0, coeffs)
    coeffs <- coeffs / sum(coeffs) 
    
    weights <- apply(weights, 2, pmax, 0)
    weights <- weights / rowSums(weights)
  }
  
  object <- list(model = "FTRL", loss.type = loss.type, loss.gradient = loss.gradient, 
                 coefficients = coeffs)
  class(object) <- "mixture"
  
  object$parameters <- list("eta" = eta, 
                            "w0" = w0,
                            "reg" = list("fun_reg" = fun_reg, "fun_reg_grad" = fun_reg_grad),
                            "constr" = list("eq" = constr_eq, "ineq" = constr_ineq, "eq_jac" = constr_eq_jac, "ineq_jac" = constr_ineq_jac),
                            "max_iter" = max_iter,
                            "default" = default)
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