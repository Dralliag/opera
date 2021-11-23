#' Implementation of FTRL (Follow The Regularized Leader)
#'
#'  FTRL \insertCite{shalev2007primal}{opera} and Chap. 5 of \insertCite{hazan2019introduction}{opera} is the online counterpart of empirical risk minimization. 
#'  It is a family of aggregation rules (including OGD) that uses at any time the empirical risk
#'  minimizer so far with an additional regularization. The online optimization can be performed
#'  on any bounded convex set that can be expressed with equality or inequality constraints. 
#'  Note that this method is still under development and a beta version. 
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
#'      \item{list}{ List with field \code{name} equal to the loss name. If using pinball loss, field \code{tau} equal to the required quantile in [0,1];}
#'      \item{function}{ A custom loss as a function of two parameters (prediction, label).}
#' }
#' @param loss.gradient \code{boolean, function} (TRUE). 
#' \describe{
#'      \item{boolean}{ If TRUE, the aggregation rule will not be directly applied to the loss function at hand,
#'      but to a gradient version of it. The aggregation rule is then similar to gradient descent aggregation rule. }
#'      \item{function}{ If loss.type is a function, the derivative of the loss in its first component should be provided to be used (it is not automatically 
#'      computed).}
#' }
#' @param w0 \code{numeric} (NULL). Vector of initialization for the weights.
#' @param max_iter \code{integer} (50). Maximum number of iterations of the optimization algorithm per round.
#' @param obj_tol \code{numeric} (1e-2). Tolerance over objective function between two iterations of the optimization.
#' @param training \code{list} (NULL). List of previous parameters.
#' @param default \code{boolean} (FALSE). Whether or not to use default parameters for fun_reg, constr_eq, constr_ineq and their grad/jac, 
#' which values are ALL ignored when TRUE.
#' @param quiet \code{boolean} (FALSE). Whether or not to display progress bars.
#'
#' @return object of class mixture.
#'
#' @references
#'   \insertAllCited{}
#'
#' @import alabama
#' @export FTRL
#'  
FTRL <- function(y, 
                 experts, 
                 eta = NULL, 
                 fun_reg = NULL, fun_reg_grad = NULL,
                 constr_eq = NULL, constr_eq_jac = NULL,
                 constr_ineq = NULL, constr_ineq_jac = NULL,
                 loss.type = list(name = "square"),
                 loss.gradient = TRUE, 
                 w0 = NULL,
                 max_iter = 50,
                 obj_tol = 1e-2,
                 training = NULL, 
                 default = FALSE,
                 quiet = TRUE) {
  
  # retrieve training values
  if (! is.null(training)) {
    eta <- training$eta
    default_eta <- training$default_eta
    fun_reg <- training$fun_reg ; fun_reg_grad <- training$fun_reg_grad ;
    constr_eq <- training$constr_eq ; constr_eq_jac <- training$constr_eq_jac ; 
    constr_ineq <- training$constr_ineq ; constr_ineq_jac <- training$constr_ineq_jac ; 
    loss.type <- training$loss.type
    loss.gradient <- training$loss.gradient
    w0 <- training$w0
    G <- training$G
    max_iter <- training$max_iter
    obj_tol <- training$obj_tol
  }
  
  # checks
  if (! is.null(loss.gradient) && ! is.function(loss.gradient) && loss.gradient == FALSE) {
    stop("loss.gradient must be provided to use the FTRL algorithm.")
  }
  if (is.null(eta)) {
    default_eta <- TRUE
    eta = Inf
  }
  else if (! exists("default_eta")) {
    default_eta <- FALSE
  }
  if (default == FALSE && (is.null(fun_reg) || ! is.function(fun_reg))) {
    stop("fun_reg cannot be missing when other optimization parameters are provided (see ?auglag... fn).")
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
  
  
  N <- ncol(experts)  # Number of experts
  T <- nrow(experts)  # Number of instants
  
  if (default) {
    if (is.null(w0)){
      w0 = rep(1/N, N)
    }
    
    fun_reg <- function(x) {sum(x * log(x/w0))}
    fun_reg_grad <- function(x) {log(x/w0) + 1}
    constr_eq <- function(x) {sum(x) - 1}
    constr_eq_jac <- function(x) matrix(1, ncol = N)
    constr_ineq <- function(x) x
    constr_ineq_jac <- function(x) diag(N)
  }
  
  parms = list()
  
  # inits
  if (is.null(training) && (is.function(loss.gradient) || loss.gradient)) {
    G <- numeric(N) 
  }
  weights <- matrix(0, nrow = T, ncol = N)
  prediction <- rep(0, T)
  if (is.null(w0) && is.null(training)) {
    parms = list(
      "par" = rep(1/N, N),
      "fn" = fun_reg, 
      "heq" = constr_eq, 
      "hin" = constr_ineq,
      "control.outer" = list(trace=FALSE, kkt2.check = FALSE, itmax = max_iter, eps = obj_tol)
    )
    parms <- parms[! sapply(parms, is.null)]
    
    # if no constraints: use optim instead of alabama
    if (is.null(parms$heq) && is.null(parms$hin)) {
      parms <- c(parms[intersect(c("par", "fn", "gr"), names(parms))], "control" = list(list("trace" = 0)))
      result <- do.call(stats::optim, parms)
    }
    else {
      result <- do.call(alabama::auglag, parms) 
    }
    weights[1, ] <- result$par
  } else {
    weights[1, ] <- if (is.null(training)) {w0} else {training$last_weights}
  }
  
  if (! quiet) steps <- init_progress(T)
  for (t in 1:T) {
    if (! quiet) update_progress(t, steps)
    
    if (is.function(loss.gradient) || loss.gradient) {
      # compute current prevision
      prediction[t] <- weights[t,] %*% experts[t,]
      
      # update gradient
      G_t <- loss(experts[t, ], y[t], prediction[t], loss.type = loss.type, loss.gradient = loss.gradient)
      G <- G + G_t
      if (default_eta) {
        eta <- 1/sqrt(1/eta^2 + sum(G_t^2))  
      }
      
      # update obj function
      obj <- function(x) (fun_reg(x) + eta * sum(G * x))
      
      # compute grad of obj function
      obj_grad <- if(is.null(fun_reg_grad)) {NULL} else {function(x) fun_reg_grad(x) + eta * G}
      
      # run optimization
      parms <- list("par" = weights[t, ],
                    "fn" = obj,
                    "gr" = obj_grad,
                    "heq" = constr_eq,
                    "heq.jac" = constr_eq_jac,
                    "hin" = constr_ineq,
                    "hin.jac" = constr_ineq_jac,
                    "control.outer" = list(trace = FALSE, itmax = max_iter, eps = obj_tol, kkt2.check = FALSE))
      
      parms <- parms[! sapply(parms, is.null)]
      # 
      
      if (is.null(parms$heq) && is.null(parms$hin)) {
        parms <- c(parms[intersect(c("par", "fn", "gr"), names(parms))], "control" = list(list("trace" = 0)))
        res_optim <- do.call(stats::optim, parms)
      }
      else {
        res_optim <- suppressWarnings(do.call(alabama::auglag, parms))
      }
      
      # update weights
      if (t < T) {
        weights[t + 1, ] <- res_optim$par
      } else {
        coeffs <- res_optim$par
      }
      # check convergence
      if (! res_optim$convergence == 0){
        warning(paste0("Optimization didn't converge at step ", t, ". Your decision space might not be compact, 
                       or your regularisation function not convex."))
      }
    }
  }
  if (! quiet) end_progress()
  
  
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
  
  object$training <- list("eta" = eta, 
                          "default_eta" = default_eta, 
                          "fun_reg" = fun_reg, "fun_reg_grad" = fun_reg_grad,
                          "constr_eq" = constr_eq, "constr_eq_jac" = constr_eq_jac,
                          "constr_ineq" = constr_ineq, "constr_ineq_jac" = constr_ineq_jac,
                          "loss.type" = loss.type,
                          "loss.gradient" = loss.gradient,
                          "w0" = w0,
                          "last_weights" = res_optim$par,
                          "G" = G,
                          "max_iter" = max_iter,
                          "obj_tol" = obj_tol)
  
  return(object)
}
