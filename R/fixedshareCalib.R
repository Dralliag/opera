fixedshareCalib <- function(y, experts, grid.eta = 1, grid.alpha = 10^(-4:-1), awake = NULL, 
  loss.type = "square", loss.gradient = TRUE, w0 = NULL, trace = F, gamma = 2, 
  training = NULL) {
  experts <- as.matrix(experts)
  
  N <- ncol(experts)  # Number of experts
  T <- nrow(experts)  # Number of instants
  
  # Uniform initial weight vector if unspecified
  if (is.null(w0)) {
    w0 <- rep(1, N)
  }
  
  if (is.null(gamma)) {
    gamma <- 2
  }
  
  awake <- as.matrix(awake)
  idx.na <- which(is.na(experts))
  awake[idx.na] <- 0
  experts[idx.na] <- 0
  
  T0 <- 0  # number of observations in previous runs
  neta <- length(grid.eta)  # Number of learning rates in the grid
  nalpha <- length(grid.alpha)  # Number of mixing rates in the grid
  bestpar <- c(floor(neta)/2, floor(nalpha)/2) + 1  # We start with the parameters in the middle of the grids
  par <- NULL
  
  cumulativeLoss <- array(0, dim = c(neta, nalpha))
  
  
  wpar <- array(w0, dim = c(N, neta, nalpha))  # Weight array (in 3 dimensions) assigned by each algorithm FixedShare(eta,alpha) 
  weights <- matrix(0, ncol = N, nrow = T)  # Matrix of weights formed by the mixture
  prediction <- rep(0, T)  # Predictions formed by the mixture
  
  if (!is.null(training)) {
    bestpar <- training$bestpar
    wpar <- training$wpar
    w0 <- training$w0
    R <- training$R
    cumulativeLoss <- training$cumulativeLoss
    T0 <- training$T
  }
  
  for (t in 1:T) {
    # Display the state of progress of the algorithm
    if (!(t%%floor(T/10)) && trace) 
      cat(floor(10 * t/T) * 10, "% -- ")
    
    # Weights, prediction forme by FixedShare(eta[t],alpha[t]) where par[t,] =
    # c(eta[t],alpha[t]) are the parameters calibrated online
    weights[t, ] <- wpar[, bestpar[1], bestpar[2]] * awake[t, ]/sum(wpar[, bestpar[1], 
      bestpar[2]] * awake[t, ])
    prediction[t] <- experts[t, ] %*% weights[t, ]
    par <- rbind(par, data.frame(eta = grid.eta[bestpar[1]], alpha = grid.alpha[bestpar[2]]))
    
    # Loop over the mixing rates alpha in the grid 'grid.alpha'
    for (k in 1:nalpha) {
      
      # Weights, prediction, and losses formed by FixedShare(eta,alpha) for (eta,alpha)
      # in the grid
      waux <- t(t(wpar[, , k] * awake[t, ])/apply(as.matrix(wpar[, , k] * awake[t, 
        ]), 2, sum))
      pred <- experts[t, ] %*% waux
      
      cumulativeLoss[, k] <- cumulativeLoss[, k] + loss(pred, y[t], loss.type)  # Non gradient cumulative losses
      lpred <- diag(lossPred(pred, y[t], pred, loss.type, loss.gradient))
      lexp <- lossPred(experts[t, ], y[t], pred, loss.type, loss.gradient)
      
      # Regret update
      R <- t(t(log(wpar[, , k]))/grid.eta) + awake[t, ] * t(lpred - t(lexp))
      
      # Weight update
      v <- truncate1(exp(t(t(matrix(R, ncol = neta)) * grid.eta)))
      v <- t(t(v)/apply(v, 2, sum))  # Renormalization of each column
      wpar[, , k] <- grid.alpha[k]/N + (1 - grid.alpha[k]) * v
    }
    
    # Grid update *************** find the best parameter in the grid
    bestpar <- which(cumulativeLoss == min(cumulativeLoss), arr.ind = TRUE)[1, 
      ]
    
    # Expand the grid if the best parameter lies on an extremity (only for the first
    # component)
    if (bestpar[1] == neta) {
      if (trace) 
        cat(" + ")
      neweta <- grid.eta[neta] * gamma^(1:3)
      grid.eta <- c(grid.eta, neweta)
      neta <- neta + length(neweta)
      wparaux <- wpar
      wpar <- array(dim = c(N, neta, nalpha))
      wpar[, 1:bestpar[1], ] <- wparaux
      for (j in 1:length(neweta)) {
        cumulativeLoss <- rbind(cumulativeLoss, 0)
        for (k in 1:nalpha) {
          perfnewpar <- fixedshare(c(training$oldY, y[1:t]), rbind(training$oldexperts, 
          matrix(experts[1:t, ], ncol = N)), neweta[j], grid.alpha[k], 
          awake = rbind(training$oldawake, matrix(awake[1:t, ], ncol = N)), 
          loss.type = loss.type, loss.gradient = loss.gradient, w0 = w0)
          cumulativeLoss[bestpar[1] + j, k] <- perfnewpar$training$cumulativeLoss
          wpar[, bestpar[1] + j, k] <- perfnewpar$coefficients
        }
      }
    }
    if (bestpar[1] == 1) {
      if (trace) 
        cat(" - ")
      neweta <- grid.eta[1]/gamma^(1:3)
      neta <- neta + length(neweta)
      bestpar[1] <- bestpar[1] + length(neweta)
      wparaux <- wpar
      wpar <- array(dim = c(N, neta, nalpha))
      wpar[, bestpar[1]:neta, ] <- wparaux
      for (j in 1:length(neweta)) {
        grid.eta <- c(neweta[j], grid.eta)
        cumulativeLoss <- rbind(0, cumulativeLoss)
        for (k in 1:nalpha) {
          perfnewpar <- fixedshare(c(training$oldY, y[1:t]), rbind(training$oldexperts, 
          matrix(experts[1:t, ], ncol = N)), neweta[j], grid.alpha[k], 
          awake = rbind(training$oldawake, matrix(awake[1:t, ], ncol = N)), 
          loss.type = loss.type, loss.gradient = loss.gradient, w0 = w0)
          cumulativeLoss[1, k] <- perfnewpar$training$cumulativeLoss
          wpar[, bestpar[1] - j, k] <- perfnewpar$coefficients
        }
      }
    }
  }
  
  # Next weights
  w <- wpar[, bestpar[1], bestpar[2]]/sum(wpar[, bestpar[1], bestpar[2]])
  
  object <- list(model = "FS", loss.type = loss.type, loss.gradient = loss.gradient, 
    coefficients = w)
  
  
  object$parameters <- list(eta = par$eta[1:T], alpha = par$alpha[1:T], grid.eta = grid.eta, 
    grid.alpha = grid.alpha)
  object$weights <- weights
  object$prediction <- prediction
  
  object$training <- list(T = T0 + T, wpar = wpar, w0 = w0, R = R, cumulativeLoss = cumulativeLoss, 
    bestpar = bestpar, grid.loss = cumulativeLoss/(T0 + T), oldexperts = rbind(training$oldexperts, 
      experts), oldY = c(training$oldY, y), oldawake = rbind(training$oldawake, 
      awake))
  
  
  rownames(object$training$grid.loss) <- grid.eta
  colnames(object$training$grid.loss) <- grid.alpha
  class(object) <- "mixture"
  
  if (trace) 
    cat("\n")
  return(object)
} 
