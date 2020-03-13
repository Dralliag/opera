
fixedshare <- function(y, experts, eta, alpha, awake = NULL, loss.type = "square", 
  loss.gradient = TRUE, w0 = NULL, training = NULL) {
  experts <- as.matrix(experts)
  
  N <- ncol(experts)  # Number of experts
  T <- nrow(experts)  # Number of instants
  
  # Uniform initial weight vector if unspecified
  if (is.null(w0)) {
    w0 <- rep(1, N)
  }
  
  awake <- as.matrix(awake)
  idx.na <- which(is.na(experts))
  awake[idx.na] <- 0
  experts[idx.na] <- 0
  
  R <- log(w0)/eta
  pred <- rep(0, T)  # Prediction vector
  cumulativeLoss <- 0  # Cumulative loss of the mixture
  weights <- matrix(0, ncol = N, nrow = T)
  
  if (!is.null(training)) {
    R <- training$R
    cumulativeLoss <- training$cumulativeLoss
  }
  
  for (t in 1:T) {
    # Weight update
    weights[t, ] <- t(truncate1(exp(eta * R)) * t(awake[t, ]))
    weights[t, ] <- weights[t, ]/sum(weights[t, ])
    
    # Prediction and loss
    pred[t] <- experts[t, ] %*% weights[t, ]
    cumulativeLoss <- cumulativeLoss + loss(pred[t], y[t], loss.type)
    lpred <- lossPred(pred[t], y[t], pred[t], loss.type, loss.gradient)
    lexp <- lossPred(experts[t, ], y[t], pred[t], loss.type, loss.gradient)
    
    # Regret and weight update
    R <- R + awake[t, ] * (c(c(lpred) - lexp))
    v <- truncate1(exp(eta * R))/sum(truncate1(exp(eta * R)))
    R <- log(alpha/N + (1 - alpha) * v)/eta
  }
  w <- t(truncate1(exp(eta * R)))
  w <- w/sum(w)
  
  object <- list(model = "FS", loss.type = loss.type, loss.gradient = loss.gradient, 
    coefficients = w)
  
  object$parameters <- list(eta = eta, alpha = alpha)
  object$weights <- weights
  object$prediction <- pred
  
  object$training <- list(R = R, cumulativeLoss = cumulativeLoss)
  
  return(object)
} 
