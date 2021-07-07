MLpol <- function(y, experts, awake = NULL, loss.type = "square", loss.gradient = TRUE, 
  training = NULL, use_cpp = getOption("opera_use_cpp", default = TRUE)) {
  
  experts <- as.matrix(experts)
  N <- ncol(experts)
  T <- nrow(experts)
  
  awake <- as.matrix(awake)
  idx.na <- which(is.na(experts))
  awake[idx.na] <- 0
  experts[idx.na] <- 0
  
  # weight assigned to each expert
  weights <- matrix(0, ncol = N, nrow = T)
  prediction <- rep(0, T)
  
  # Initialization or the learning parameter
  B <- 0
  eta <- matrix(exp(700), ncol = N, nrow = T + 1)
  # regret suffered by each expert
  R <- rep(0, N)
  
  if (!is.null(training)) {
    eta[1, ] <- training$eta
    R <- training$R
    B <- training$B
  } else {
    training <- list(eta = eta[1, ])
  }
  
  w <- rep(0, N)
  
  if (use_cpp){
    loss_tau <- ifelse(! is.null(loss.type$tau), loss.type$tau, 0)
    loss_name <- loss.type$name
    B <- computeMLPolEigen(awake,eta,experts,weights,y,prediction,
                           R,w,B,loss_name,loss_tau,loss.gradient)
  }
  else{
    steps <- init_progress(T)
    
    for (t in 1:T) {
      update_progress(t, steps)
      
      # We check if there is at least one expert with positive weight
      if (max(awake[t, ] * R) > 0) {
        w <- eta[t, ] * pmax(R, 0)/sum(eta[t, ] * pmax(R, 0))
      } else {
        w <- rep(1, N)
      }
      
      # form the mixture and the prediction
      p <- awake[t, ] * w/sum(awake[t, ] * w)
      pred <- experts[t, ] %*% p
      
      # save the mixture and the prediction
      weights[t, ] <- p
      prediction[t] <- pred
      
      # Observe losses
      lpred <- lossPred(pred, y[t], pred, loss.type = loss.type, loss.gradient = loss.gradient)
      lexp <- lossPred(experts[t, ], y[t], pred, loss.type = loss.type, loss.gradient = loss.gradient)
      
      # Update the regret and the weight
      r <- awake[t, ] * c(c(lpred) - lexp)
      R <- R + r
      
      # Update the learning rate
      newB <- max(B, max(r^2))
      eta[t + 1, ] <- 1/(1/eta[t, ] + r^2 + newB - B)
      B <- newB
    }
    end_progress()
    
    # We check if there is at least one expert with positive weight
    if (max(R) > 0) {
      w <- eta[T + 1, ] * pmax(R, 0)/sum(eta[T + 1, ] * pmax(R, 0))
    } else {
      w <- rep(1/N, N)
    }
  }
  object <- list(model = "MLpol", loss.type = loss.type, loss.gradient = loss.gradient, 
                 coefficients = w)
  
  object$parameters <- list(eta = eta[1:T, ])
  object$weights <- weights
  object$prediction <- prediction
  
  object$training <- list(eta = eta[T + 1, ], R = R, B = B)
  class(object) <- "mixture"
  return(object)
} 
