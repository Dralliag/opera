MLprod <- function(y, experts, awake = NULL, loss.type = "square", loss.gradient = TRUE, 
                   w0 = NULL, training = NULL, use_cpp = getOption("opera_use_cpp", default = TRUE)) {
  
  experts <- as.matrix(experts)
  N <- ncol(experts)
  T <- nrow(experts)
  
  # Uniform initial weight vector if unspecified
  if (is.null(w0)) {
    w0 <- rep(1, N)
  }
  
  awake <- as.matrix(awake)
  idx.na <- which(is.na(experts))
  awake[idx.na] <- 0
  experts[idx.na] <- 0
  
  weights <- matrix(0, ncol = N, nrow = T)
  R <- log(w0)
  
  L <- rep(1, N)
  maxloss <- 0
  
  eta <- matrix(exp(700), ncol = N, nrow = T + 1)
  prediction <- rep(0, T)
  
  if (!is.null(training)) {
    eta[1, ] <- training$eta
    R <- training$R
    L <- training$L
    maxloss <- training$maxloss
  }
  
  if (use_cpp){
    loss_tau <- ifelse(! is.null(loss.type$tau), loss.type$tau, 0)
    loss_name <- loss.type$name
    B <- computeMLProdEigen(awake,eta,experts,weights,y,prediction,
                            R,L,maxloss,loss_name,loss_tau,loss.gradient);
  }
  else{
    steps <- init_progress(T)
    
    for (t in 1:T) {
      update_progress(t, steps)
      
      # Update weights
      idx <- awake[t,] > 0
      R.max <- max(R[idx])
      w <- numeric(N)
      w[idx] <- exp(R[idx]-R.max)
      w[idx] <- eta[t, idx] * w[idx]/sum(eta[t, idx] * w[idx])
      
      p <- awake[t, ] * w/sum(awake[t, ] * w)
      pred <- experts[t, ] %*% p  # Predict
      
      weights[t, ] <- p
      prediction[t] <- pred
      
      # Observe losses
      lpred <- lossPred(pred, y[t], pred, loss.type = loss.type, loss.gradient = loss.gradient)
      lexp <- lossPred(experts[t, ], y[t], pred, loss.type = loss.type, loss.gradient = loss.gradient)
      
      r <- awake[t, ] * (c(c(lpred) - lexp))
      L <- L + r^2
      maxloss <- pmax(maxloss, abs(r))
      neweta <- pmin(1/(2 * maxloss), sqrt(log(N)/L))
      
      # Update regret and learning parameter
      R <- neweta/eta[t, ] * R + log(1 + awake[t, ] * neweta * (c(c(lpred) - lexp)))
      eta[t + 1, ] <- neweta
      
      if (is.na(sum(R))) {
        browser("Nan in R")
      }
    }
    end_progress()
  }
  
  R.max <- max(R)
  w <- exp(R-R.max)
  w <- eta[T + 1, ] * w/sum(eta[T + 1, ] * w)
  
  object <- list(model = "MLprod", loss.type = loss.type, loss.gradient = loss.gradient, 
                 coefficients = w)
  
  object$parameters <- list(eta = eta[1:T, ])
  object$weights <- weights
  object$prediction <- prediction
  
  object$training <- list(eta = eta[T + 1, ], R = R, L = L, maxloss = maxloss)
  class(object) <- "mixture"
  
  return(object)
} 
