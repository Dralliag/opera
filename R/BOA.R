
BOA <- function(y, experts, awake = NULL, loss.type = "square", loss.gradient = TRUE, 
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
  
  R <- rep(0, N)
  R.reg <- rep(0, N)
  weights <- matrix(0, ncol = N, nrow = T)
  prediction <- rep(0, T)
  eta <- matrix(exp(350), ncol = N, nrow = T + 1)
  V <- rep(0, N)
  B <- rep(0, N)
  
  if (!is.null(training)) {
    w0 <- training$w0
    R <- training$R
    R.reg <- training$R.reg
    eta[1, ] <- training$eta
    B <- training$B
    V <- training$V
  }
  
  if (use_cpp){
    loss_tau <- ifelse(! is.null(loss.type$tau), loss.type$tau, 0) 
    loss_name <- loss.type$name
    computeBOAEigen(awake,eta,experts,weights,y,prediction,
                    w,w0,R,R.reg,B,V,loss_name,loss_tau,loss.gradient);
  }
  else{
    for (t in 1:T) {
      idx <- awake[t,] > 0
      R.aux <- log(w0) + eta[t, ] * R.reg
      R.max <- max(R.aux[idx])
      w <- numeric(N)
      w[idx] <- exp(R.aux[idx] - R.max)
      
      p <- awake[t, ] * w/sum(awake[t, ] * w)
      pred <- experts[t, ] %*% p
      
      weights[t, ] <- p
      prediction[t] <- pred
      
      lpred <- lossPred(pred, y[t], pred, loss.type = loss.type, loss.gradient = loss.gradient)
      lexp <- lossPred(experts[t, ], y[t], pred, loss.type = loss.type, loss.gradient = loss.gradient)
      
      # Instantaneous regret
      r <-  awake[t, ] * c(c(lpred) - lexp)
      
      # Update the learning rates
      B <- pmax(B, abs(r))
      V <- V + r^2
      eta[t + 1, ] <- pmin(pmin(1/(2 * B), sqrt(log(N)/V)),exp(350))
      
      if (max(eta[t+1, ]) > exp(300)) {
        # if some losses still have not been observed
        r.reg <- r
      } else {
        r.reg <- r - eta[t+1, ] * r^2
      }
      
      # Update the regret and the regularized regret used by BOA
      R <- R + r
      R.reg <- R.reg + r.reg
    }
  }
  
  R.aux <- log(w0) + eta[T + 1, ] * R.reg
  R.max <- max(R.aux)
  w <- exp(R.aux - R.max)
  
  object <- list(model = "BOA", loss.type = loss.type, loss.gradient = loss.gradient, 
                 coefficients = w/sum(w))
  
  object$parameters <- list(eta = eta[1:T, ])
  object$weights <- weights
  object$prediction <- prediction
  
  object$training <- list(eta = eta[T + 1, ], R = R, w0 = w0, R.reg = R.reg, V= V, B=B)
  class(object) <- "mixture"
  return(object)
} 
