# best sequence of experts oracle
bestShifts <- function(y, experts, awake = NULL, loss.type = list(name = "square"), 
  trace = FALSE) {
  N <- ncol(experts)
  T <- nrow(experts)
  INF <- exp(700)
  # m-1 shifts, expert
  
  # Full activation if unspecified
  if (is.null(awake)) {
    awake <- matrix(1, nrow = T, ncol = N)
  }
  
  awake <- as.matrix(awake)
  idx.na <- which(is.na(experts))
  awake[idx.na] <- 0
  experts[idx.na] <- 0
  
  loss.name <- c("square", "absolute", "percentage")
  loss.number <- which(loss.name == loss.type$name)
  
  L <- array(INF, dim = c(T, N, 3))
  L[1, , ] <- 0
  instanceLoss <- array(0, dim = c(N, 3))
  
  # browser()
  for (t in 1:T) {
    if (!(t%%100) && trace) {
      cat(floor(t^2/T^2 * 10000)/100, "% -- ")
      cat("t = ", t - 1, ": average minimal loss with t shifts : ", min(sqrt(L[t - 
        1, , loss.number]/(t - 1))), " without shifts :", min(sqrt(L[1, , 
        loss.number]/(t - 1))), "\n")
    }
    
    Et1 <- which(awake[t - 1, ] > 0)
    Et <- which(awake[t, ] > 0)
    for (l in 1:3) {
      instanceLoss[, l] <- loss(experts[t, ], y[t], loss.name[l]) * awake[t, 
        ]
    }
    L[1:t, -Et, ] <- INF
    if (t > 1) {
      L1 <- L[1:t, Et1, ]
      idx_min <- apply(L1[, , loss.number], 1, order)[1:2, ]
      for (m in t:2) {
        for (i in Et) {
          if (idx_min[1, m - 1] == i) 
          aux <- idx_min[2, m - 1] else aux <- idx_min[1, m - 1]
          
          if (L[m, i, loss.number] < L1[m - 1, aux, loss.number]) 
          L[m, i, ] <- L[m, i, ] + instanceLoss[i, ] else L[m, i, ] <- L1[m - 1, aux, ] + instanceLoss[i, ]
        }
      }
    }
    L[1, , ] <- L[1, , ] + instanceLoss
  }
  loss.experts <- L[, , loss.number]/T
  loss <- apply(loss.experts, 1, min)
  res <- list(loss = loss)
  if (loss.type == "square") {
    res <- list(loss = loss, rmse = sqrt(loss))
  }
  return(res)
} 
