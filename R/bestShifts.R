# best sequence of experts oracle
bestShifts <- function(y, experts, awake = NULL, loss.type = list(name = "square")) {
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
  
  L <- array(INF, dim = c(T, N))
  L[1, ] <- 0
  
  for (t in 1:T) {
    Et1 <- which(awake[t - 1, ] > 0)
    Et <- which(awake[t, ] > 0)
    for (l in 1:3) {
      instanceLoss <- loss(experts[t, ], y[t], loss.type) * awake[t, 
        ]
    }
    L[1:t, -Et] <- INF
    if (t > 1) {
      L1 <- L[1:t, Et1]
      idx_min <- apply(L1[, ], 1, order)[1:2, ]
      for (m in t:2) {
        for (i in Et) {
          if (idx_min[1, m - 1] == i) 
            aux <- idx_min[2, m - 1] else aux <- idx_min[1, m - 1]
            
            if (L[m, i] < L1[m - 1, aux]) 
              L[m, i] <- L[m, i] + instanceLoss[i] else L[m, i] <- L1[m - 1, aux] + instanceLoss[i]
        }
      }
    }
    L[1, ] <- L[1, ] + instanceLoss
  }
  loss.experts <- L[, ]/T
  loss <- apply(loss.experts, 1, min)
  for (i in 2:T) {
    if (loss[i] > loss[i-1]) {
      loss[i] <- loss[i-1]
    }
  }
  res <- list(loss = loss)
  if (is.list(loss.type) && loss.type$name == "square") {
    res <- list(loss = loss, rmse = sqrt(loss))
  }
  return(res)
} 
