lossConv <- function(p, y, experts, awake = NULL, loss.type = list(name = "square")) {
  
  N <- ncol(experts)  # Number of experts
  T <- nrow(experts)  # Number of instants
  
  p <- matrix(as.numeric(p), nrow = N)
  
  # Experts are always active if awake is unspecified
  if (is.null(awake)) {
    awake <- matrix(1, nrow = T, ncol = N)
  }
  idx.na <- which(is.na(experts))
  awake[idx.na] <- 0
  experts[idx.na] <- 0
  
  pond <- awake %*% p
  pred <- ((experts * awake) %*% p)/pond
  l <- mean(lossPred(x = pred, y = y, loss.type = loss.type))
  return(l)
} 
