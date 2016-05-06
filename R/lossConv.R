lossConv <- function(p, y, experts, awake = NULL, loss.type = "square") {
  
  experts <- matrix(as.numeric(as.matrix(experts)), nrow = length(y))
  N <- ncol(experts)  # Number of experts
  T <- nrow(experts)  # Number of instants
  
  p <- matrix(as.numeric(as.matrix(p)), nrow = N)
  
  # Experts are always active if awake is unspecified
  if (is.null(awake)) {
    awake <- matrix(1, nrow = T, ncol = N)
  }
  awake <- matrix(as.numeric(as.matrix(awake)), nrow = length(y))
  idx.na <- which(is.na(experts))
  awake[idx.na] <- 0
  experts[idx.na] <- 0
  
  pond <- awake %*% p
  pred <- ((experts * awake) %*% p)/pond
  l <- mean(loss(pred, y, loss.type = loss.type))
  return(l)
} 
