convexLoss <-
function(weights,y,experts,awake=NULL) {
  N <- ncol(experts)
  T <- nrow(experts)

  # Experts are always active if awake is unspecified
  if (is.null(awake)) {awake = matrix(1, nrow = T, ncol = N)} 
  awake = as.matrix(awake)
  idx.na <- which(is.na(experts))
  awake[idx.na] <- 0
  experts[idx.na] <- 0

  pond <- awake %*% weights
  pred <- (experts * awake) %*% weights / pond

  err1 <- loss(pred,y,'squareloss') * pond
  err2 <- loss(pred,y,'mae') * pond 
  err3 <- loss(pred,y,'mape') * pond
  
  mpred <- mean(pred)
  mY <- mean(y)
  merr1 <- sum(err1)  / sum(pond)
  merr2 <- sum(err2)  / sum(pond)
  merr3 <- sum(err3)  / sum(pond)
  
  corr <- sum( (pred - mpred) * (y - mY)) / sqrt(sum((pred - mpred)^2) * sum((y - mY)^2))
  disp1 <- 1.96 / sqrt(T) * sqrt(mean((err1 - merr1)^2) / (4 * merr1))
  disp2 <- 1.96 / sqrt(T) * sd(err2)
  disp3 <- 1.96 / sqrt(T) * sd(err3)
  return(data.frame(rmse = sqrt(merr1), mae = merr2, mape = merr3, disp_rmse = disp1, disp_mae = disp2, disp_mape = disp3, corr = corr))
}
