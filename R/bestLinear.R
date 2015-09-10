# best linear oracle
bestLinear <-
function(y,experts, lambda = 0)
{
  if (sum(is.na(experts)>0)) {warning("NA not allowed in expert advice for linear oracle")}
  
  experts <- as.matrix(experts)
  u <- solve(lambda * diag(1,ncol(experts)) + t(experts) %*% experts,t(experts)%*%y)
  prev <- experts %*% u
  loss <- mean((prev-y)^2)
  return(list(loss=loss,weights=u,prediction=prev,rmse=sqrt(loss)))
}
