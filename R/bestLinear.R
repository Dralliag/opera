bestLinear <-
function(y,experts, lambda = 1, awake=NULL, loss.type='squareloss')
{
  if (!is.null(awake)) {stop('Sleeping not allowed here!')}
  if (sum(is.na(experts)>0)) {warning("There are NA's in expert advice")}
  
  experts <- as.matrix(experts)
  u <- solve(lambda * diag(1,ncol(experts)) + t(experts) %*% experts,t(experts)%*%y)
  prev <- experts %*% u
  loss <- sqrt(mean((prev-y)^2))
  return(list(loss=loss,prediction=prev,u=u))
}
