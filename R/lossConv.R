# Error of a convex combination The function \code{lossConv} computes the loss suffered by the
# constant convex mixture \code{p}.  @param p A non-negative vector that sum to one. It is of length
# \code{N}, the number of experts.  @param y A vector that contains the observations to be
# predicted.  @param experts A matrix containing the experts forecasts. Each column corresponds to
# the predictions proposed by an expert to predict \code{Y}. It has as many columns as there are
# experts.  @param awake A matrix specifying the activation coefficients of the experts. Its entries
# lie in \code{[0,1]}.  Needed if some experts are specialists and do not always form and suggest
# prediction.  If the expert number \code{k} at instance \code{t} does not form any prediction of
# observation \code{Y_t}, we can put \code{awake[t,k]=0} so that the mixture does not consider
# expert \code{k} in the mixture to predict \code{Y_t}.  @param loss.type A string specifying the
# loss function considered to evaluate the performance.  It can be "square", "absolute",
# "percentage", or "pinball". See \code{\link{loss}} for more details.  @param tau Quantile to be
# predicted if loss.type = "pinball" (default value 0.5 to predict the median).  @return The average
# errors suffered by the mixture. Note that if \code{loss.type = "square"}, the \code{rmse} is
# returned.  Note also that the instance are weighted according to the number of activated experts.
# @note The function \code{lossConv} is for instance used to compute the best convex combination in
# hindsight.  @seealso \code{\link{loss}}, @keywords ~kwd1 ~kwd2
lossConv <- function(p, y, experts, awake = NULL, loss.type = "square") {
  
  experts <- as.matrix(experts)
  N <- ncol(experts)  # Number of experts
  T <- nrow(experts)  # Number of instants
  
  p <- matrix(p, nrow = N)
  # Experts are always active if awake is unspecified
  if (is.null(awake)) {
    awake <- matrix(1, nrow = T, ncol = N)
  }
  awake <- as.matrix(awake)
  idx.na <- which(is.na(experts))
  awake[idx.na] <- 0
  experts[idx.na] <- 0
  
  pond <- awake %*% p
  pred <- ((experts * awake) %*% p)/pond
  l <- mean(loss(pred, y, loss.type = loss.type))
  return(l)
} 
