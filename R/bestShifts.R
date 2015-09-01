#' best sequence of experts oracle
#' 
#'  The
#' function \code{bestShifts} computes for all number m of stwitches the
#' sequence of experts with at most $m$ shifts that would have performed the
#' best to predict the sequence of observations in \code{y}.
#' 
#' 
#' @param y  vector that contains the observations
#' to be predicted.
#' @param experts A matrix containing the
#' experts forecasts. Each column corresponds to the predictions proposed by an
#' expert to predict \code{Y}. It has as many columns as there are experts.
#' @param awake A matrix specifying the
#' activation coefficients of the experts. Its entries lie in \code{[0,1]}.
#' Needed if some experts are specialists and do not always form and suggest
#' prediction.  If the expert number \code{k} at instance \code{t} does not
#' form any prediction of observation \code{Y_t}, we can put
#' \code{awake[t,k]=0} so that the mixture does not consider expert \code{k} in
#' the mixture to predict \code{Y_t}.
#' @param loss.type A string specifying
#' the loss function considered to evaluate the performance. It can be
#' "squareloss", "mae", "mape", or "pinballloss". See \code{\link{loss}} for
#' more details.
#' @return  %% Returns a matrix of dimension \code{c(T,N,3)} where
#' \code{T} is the number of instance to be predicted (i.e., the length of the
#' sequence \code{y}) and \code{N} is the number of experts.  A matrix \code{L}
#' of dimension \code{c(T,N,3)} where the third dimension is the type of loss
#' (1:squareloss, 2:mae, 3:mape), and the value of $L(m,k,l)$ is the loss
#' (determined by \code{l}) suffered by the best sequence of expert with at
#' most $m-1$ shifts and finishing with expert number $k$.
#' @author Pierre Gaillard <pierre-p.gaillard@@edf.fr>
#' @seealso 
#' \code{\link{bestConvex}}, \code{\link{bestShiftsDay}}, \code{\link{loss}}
#' @keywords ~kwd1 ~kwd2
#' @export bestShifts
bestShifts <-
function(y, experts, awake=NULL, loss.type = 'squareloss')
{
    N <- ncol(experts)
    T <- nrow(experts)
    INF <- exp(700)
    # m-1 shifts, expert
    
    if (is.null(awake)) {awake = matrix(1, nrow = T, ncol = N)} # Full activation if unspecified

    awake <- as.matrix(awake)
    idx.na <- which(is.na(experts))
    awake[idx.na] <- 0
    experts[idx.na] <- 0

    loss.name <- c('squareloss', 'mae','mape')
    loss.number <- which(loss.name == loss.type)
    
    L <- array(INF, dim = c(T,N,3))
    L[1,,] <- 0
    instanceLoss <- array(0,dim=c(N,3))
    for (t in 1:T)
    {
      if (!(t %% 100)) {
        cat(floor(t^2/T^2 * 10000)/100, '% -- ')
        cat('actual minimal loss with t shifts : ', min(sqrt(L[t-1,,1]/(t-1))), ' without shifts :', min(sqrt(L[1,,1]/(t-1))), '\n')
        save(L,file = "L.rdata")
      }
        
      Et1 <- which(awake[t-1,] > 0)
      Et <- which(awake[t,] > 0)
      for (l in 1:3) {
        instanceLoss[,l] <- loss(experts[t,],y[t],loss.name[l]) * awake[t,]
      }
      L[1:t,-Et,] <- INF
      if (t > 1) {
        L1 <- L[1:t, Et1,]
        idx_min <- apply(L1[,,loss.number], 1, order)[1:2,]
        for (m in t:2) {
          for (i in Et) {
              if (idx_min[1,m-1] == i)
                aux = idx_min[2,m-1]
              else
                aux = idx_min[1,m-1]
                
              if (L[m,i,loss.number] < L1[m-1,aux,loss.number])
                L[m,i,] <- L[m,i,] + instanceLoss[i,]
              else
                L[m,i,] <- L1[m-1,aux,] + instanceLoss[i,]
          }
        }
      }
      L[1,,] <- L[1,,] + instanceLoss
    }
  save(L,file = "bestShifts.rdata")
  return(L)
}
