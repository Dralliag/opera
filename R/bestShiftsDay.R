#' best sequence of experts switching at midnights
#' 
#' 
#' Similarly as \code{\link{bestShifts}}, however no switch between experts are
#' allowed within a day. Switches only occur at midnight (every 48 instances).
#' 
#' 
#' @param y  A vector that contains the
#' observations to be predicted.
#' @param experts matrix that contains the
#' experts forecasts. Each column corresponds to the predictions proposed by an
#' expert to predict \code{y}. It has as many columns as there are experts.
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
#' @return A matrix \code{L} of dimension
#' \code{c(numberofdays,numberofexperts,3)} where the third dimension is the
#' type of loss (1:squareloss, 2:mae, 3:mape), and the value of $L(m,k,l)$ is
#' the loss (determined by \code{l}) suffered by the best sequence of expert
#' with at most $m-1$ shifts and finishing with expert number $k$.
#' @author Pierre Gaillard <pierre-p.gaillard@@edf.fr>
#' @seealso 
#' \code{\link{bestShifts}}
#' @export bestShiftsDay
bestShiftsDay <-
function(y,experts, awake=NULL, loss.type = 'squareloss')
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
    
    nbDays <- T/48
    L <- array(INF, dim = c(nbDays,N,3))
    dayLoss <- array(0,dim=c(N,3))
    L[1,,] <- 0
    for (j in 1:nbDays)
    {
      if (!(j %% floor(nbDays/10))) {
        cat(100 * (j / nbDays),'% -- ')
        print(min(L[(j-1),,1]))      
      }
      t <- 48 * (j-1)+1
      Et1 <- which(awake[t-48,] > 0)
      Et <- which(awake[t,] > 0)
      for (l in 1:3) {
        dayLoss[,l] <- apply(loss(experts[t:(t+47),],y[t:(t+47)],loss.name[l]) * awake[t:(t+47),],2,sum)
      }
      L[1:j,-Et,] <- INF  
      if (j > 1) {
        L1 <- L[1:j, Et1,]
        idx_min <- apply(L1[,,loss.number], 1, order)[1:2,]
        
        for (m in j:2) {
          for (i in Et) {
            if (idx_min[1,m-1] == i)
              aux = idx_min[2,m-1]
            else
              aux = idx_min[1,m-1]
                
            if (L[m,i,loss.number] < L1[m-1,aux,loss.number])
            {
              L[m,i,] <- L[m,i,] + dayLoss[i,]
            }
            else {
              L[m,i,] <- L1[m-1,aux,] + dayLoss[i,]
            }
          }
        }
      }
      L[1,,] <- L[1,,] + dayLoss
    }
    return(L)
}
