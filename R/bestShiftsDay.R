bestShiftsDay <-
function(y,experts, awake=NULL, loss.type = 'squareloss')
{
    N <- ncol(experts)
    T <- nrow(experts)
    INF <- exp(700)
    # m-1 shifts, expert    

    if (is.null(awake)) {awake = matrix(1, nrow = T, ncol = N)} # Activation 1 si non spécifiée
    awake = as.matrix(awake)
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
