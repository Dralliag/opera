lossConv <-
function(p,y,experts,awake=NULL,loss.type = 'squareloss') {

   experts <- as.matrix(experts)
   N <- ncol(experts)  # Nombre d'experts
   T <- nrow(experts)  # Nombre d'instants

   # Experts are always active if awake is unspecified
   if (is.null(awake)) {awake = matrix(1, nrow = T, ncol = N)} 
   awake = as.matrix(awake)
   idx.na <- which(is.na(experts))
   awake[idx.na] <- 0
   experts[idx.na] <- 0

   pond <- awake %*% p
   pred <- ((experts* awake) %*% p) / pond
   if (loss.type == 'squareloss')
      l = sqrt(sum((pred-y)^2 * abs(pond)) / sum(abs(pond)))
   else if (loss.type == 'mae')
      l = sum(abs(pred-y) * abs(pond)) / sum(abs(pond))
   else 
      l = sum(abs(pred-y)/y * abs(pond)) / sum(abs(pond))
   return(l)
}
