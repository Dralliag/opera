loss <-
function(x,y,loss.type = 'squareloss', tau = 0.1) {
   if (loss.type == 'squareloss')
      l <- (x-y)^2
   else if (loss.type == 'mae')
      l <- abs(x-y)
   else if (loss.type == 'mape')
      l <- abs(x-y)/y
   else if (loss.type == 'pinballloss')
      l <- (tau - (y>x)) * (x - y)
   return(l)
}
