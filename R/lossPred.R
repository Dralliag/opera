lossPred <-
function(x, y, pred = NULL, loss.type = 'squareloss', loss.gradient = FALSE, tau = 0.1) {
   npred <- length(pred)
   nx <- length(x)
   if (npred > 1 && nx > 1) {
      if (!loss.gradient) {
         if (loss.type == 'squareloss')
            l = matrix(rep((x-y)^2, npred), ncol = npred)
         else if (loss.type == 'mae')
            l = matrix(rep(abs(x - y), npred), ncol = npred)
         else if (loss.type == 'mape')
            l = matrix(rep(abs(x - y) / y, npred), ncol = npred)
         else if (loss.type == 'pinballloss')
            l = matrix(rep(((y<x)-tau) * (x - y), npred), ncol = npred)
      } else {
         if (loss.type == 'squareloss') 
            l = 2 * t(matrix(rep(pred - y, nx), ncol = nx)) * matrix(rep(x, npred), ncol = npred)
         else if (loss.type == 'mae')
            l = t(matrix(rep(sign(pred - y), nx), ncol = nx)) * matrix(rep(x, npred), ncol = npred)
         else if (loss.type == 'mape')
            l =  matrix(rep(x, npred), ncol = npred) / y * t(matrix(rep(sign(pred - y), nx), ncol = nx))
         else if (loss.type == 'pinballloss')
            l = t(matrix(rep((y < pred)-tau, nx), ncol = nx)) * matrix(rep(x, npred), ncol = npred)
      }
   } else {
      if (!loss.gradient) {
         if (loss.type == 'squareloss')
            l = (x-y)^2
         else if (loss.type == 'mae')
            l = abs(x - y)
         else if (loss.type == 'mape')
            l = abs(x - y) / y
         else if (loss.type == 'pinballloss')
            l = ((y<x)-tau) * (x - y)
      } else {
         if (loss.type == 'squareloss') 
            l = 2 * (pred - y) * x
         else if (loss.type == 'mae')
            l = sign(pred - y) * x
         else if (loss.type == 'mape')
            l = x / y * sign(pred - y)  
         else if (loss.type == 'pinballloss') 
            l = ((y < pred)-tau) * x
      }
   }
   return(l)
}
