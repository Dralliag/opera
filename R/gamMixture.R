
gamMixture <-
function(y, experts, z, lambda, nknots = 5, degree = 3, loss.type = 'squareloss',
                       href = 1, period = 1, uniform = F, knots = NULL, tau = 0.5) {
  # création de la matrice de design
  
  if (is.null(knots)) {
    if (!uniform) {
      knots = quantile(z,probs = seq(0,1,length.out=nknots+2)[-c(1,nknots+2)])
    } else {
      knots = seq(min(z),max(z),length=nknots+2)[-c(1,nknots+2)]
    }
  } else {
    nknots = length(knots)
  }
  B = bs(z, knots=knots, degree = degree, intercept=T)
  X = matrix(apply(experts,2,function(a){a*B}),nrow=length(z))
  
  N = ncol(experts)
  # On applique un ridge normal sur cette nouvelle matrice d'experts
  #u <- ridgeHour(y, X, lambda, href = Hour, period = 48, w0 = rep(c(1,rep(0,nknots+degree)),N)/N)
  u <- pinballHour(y, X, lambda, href = href, period = period, tau = tau)
  prediction = apply(X*u$weights,1,sum)
  prediction[1:(href+period)] = apply(experts[1:(href+period),],1,mean)
  
  return(list(weight=u$weights, design = B, knots = knots, prediction=prediction))
}
