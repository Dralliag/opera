
gamMixture <-
function(y, experts, z, lambda, nknots = 5, 
  degree = 3, loss.type = 'squareloss', 
  uniform = FALSE, knots = NULL) {

  # building the design matrix with splines  (experts that may smoothly depend on z)  
  if (is.null(knots)) {
    if (!uniform) {
      knots = quantile(z,probs = seq(0,1,length.out=nknots+2)[-c(1,nknots+2)])
    } else {
      knots = seq(min(z),max(z),length=nknots+2)[-c(1,nknots+2)]
    }
  } else {
    nknots = length(knots)
  }
  B = splines::bs(z, knots=knots, degree = degree, intercept=T)
  X = matrix(apply(experts,2,function(a){a*B}),nrow=length(z))
  
  N = ncol(experts)

  # We apply Ridge to this new matrix of experts
  u <- ridge(y, X, lambda, w0 = rep(c(1,rep(0,nknots+degree)),N)/N)
  prediction = apply(X*u$weights,1,sum)
  prediction[1] = mean(experts[1,])
  
  return(list(weight=u$weights, design = B, knots = knots, prediction=prediction))
}
