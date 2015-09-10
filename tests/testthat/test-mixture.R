# Unit tests of opera package using testhat package

context("Testing mixture function")

# load some basic data to perform tests
n = 50
X = cbind(rep(0,n),rep(1,n))
Y = rep(0.4,n)
X[n,] = c(1,1)
Y[n] = 1
awake = cbind(rep(c(0,1),n/2),1)


# Test of mixture functions
test_that("EWA functions are ok", {
  m.ewa = mixture(y = Y, experts = X, aggregationRule = "EWA")
  expect_true(abs(m.ewa$weights.forecast[1]-0.6)<1e-10)
  
  eta = 2
  m.ewa.fixed = mixture(y = Y, experts = X, aggregationRule = list(name = "EWA", eta = eta))
  idx.eta  = which(m.ewa$grid == eta)
  expect_equal(m.ewa$gridloss[idx.eta], sum(loss(m.ewa.fixed$prediction,Y)))
})