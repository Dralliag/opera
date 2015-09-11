# Unit tests of opera package using testhat package

context("Testing mixture function")

# load some basic data to perform tests
n = 50
X = cbind(rep(0,n),rep(1,n))
Y = rep(0.4,n)
X[n,] = c(1,1)
Y[n] = 1
awake = cbind(rep(c(0,1),n/2),1)


# Test of EWA
test_that("EWA functions are ok", {
  w0 = c(0.3,0.7)
  possible_loss_type = c("mape", "mae", "squareloss", "pinballloss")
  i.loss = sample(1:4,1)
  tau = runif(1)
  m.ewa = mixture(y = Y, experts = X, aggregationRule = list(name = "EWA", loss.type = possible_loss_type[i.loss], tau = tau),
                  w0 = w0)
  expect_true(abs(m.ewa$weights.forecast[1]-0.6)<1e-1)
  expect_equal(m.ewa$loss, mean(loss(m.ewa$prediction,Y,loss.type = possible_loss_type[i.loss], tau = tau)))
  expect_identical(m.ewa$weights[1,],w0)
  
  eta = 2
  m.ewa.fixed = mixture(y = Y, experts = X, aggregationRule = list(name = "EWA", eta = eta, loss.type = possible_loss_type[i.loss], tau = tau), w0 = w0)
  idx.eta  = which(m.ewa$grid.eta == eta)
  expect_equal(m.ewa$grid.loss[idx.eta], mean(loss(m.ewa.fixed$prediction,Y, loss.type = possible_loss_type[i.loss], tau =tau)))
  expect_equal(m.ewa.fixed$loss, m.ewa$grid.loss[idx.eta])
  expect_identical(m.ewa.fixed$weights[1,],w0)
  
  m.ewa = mixture(y = Y, experts = X, aggregationRule = "EWA", awake = awake)
  expect_true(abs(m.ewa$weights.forecast[1]-0.6)<1e-1)
  expect_equal(m.ewa$loss, mean(loss(m.ewa$prediction,Y)))
  
})

test_that("Quantile mixture are ok", {
  # test de la partie quantile
  n = 200
  quantiles = seq(0.1,0.9,by = 0.1)
  K = length(quantiles)
  Y = rnorm(n, mean = 0, sd = 1)
  X = t(matrix(rep(quantile(Y, probs = quantiles),n),nrow = K))
  i = sample(1:K,1) 
  
  m.ewa = mixture(y = Y, experts = X, 
                  aggregationRule = list(name = "EWA", loss.type = "pinballloss", loss.gradient = FALSE, eta = 10e-1, tau = quantiles[i], gamma = 100))
  expect_equal(m.ewa$loss,mean(loss(m.ewa$prediction,Y,loss.type = "pinballloss", tau = quantiles[i])))
  expect_less_than(abs(sum(X[1,] * m.ewa$weights.forecast) - X[1,i]), 0.4)
  
  m.ewa = mixture(y = Y, experts = X[,c(1,K)], 
                  aggregationRule = list(name = "EWA", loss.type = "pinballloss", tau = quantiles[i], gamma = 100))
  expect_equal(m.ewa$loss,mean(loss(m.ewa$prediction,Y,loss.type = "pinballloss", tau = quantiles[i])))
  expect_less_than(abs(sum(X[1,c(1,K)] * m.ewa$weights.forecast) - X[1,i]), 0.4)
})
