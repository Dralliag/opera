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
test_that("EWA is ok", {
  w0 = c(0.3,0.7)
  possible_loss_type = c("percentage", "absolute", "square", "pinball")
  i.loss = sample(1:4,1)
  tau = runif(1)
  m = mixture(y = Y, experts = X, aggregationRule = list(name = "EWA", loss.type = possible_loss_type[i.loss], tau = tau),
              w0 = w0)
  expect_true(abs(m$coefficients[1]-0.6)<1e-1)
  expect_equal(m$loss, mean(loss(m$prediction,Y,loss.type = possible_loss_type[i.loss], tau = tau)))
  expect_identical(m$weights[1,],w0)
  
  eta = 0.5
  m.fixed = mixture(y = Y, experts = X, aggregationRule = list(name = "EWA", eta = eta, loss.type = possible_loss_type[i.loss], tau = tau), w0 = w0)
  idx.eta  = which(m$grid.eta == eta)
  expect_equal(m$grid.loss[idx.eta], mean(loss(m.fixed$prediction,Y, loss.type = possible_loss_type[i.loss], tau =tau)))
  expect_equal(m.fixed$loss, m$grid.loss[idx.eta])
  expect_identical(m.fixed$weights[1,],w0)
  
  m = mixture(y = Y, experts = X, aggregationRule = "EWA", awake = awake)
  expect_true(abs(m$coefficients[1]-0.6)<1e-1)
  expect_equal(m$loss, mean(loss(m$prediction,Y)))
  
  grid.eta = runif(5)  
  m = mixture(y = Y, experts = X, aggregationRule = list(name = "EWA", grid.eta = grid.eta), awake = awake)
  expect_equal(sum(!(grid.eta %in% m$grid.eta)),0)
})

# Test of Fixed-share
test_that("Fixed-share is ok", {
  w0 = c(0.3,0.7)
  possible_loss_type = c("percentage", "absolute", "square", "pinball")
  i.loss = sample(1:4,1)
  tau = runif(1)
  m = mixture(y = Y, experts = X, aggregationRule = list(name = "FS", loss.type = possible_loss_type[i.loss], tau = tau),
              w0 = w0)
  expect_true(abs(m$coefficients[1]-0.6)<1e-1)
  expect_equal(m$loss, mean(loss(m$prediction,Y,loss.type = possible_loss_type[i.loss], tau = tau)))
  expect_identical(m$weights[1,],w0)
  
  eta = 2
  alpha = 0.01
  m.fixed = mixture(y = Y, experts = X, 
                    aggregationRule = list(name = "FS", eta = eta, alpha = alpha, loss.type = possible_loss_type[i.loss], tau = tau), 
                    w0 = w0)
  idx.eta  = which(m$grid.eta == eta)
  idx.alpha = which(m$grid.alpha == alpha)
  expect_equal(m$grid.loss[idx.eta,idx.alpha], mean(loss(m.fixed$prediction,Y, loss.type = possible_loss_type[i.loss], tau =tau)))
  expect_equal(m.fixed$loss, m$grid.loss[idx.eta,idx.alpha])
  expect_identical(m.fixed$weights[1,],w0)
  
  m = mixture(y = Y, experts = X, aggregationRule = "FS", awake = awake)
  expect_true(abs(m$coefficients[1]-0.6)<1e-1)
  expect_equal(m$loss, mean(loss(m$prediction,Y)))
  
  grid.eta = runif(5)
  grid.alpha = runif(3)
  m = mixture(y = Y, experts = X, aggregationRule = list(name = "FS", grid.eta = grid.eta, grid.alpha = grid.alpha), awake = awake)
  expect_equal(sum(!(grid.eta %in% m$grid.eta)),0)
  expect_identical(grid.alpha, m$grid.alpha)
  
})

# Test of Ridge
test_that("Ridge is ok", {
  w0 = c(0.5,0.5)
  m = mixture(y = Y, experts = X, aggregationRule = list(name = "Ridge"), w0 = w0)
  expect_equal(m$loss, mean(loss(m$prediction,Y)))
  expect_identical(m$weights[1,],w0)
  expect_true(!is.na(sum(m$weights)))
  
  lambda = 2
  m.fixed = mixture(y = Y, experts = X, aggregationRule = list(name = "Ridge", lambda = lambda), w0 = w0)
  idx.lambda  = which(m$grid.lambda == lambda)
  expect_equal(m$grid.loss[idx.lambda], mean(loss(m.fixed$prediction,Y)))
  expect_equal(m.fixed$loss, m$grid.loss[idx.lambda])
  expect_identical(m.fixed$weights[1,],w0)
  
  grid.lambda = runif(3)
  m = mixture(y = Y, experts = X, aggregationRule = list(name = "Ridge", grid.lambda = grid.lambda, gamma = 100), awake = awake)
  expect_equal(sum(!(grid.lambda %in% m$grid.lambda)),0)
  
})


# test of MLPol
test_that("MLpol, MLprod, MLewa, and BOA are ok", {
  possible_loss_type = c("percentage", "absolute", "square", "pinball")
  i.loss = sample(1:4,1)
  tau = runif(1)
  m = mixture(y = Y, experts = X, 
              aggregationRule = list(name = "MLpol", loss.type = possible_loss_type[i.loss], tau = tau))
  expect_true(abs(m$coefficients[1]-0.6)<2e-1)
  expect_equal(m$loss, mean(loss(m$prediction,Y,loss.type = possible_loss_type[i.loss], tau = tau)))
  
  w0 = c(0.3, 0.7)
  m = mixture(y = Y, experts = X, 
              aggregationRule = list(name = "MLprod", loss.type = possible_loss_type[i.loss], tau = tau),
              w0 = w0)
  expect_true(abs(m$coefficients[1]-0.6)<2e-1)
  expect_equal(m$loss, mean(loss(m$prediction,Y,loss.type = possible_loss_type[i.loss], tau = tau)))
  expect_identical(m$weights[1,],w0)
  
  m = mixture(y = Y, experts = X, 
              aggregationRule = list(name = "MLewa", loss.type = possible_loss_type[i.loss], tau = tau),
              w0 = w0)
  expect_true(abs(m$coefficients[1]-0.6)<2e-1)
  expect_equal(m$loss, mean(loss(m$prediction,Y,loss.type = possible_loss_type[i.loss], tau = tau)))
  expect_identical(m$weights[1,],w0)
  
  m = mixture(y = Y, experts = X, 
              aggregationRule = list(name = "BOA", loss.type = possible_loss_type[i.loss], tau = tau),
              w0 = w0)
  expect_true(abs(m$coefficients[1]-0.6)<2e-1)
  expect_equal(m$loss, mean(loss(m$prediction,Y,loss.type = possible_loss_type[i.loss], tau = tau)))
  expect_identical(m$weights[1,],w0)
  
  
  m = mixture(y = Y, experts = X, aggregationRule = "BOA", awake = awake)
  expect_true(abs(m$coefficients[1]-0.6)<2e-1)
  expect_equal(m$loss, mean(loss(m$prediction,Y)))
})



test_that("Quantile mixture are ok", {
  # test de la partie quantile
  n = 200
  quantiles = seq(0.1,0.9,by = 0.1)
  K = length(quantiles)
  Y = rnorm(n, mean = 0, sd = 1)
  X = t(matrix(rep(quantile(Y, probs = quantiles),n),nrow = K))
  i = sample(1:K,1) 
  
  m = mixture(y = Y, experts = X, 
              aggregationRule = list(name = "EWA", loss.type = "pinball", loss.gradient = FALSE, eta = 10e-1, tau = quantiles[i], gamma = 100))
  expect_equal(m$loss,mean(loss(m$prediction,Y,loss.type = "pinball", tau = quantiles[i])))
  expect_less_than(abs(sum(X[1,] * m$coefficients) - X[1,i]), 0.4)
  
  m = mixture(y = Y, experts = X[,c(1,K)], 
              aggregationRule = list(name = "EWA", loss.type = "pinball", tau = quantiles[i], gamma = 100))
  expect_equal(m$loss,mean(loss(m$prediction,Y,loss.type = "pinball", tau = quantiles[i])))
  expect_less_than(abs(sum(X[1,c(1,K)] * m$coefficients) - X[1,i]), 0.4)
  
  # Fixed share
  m = mixture(y = Y, experts = X, 
              aggregationRule = list(name = "FS", loss.type = "pinball", loss.gradient = FALSE, eta = 10e-1, alpha = 0.01, tau = quantiles[i]))
  expect_equal(m$loss,mean(loss(m$prediction,Y,loss.type = "pinball", tau = quantiles[i])))
  expect_less_than(abs(sum(X[1,] * m$coefficients) - X[1,i]), 0.4)
  
  m = mixture(y = Y, experts = X[,c(1,K)], 
              aggregationRule = list(name = "FS", loss.type = "pinball", tau = quantiles[i], gamma = 10))
  expect_equal(m$loss,mean(loss(m$prediction,Y,loss.type = "pinball", tau = quantiles[i])))
  expect_less_than(abs(sum(X[1,c(1,K)] * m$coefficients) - X[1,i]), 0.4)
  
  m = mixture(y = Y, experts = X, 
              aggregationRule = list(name = "MLpol", loss.type = "pinball", loss.gradient = FALSE,
                                     tau = quantiles[i]))
  expect_equal(m$loss,mean(loss(m$prediction,Y,loss.type = "pinball", tau = quantiles[i])))
  expect_less_than(abs(sum(X[1,] * m$coefficients) - X[1,i]), 0.4)
  
  m = mixture(y = Y, experts = X[,c(1,K)], 
              aggregationRule = list(name = "MLpol", loss.type = "pinball", tau = quantiles[i]))
  expect_equal(m$loss,mean(loss(m$prediction,Y,loss.type = "pinball", tau = quantiles[i])))
  expect_less_than(abs(sum(X[1,c(1,K)] * m$coefficients) - X[1,i]), 0.8)
  
  m = mixture(y = Y, experts = X, 
              aggregationRule = list(name = "MLprod", loss.type = "pinball", loss.gradient = FALSE,
                                     tau = quantiles[i]))
  expect_equal(m$loss,mean(loss(m$prediction,Y,loss.type = "pinball", tau = quantiles[i])))
  expect_less_than(abs(sum(X[1,] * m$coefficients) - X[1,i]), 0.4)
  
  m = mixture(y = Y, experts = X[,c(1,K)], 
              aggregationRule = list(name = "MLprod", loss.type = "pinball", tau = quantiles[i]))
  expect_equal(m$loss,mean(loss(m$prediction,Y,loss.type = "pinball", tau = quantiles[i])))
  expect_less_than(abs(sum(X[1,c(1,K)] * m$coefficients) - X[1,i]), 0.8)
  
  m = mixture(y = Y, experts = X, 
              aggregationRule = list(name = "MLewa", loss.type = "pinball", loss.gradient = FALSE,
                                     tau = quantiles[i]))
  expect_equal(m$loss,mean(loss(m$prediction,Y,loss.type = "pinball", tau = quantiles[i])))
  expect_less_than(abs(sum(X[1,] * m$coefficients) - X[1,i]), 0.4)
  
  m = mixture(y = Y, experts = X[,c(1,K)], 
              aggregationRule = list(name = "MLewa", loss.type = "pinball", tau = quantiles[i]))
  expect_equal(m$loss,mean(loss(m$prediction,Y,loss.type = "pinball", tau = quantiles[i])))
  expect_less_than(abs(sum(X[1,c(1,K)] * m$coefficients) - X[1,i]), 0.8)
  
  m = mixture(y = Y, experts = X, 
              aggregationRule = list(name = "BOA", loss.type = "pinball", loss.gradient = FALSE,
                                     tau = quantiles[i]))
  expect_equal(m$loss,mean(loss(m$prediction,Y,loss.type = "pinball", tau = quantiles[i])))
  expect_less_than(abs(sum(X[1,] * m$coefficients) - X[1,i]), 0.4)
  
  m = mixture(y = Y, experts = X[,c(1,K)], 
              aggregationRule = list(name = "BOA", loss.type = "pinball", tau = quantiles[i]))
  expect_equal(m$loss,mean(loss(m$prediction,Y,loss.type = "pinball", tau = quantiles[i])))
  expect_less_than(abs(sum(X[1,c(1,K)] * m$coefficients) - X[1,i]), 0.8)
})
