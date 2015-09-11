# Unit tests of opera package using testhat package

context("Testing oracle function")

# load some basic data to perform tests
n = 50
X = cbind(rep(0,n),rep(1,n))
Y = rep(0.4,n)
X[n,] = c(1,1)
Y[n] = 1
awake = cbind(rep(c(0,1),n/2),1)
  
  
# Test of loss functions
test_that("loss functions return correct values", {
  expect_that(dim(loss(X,Y, loss.type="squareloss"))[1], equals(n))
})

# Test of oracle functions
test_that("Best expert oracle is ok", {
  m.best_expert = oracle(y = Y,experts = X, oracle = "expert")
  expect_that(m.best_expert$weights[1], equals(1))
  expect_that(m.best_expert$loss, equals(mean((X[,1]-Y)^2)))
  expect_that(sum(m.best_expert$prediction), equals(sum(X[,1])))
  expect_that(m.best_expert$rmse, equals(rmse(X[,1],Y)))
  
  expect_error(oracle(y = Y,experts = X,oracle = "expert", awake = awake),"Sleeping not allowed for best expert oracle.")
  expect_warning(oracle(y = Y,experts = X,oracle = list(name = "expert", lambda = 3)), 
                        "unused lambda parameter")  
  expect_warning(oracle(y = Y,experts = X,oracle = list(name = "expert", niter = 3)), 
                                                                   "unused niter parameter")
  
  # test de la partie quantile
  quantiles = seq(0.1,0.9,by = 0.1)
  K = length(quantiles)
  Y = rnorm(n, mean = 0, sd = 1)
  X = t(matrix(rep(quantile(Y, probs = quantiles),n),nrow = K))
  i = sample(1:K,1) 
  m.best_expert = oracle(y = Y,experts = X, oracle = list(name = "expert", loss.type="pinballloss", tau = quantiles[i]))
  expect_equal(which(m.best_expert$weights == 1),i)
})

test_that("Best convex oracle is ok", {
  m = oracle(y = Y,experts = X, oracle = "convex")
  expect_equal(m$weights[1], 0.6)
  expect_equal(sum(m$weights),1)
  expect_equal(m$loss, 0)
  expect_true(sum(abs(m$prediction - Y))<1e-10)
  expect_equal(m$rmse, 0)
  
  
  m = oracle(y = Y,experts = X, oracle = list(name = "convex",loss.type="mape"))
  expect_true(abs(m$weights[1] - 0.6)<1e-4)
  expect_true(m$loss<1e-4)
  expect_true(sum(abs(m$prediction - Y))<1e-4)
  expect_equal(m$rmse, NULL)
  
  m = oracle(y = Y,experts = X, oracle = list(name = "convex",loss.type="mae"), awake=awake)
  expect_true(abs(m$weights[1] - 0.6)<1e-1)
  expect_equal(m$rmse, NULL)
  expect_equal(mean(loss(m$prediction,Y,"mae")), lossConv(m$weights,Y,X,awake,"mae"))
  expect_equal(m$loss, mean(loss(m$prediction,Y,"mae")))
})

test_that("Best linear oracle is ok", {
  m = oracle(y = Y, experts = X, oracle = "linear")
  expect_equal(m$weights[1],0.6)
  expect_equal(sum(m$weights),1)
  expect_equal(m$loss, 0)
  expect_true(sum(abs(m$prediction - Y))<1e-10)
  expect_equal(m$rmse, 0)
  expect_error(oracle(y = Y, experts = X, oracle = "linear", awake = awake), "Sleeping not allowed for linear oracle!")
})

test_that("Best shifting oracle is ok",{
  m = oracle(y = Y,  experts = X, oracle = "shifting")
  expect_equal(m$loss[1], min(mean(loss(X[,1],Y)), mean(loss(X[,2],Y))))
})


