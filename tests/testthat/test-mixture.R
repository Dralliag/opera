# Unit tests of opera package using testhat package

context("Testing mixture function")

# load some basic data to perform tests
n <- 50
X <- cbind(rep(0, n), rep(1, n))
Y <- rep(0.4, n)
X[n, ] <- c(1, 1)
Y[n] <- 1
awake <- cbind(rep(c(0, 1), n/2), 1)


# Test of EWA
test_that("EWA is ok", {
  w0 <- c(0.3, 0.7)
  possible_loss_type <- c("percentage", "absolute", "square", "pinball")
  i.loss <- sample(1:4, 1)
  m <- mixture(Y = Y, experts = X, model = "EWA", loss.type = possible_loss_type[i.loss], coefficients = w0)
  expect_true(abs(m$coefficients[1] - 0.6) < 0.1)
  expect_equal(m$loss, mean(loss(m$prediction, Y, loss.type = possible_loss_type[i.loss])))
  expect_identical(m$weights[1, ], w0)
  # expect_null(predict(m)) e <- c(0.3,0.5) expect_equal(c(predict(m,e)), sum(c(e)*c(m$coefficients)))
  
  
  m <- mixture(Y = Y, experts = X, model = "EWA", loss.type = possible_loss_type[i.loss], coefficients = w0)
  
  
  eta <- 0.5
  m.fixed <- mixture(Y = Y, experts = X, model = "EWA", parameters = list(eta = eta), loss.type = possible_loss_type[i.loss], 
    coefficients = w0)
  idx.eta <- which(m$parameters$grid.eta == eta)
  expect_equal(m$training$grid.loss[idx.eta], mean(loss(m.fixed$prediction, Y, loss.type = possible_loss_type[i.loss])))
  expect_equal(m.fixed$loss, m$training$grid.loss[idx.eta])
  expect_identical(m.fixed$weights[1, ], w0)
  
  m <- mixture(Y = Y, experts = X, model = "EWA", awake = awake)
  expect_true(abs(m$coefficients[1] - 0.6) < 0.1)
  expect_equal(m$loss, mean(loss(m$prediction, Y)))
  # e <- c(0.3,0.5) expect_equal(c(predict(m,e)), sum(c(e)*c(m$coefficients)))
  
  grid.eta <- runif(5)
  m <- mixture(Y = Y, experts = X, model = "EWA", parameters = list(grid.eta = grid.eta), awake = awake)
  expect_equal(sum(!(grid.eta %in% m$parameters$grid.eta)), 0)
  
  m1 <- mixture(Y = Y[1:10], experts = X[1:10, ], model = "EWA", parameters = list(grid.eta = grid.eta), 
    awake = awake[1:10, ])
  m1 <- predict(object = m1, newexperts = X[-c(1:10), ], newY = Y[-c(1:10)], awake = awake[-c(1:10), 
    ], online = TRUE, type = "model")
  expect_equal(m, m1)
  
})

# Test of Fixed-share
test_that("Fixed-share is ok", {
  w0 <- c(0.3, 0.7)
  possible_loss_type <- c("percentage", "absolute", "square", "pinball")
  i.loss <- sample(1:4, 1)
  m <- mixture(Y = Y, experts = X, model = "FS", loss.type = possible_loss_type[i.loss], coefficients = w0)
  expect_true(abs(m$coefficients[1] - 0.6) < 0.1)
  expect_equal(m$loss, mean(loss(m$prediction, Y, loss.type = possible_loss_type[i.loss])))
  expect_identical(m$weights[1, ], w0)
  
  # e <- c(0.3,0.5) expect_equal(c(predict(m,e)), sum(c(e)*c(m$coefficients)))
  eta <- 2
  alpha <- 0.01
  m.fixed <- mixture(Y = Y, experts = X, model = "FS", parameters = list(eta = eta, alpha = alpha), 
    loss.type = possible_loss_type[i.loss], coefficients = w0)
  idx.eta <- which(m$parameters$grid.eta == eta)
  idx.alpha <- which(m$parameters$grid.alpha == alpha)
  expect_equal(m$training$grid.loss[idx.eta, idx.alpha], mean(loss(m.fixed$prediction, Y, loss.type = possible_loss_type[i.loss])))
  expect_equal(m.fixed$loss, m$training$grid.loss[idx.eta, idx.alpha])
  expect_identical(m.fixed$weights[1, ], w0)
  
  m <- mixture(Y = Y, experts = X, model = "FS", awake = awake)
  expect_true(abs(m$coefficients[1] - 0.6) < 0.1)
  expect_equal(m$loss, mean(loss(m$prediction, Y)))
  
  grid.eta <- runif(5)
  grid.alpha <- runif(3)
  m <- mixture(Y = Y, experts = X, model = "FS", parameters = list(grid.eta = grid.eta, grid.alpha = grid.alpha), 
    awake = awake)
  expect_equal(sum(!(grid.eta %in% m$parameters$grid.eta)), 0)
  expect_identical(grid.alpha, m$parameters$grid.alpha)
})

# Test of Ridge
test_that("Ridge is ok", {
  w0 <- c(0.5, 0.5)
  m <- mixture(Y = Y, experts = X, model = "Ridge", coefficients = w0)
  expect_equal(m$loss, mean(loss(m$prediction, Y)))
  expect_identical(m$weights[1, ], w0)
  expect_true(!is.na(sum(m$weights)))
  
  # e <- c(0.3,0.5) expect_equal(c(predict(m,e)), sum(c(e)*c(m$coefficients)))
  
  lambda <- 2
  m.fixed <- mixture(Y = Y, experts = X, model = "Ridge", parameters = list(lambda = lambda), coefficients = w0)
  idx.lambda <- which(m$parameters$grid.lambda == lambda)
  expect_equal(m$training$grid.loss[idx.lambda], mean(loss(m.fixed$prediction, Y)))
  expect_equal(m.fixed$loss, m$training$grid.loss[idx.lambda])
  expect_identical(m.fixed$weights[1, ], w0)
  
  grid.lambda <- runif(3)
  m <- mixture(Y = Y, experts = X, model = "Ridge", parameters = list(grid.lambda = grid.lambda, gamma = 100), 
    awake = awake)
  expect_equal(sum(!(grid.lambda %in% m$parameters$grid.lambda)), 0)
})


# test of MLPol,...
test_that("MLpol, MLprod, MLewa, and BOA are ok", {
  possible_loss_type <- c("percentage", "absolute", "square", "pinball")
  i.loss <- sample(1:4, 1)
  m <- mixture(Y = Y, experts = X, loss.type = possible_loss_type[i.loss])
  expect_true(abs(m$coefficients[1] - 0.6) < 0.2)
  expect_equal(m$loss, mean(loss(m$prediction, Y, loss.type = possible_loss_type[i.loss])))
  
  m1 <- mixture(loss.type = possible_loss_type[i.loss])
  m1 <- predict(object = m1, newexperts = X, newY = Y, online = TRUE, type = "model")
  expect_equal(m, m1)
  
  m1 <- mixture(Y = Y[1:10], experts = X[1:10, ], loss.type = possible_loss_type[i.loss])
  m1 <- predict(object = m1, newexperts = X[-c(1:10), ], newY = Y[-c(1:10)], online = TRUE, type = "model")
  expect_equal(m, m1)
  
  w0 <- c(0.3, 0.7)
  m <- mixture(Y = Y, experts = X, model = "MLprod", coefficients = w0, loss.type = possible_loss_type[i.loss])
  expect_true(abs(m$coefficients[1] - 0.6) < 0.2)
  expect_equal(m$loss, mean(loss(m$prediction, Y, loss.type = possible_loss_type[i.loss])))
  expect_identical(m$weights[1, ], w0)
  
  m <- mixture(Y = Y, experts = X, model = "MLewa", coefficients = w0, loss.type = possible_loss_type[i.loss])
  expect_true(abs(m$coefficients[1] - 0.6) < 0.2)
  expect_equal(m$loss, mean(loss(m$prediction, Y, loss.type = possible_loss_type[i.loss])))
  expect_identical(m$weights[1, ], w0)
  e <- c(0.3, 0.5)
  # expect_equal(c(predict(m,e)), sum(c(e)*c(m$coefficients))) à vérifier plus tard
  
  w0 <- c(0.3, 0.7)
  m <- mixture(Y = Y, experts = X, model = "BOA", coefficients = w0, loss.type = possible_loss_type[i.loss])
  expect_true(abs(m$coefficients[1] - 0.6) < 0.2)
  expect_equal(m$loss, mean(loss(m$prediction, Y, loss.type = possible_loss_type[i.loss])))
  expect_identical(m$weights[1, ], w0)
  
  
  m <- mixture(Y = Y[1:5], experts = X[1:5, ], model = "MLewa", awake = awake[1:5, ])
  m <- predict(m, newexperts = X[-c(1:5), ], newY = Y[-c(1:5)], awake = awake[-c(1:5), ])
  expect_true(abs(m$coefficients[1] - 0.6) < 0.2)
  expect_equal(m$loss, mean(loss(m$prediction, Y)))
  
  m1 <- mixture(Y = Y, experts = X, model = "MLewa", awake = awake)
  expect_equal(m, m1)
})


test_that("Quantile mixture are ok", {
  # test de la partie quantile
  n <- 200
  quantiles <- seq(0.1, 0.9, by = 0.1)
  K <- length(quantiles)
  Y <- rnorm(n, mean = 0, sd = 1)
  X <- t(matrix(rep(quantile(Y, probs = quantiles), n), nrow = K))
  i <- sample(1:K, 1)
  l <- list(name = "pinball", tau = quantiles[i])
  m <- mixture(Y = Y, experts = X, model = "EWA", loss.type = l, loss.gradient = FALSE, parameters = list(eta = 1, 
    gamma = 100))
  expect_equal(m$loss, mean(loss(m$prediction, Y, loss.type = l)))
  expect_less_than(abs(sum(X[1, ] * m$coefficients) - X[1, i]), 0.4)
  # e <- rnorm(K) expect_equal(c(predict(m,e)), sum(c(e)*c(m$coefficients)))
  
  m <- mixture(Y = Y, experts = X[, c(1, K)], model = "EWA", loss.type = l, parameters = list(gamma = 100))
  expect_equal(m$loss, mean(loss(m$prediction, Y, loss.type = l)))
  expect_less_than(abs(sum(X[1, c(1, K)] * m$coefficients) - X[1, i]), 0.4)
  
  # Fixed share
  m <- mixture(Y = Y, experts = X, model = "FS", loss.type = l, loss.gradient = FALSE, parameters = list(eta = 1, 
    alpha = 0.01))
  expect_equal(m$loss, mean(loss(m$prediction, Y, loss.type = l)))
  expect_less_than(abs(sum(X[1, ] * m$coefficients) - X[1, i]), 0.4)
  
  m <- mixture(Y = Y, experts = X[, c(1, K)], model = "FS", loss.type = l, parameters = list(gamma = 10))
  expect_equal(m$loss, mean(loss(m$prediction, Y, loss.type = l)))
  expect_less_than(abs(sum(X[1, c(1, K)] * m$coefficients) - X[1, i]), 0.4)
  
  m <- mixture(Y = Y, experts = X, model = "MLpol", loss.type = l, loss.gradient = FALSE)
  expect_equal(m$loss, mean(loss(m$prediction, Y, loss.type = l)))
  expect_less_than(abs(sum(X[1, ] * m$coefficients) - X[1, i]), 0.4)
  
  m <- mixture(Y = Y, experts = X[, c(1, K)], model = "MLpol", loss.type = l)
  expect_equal(m$loss, mean(loss(m$prediction, Y, loss.type = l)))
  expect_less_than(abs(sum(X[1, c(1, K)] * m$coefficients) - X[1, i]), 0.8)
  
  m <- mixture(Y = Y, experts = X, model = "MLprod", loss.type = l, loss.gradient = FALSE)
  expect_equal(m$loss, mean(loss(m$prediction, Y, loss.type = l)))
  expect_less_than(abs(sum(X[1, ] * m$coefficients) - X[1, i]), 0.4)
  
  m <- mixture(Y = Y, experts = X[, c(1, K)], model = "MLprod", loss.type = l)
  expect_equal(m$loss, mean(loss(m$prediction, Y, loss.type = l)))
  expect_less_than(abs(sum(X[1, c(1, K)] * m$coefficients) - X[1, i]), 0.8)
  
  # expect_equal(c(predict(m,e[c(1,K)])), sum(c(e[c(1,K)])*c(m$coefficients)))
  
  m <- mixture(Y = Y, experts = X, model = "MLewa", loss.type = l, loss.gradient = FALSE)
  expect_equal(m$loss, mean(loss(m$prediction, Y, loss.type = l)))
  expect_less_than(abs(sum(X[1, ] * m$coefficients) - X[1, i]), 0.4)
  
  m <- mixture(Y = Y, experts = X[, c(1, K)], model = "MLewa", loss.type = l)
  expect_equal(m$loss, mean(loss(m$prediction, Y, loss.type = l)))
  expect_less_than(abs(sum(X[1, c(1, K)] * m$coefficients) - X[1, i]), 0.8)
  
  m <- mixture(Y = Y, experts = X, model = "BOA", loss.type = l, loss.gradient = FALSE)
  expect_equal(m$loss, mean(loss(m$prediction, Y, loss.type = l)))
  expect_less_than(abs(sum(X[1, ] * m$coefficients) - X[1, i]), 0.4)
  
  m <- mixture(Y = Y, experts = X[, c(1, K)], model = "BOA", loss.type = l)
  expect_equal(m$loss, mean(loss(m$prediction, Y, loss.type = l)))
  expect_less_than(abs(sum(X[1, c(1, K)] * m$coefficients) - X[1, i]), 0.8)
}) 
