# Unit tests of opera package using testhat package

context("Testing mixture function")

# load some basic data to perform tests
set.seed(3)
n <- 50
X <- cbind(rep(0, n), rep(1, n))
colnames(X) <- c("Exp1","Exp2")
Y <- rep(0.4, n)
X[n, ] <- c(1, 1)
Y[n] <- 1
awake <- cbind(rep(c(0, 1), n/2), 1)

# Test warnings and errors
test_that("Errors are explained", {
  expect_error(mixture(loss.type = "plop"), "loss.type")
  expect_error(oracle(Y, X, loss.type = "plop"), "loss.type")
  
  expect_error(mixture(Y = Y), "Bad dimensions")
  expect_error(mixture(experts = X), "Bad dimensions")
  expect_error(mixture(Y = Y[1:2], experts = X), "Bad dimensions")
  expect_error(oracle(Y = Y[1:2], experts = X), "Bad dimensions")
  
})

# Test of EWA
test_that("EWA is ok", {
  w0 <- c(0.3, 0.7)
  possible_loss_type <- c("percentage", "absolute", "square", "pinball")
  i.loss <- sample(1:4, 1)
  
  expect_silent(
    m <- mixture(Y = Y, experts = X, model = "EWA", loss.type = possible_loss_type[i.loss],
                                        coefficients = w0)
    )

  expect_true(abs(m$coefficients[1] - 0.6) < 0.1)
  expect_equal(m$loss, mean(loss(m$prediction, Y, loss.type = possible_loss_type[i.loss])))
  expect_identical(as.numeric(m$weights[1, ]), w0)


  expect_failure(expect_warning(
    m <- mixture(Y = Y, experts = X, model = "EWA", loss.type = possible_loss_type[i.loss],
                                        coefficients = w0)))


  eta <- 0.5
  expect_failure(expect_warning(
    m.fixed <- mixture(Y = Y, experts = X, model = "EWA", parameters = list(eta = eta),
                                        loss.type = possible_loss_type[i.loss], coefficients = w0)))
  idx.eta <- which(m$parameters$grid.eta == eta)
  expect_equal(m$training$grid.loss[idx.eta], mean(loss(m.fixed$prediction, Y,
                                                        loss.type = possible_loss_type[i.loss])))
  expect_equal(m.fixed$loss, m$training$grid.loss[idx.eta])
  expect_identical(as.numeric(m.fixed$weights[1, ]), w0)

  expect_failure(expect_warning(
    m <- mixture(Y = Y, experts = X, model = "EWA", awake = awake)))
  expect_true(abs(m$coefficients[1] - 0.6) < 0.1)
  expect_equal(m$loss, mean(loss(m$prediction, Y)))
  # e <- c(0.3,0.5) expect_equal(c(predict(m,e)), sum(c(e)*c(m$coefficients)))

  grid.eta <- runif(5)
  expect_failure(expect_warning(
    m <- mixture(Y = Y, experts = X, model = "EWA", parameters = list(grid.eta = grid.eta),
               awake = awake)))
  expect_equal(sum(!(grid.eta %in% m$parameters$grid.eta)), 0)

  expect_failure(expect_warning(
    m1 <- mixture(Y = Y[1:10], experts = X[1:10, ], model = "EWA", parameters = list(grid.eta = grid.eta),
                awake = awake[1:10, ])))
  m1 <- predict(object = m1, newexperts = X[-c(1:10), ], newY = Y[-c(1:10)], awake = awake[-c(1:10),], online = TRUE, type = "model")
  expect_equal(m, m1)

})

# Test of Fixed-share
test_that("Fixed-share is ok", {
  w0 <- c(0.3, 0.7)
  possible_loss_type <- c("percentage", "absolute", "square", "pinball")
  i.loss <- sample(1:4, 1)
  m <- mixture(Y = Y, experts = X, model = "FS", loss.type = possible_loss_type[i.loss],
               coefficients = w0)
  expect_true(abs(m$coefficients[1] - 0.6) < 0.1)
  expect_equal(m$loss, mean(loss(m$prediction, Y, loss.type = possible_loss_type[i.loss])))
  expect_identical(as.numeric(m$weights[1, ]), w0)

  # e <- c(0.3,0.5) expect_equal(c(predict(m,e)), sum(c(e)*c(m$coefficients)))
  eta <- 2
  alpha <- 0.01
  m.fixed <- mixture(Y = Y, experts = X, model = "FS", parameters = list(eta = eta,
                                                                         alpha = alpha), loss.type = possible_loss_type[i.loss], coefficients = w0)
  idx.eta <- which(m$parameters$grid.eta == eta)
  idx.alpha <- which(m$parameters$grid.alpha == alpha)
  expect_equal(m$training$grid.loss[idx.eta, idx.alpha], mean(loss(m.fixed$prediction,
                                                                   Y, loss.type = possible_loss_type[i.loss])))
  expect_equal(m.fixed$loss, m$training$grid.loss[idx.eta, idx.alpha])
  expect_identical(as.numeric(m.fixed$weights[1, ]), w0)

  m <- mixture(Y = Y, experts = X, model = "FS", awake = awake)
  expect_true(abs(m$coefficients[1] - 0.6) < 0.1)
  expect_equal(m$loss, mean(loss(m$prediction, Y)))

  grid.eta <- runif(5)
  grid.alpha <- runif(3)
  m <- mixture(Y = Y, experts = X, model = "FS", parameters = list(grid.eta = grid.eta,
                                                                   grid.alpha = grid.alpha), awake = awake)
  expect_equal(sum(!(grid.eta %in% m$parameters$grid.eta)), 0)
  expect_identical(grid.alpha, m$parameters$grid.alpha)
})

# Test of Ridge
test_that("Ridge is ok", {
  w0 <- c(0.5, 0.5)
  m <- mixture(Y = Y, experts = X, model = "Ridge", coefficients = w0)
  expect_equal(m$loss, mean(loss(m$prediction, Y)))
  expect_identical(as.numeric(m$weights[1, ]), w0)
  expect_true(!is.na(sum(m$weights)))

  # e <- c(0.3,0.5) expect_equal(c(predict(m,e)), sum(c(e)*c(m$coefficients)))

  lambda <- 2
  m.fixed <- mixture(Y = Y, experts = X, model = "Ridge", parameters = list(lambda = lambda),
                     coefficients = w0)
  idx.lambda <- which(m$parameters$grid.lambda == lambda)
  expect_equal(m$training$grid.loss[idx.lambda], mean(loss(m.fixed$prediction,
                                                           Y)))
  expect_equal(m.fixed$loss, m$training$grid.loss[idx.lambda])
  expect_identical(as.numeric(m.fixed$weights[1, ]), w0)

  grid.lambda <- runif(3)
  m <- mixture(Y = Y, experts = X, model = "Ridge", parameters = list(grid.lambda = grid.lambda,
                                                                      gamma = 100))
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
  m1 <- predict(object = m1, newexperts = X[-c(1:10), ], newY = Y[-c(1:10)], online = TRUE,
                type = "model")
  expect_equal(m, m1)

  w0 <- c(0.3, 0.7)
  m <- mixture(Y = Y, experts = X, model = "MLprod", coefficients = w0, loss.type = possible_loss_type[i.loss])
  expect_true(abs(m$coefficients[1] - 0.6) < 0.2)
  expect_equal(m$loss, mean(loss(m$prediction, Y, loss.type = possible_loss_type[i.loss])))
  expect_identical(as.numeric(m$weights[1, ]), w0)

  m <- mixture(Y = Y, experts = X, model = "MLewa", coefficients = w0, loss.type = possible_loss_type[i.loss])
  expect_true(abs(m$coefficients[1] - 0.6) < 0.2)
  expect_equal(m$loss, mean(loss(m$prediction, Y, loss.type = possible_loss_type[i.loss])))
  expect_identical(as.numeric(m$weights[1, ]), w0)
  e <- c(0.3, 0.5)
  # expect_equal(c(predict(m,e)), sum(c(e)*c(m$coefficients))) à vérifier plus tard

  w0 <- c(0.3, 0.7)
  m <- mixture(Y = Y, experts = X, model = "BOA", coefficients = w0, loss.type = possible_loss_type[i.loss])
  expect_true(abs(m$coefficients[1] - 0.6) < 0.2)
  expect_equal(m$loss, mean(loss(m$prediction, Y, loss.type = possible_loss_type[i.loss])))
  expect_identical(as.numeric(m$weights[1, ]), w0)


  m <- mixture(Y = Y[1:5], experts = X[1:5, ], model = "MLewa", awake = awake[1:5,
                                                                              ])
  m <- predict(m, newexperts = X[-c(1:5), ], newY = Y[-c(1:5)], awake = awake[-c(1:5),
                                                                              ])
  expect_true(abs(m$coefficients[1] - 0.6) < 0.2)
  expect_equal(m$loss, mean(loss(m$prediction, Y)))

  m1 <- mixture(Y = Y, experts = X, model = "MLewa", awake = awake)
  expect_equal(m, m1)
})

# test of quantile prediction
test_that("Quantile mixture are ok", {
  # test de la partie quantile
  n <- 200
  quantiles <- seq(0.1, 0.9, by = 0.1)
  K <- length(quantiles)
  Y <- rnorm(n, mean = 0, sd = 1)
  X <- t(matrix(rep(quantile(Y, probs = quantiles), n), nrow = K))
  i <- sample(1:K, 1)
  l <- list(name = "pinball", tau = quantiles[i])
  m <- mixture(Y = Y, experts = X, model = "EWA", loss.type = l, loss.gradient = FALSE,
               parameters = list(eta = 1, gamma = 100))
  expect_equal(m$loss, mean(loss(m$prediction, Y, loss.type = l)))
  expect_lt(abs(sum(X[1, ] * m$coefficients) - X[1, i]), 0.4)
  # e <- rnorm(K) expect_equal(c(predict(m,e)), sum(c(e)*c(m$coefficients)))

  m <- mixture(Y = Y, experts = X[, c(1, K)], model = "EWA", loss.type = l, parameters = list(gamma = 100))
  expect_equal(m$loss, mean(loss(m$prediction, Y, loss.type = l)))
  expect_lt(abs(sum(X[1, c(1, K)] * m$coefficients) - X[1, i]), 0.4)

  # Fixed share
  m <- mixture(Y = Y, experts = X, model = "FS", loss.type = l, loss.gradient = FALSE,
               parameters = list(eta = 1, alpha = 0.01))
  expect_equal(m$loss, mean(loss(m$prediction, Y, loss.type = l)))
  expect_lt(abs(sum(X[1, ] * m$coefficients) - X[1, i]), 0.4)

  m <- mixture(Y = Y, experts = X[, c(1, K)], model = "FS", loss.type = l, parameters = list(gamma = 10))
  expect_equal(m$loss, mean(loss(m$prediction, Y, loss.type = l)))
  expect_lt(abs(sum(X[1, c(1, K)] * m$coefficients) - X[1, i]), 0.4)

  m <- mixture(Y = Y, experts = X, model = "MLpol", loss.type = l, loss.gradient = FALSE)
  expect_equal(m$loss, mean(loss(m$prediction, Y, loss.type = l)))
  expect_lt(abs(sum(X[1, ] * m$coefficients) - X[1, i]), 0.4)

  m <- mixture(Y = Y, experts = X[, c(1, K)], model = "MLpol", loss.type = l)
  expect_equal(m$loss, mean(loss(m$prediction, Y, loss.type = l)))
  expect_lt(abs(sum(X[1, c(1, K)] * m$coefficients) - X[1, i]), 0.8)

  m <- mixture(Y = Y, experts = X, model = "MLprod", loss.type = l, loss.gradient = FALSE)
  expect_equal(m$loss, mean(loss(m$prediction, Y, loss.type = l)))
  expect_lt(abs(sum(X[1, ] * m$coefficients) - X[1, i]), 0.4)

  m <- mixture(Y = Y, experts = X[, c(1, K)], model = "MLprod", loss.type = l)
  expect_equal(m$loss, mean(loss(m$prediction, Y, loss.type = l)))
  expect_lt(abs(sum(X[1, c(1, K)] * m$coefficients) - X[1, i]), 0.8)

  # expect_equal(c(predict(m,e[c(1,K)])), sum(c(e[c(1,K)])*c(m$coefficients)))

  m <- mixture(Y = Y, experts = X, model = "MLewa", loss.type = l, loss.gradient = FALSE)
  expect_equal(m$loss, mean(loss(m$prediction, Y, loss.type = l)))
  expect_lt(abs(sum(X[1, ] * m$coefficients) - X[1, i]), 0.4)

  m <- mixture(Y = Y, experts = X[, c(1, K)], model = "MLewa", loss.type = l)
  expect_equal(m$loss, mean(loss(m$prediction, Y, loss.type = l)))
  expect_lt(abs(sum(X[1, c(1, K)] * m$coefficients) - X[1, i]), 0.8)

  m <- mixture(Y = Y, experts = X, model = "BOA", loss.type = l, loss.gradient = FALSE)
  expect_equal(m$loss, mean(loss(m$prediction, Y, loss.type = l)))
  expect_lt(abs(sum(X[1, ] * m$coefficients) - X[1, i]), 0.4)

  m <- mixture(Y = Y, experts = X[, c(1, K)], model = "BOA", loss.type = l)
  expect_equal(m$loss, mean(loss(m$prediction, Y, loss.type = l)))
  expect_lt(abs(sum(X[1, c(1, K)] * m$coefficients) - X[1, i]), 0.8)
})

# test of predict function
test_that("Predict method is ok", {
  for (model in c("BOA", "MLpol", "MLprod", "MLewa", "FS", "Ridge")) {
    l <- sample(c("square", "pinball", "percentage", "absolute"), 1)
    if (model == "Ridge") {
      l <- "square"
      awake <- NULL
    }
    m <- mixture(model = model, loss.type = l)
    expect_warning(predict(m))

    # single online prediction and sequential prediction return similar models
    m1 <- predict(m, newY = Y, newexperts = X, type = "m", awake = awake)
    m2 <- m
    for (t in 1:n) {
      m2 <- predict(m2, newY = Y[t], newexperts = X[t, ], online = TRUE, type = "model",
                    awake = awake[t, ])
    }
    expect_equal(m1, m2)

    # batch prediction is ok
    m1 <- predict(m, newY = Y, newexperts = X, type = "m", online = FALSE, awake = awake)
    expect_equal(m1$coefficients, m2$coefficients)

    expect_warning(m2 <- predict(m, newexperts = X, type = "r", online = FALSE, awake = awake))
    expect_equal(m1$prediction, m2)
  }
})

# test that regret and cumulative loss of the expert are coherent
test_that("Regrets and Losses are coherent", {
  for (model in c("BOA", "MLpol", "MLprod", "MLewa", "FS", "EWA")) {
    l <- sample(c("square", "pinball", "percentage", "absolute"), 1)

    m <- mixture(Y = Y, experts = X, model = "EWA", loss.type = l, loss.gradient = FALSE,parameters = list(eta = 1, alpha = 0.1))
    o <- oracle(Y = Y, experts = X, model = "expert", loss.type = l)
    l1 <- m$loss*n - m$training$R
    l2 <- o$loss.experts*n
    l3 = apply(loss(X,Y,loss.type=l),2,sum)
    expect_equal(as.numeric(l1),as.numeric(l2))
    expect_equal(l2,l3)
  }
})

# test multi-dimensional data
test_that("Dimension d>1 is ok",{
  # load some basic data to perform tests
  n <- 10
  d <- 3
  for (algo in c("BOA", "MLpol", "MLprod", "MLewa", "FS", "Ridge","OGD")) {
    l <- sample(c("square", "pinball", "percentage", "absolute"), 1)
    if (algo == "Ridge") {
      l <- "square"
      awake <- NULL
    }
    # Une petite fonction pour creer les prévisions de la base canonique
    base_predictions = function(d,n) {
      decimals <- c(0:(2^d-1))
      m <- cbind(diag(d),-diag(d))
      return(t(matrix(rep(t(m),n),nrow = 2*d)))
    }
    X <- base_predictions(d,n) # X is the canonical basis
    theta.star <- sign(rnorm(d)) * runif(d) # point to be predicted
    theta.star <- runif(1) * theta.star / sum(abs(theta.star))  # the target point is in the L1 unit ball
    if (l == "percentage") {
      X <- abs(X)
      theta.star <- abs(theta.star)
    }
    Y <- rep(theta.star, n)

    m <- mixture(model = algo, loss.type = l)
    for (i in seq(1,n*d,by=d)) {
      idx = i + 0:(d-1)
      m <- predict(object = m, newY = Y[idx],newexperts = X[idx,], online = FALSE)
    }
    # theta.star - predict(m,newexperts = X[1:3,], online = FALSE, type = "response")
    m$d <- d
    m$T <- m$T / d
    m$prediction <- t(matrix(m$prediction,nrow = d))
    m$Y <- t(matrix(m$Y,nrow = d))
    m$weights <- m$weights[seq(1,n*d,by=d),]

    summary(m)
    plot(m)
    X1 <- seriesToBlock(X, d = d)
    Y1 <- seriesToBlock(Y, d = d)
    m1 <- mixture(Y = Y1, experts= X1, model = algo, loss.type = l)
    expect_equal(m,m1)
  }
})

test_that("Names are correctly rendered",{
  colnames(X) <- c("bob","alice")
  m <- mixture(Y = Y, experts=X)
  expect_equal(m$names.experts,colnames(X))

  m <- mixture()
  m <- predict(m,newY = Y[1],newexperts = X[1,])
  expect_equal(m$names.experts,colnames(X))
})

