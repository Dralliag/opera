# Unit tests of opera package using testhat package

context("Testing oracle function")

# load some basic data to perform tests
n <- 50
X <- cbind(rep(0, n), rep(1, n))
Y <- rep(0.4, n)
X[n, ] <- c(1, 1)
colnames(X) <- c("Exp1","Exp2")
Y[n] <- 1
awake <- cbind(rep(c(0, 1), n/2), 1)


# Test of loss functions
test_that("loss functions return correct values", {
  expect_that(dim(loss(X, Y, loss.type = "square"))[1], equals(n))
})

# Test of oracle functions
test_that("Best expert oracle is ok", {
  m <- oracle(Y = Y, experts = X, model = "expert")
  expect_that(as.numeric(m$coefficients[1]), equals(1))
  expect_that(m$loss, equals(mean((X[, 1] - Y)^2)))
  expect_that(sum(m$prediction), equals(sum(X[, 1])))
  expect_that(m$rmse, equals(sqrt(mean((X[, 1]- Y)^2))))
  
  expect_warning(oracle(Y = Y, experts = X, model = "expert", awake = awake), "When experts are unactive")
  expect_warning(oracle(Y = Y, experts = X, model = "expert", lambda = 3), "Unused lambda parameter")
  expect_warning(oracle(Y = Y, experts = X, model = "expert", niter = 3), "Unused niter parameter")
})

test_that("Best convex oracle is ok", {
  m <- oracle(Y = Y, experts = X, model = "convex")
  expect_equal(m$coefficients[1], 0.6)
  expect_equal(sum(m$coefficients), 1)
  expect_equal(m$loss, 0)
  expect_true(sum(abs(m$prediction - Y)) < 1e-10)
  expect_equal(m$rmse, 0)
  
  expect_warning(m <- oracle(Y = Y, experts = X, model = "convex", loss.type = "percentage"))
  expect_true(abs(m$coefficients[1] - 0.6) < 1e-04)
  expect_true(m$loss < 1e-04)
  expect_true(sum(abs(m$prediction - Y)) < 1e-04)
  
  expect_warning(m <- oracle(Y = Y, experts = X, model = "convex", loss.type = "absolute", awake = awake))
  expect_true(abs(m$coefficients[1] - 0.6) < 0.1)
  l <- getAnywhere(lossConv)$objs[[1]]
  expect_equal(mean(loss(m$prediction, Y, "absolute")), l(m$coefficients, Y, X, 
    awake, "absolute"))
  expect_equal(m$loss, mean(loss(m$prediction, Y, "absolute")))
})

test_that("Best linear oracle is ok", {
  m <- oracle(Y = Y, experts = X, model = "linear")
  expect_equal(m$coefficients[1], 0.6)
  expect_equal(sum(m$coefficients), 1)
  expect_equal(m$loss, 0)
  expect_true(sum(abs(m$prediction - Y)) < 1e-10)
  expect_equal(m$rmse, 0)
  expect_error(oracle(Y = Y, experts = X, model = "linear", awake = awake), "Sleeping or missing values not allowed")
  
  expect_warning(m <- oracle(Y = Y, experts = X, model = "linear", loss.type = "percentage"))
  expect_equal(m$loss, mean(loss(m$prediction, Y, loss.type = "percentage")))
})

test_that("Quantile oracles are ok", {
  set.seed(1)
  
  # test of quantile oracles
  quantiles <- seq(0.1, 0.9, by = 0.1)
  K <- length(quantiles)
  Y <- rnorm(n, mean = 0, sd = 1)
  X <- t(matrix(rep(quantile(Y, probs = quantiles), n), nrow = K))
  i <- sample(1:K, 1)
  
  l <- list(name = "pinball", tau = quantiles[i])
  # best expert oracle
  m.best_expert <- oracle(Y = Y, experts = X, model = "expert", loss.type = l)
  expect_equal(which(m.best_expert$coefficients == 1), i)
  expect_equal(m.best_expert$loss, mean(loss(m.best_expert$prediction, Y, loss.type = l)))
  
  # best convex oracle
  expect_warning(m <- oracle(Y = Y, experts = X[, c(1, K)], model = "convex", loss.type = l))
  expect_lt(abs(sum(X[1, c(1, K)] * m$coefficients) - X[1, i]), 0.1)
  expect_equal(m$loss, mean(loss(m$prediction, Y, loss.type = l)))
  expect_warning(oracle(Y = Y, experts = X[, c(1, K)], model = "convex", loss.type = l))
  
  # best linear oracle (with singular matrix)
  expect_warning(m <- oracle(Y = Y, experts = X[, c(1, K)], model = "linear", loss.type = l, niter = 10))
  expect_lt(abs(sum(X[1, c(1, K)] * m$coefficients) - X[1, i]), 0.1)
  expect_equal(m$loss, mean(loss(m$prediction, Y, loss.type = l)))
  expect_warning(oracle(Y = Y, experts = X[, c(1, K)], model = "linear", loss.type = l))
  
  
  # best linear oracle (with direct computation using rq)
  X[n, ] <- 1
  Y[n] <- 1
  m <- oracle(Y = Y, experts = X[, c(1, K)], model = "linear", loss.type = l)
  expect_lt(abs(sum(X[1, c(1, K)] * m$coefficients) - X[1, i]), 0.1)
  expect_equal(m$loss, mean(loss(m$prediction, Y, loss.type = l)))
})

test_that("Best shifting oracle is ok", {
  m <- oracle(Y = Y, experts = X, model = "shifting", loss.type = "square")
  expect_equal(m$loss[1], min(mean(loss(X[, 1], Y)), mean(loss(X[, 2], Y))))
  expect_equal(class(m), "oracle")
  expect_equal(class(summary(m)), "summary.oracle")
}) 

# test multi-dimensional data

test_that("Dimension d>1 is ok",{
  set.seed(1)
  # load some basic data to perform tests
  n <- 10
  d <- 3
  for (model in c("expert", "convex", "linear")) {
    for (l in c("square", "pinball", "percentage", "absolute")) {
    
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
    
    cat(model, l, "\n")
    m <- oracle(Y = Y,experts = X, model = model, loss.type = l)
    m$d <- d
    m$prediction <- seriesToBlock(m$prediction,d)
    m$Y <- seriesToBlock(m$Y,d)
    m$residuals <- seriesToBlock(m$residuals,d)
    m$experts <- seriesToBlock(m$experts,d)
    summary(m)
    plot(m)
    
    X <- seriesToBlock(X, d = d)
    Y <- seriesToBlock(Y, d = d)
    m1 <- oracle(Y = Y, experts= X, model = model, loss.type = l)
    expect_equal(m$experts,m1$experts)
    expect_true(mean(abs(m$prediction - m1$prediction)) < mean(abs(Y))/10)
    }
  }
})

