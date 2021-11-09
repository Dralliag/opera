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
  for (possible_loss in c("percentage", "absolute", "square", "pinball")) {
    cur_loss <- list("name" = possible_loss)
    if (possible_loss == "pinball") {cur_loss$tau <- 0.5}
    
    expect_silent(
      m <- mixture(Y = Y, experts = X, model = "EWA", loss.type = cur_loss,
                   coefficients = w0, quiet = TRUE)
    )
    
    expect_true(abs(m$coefficients[1] - 0.6) < 0.1)
    expect_equal(m$loss, mean(loss(m$prediction, Y, loss.type = cur_loss)))
    expect_identical(as.numeric(m$weights[1, ]), w0)
    
    
    expect_failure(expect_warning(
      m <- mixture(Y = Y, experts = X, model = "EWA", loss.type = cur_loss,
                   coefficients = w0, parameters = list("grid.eta" = 1), quiet = TRUE)))
    
    
    eta <- 0.5
    expect_failure(expect_warning(
      m.fixed <- mixture(Y = Y, experts = X, model = "EWA", parameters = list(eta = eta),
                         loss.type = cur_loss, coefficients = w0, quiet = TRUE)))
    idx.eta <- which(m$parameters$grid.eta == eta)
    expect_equal(m$training$grid.loss[idx.eta], mean(loss(m.fixed$prediction, Y,
                                                              loss.type = cur_loss)))
    expect_equal(m.fixed$loss, m$training$grid.loss[idx.eta])
    expect_identical(as.numeric(m.fixed$weights[1, ]), w0)
    
    expect_failure(expect_warning(
      m <- mixture(Y = Y, experts = X, model = "EWA", awake = awake, quiet = TRUE)))
    expect_true(abs(m$coefficients[1] - 0.6) < 0.1)
    expect_equal(m$loss, mean(loss(m$prediction, Y)))
    # e <- c(0.3,0.5) expect_equal(c(predict(m,e)), sum(c(e)*c(m$coefficients)))
    
    grid.eta <- runif(5)
    expect_failure(expect_warning(
      m <- mixture(Y = Y, experts = X, model = "EWA", parameters = list(grid.eta = grid.eta),
                   awake = awake, quiet = TRUE)))
    expect_equal(sum(!(grid.eta %in% m$parameters$grid.eta)), 0)
    
    expect_failure(expect_warning(
      m1 <- mixture(Y = Y[1:10], experts = X[1:10, ], model = "EWA", parameters = list(grid.eta = grid.eta),
                    awake = awake[1:10, ], quiet = TRUE)))
      m1 <- predict(object = m1, newexperts = X[-c(1:10), ], newY = Y[-c(1:10)], awake = awake[-c(1:10),], online = TRUE, type = "model", quiet = TRUE)
    expect_equal(m, m1)
  }
})

# Test of Fixed-share
test_that("Fixed-share is ok", {
  w0 <- c(0.3, 0.7)
  for (possible_loss in c("percentage", "absolute", "square", "pinball")) {
    cur_loss <- list("name" = possible_loss)
    if (possible_loss == "pinball") {cur_loss$tau <- 0.5}
    
    possible_loss_type <- c("percentage", "absolute", "square", "pinball")
    i.loss <- sample(1:4, 1)
    m <- mixture(Y = Y, experts = X, model = "FS", loss.type = cur_loss,
                 coefficients = w0, parameters = list(grid.eta = 1), quiet = TRUE)
    expect_true(abs(m$coefficients[1] - 0.6) < 0.1)
    expect_equal(m$loss, mean(loss(m$prediction, Y, loss.type = cur_loss)))
    expect_identical(as.numeric(m$weights[1, ]), w0)
    
    # e <- c(0.3,0.5) expect_equal(c(predict(m,e)), sum(c(e)*c(m$coefficients)))
    eta <- 2
    alpha <- 0.01
    m.fixed <- mixture(Y = Y, experts = X, model = "FS", parameters = list(eta = eta, alpha = alpha), 
                       loss.type = cur_loss, coefficients = w0, quiet = TRUE)
    idx.eta <- which(m$parameters$grid.eta == eta)
    idx.alpha <- which(m$parameters$grid.alpha == alpha)
    expect_equal(m$training$grid.loss[idx.eta, idx.alpha], mean(loss(m.fixed$prediction,
                                                                         Y, loss.type = cur_loss)))
    expect_equal(m.fixed$loss, m$training$grid.loss[idx.eta, idx.alpha])
    expect_identical(as.numeric(m.fixed$weights[1, ]), w0)
    
    m <- mixture(Y = Y, experts = X, model = "FS", awake = awake, quiet = TRUE)
    expect_true(abs(m$coefficients[1] - 0.6) < 0.1)
    expect_equal(m$loss, mean(loss(m$prediction, Y)))
    
    grid.eta <- runif(5)
    grid.alpha <- runif(3)
    m <- mixture(Y = Y, experts = X, model = "FS", parameters = list(grid.eta = grid.eta,
                                                                     grid.alpha = grid.alpha), awake = awake, quiet = TRUE)
    expect_equal(sum(!(grid.eta %in% m$parameters$grid.eta)), 0)
    expect_identical(grid.alpha, m$parameters$grid.alpha)
  }
})

# Test of Ridge
test_that("Ridge is ok", {
  w0 <- c(0.5, 0.5)
  m <- mixture(Y = Y, experts = X, model = "Ridge", coefficients = w0, quiet = TRUE)
  expect_equal(m$loss, mean(loss(m$prediction, Y)))
  expect_identical(as.numeric(m$weights[1, ]), w0)
  expect_true(!is.na(sum(m$weights)))
  nlambda <- length(m$parameters$grid.lambda)
  lambda <- m$parameters$grid.lambda[nlambda]
  
  m.fixed <- mixture(Y = Y, experts = X, model = "Ridge", parameters = list(lambda = lambda),
                     coefficients = w0, quiet = TRUE)
  idx.lambda <- which(m$parameters$grid.lambda == lambda)
  expect_equal(m$training$grid.loss[idx.lambda], mean(loss(m.fixed$prediction,
                                                               Y)))
  expect_equal(m.fixed$loss, m$training$grid.loss[idx.lambda])
  expect_identical(as.numeric(m.fixed$weights[1, ]), w0)
  
  grid.lambda <- runif(3)
  m <- mixture(Y = Y, experts = X, model = "Ridge", parameters = list(grid.lambda = grid.lambda,
                                                                      gamma = 100), quiet = TRUE)
  expect_equal(sum(!(grid.lambda %in% m$parameters$grid.lambda)), 0)
})

# test of MLPol,...
test_that("MLpol, MLprod, MLewa, and BOA are ok", {
  for (possible_loss in c("percentage", "absolute", "square", "pinball")) {
    cur_loss <- list("name" = possible_loss)
    if (possible_loss == "pinball") {cur_loss$tau <- 0.5}
    
    m <- mixture(Y = Y, experts = X, loss.type = cur_loss, quiet = TRUE)
    expect_true(abs(m$coefficients[1] - 0.6) < 0.2)
    expect_equal(m$loss, mean(loss(m$prediction, Y, loss.type = cur_loss)))
    
    m1 <- mixture(loss.type = cur_loss, quiet = TRUE)
    m1 <- predict(object = m1, newexperts = X, newY = Y, online = TRUE, type = "model", quiet = TRUE)
    expect_equal(m, m1)
    
    m1 <- mixture(Y = Y[1:10], experts = X[1:10, ], loss.type = cur_loss, quiet = TRUE)
    m1 <- predict(object = m1, newexperts = X[-c(1:10), ], newY = Y[-c(1:10)], online = TRUE,
                  type = "model", quiet = TRUE)
    expect_equal(m, m1)
    
    w0 <- c(0.3, 0.7)
    m <- mixture(Y = Y, experts = X, model = "MLprod", coefficients = w0, loss.type = cur_loss, quiet = TRUE)
    expect_true(abs(m$coefficients[1] - 0.6) < 0.2)
    expect_equal(m$loss, mean(loss(m$prediction, Y, loss.type = cur_loss)))
    expect_equal(as.numeric(m$weights[1, ]), w0)
    
    m <- mixture(Y = Y, experts = X, model = "MLewa", coefficients = w0, loss.type = cur_loss, quiet = TRUE)
    expect_true(abs(m$coefficients[1] - 0.6) < 0.2)
    expect_equal(m$loss, mean(loss(m$prediction, Y, loss.type = cur_loss)))
    expect_identical(as.numeric(m$weights[1, ]), w0)
    e <- c(0.3, 0.5)
    # expect_equal(c(predict(m,e)), sum(c(e)*c(m$coefficients))) à vérifier plus tard
    
    w0 <- c(0.3, 0.7)
    m <- mixture(Y = Y, experts = X, model = "BOA", coefficients = w0, loss.type = cur_loss, quiet = TRUE)
    expect_true(abs(m$coefficients[1] - 0.6) < 0.2)
    expect_equal(m$loss, mean(loss(m$prediction, Y, loss.type = cur_loss)))
    expect_identical(as.numeric(m$weights[1, ]), w0)
    
    
    m <- mixture(Y = Y[1:5], experts = X[1:5, ], model = "MLewa", awake = awake[1:5,
    ], quiet = TRUE)
    m <- predict(m, newexperts = X[-c(1:5), ], newY = Y[-c(1:5)], awake = awake[-c(1:5),
    ], quiet = TRUE)
    expect_true(abs(m$coefficients[1] - 0.6) < 0.2)
    expect_equal(m$loss, mean(loss(m$prediction, Y)))
    
    m1 <- mixture(Y = Y, experts = X, model = "MLewa", awake = awake, quiet = TRUE)
    expect_equal(m, m1)
    
  }
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
  l <- list("name" = "pinball", "tau" = quantiles[i])
  m <- mixture(Y = Y, experts = X, model = "EWA", loss.type = l, loss.gradient = FALSE,
               parameters = list(eta = 1, gamma = 100), quiet = TRUE)
  expect_equal(m$loss, mean(loss(m$prediction, Y, loss.type = l)))
  expect_lt(abs(sum(X[1, ] * m$coefficients) - X[1, i]), 0.4)
  
  m <- mixture(Y = Y, experts = X[, c(1, K)], model = "EWA", loss.type = l, parameters = list(gamma = 100), quiet = TRUE)
  expect_equal(m$loss, mean(loss(m$prediction, Y, loss.type = l)))
  expect_lt(abs(sum(X[1, c(1, K)] * m$coefficients) - X[1, i]), 0.4)
  
  # Fixed share
  m <- mixture(Y = Y, experts = X, model = "FS", loss.type = l, loss.gradient = FALSE,
               parameters = list(eta = 1, alpha = 0.01), quiet = TRUE)
  expect_equal(m$loss, mean(loss(m$prediction, Y, loss.type = l)))
  expect_lt(abs(sum(X[1, ] * m$coefficients) - X[1, i]), 0.4)
  
  m <- mixture(Y = Y, experts = X[, c(1, K)], model = "FS", loss.type = l, parameters = list(gamma = 10), quiet = TRUE)
  expect_equal(m$loss, mean(loss(m$prediction, Y, loss.type = l)))
  expect_lt(abs(sum(X[1, c(1, K)] * m$coefficients) - X[1, i]), 0.4)
  
  m <- mixture(Y = Y, experts = X, model = "MLpol", loss.type = l, loss.gradient = FALSE, quiet = TRUE)
  expect_equal(m$loss, mean(loss(m$prediction, Y, loss.type = l)))
  expect_lt(abs(sum(X[1, ] * m$coefficients) - X[1, i]), 0.4)
  
  m <- mixture(Y = Y, experts = X[, c(1, K)], model = "MLpol", loss.type = l, quiet = TRUE)
  expect_equal(m$loss, mean(loss(m$prediction, Y, loss.type = l)))
  expect_lt(abs(sum(X[1, c(1, K)] * m$coefficients) - X[1, i]), 0.8)
  
  m <- mixture(Y = Y, experts = X, model = "MLprod", loss.type = l, loss.gradient = FALSE, quiet = TRUE)
  expect_equal(m$loss, mean(loss(m$prediction, Y, loss.type = l)))
  expect_lt(abs(sum(X[1, ] * m$coefficients) - X[1, i]), 0.4)
  
  m <- mixture(Y = Y, experts = X[, c(1, K)], model = "MLprod", loss.type = l, quiet = TRUE)
  expect_equal(m$loss, mean(loss(m$prediction, Y, loss.type = l)))
  expect_lt(abs(sum(X[1, c(1, K)] * m$coefficients) - X[1, i]), 0.8)
  
  m <- mixture(Y = Y, experts = X, model = "MLewa", loss.type = l, loss.gradient = FALSE, quiet = TRUE)
  expect_equal(m$loss, mean(loss(m$prediction, Y, loss.type = l)))
  expect_lt(abs(sum(X[1, ] * m$coefficients) - X[1, i]), 0.4)
  
  m <- mixture(Y = Y, experts = X[, c(1, K)], model = "MLewa", loss.type = l, quiet = TRUE)
  expect_equal(m$loss, mean(loss(m$prediction, Y, loss.type = l)))
  expect_lt(abs(sum(X[1, c(1, K)] * m$coefficients) - X[1, i]), 0.8)
  
  m <- mixture(Y = Y, experts = X, model = "BOA", loss.type = l, loss.gradient = FALSE, quiet = TRUE)
  expect_equal(m$loss, mean(loss(m$prediction, Y, loss.type = l)))
  expect_lt(abs(sum(X[1, ] * m$coefficients) - X[1, i]), 0.4)
  
  m <- mixture(Y = Y, experts = X[, c(1, K)], model = "BOA", loss.type = l, quiet = TRUE)
  expect_equal(m$loss, mean(loss(m$prediction, Y, loss.type = l)))
  expect_lt(abs(sum(X[1, c(1, K)] * m$coefficients) - X[1, i]), 0.8)
})

# test of predict function
test_that("Predict method is ok, with and without awake, use_cpp or not", {
  for (model in c("MLpol", "MLprod", "MLewa", "FS", "Ridge", "BOA", "EWA")) {
    for (possible_loss in c("percentage", "absolute", "square", "pinball")) {
      cur_loss <- list("name" = possible_loss)
      
      if (model == "Ridge") {
        cur_loss <- list("name" = "square")
        awake <- NULL
      }
      m <- mixture(model = model, loss.type = cur_loss, quiet = TRUE)
      expect_warning(predict(m, quiet = TRUE))
      
      # with awake
      # single online prediction and sequential prediction return similar models
      m1_cpp <- predict(m, newY = Y, newexperts = X, type = "m", awake = awake, quiet = TRUE, use_cpp = TRUE)
      m2_cpp <- m
      for (t in 1:n) {
        m2_cpp <- predict(m2_cpp, newY = Y[t], newexperts = X[t, ], online = TRUE, type = "model",
                      awake = awake[t, ], quiet = TRUE, use_cpp = TRUE)
      }
      expect_equal(m1_cpp, m2_cpp)
      
      m1_r <- predict(m, newY = Y, newexperts = X, type = "m", awake = awake, quiet = TRUE, use_cpp = FALSE)
      m2_r <- m
      for (t in 1:n) {
        m2_r <- predict(m2_r, newY = Y[t], newexperts = X[t, ], online = TRUE, type = "model",
                      awake = awake[t, ], quiet = TRUE, use_cpp = FALSE)
      }
      expect_equal(m1_r, m2_r)
      
      expect_equal(m1_r, m1_cpp) 
      
      # batch prediction is ok
      m1 <- predict(m, newY = Y, newexperts = X, type = "m", online = FALSE, awake = awake, quiet = TRUE, use_cpp = TRUE)
      expect_equal(m1$coefficients, m2_cpp$coefficients)
      
      expect_warning(m2 <- predict(m, newexperts = X, type = "r", online = FALSE, awake = awake, quiet = TRUE, use_cpp = TRUE))
      expect_equal(m1$prediction, m2)
      
      
      # without awake
      # single online prediction and sequential prediction return similar models
      m1_cpp <- predict(m, newY = Y, newexperts = X, type = "m", awake = NULL, quiet = TRUE, use_cpp = TRUE)
      m2_cpp <- m
      for (t in 1:n) {
        m2_cpp <- predict(m2_cpp, newY = Y[t], newexperts = X[t, ], online = TRUE, type = "model",
                          awake = NULL, quiet = TRUE, use_cpp = TRUE)
      }
      expect_equal(m1_cpp, m2_cpp)
      
      m1_r <- predict(m, newY = Y, newexperts = X, type = "m", awake = NULL, quiet = TRUE, use_cpp = FALSE)
      m2_r <- m
      for (t in 1:n) {
        m2_r <- predict(m2_r, newY = Y[t], newexperts = X[t, ], online = TRUE, type = "model",
                        awake = NULL, quiet = TRUE, use_cpp = FALSE)
      }
      expect_equal(m1_r, m2_r)
      
      expect_equal(m1_r, m1_cpp) 
      
      # batch prediction is ok
      m1 <- predict(m, newY = Y, newexperts = X, type = "m", online = FALSE, awake = NULL, quiet = TRUE, use_cpp = TRUE)
      expect_equal(m1$coefficients, m2_cpp$coefficients)
      
      expect_warning(m2 <- predict(m, newexperts = X, type = "r", online = FALSE, awake = NULL, quiet = TRUE, use_cpp = TRUE))
      expect_equal(m1$prediction, m2)
    }
  }
})

# test of predict function on FTRL
test_that("Predict FTRL is ok, use_cpp or not", {
  model <- "FTRL"
  
  for (possible_loss in c("percentage", "absolute", "square", "pinball")) {
    cur_loss <- list("name" = possible_loss)
    # if (possible_loss == "pinball") {cur_loss$tau <- 0.5}
    if (model == "Ridge") {
      cur_loss <- list("name" = "square")
      awake <- NULL
    }
    m <- mixture(model = model, loss.type = cur_loss, quiet = TRUE)
    expect_warning(predict(m, quiet = TRUE))
    
    # without awake
    # single online prediction and sequential prediction return similar models
    m1_cpp <- predict(m, newY = Y, newexperts = X, type = "m", quiet = TRUE, use_cpp = TRUE)
    m2_cpp <- m
    for (t in 1:n) {
      m2_cpp <- predict(m2_cpp, newY = Y[t], newexperts = X[t, ], online = TRUE, type = "model",
                        quiet = TRUE, use_cpp = TRUE)
    }
    expect_equal(m1_cpp[! which(names(m1_cpp) %in% c("training", "parameters"))], m2_cpp[! which(names(m2_cpp) %in% c("training", "parameters"))])
    expect_equal(m1_cpp[[which(names(m1_cpp) == "training")]][c(1:2, 9:15)], m2_cpp[[which(names(m2_cpp) == "training")]][c(1:2, 9:15)])
    expect_equal(m1_cpp[[which(names(m1_cpp) == "parameters")]][c(1:2, 5, 7)], m2_cpp[[which(names(m2_cpp) == "parameters")]][c(1:2, 5, 7)])
    # /!\ environment is kept with function declaration
    
    m1_r <- predict(m, newY = Y, newexperts = X, type = "m", quiet = TRUE, use_cpp = FALSE)
    m2_r <- m
    for (t in 1:n) {
      m2_r <- predict(m2_r, newY = Y[t], newexperts = X[t, ], online = TRUE, type = "model",
                      quiet = TRUE, use_cpp = FALSE)
    }
    expect_equal(m1_r[! which(names(m1_r) %in% c("training", "parameters"))], m2_r[! which(names(m2_r) %in% c("training", "parameters"))])
    expect_equal(m1_r[[which(names(m1_r) == "training")]][c(1:2, 9:15)], m2_r[[which(names(m2_r) == "training")]][c(1:2, 9:15)])
    expect_equal(m1_r[[which(names(m1_r) == "parameters")]][c(1:2, 5, 7)], m2_r[[which(names(m2_r) == "parameters")]][c(1:2, 5, 7)])
    # /!\ environment is kept with function declaration
    
    expect_equal(m1_r[! which(names(m1_r) %in% c("training", "parameters"))], m1_cpp[! which(names(m1_cpp) %in% c("training", "parameters"))])
    expect_equal(m1_r[[which(names(m1_r) == "training")]][c(1:2, 9:15)], m1_cpp[[which(names(m1_cpp) == "training")]][c(1:2, 9:15)])
    expect_equal(m1_r[[which(names(m1_r) == "parameters")]][c(1:2, 5, 7)], m1_cpp[[which(names(m1_cpp) == "parameters")]][c(1:2, 5, 7)])
    # /!\ environment is kept with function declaration
    
    # batch prediction is ok
    m1 <- predict(m, newY = Y, newexperts = X, type = "m", online = FALSE, quiet = TRUE, use_cpp = TRUE)
    expect_equal(m1$coefficients, m2_cpp$coefficients)
    
    expect_warning(m2 <- predict(m, newexperts = X, type = "r", online = FALSE, quiet = TRUE, use_cpp = TRUE))
    expect_equal(m1$prediction, m2)
  }
})

# test that regret and cumulative loss of the expert are coherent
test_that("Regrets and Losses are coherent", {
  for (model in c("BOA", "EWA", "MLpol", "MLprod", "MLewa", "FS", "FTRL")) {
    for (possible_loss in c("percentage", "absolute", "square", "pinball")) {
      cur_loss <- list("name" = possible_loss)
      if (possible_loss == "pinball") {cur_loss$tau <- 0.5}
      
      m <- mixture(Y = Y, experts = X, model = "EWA", loss.type = cur_loss, loss.gradient = FALSE,parameters = list(eta = 1, alpha = 0.1), quiet = TRUE)
      o <- oracle(Y = Y, experts = X, model = "expert", loss.type = cur_loss)
      l1 <- m$loss*n - m$training$R
      l2 <- o$loss.experts*n
      l3 = colSums(apply(X, 2, function(x) loss(x,Y,loss.type=cur_loss)))
      expect_equal(as.numeric(l1),as.numeric(l2))
      expect_equal(l2,l3)
    }
  }
})

# test multi-dimensional data
test_that("Dimension d>1 is ok",{
  # load some basic data to perform tests
  n <- 10
  d <- 3
  for (algo in c("BOA", "EWA", "MLpol", "MLprod", "MLewa", "FS", "Ridge","OGD", "FTRL")) {
    for (possible_loss in c("percentage", "absolute", "square", "pinball")) {
      cur_loss <- list("name" = possible_loss)
      if (possible_loss == "pinball") {cur_loss$tau <- 0.5}
      if (algo == "Ridge") {
        cur_loss <- list("name" = "square")
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
      if (cur_loss$name == "percentage") {
        X <- abs(X)
        theta.star <- abs(theta.star)
      }
      Y <- rep(theta.star, n)
      
      m <- mixture(model = algo, loss.type = cur_loss, quiet = TRUE)
      for (i in seq(1,n*d,by=d)) {
        idx = i + 0:(d-1)
        m <- predict(object = m, newY = Y[idx], newexperts = X[idx,], online = FALSE, quiet = TRUE)
      }
      
      expect_output(summary(m), NA)
      expect_error(plot(m), NA)
      expect_output(print(m))
      
      m$d <- d
      m$T <- m$T / d
      m$prediction <- t(matrix(m$prediction,nrow = d))
      m$Y <- t(matrix(m$Y,nrow = d))
      m$weights <- m$weights[seq(1,n*d,by=d),]
      
      X1 <- seriesToBlock(X, d = d)
      Y1 <- seriesToBlock(Y, d = d)
      m1 <- mixture(Y = Y1, experts= X1, model = algo, loss.type = cur_loss, quiet = TRUE)
      expect_equal(m,m1)
    }
  }
})

test_that("Names are correctly rendered",{
  colnames(X) <- c("bob","alice")
  m <- mixture(Y = Y, experts=X, quiet = TRUE)
  expect_equal(m$names.experts,colnames(X))
  
  m <- mixture()
  m <- predict(m,newY = Y[1],newexperts = X[1,], quiet = TRUE)
  expect_equal(m$names.experts,colnames(X))
})
