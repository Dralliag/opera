context("Testing seriesToBlock and blockToSeries")

# load some basic data to perform tests
n <- 10
d <- 3
K <- 23
X <- matrix(runif(K*n*d), ncol = K, nrow = K*d)
Y <- rnorm(d*n)

test_that("block functions are OK", {
  expect_equal(blockToSeries(seriesToBlock(X,d)),X)
  expect_equal(seriesToBlock(blockToSeries(X),dim(X)[2]),X)
  expect_equal(blockToSeries(seriesToBlock(Y,d)),Y)
})

