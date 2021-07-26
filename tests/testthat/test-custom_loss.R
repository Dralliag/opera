# Unit tests of opera package using testhat package

context("Testing custom losses")


# Test of loss functions
test_that("custom loss functions return correct values", {
  x <- runif(10)
  y <- runif(1)
  pred <- runif(10)
  
  # SQUARE
  lossp_custom <- loss(x, y, pred, loss.type = function(x, y) (x-y)**2, loss.gradient = function(x, y) 2 * c(x - y))
  lossp_opera <- loss(x, y, pred, loss.type = list(name = "square"), loss.gradient = T)
  expect_equal(lossp_custom, lossp_opera)
  
  # PERCENTAGE
  lossp_custom <- loss(x, y, pred, loss.type = function(x, y) abs(x - y)/y, loss.gradient = function(x, y) 1/y * sign(c(x - y)))
  lossp_opera <- loss(x, y, pred, loss.type = list(name = "percentage"), loss.gradient = T)
  expect_equal(lossp_custom, lossp_opera)
  
  # ABSOLUTE
  lossp_custom <- loss(x, y, pred, loss.type = function(x, y) abs(x - y), loss.gradient = function(x, y) sign(c(x - y)))
  lossp_opera <- loss(x, y, pred, loss.type = list(name = "absolute"), loss.gradient = T)
  expect_equal(lossp_custom, lossp_opera)
  
  # PINBALL
  lossp_custom <- loss(x, y, pred, loss.type = function(x, y) (0.5 - (y < x)) * (y - x), loss.gradient = function(x, y) c((y < x) - 0.5))
  lossp_opera <- loss(x, y, pred, loss.type = list(name = "pinball", tau = 0.5), loss.gradient = T)
  expect_equal(lossp_custom, lossp_opera)
})
