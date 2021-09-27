## ----echo=FALSE---------------------------------------------------------------
suppressPackageStartupMessages({
  suppressMessages({
    suppressWarnings({
      library("knitr")
      library("htmltools")
    })
  })
})


knitr::opts_chunk$set(
  comment = "#>",
  collapse = TRUE,
  warning = FALSE,
  message = FALSE,
  fig.width = 4.5,
  fig.height = 3,
  fig.path = "../inst/img/",
  out.width = 300,
  cache.path = "../inst/cache/"
)
knitr::knit_hooks$set(imgcenter = function(before, options, envir){
  if (before) {
    HTML("<p align='center'>")
  } else {
    HTML("</p>")
  }
})

## ----eval=FALSE---------------------------------------------------------------
#  install.packages("opera")

## ----eval=FALSE---------------------------------------------------------------
#  install.packages("devtools")
#  devtools::install_github("dralliag/opera")

## ---- echo=TRUE---------------------------------------------------------------
library(opera)
set.seed(1)

## ---- echo=FALSE--------------------------------------------------------------
dynamic <- TRUE
if(exists("readme") && readme){
  dynamic <- FALSE
}

## -----------------------------------------------------------------------------
data(electric_load)
attach(electric_load)
idx_data_test <- 620:nrow(electric_load)
data_train <- electric_load[-idx_data_test, ] 
data_test <- electric_load[idx_data_test, ]  


## ----Load, echo = TRUE, eval = FALSE------------------------------------------
#  plot(Load, type = "l", main = "The electric Load")

## ----Load, echo = FALSE, imgcenter = TRUE-------------------------------------
plot(Load, type = "l", main = "The electric Load")

## ----Temp, echo = TRUE, eval = FALSE------------------------------------------
#  plot(Temp, Load, pch = 16, cex = 0.5, main = "Temperature vs Load")

## ----Temp, echo = FALSE, imgcenter = TRUE-------------------------------------
plot(Temp, Load, pch = 16, cex = 0.5, main = "Temperature vs Load")

## ----NumWeek, echo = TRUE, eval = FALSE---------------------------------------
#  plot(NumWeek, Load, pch = 16, cex = 0.5, main = "Annual seasonality")

## ----NumWeek, echo = FALSE, imgcenter = TRUE----------------------------------
plot(NumWeek, Load, pch = 16, cex = 0.5, main = "Annual seasonality")

## ----gam, results="hide",message=F, warning=F---------------------------------
library(mgcv)
gam.fit <- gam(Load ~ s(IPI) + s(Temp) + s(Time, k=3) + 
                 s(Load1) + as.factor(NumWeek), data = data_train)
gam.forecast <- predict(gam.fit, newdata = data_test)

## ----gamar,results="hide"-----------------------------------------------------
# medium term model
medium.fit <- gam(Load ~ s(Time,k=3) + s(NumWeek) + s(Temp) + s(IPI), data = data_train)
electric_load$Medium <- c(predict(medium.fit), predict(medium.fit, newdata = data_test))
electric_load$Residuals <- electric_load$Load - electric_load$Medium

# autoregressive correction
ar.forecast <- numeric(length(idx_data_test))
for (i in seq(idx_data_test)) {
  ar.fit <- ar(electric_load$Residuals[1:(idx_data_test[i] - 1)])
  ar.forecast[i] <- as.numeric(predict(ar.fit)$pred) + electric_load$Medium[idx_data_test[i]]
}

## ----gbm,results="hide",message=FALSE,warning=FALSE---------------------------
library(caret)
gbm.fit <- train(Load ~ IPI + IPI_CVS + Temp + Temp1 + Time + Load1 + NumWeek, data = data_train, method = "gbm")
gbm.forecast <- predict(gbm.fit, newdata = data_test)

## ----loadAndForecasts, echo = TRUE, eval = FALSE------------------------------
#  Y <- data_test$Load
#  X <- cbind(gam.forecast, ar.forecast, gbm.forecast)
#  matplot(cbind(Y, X), type = "l", col = 1:6, ylab = "Weekly load",
#          xlab = "Week", main = "Expert forecasts and observations")

## ----loadAndForecasts, echo = FALSE, imgcenter = TRUE-------------------------
Y <- data_test$Load
X <- cbind(gam.forecast, ar.forecast, gbm.forecast)
matplot(cbind(Y, X), type = "l", col = 1:6, ylab = "Weekly load", 
        xlab = "Week", main = "Expert forecasts and observations")

## ---- echo=FALSE--------------------------------------------------------------
Y <- data_test$Load
X <- cbind(gam.forecast, ar.forecast, gbm.forecast)
colnames(X) <- c("gam", "ar", "gbm")

## ---- echo = TRUE, eval = FALSE, warning=FALSE--------------------------------
#  oracle.convex <- oracle(Y = Y, experts = X, loss.type = "square", model = "convex")
#  print(oracle.convex)
#  plot(oracle.convex)

## ----oracle, echo = FALSE, eval = FALSE, warning=FALSE------------------------
#  oracle.convex <- oracle(Y = Y, experts = X, loss.type = "square", model = "convex")
#  print(oracle.convex)
#  plot(oracle.convex, dynamic = dynamic)

## ----oracle, echo = FALSE, imgcenter = TRUE-----------------------------------
oracle.convex <- oracle(Y = Y, experts = X, loss.type = "square", model = "convex")
print(oracle.convex)
plot(oracle.convex, dynamic = dynamic)

## ----MLpolInit----------------------------------------------------------------
MLpol0 <- mixture(model = "MLpol", loss.type = "square")

## ----MLpolFit-----------------------------------------------------------------
MLpol <- MLpol0
for (i in 1:length(Y)) {
  MLpol <- predict(MLpol, newexperts = X[i, ], newY = Y[i])
}

## ----MLpolsummary-------------------------------------------------------------
summary(MLpol)

## ---- echo = TRUE, eval = FALSE-----------------------------------------------
#  plot(MLpol, pause = TRUE)

## ----MLpol, echo = FALSE, eval = FALSE----------------------------------------
#  plot(MLpol, dynamic = dynamic, pause = TRUE)

## ----MLpol, echo = FALSE, imgcenter = TRUE------------------------------------
plot(MLpol, dynamic = dynamic, pause = TRUE)

## ----MLpolPredict, eval = FALSE-----------------------------------------------
#  MLpol <- predict(MLpol0, newexpert = X, newY = Y, online = TRUE)

## ----MLpolDirect, eval = FALSE------------------------------------------------
#  MLpol <- mixture(Y = Y, experts = X, model = "MLpol", loss.type = "square")

## ---- eval=FALSE--------------------------------------------------------------
#  predictions <- MLpol$prediction

## ---- eval=FALSE--------------------------------------------------------------
#  newexperts <- X[1:3, ] # Experts forecasts to predict 3 new points
#  pred <- predict(MLpol, newexperts = newexperts, online = FALSE, type = 'response')

## ---- eval=TRUE, echo=FALSE---------------------------------------------------
newexperts <- X[1:3, ] # Experts forecasts to predict 3 new points  
pred <- predict(MLpol, newexperts = newexperts, online = FALSE, type = 'response')
print(c(pred))

## ---- eval=FALSE--------------------------------------------------------------
#  pred = newexperts %*% MLpol$coefficients

## ---- eval=TRUE, echo=FALSE---------------------------------------------------
newexperts <- X[1:3, ] # Experts forecasts to predict 3 new points  
pred = c(newexperts %*% MLpol$coefficients)
print(pred)

## ----block_transform, eval = TRUE---------------------------------------------
YBlock <- seriesToBlock(X = Y, d = 4)
XBlock <- seriesToBlock(X = X, d = 4)

## ----block_mixture, eval = TRUE-----------------------------------------------
MLpolBlock <- mixture(Y = YBlock, experts = XBlock, model = "MLpol", loss.type = "square")

## ----block_tranformback, eval = TRUE------------------------------------------
prediction <- blockToSeries(MLpolBlock$prediction)

## ----blockbyblock, eval = TRUE------------------------------------------------
MLpolBlock <- MLpol0
d = 4
n <- length(Y)/d
for (i in 0:(n-1)) { 
  idx <- 4*i + 1:4 # next four weeks to be predicted
  MLpolBlock <- predict(MLpolBlock, newexperts = X[idx, ], newY = Y[idx], online = FALSE)
}

## ---- echo = FALSE, eval = TRUE-----------------------------------------------
  print(head(MLpolBlock$weights))

