## ----echo=FALSE----------------------------------------------------------
library("knitr")
library("htmltools")

knitr::opts_chunk$set(
  comment = "#>",
  collapse = TRUE,
  warning = TRUE,
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

## ----eval=FALSE----------------------------------------------------------
#  install.packages("opera")

## ----eval=FALSE----------------------------------------------------------
#  install.packages("devtools")
#  devtools::install_github("dralliag/opera")

## ---- echo=FALSE---------------------------------------------------------
library(opera)
set.seed(1)

## ------------------------------------------------------------------------
data(electric_load)
attach(electric_load)
idx_data_test <- 620:nrow(electric_load)
data_train <- electric_load[-idx_data_test, ] 
data_test <- electric_load[idx_data_test, ]  

## ----Load, echo = TRUE, eval = FALSE-------------------------------------
#  plot(Load, type = "l", main = "The electric Load")

## ----Load, echo = FALSE, imgcenter = TRUE--------------------------------
plot(Load, type = "l", main = "The electric Load")

## ----Temp, echo = TRUE, eval = FALSE-------------------------------------
#  plot(Temp, Load, pch = 16, cex = 0.5, main = "Temperature vs Load")

## ----Temp, echo = FALSE, imgcenter = TRUE--------------------------------
plot(Temp, Load, pch = 16, cex = 0.5, main = "Temperature vs Load")

## ----NumWeek, echo = TRUE, eval = FALSE----------------------------------
#  plot(NumWeek, Load, pch = 16, cex = 0.5, main = "Annual seasonality")

## ----NumWeek, echo = FALSE, imgcenter = TRUE-----------------------------
plot(NumWeek, Load, pch = 16, cex = 0.5, main = "Annual seasonality")

## ----gam, results="hide",message=F, warning=F----------------------------
library(mgcv)
gam.fit <- gam(Load ~ s(IPI) + s(Temp) + s(Time, k=3) + 
                 s(Load1) + as.factor(NumWeek), data = data_train)
gam.forecast <- predict(gam.fit, newdata = data_test)

## ----gamar,results="hide"------------------------------------------------
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

## ----gbm,results="hide",message=FALSE,warning=FALSE----------------------
library(caret)
gbm.fit <- train(Load ~ IPI + IPI_CVS + Temp + Temp1 + Time + Load1 + NumWeek, 
                 data = data_train, method = "gbm")
gbm.forecast <- predict(gbm.fit, newdata = data_test)

## ----loadAndForecasts, echo = TRUE, eval = FALSE-------------------------
#  Y <- data_test$Load
#  X <- cbind(gam.forecast, ar.forecast, gbm.forecast)
#  matplot(cbind(Y, X), type = "l", col = 1:6, ylab = "Weekly load",
#          xlab = "Week", main = "Expert forecasts and observations")

## ----loadAndForecasts, echo = FALSE, imgcenter = TRUE--------------------
Y <- data_test$Load
X <- cbind(gam.forecast, ar.forecast, gbm.forecast)
matplot(cbind(Y, X), type = "l", col = 1:6, ylab = "Weekly load", 
        xlab = "Week", main = "Expert forecasts and observations")

## ---- echo=FALSE---------------------------------------------------------
colnames(X) <- c("gam", "ar", "gbm")

## ----oracle, echo = TRUE, eval = FALSE-----------------------------------
#  oracle.convex <- oracle(Y = Y, experts = X, loss.type = "square", model = "convex")
#  plot(oracle.convex)

## ----oracle, echo = FALSE, imgcenter = TRUE------------------------------
oracle.convex <- oracle(Y = Y, experts = X, loss.type = "square", model = "convex")
plot(oracle.convex)

## ----printoracle---------------------------------------------------------
print(oracle.convex)

## ----MLpolInit-----------------------------------------------------------
MLpol0 <- mixture(model = "MLpol", loss.type = "square")

## ----MLpolFit------------------------------------------------------------
MLpol <- MLpol0
for (i in 1:length(Y)) {
  MLpol <- predict(MLpol, newexperts = X[i, ], newY = Y[i])
}

## ----MLpolsummary--------------------------------------------------------
summary(MLpol)

## ----MLpol, echo = TRUE, eval = FALSE------------------------------------
#  plot(MLpol, pause = TRUE, col = brewer.pal(3,name = "Set1"))

## ----MLpol, echo = FALSE, imgcenter = TRUE-------------------------------
plot(MLpol, pause = TRUE)

## ----MLpolPredict, eval = FALSE------------------------------------------
#  MLpol <- predict(MLpol0, newexpert = X, newY = Y, online = TRUE)

## ----MLpolDirect, eval = FALSE-------------------------------------------
#  MLpol <- mixture(Y = Y, experts = X, model = "MLpol", loss.type = "square")

