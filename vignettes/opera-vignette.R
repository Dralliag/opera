## ---- echo=FALSE---------------------------------------------------------
library(opera)
set.seed(1)

## ------------------------------------------------------------------------
data(electric_load)
idx_data_test <- 620:nrow(electric_load)
data_train <- electric_load[-idx_data_test, ] 
data_test <- electric_load[idx_data_test, ]  

## ---- fig.height=3,fig.width=4.5-----------------------------------------
attach(electric_load)
plot(Load, type = "l", main = "The electric Load")
plot(Temp, Load, pch = 16, cex = 0.5, main = "Temperature vs Load")
plot(NumWeek, Load, pch = 16, cex = 0.5, main = "Annual seasonality")

## ----chunk_name, results="hide",message=F, warning=F---------------------
library(mgcv)
gam.fit <- gam(Load ~ s(IPI) + s(Temp) + s(Time, k=3) + 
                s(Load1) + as.factor(NumWeek), data = data_train)
gam.forecast <- predict(gam.fit, newdata = data_test)

## ----results="hide"------------------------------------------------------
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

## ----results="hide",message=FALSE,warning=FALSE--------------------------
library(caret)
gbm.fit <- train(Load ~ IPI + IPI_CVS + Temp + Temp1 + Time + Load1 + NumWeek, 
                  data = data_train, method = "gbm")
gbm.forecast <- predict(gbm.fit, newdata = data_test)

## ----fig.height=3,fig.width=4.5------------------------------------------
Y <- data_test$Load
X <- cbind(gam.forecast, ar.forecast, gbm.forecast)
matplot(cbind(Y, X), type = "l", col = 1:6, ylab = "Weekly load", xlab = "Week", main = "Expert forecasts and observations")

## ---- echo=FALSE---------------------------------------------------------
colnames(X) <- c("gam", "ar", "gbm")

## ----fig.height=3,fig.width=4.5------------------------------------------
oracle.convex <- oracle(Y = Y, experts = X, loss.type = "square", model = "convex")
plot(oracle.convex)
print(oracle.convex)

## ------------------------------------------------------------------------
MLpol0 <- mixture(model = "MLpol", loss.type = "square")

## ------------------------------------------------------------------------
MLpol <- MLpol0
for (i in 1:length(Y)) {
  MLpol <- predict(MLpol, newexperts = X[i, ], newY = Y[i])
}

## ----fig.height=3,fig.width=4.5------------------------------------------
summary(MLpol)
plot(MLpol, pause = TRUE)

## ---- eval = FALSE-------------------------------------------------------
#  MLpol <- predict(MLpol0, newexpert = X, newY = Y, online = TRUE)

## ---- eval = FALSE-------------------------------------------------------
#  MLpol <- mixture(Y = Y, experts = X, model = "MLpol", loss.type = "square")

