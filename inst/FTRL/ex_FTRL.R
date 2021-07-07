library('opera')  # load the package
set.seed(1)

# Example: find the best one week ahead forecasting strategy (weekly data)
# packages
library(mgcv)

# import data
data(electric_load)
idx_data_test <- 620:nrow(electric_load)
data_train <- electric_load[-idx_data_test, ]
data_test <- electric_load[idx_data_test, ]

# Build the expert forecasts 
# ##########################
# 1) A generalized additive model
gam.fit <- gam(Load ~ s(IPI) + s(Temp) + s(Time, k=3) + 
                 s(Load1) + as.factor(NumWeek), data = data_train)
gam.forecast <- predict(gam.fit, newdata = data_test)

# 2) An online autoregressive model on the residuals of a medium term model
# Medium term model to remove trend and seasonality (using generalized additive model)
detrend.fit <- gam(Load ~ s(Time,k=3) + s(NumWeek) + s(Temp) + s(IPI), data = data_train)
electric_load$Trend <- c(predict(detrend.fit), predict(detrend.fit,newdata = data_test))
electric_load$Load.detrend <- electric_load$Load - electric_load$Trend


# Residual analysis
ar.forecast <- numeric(length(idx_data_test))
for (i in seq(idx_data_test)) {
  ar.fit <- ar(electric_load$Load.detrend[1:(idx_data_test[i] - 1)])
  ar.forecast[i] <- as.numeric(predict(ar.fit)$pred) + electric_load$Trend[idx_data_test[i]]
}

# Aggregation of experts
###########################
X <- cbind(gam.forecast, ar.forecast, gam.forecast + 1, ar.forecast - 1, gam.forecast + 2, 
           ar.forecast - 2, gam.forecast + 3, ar.forecast - 3, gam.forecast + 4, 
           ar.forecast - 4, gam.forecast + 5, ar.forecast - 5, gam.forecast + 6, 
           ar.forecast - 6, gam.forecast + 7, ar.forecast - 7, gam.forecast + 8, 
           ar.forecast - 8, gam.forecast + 9, ar.forecast - 9)
names(X) <- c("gam.forecast", "ar.forecast", "gam.forecast1", "ar.forecast1", "ar.forecast1", 
              "ar.forecast2", "gam.forecast3", "ar.forecast3", "gam.forecast4", 
              "ar.forecast4", "gam.forecast5", "ar.forecast5", "gam.forecast6", 
              "ar.forecast6", "gam.forecast7", "ar.forecast7", "gam.forecast8", 
              "ar.forecast8", "gam.forecast9", "ar.forecast9")
colnames(X) <- c('gam', 'ar', 'gam1', 'ar1', 'gam2', 'ar2', 'gam3', 'ar3', 'gam4', 
                 'ar4', 'gam5', 'ar5', 'gam6', 'ar6', 'gam7', 'ar7', 'gam8', 'ar8', 'gam9', 'ar9')
Y <- data_test$Load

matplot(cbind(Y, X), type = 'l', col = 1:6, ylab = 'Weekly load', xlab = 'Week')


res <- mixture(Y = Y, experts = X, 
               model = 'FTRL', 
               loss.gradient = TRUE, 
               use_cpp = FALSE, 
               parameters = list("eta" = 0.1, 
                                 "fn" = function(x) sqrt(sum((x)**2)),
                                 "heq" = function(x) sum(x) - 1, 
                                 "hin" = function(x) x))



y = Y
experts = X 
eta = 0.1
reg = function(x) sqrt(sum(x^2))
heq = function(x) sum(x) - 1
heq_jac = NULL
hin = function(x) x
hin_jac = NULL
loss.type = "square"
loss.gradient = TRUE
w0 = 0
itmax = 50
