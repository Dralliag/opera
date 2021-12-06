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

X <- cbind(gam.forecast, ar.forecast)
colnames(X) <- c('gam', 'ar')
Y <- data_test$Load

matplot(cbind(Y, X), type = 'l', col = 1:6, ylab = 'Weekly load', xlab = 'Week')


# How good are the expert? Look at the oracles
oracle.convex <- oracle(Y = Y, experts = X, loss.type = 'square', model = 'convex')

if(interactive()){
  plot(oracle.convex)
}

oracle.convex

# Is a single expert the best over time ? Are there breaks ?
oracle.shift <- oracle(Y = Y, experts = X, loss.type = 'percentage', model = 'shifting')
if(interactive()){
  plot(oracle.shift)
}
oracle.shift

# Online aggregation of the experts with BOA
#############################################

# Initialize the aggregation rule
m0.BOA <- mixture(model = 'BOA', loss.type = 'square')

# Perform online prediction using BOA There are 3 equivalent possibilities 1)
# start with an empty model and update the model sequentially
m1.BOA <- m0.BOA
for (i in 1:length(Y)) {
  m1.BOA <- predict(m1.BOA, newexperts = X[i, ], newY = Y[i], quiet = TRUE)
}

# 2) perform online prediction directly from the empty model
m2.BOA <- predict(m0.BOA, newexpert = X, newY = Y, online = TRUE, quiet = TRUE)

# 3) perform the online aggregation directly
m3.BOA <- mixture(Y = Y, experts = X, model = 'BOA', loss.type = 'square', quiet = TRUE)

# These predictions are equivalent:
identical(m1.BOA, m2.BOA)  # TRUE
identical(m1.BOA, m3.BOA)  # TRUE

# Display the results
summary(m3.BOA)
if(interactive()){
  plot(m1.BOA)
}

# Using d-dimensional time-series
##################################

# Consider the above exemple of electricity consumption 
#  to be predicted every four weeks
YBlock <- seriesToBlock(X = Y, d = 4)
XBlock <- seriesToBlock(X = X, d = 4)

# The four-week-by-four-week predictions can then be obtained 
# by directly using the `mixture` function as we did earlier. 

MLpolBlock <- mixture(Y = YBlock, experts = XBlock, model = "MLpol", loss.type = "square", 
                      quiet = TRUE)


# The predictions can finally be transformed back to a 
# regular one dimensional time-series by using the function `blockToSeries`.

prediction <- blockToSeries(MLpolBlock$prediction)

#### Using the `online = FALSE` option

# Equivalent solution is to use the `online = FALSE` option in the predict function. 
# The latter ensures that the model coefficients are not 
# updated between the next four weeks to forecast.
MLpolBlock <- mixture(model = "BOA", loss.type = "square")
d = 4
n <- length(Y)/d
for (i in 0:(n-1)) { 
  idx <- 4*i + 1:4 # next four weeks to be predicted
  MLpolBlock <- predict(MLpolBlock, newexperts = X[idx, ], newY = Y[idx], online = FALSE, 
                        quiet = TRUE)
}

print(head(MLpolBlock$weights))

