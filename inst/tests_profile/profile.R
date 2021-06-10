require(profvis)
require(opera)
require(mgcv)
require(caret)
require(gbm)

data(electric_load)

n <- 5
electric_load_more <- electric_load
for(yy in 1:n){
  tmp <- electric_load
  step <- max(electric_load_more$Year) - min(tmp$Year)
  tmp$Year <- tmp$Year + step + 1
  electric_load_more <- rbind(electric_load_more, tmp)
}

idx_data_train <- round(nrow(electric_load_more) * 0.50)


data_train <- electric_load_more[1:idx_data_train, ] 
data_test <- electric_load_more[-c(1:idx_data_train), ]  


gam.fit <- gam(Load ~ s(IPI) + s(Temp) + s(Time, k=3) + 
                 s(Load1) + as.factor(NumWeek), data = data_train)
gam.forecast <- predict(gam.fit, newdata = data_test)

# medium term model
medium.fit <- gam(Load ~ s(Time,k=3) + s(NumWeek) + s(Temp) + s(IPI), data = data_train)
electric_load_more$Medium <- c(predict(medium.fit), predict(medium.fit, newdata = data_test))
electric_load_more$Residuals <- electric_load_more$Load - electric_load_more$Medium



gbm.fit <- train(Load ~ IPI + IPI_CVS + Temp + Temp1 + Time + Load1 + NumWeek, 
                 data = data_train, method = "gbm")
gbm.forecast <- predict(gbm.fit, newdata = data_test)

Y <- data_test$Load
X <- cbind(gam.forecast, gbm.forecast)

res <- X
for(n in 1:100){
  tmp <- X + runif(min = -10000, max = 10000, n = nrow(X) * 2)
  colnames(tmp) <- paste0(colnames(X), "_", n)
  res <- cbind(res, tmp)
}

n <- 100
prof_oracl_convex <- profvis({
  oracle.convex <- oracle(Y = Y[1:n], experts = res[1:n, ], loss.type = "square", model = "convex")
}, interval = 0.005)
htmlwidgets::saveWidget(prof_oracl_convex, file = "prof_oracl_convex.html")

# prof_oracl_linear <- profvis({
#   oracle.convex <- oracle(Y = Y[1:n], experts = res[1:n, ], loss.type = "square", model = "linear")
# }, interval = 0.005)
# htmlwidgets::saveWidget(prof_oracl_linear, file = "prof_oracl_linear.html")

n <- 500
prof_oracl_shifting <- profvis({
  oracle.convex <- oracle(Y = Y[1:n], experts = res[1:n, ], loss.type = "square", model = "shifting")
}, interval = 0.005)
htmlwidgets::saveWidget(prof_oracl_shifting, file = "prof_oracl_shifting.html")

n <- 500
prof_oracl_expert <- profvis({
  oracle.convex <- oracle(Y = Y, experts = res, loss.type = "square", model = "expert")
}, interval = 0.005)
htmlwidgets::saveWidget(prof_oracl_expert, file = "prof_oracl_expert.html")

n <- 500
profvis({
  oracle.convex <- oracle(Y = Y[1:n], experts = res[1:n, ], loss.type = "square", model = "convex")
}, interval = 0.01)

