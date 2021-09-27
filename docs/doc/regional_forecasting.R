## -----------------------------------------------------------------------------
rm(list=objects())
library(magrittr)
library(opera)

regional_load <- readRDS(system.file("extdata", "regional_load.RDS", package = "opera", mustWork = TRUE))

n0 <- c(1:2011)
Data0 <- regional_load[n0,]
Data1 <- regional_load[-n0,]


name_reg <- c("Nouvelle_A", "Auvergne_R", "Bourgogne", "Occitanie", 
              "Hauts_de_F", "Normandie", "Bretagne", 
              "Centre_Val", "Ile_de_Fra", "Pays_de_la_Loire", 
              "Provence_A", "Grand_Est")

nb.cols <- length(name_reg)
col <- colorRampPalette( RColorBrewer::brewer.pal(8, "Set3"))(nb.cols)

plot(Data0$Date, Data0[,name_reg[1]], type = 'l', xlim = range(Data0$Date, Data1$Date), 
     ylab = 'Regional load (MW)', xlab = "Date", col = col[1], ylim = range(Data0[,name_reg]))
lines(Data1$Date, Data1[,name_reg[1]], col = col[1], lty = 'dashed')
for(i in c(2:length(name_reg))){
  lines(Data0$Date, Data0[,name_reg[i]], col = col[i])
  lines(Data1$Date, Data1[,name_reg[i]], col = col[i], type = 'l', lty = 'dashed')
}
abline(v = Data1$Date[1], lty = "dashed")

## -----------------------------------------------------------------------------

plot(Data0$Date, Data0$Load, type = 'l', 
     xlim = range(Data0$Date, Data1$Date), 
     ylab = 'National load (MW)', xlab = "Date",)
lines(Data1$Date, Data1$Load, col = 'red', type = 'l')
legend("topleft", c("learning set", "test set"), col = c("black", "red"), 
       lty = 1, bty = 'n', ncol = 2)
abline(v = Data1$Date[1], lty = "dashed")


## ---- eval=T------------------------------------------------------------------
eq <- Load ~ s(as.numeric(Date), k=3)  + WeekDays:DLS  + s(toy, k = 20, bs = 'cc') + 
  s(as.numeric(Date),T_national, k=4) + s(T_national_s95, k=5) + s(T_national_s99, k=5) + 
  s(T_national_s99_min, T_national_s99_max) + Load.48:WeekDays + Load.336 

g <- mgcv::gam(eq , data = Data0, family = "gaussian")
gam_nat.forecast <- predict(g, newdata = Data1)

## ---- eval=FALSE--------------------------------------------------------------
#  ### extract the GAM terms
#  terms0 <- predict(g, newdata = Data0, type = 'terms')
#  terms1 <- predict(g, newdata = Data1, type = 'terms')
#  colnames(terms0) <- paste0("gterms_", c(1:ncol(terms0)))
#  colnames(terms1) <- paste0("gterms_", c(1:ncol(terms1)))
#  
#  ### include it in the stacked RF
#  #### first build a data frame
#  data0_rf <- data.frame(Data0, terms0)
#  data0_rf$res <- g$residuals  # fitted residuals of the GAM
#  data0_rf$res.48 <- c(data0_rf$res[1], data0_rf$res[1:(length(data0_rf$res)-1)])  # lags of the fitted residuals
#  data0_rf$res.336 <- c(data0_rf$res[1:7], data0_rf$res[1:(length(data0_rf$res)-7)])
#  
#  data1_rf <- data.frame(Data1, terms1)
#  data1_rf$g.forecast <- predict(g, newdata=Data1)
#  residuals <- data1_rf$Load -  data1_rf$g.forecast #forecast  residuals
#  data1_rf$res.48 <- c(residuals[1], residuals[1:(length(residuals)-1)]) #lags of the forecast residuals
#  data1_rf$res.336 <- c(residuals[1:7], residuals[1:(length(residuals)-7)])
#  
#  ### equation of the RF
#  cov <- c(paste0('Load', c(".48", ".336")),
#           paste0('T_', "national"),
#           paste0('T_', "national", '_s99'),
#           paste0('T_', "national", '_s95'),
#           paste0('T_', "national", '_s95_min'),
#           paste0('T_', "national", '_s95_max'),
#           'toy', 'BH', 'tod', 'DLS', 'Summer_break', 'Christmas_break', 'res.48', 'res.336')
#  
#  gterm <- paste0("gterms_", c(1:ncol(terms0)))
#  cov <- paste0(c(cov, gterm), collapse = '+')
#  formule_rf <- paste0("res", "~", cov)
#  
#  gamrf_nat <- ranger::ranger(formule_rf, data = data0_rf, quantreg = T)

## ---- eval=FALSE--------------------------------------------------------------
#  qu <- c(0.05, 0.1,0.5,0.9, 0.95)
#  gamrf_nat.forecast <- predict(gamrf_nat, data = data1_rf, quantiles = qu, type='quantiles')$predictions +
#    as.vector(gam_nat.forecast)

## -----------------------------------------------------------------------------
experts <- readRDS(system.file("extdata", "regional_experts.RDS", package = "opera", mustWork = TRUE))
experts_raw <- readRDS(system.file("extdata", "regional_experts_raw.RDS", package = "opera", mustWork = TRUE))

dim(experts)
agg <- opera::mixture(Y = as.numeric(Data1$Load), experts=experts)

## Run if you want to see dynamic ploy
# plot(agg)

mape <- function(y, ychap, digits = 3){
  return(signif(100*mean(abs(y-ychap)/abs(y),na.rm=TRUE),digits=digits))
}

## RTE national forecast, dayahead
mape(y = Data1$Load, ychap = Data1$Forecast_RTE_dayahead) 
## National forecast
mape(y = Data1$Load, ychap = experts[,  "nat0.5"]) 
## Aggregation 
mape(y = Data1$Load, ychap = agg$prediction) 
## sum of regional forecast
mape(y = Data1$Load, ychap = rowSums(experts_raw[, grep("0.5", colnames(experts_raw))])) 

## Pre_covid
pre_covid <- which(Data1$Date >= as.POSIXct(strptime("2019-09-01 00:00:00", "%Y-%m-%d %H:%M:%S"), tz="UTC") &
                 Data1$Date < as.POSIXct(strptime("2020-03-16 00:00:00", "%Y-%m-%d %H:%M:%S"), tz="UTC"))
mape(y = Data1$Load[pre_covid], ychap = agg$prediction[pre_covid])

## hard lockdown
hard_lockdown <- which(Data1$Date >= as.POSIXct(strptime("2020-03-16 00:00:00", "%Y-%m-%d %H:%M:%S"), tz="UTC") &
                 Data1$Date < as.POSIXct(strptime("2020-04-16 00:00:00", "%Y-%m-%d %H:%M:%S"), tz="UTC"))
mape(y = Data1$Load[hard_lockdown], ychap = agg$prediction[hard_lockdown])

## "post" lockdown
post_lockdown <- which(Data1$Date >= as.POSIXct(strptime("2020-04-16 00:00:00", "%Y-%m-%d %H:%M:%S"), tz="UTC"))
mape(y = Data1$Load[post_lockdown], ychap = agg$prediction[post_lockdown])


## -----------------------------------------------------------------------------

reg_scale <- Data1[,  name_reg]

m <- matrix(c(1,mean(Data0$Load)/colMeans(Data0[,name_reg])), 
            nrow = nrow(Data1), ncol = length(name_reg) + 1, byrow = T)
Y <- as.matrix(Data1[, c("Load", name_reg)]*m)
experts2 <- array(0, dim = c(nrow(Data1), length(name_reg) + 1, 5))
name <-  c("nat", name_reg)
for(k in c(1:(length(name_reg)+1))){
  experts2[,k,] <- experts[, grep(name[k], colnames(experts))]
}

agg_vector<- opera::mixture(Y = Y, experts = experts2)

plot(agg_vector, dynamic = FALSE)
# plot(agg_vector)

mape(y = Data1$Load, ychap = agg_vector$prediction[,1]) 
mape(y = Data1$Load[pre_covid], ychap = agg_vector$prediction[pre_covid,1])
mape(y = Data1$Load[hard_lockdown], ychap = agg_vector$prediction[hard_lockdown,1])
mape(y = Data1$Load[post_lockdown], ychap = agg_vector$prediction[post_lockdown,1])



## -----------------------------------------------------------------------------
X <- experts
Y <- Data1$Load
agg <- opera::mixture(Y=as.numeric(Data1$Load), experts=experts)


change_horizon <- function(h, agg){
  ## the h first weights stay equal to inital values of the weights
  rep(agg$weights[1,], h) %>% matrix(ncol = ncol(agg$weights)) %>% dim 
  w <- rbind(rep(agg$weights[1,], h) %>% 
               matrix(ncol = ncol(agg$weights)), agg$weights[1:(nrow(agg$weights)-h),])
  
  agg_horizon <- agg
  agg_horizon$weights <- w
  agg_horizon$prediction <-  rowSums(w*agg$experts)
  return(agg_horizon)
}


map_horizon <- mape(y = Y, ychap = agg$prediction)
for(h in c(1:14)){
  agg_horizon <- change_horizon(h = h, agg)
  map_horizon <- c(map_horizon, mape(y = Y, ychap = agg_horizon$prediction))
}
plot(map_horizon, type='l')


