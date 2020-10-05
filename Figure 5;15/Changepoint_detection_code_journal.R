#####Incidence Data
#####Cyprus Data from 01/03/2020 until 31/05/2020
############################################################
data_used <- c(0,0,0,0,0,0,1,0,1,4,2,6,12,8,12,5,8,16,13,17,10,6,13,18,14,32,31,
    26,32,37,39,29,46,21,25,9,38,39,23,17,18,37,21,26,24,15,9,5,6,0,11,6,
    7,13,2,8,4,16,5,9,4,7,3,2,4,5,6,2,3,8,1,0,3,2,3,2,4,0,1,5,3,8,1,2,0,2,
    0,3,0,2,3,1)

l <- length(data_used)

##### Code and plots for continuous, piecewise-linear scenario

# Change-point detection
library(IDetect)
cpt_all_info <- cplm_ic(data_used, points = 7, th_const = 1.2)
cpt <- cpt_all_info$cpt_ic$ssic_pen

# Forecasting using only the last segment
library(forecast)

fcast_linear <- forecast(data_used[49:l], 7, level = 0.95) ## based on the last detected change-point

# Deriving the point estimators and the lower and upper bounds for forecasting
lower_vector<- numeric()
for (i in 1:length(fcast_linear$lower)){
  lower_vector[i] <- ifelse(fcast_linear$lower[i] < 0, 0, floor(fcast_linear$lower[i]))
}

mean_vector<- numeric()
for (i in 1:length(fcast_linear$mean)){
  mean_vector[i] <- ifelse(fcast_linear$mean[i] < 0, 0, round(fcast_linear$mean[i]))
}

upper_vector <- numeric()
for (i in 1:length(fcast_linear$upper)){
  upper_vector[i] <- ifelse(fcast_linear$upper[i] < 0, 0, round(fcast_linear$upper[i]))
}

# Creating the plot
par(mfrow=c(1,1),mar=c(5,6,4,1))
par(cex.lab=2.35,cex.main=2.35,cex.axis=2.35,bg="white")
plot.ts(data_used, main = "Cases, estimated change-points, and the signal", ylab = "Number of cases", xlab = "Date", xaxt="n", xlim = c(0,99))
axis(1, at = seq(1,99,7), labels=c("01/03", "08/03", "15/03", "22/03", "29/03", "05/04", "12/04", "19/04", "26/04", "03/05", "10/05", "17/05", "24/05", "31/05", "07/06"))
lines(est_signal(data_used, cpt, type = "slope"), col = "red", lwd = 3)
abline(v = cpt, col = "blue", lty = 2, lwd = 2)
points(93:99, mean_vector, col = "black", pch = 16, cex = 1.2)
points(93:99, lower_vector, col = "green4", pch = 16, cex = 1.2)
points(93:99, upper_vector, col = "red", pch = 16, cex = 1.2)
for(i in 1:7)
{lines(c(92+i,92+i), c(lower_vector[i], upper_vector[i]), lwd = 2)}


##### Code and plots for piecewise-constant scenario

# Change-point detection
library(IDetect)
cpt_all_info_constant <- pcm_ic(data_used, points = 7, th_const = 1.6)
cpt_constant <- cpt_all_info_constant$cpt_ic$ssic_pen

# Forecasting using only the last segment
library(forecast)

fcast_constant <- forecast(data_used[63:l], 7, level = 0.95) ## based on the last detected change-point

# Deriving the point estimators and the lower and upper bounds for forecasting
lower_vector_const<- numeric()
for (i in 1:length(fcast_constant$lower)){
  lower_vector_const[i] <- ifelse(fcast_constant$lower[i] < 0, 0, floor(fcast_constant$lower[i]))
}

mean_vector_const<- numeric()
for (i in 1:length(fcast_constant$mean)){
  mean_vector_const[i] <- ifelse(fcast_constant$mean[i] < 0, 0, round(fcast_constant$mean[i]))
}

upper_vector_const <- numeric()
for (i in 1:length(fcast_constant$upper)){
  upper_vector_const[i] <- ifelse(fcast_constant$upper[i] < 0, 0, round(fcast_constant$upper[i]))
}

# Creating the plot
par(mfrow=c(1,1),mar=c(5,6,4,1))
par(cex.lab=2.35,cex.main=2.35,cex.axis=2.35,bg="white")
plot.ts(data_used, main = "Cases, estimated change-points, and the signal", ylab = "Number of cases", xlab = "Date", xaxt="n", xlim = c(0,99))
axis(1, at = seq(1,99,7), labels=c("01/03", "08/03", "15/03", "22/03", "29/03", "05/04", "12/04", "19/04", "26/04", "03/05", "10/05", "17/05", "24/05", "31/05", "07/06"))
lines(est_signal(data_used, cpt_constant, type = "mean"), col = "red", lwd = 3)
abline(v = cpt_constant, col = "blue", lty = 2, lwd = 2)
points(93:99, mean_vector_const, col = "black", pch = 16, cex = 1.2)
points(93:99, lower_vector_const, col = "green4", pch = 16, cex = 1.2)
points(93:99, upper_vector_const, col = "red", pch = 16, cex = 1.2)
for(i in 1:7)
{lines(c(92+i,92+i), c(lower_vector_const[i], upper_vector_const[i]), lwd = 2)}
