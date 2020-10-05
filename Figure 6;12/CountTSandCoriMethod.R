#####Incidence Data
#####Cyprus Data from 01/03/2020 until 31/05/2020
############################################################
Y=c(0,0,0,0,0,0,1,0,1,4,2,6,12,8,12,5,8,16,13,17,10,6,13,18,14,32,31,
   26,32,37,39,29,46,21,25,9,38,39,23,17,18,37,21,26,24,15,9,5,6,0,11,6,
    7,13,2,8,4,16,5,9,4,7,3,2,4,5,6,2,3,8,1,0,3,2,3,2,4,0,1,5,3,8,1,2,0,2,
    0,3,0,2,3,1)


#####Define corresponding dates 
dates <- seq(as.Date("2020/3/1"), length=length(Y), by="day")
dates
length(dates)
Sys.setlocale("LC_TIME", "C")
format(Sys.Date(), "%Y-%b-%d")



#####Analysis starts on 4/03/2020
dates=dates[4:length(Y)]
Y=Y[4:length(Y)]

###############################
#################################
library(tscount)
set.seed(12345)

#####Poisson Analysis for daily number 
#######Whole data analysis without interventions based on
#######log-linear model
fit.all.log <- tsglm(Y, model=list(past_obs=1, past_mean=1), link="log")
summary(fit.all.log)
###################################################################



#####Fokianos and Fried method for detecting interventions (takes some time to run)
Covidfit_intervmultiple <- interv_multiple(fit=fit.all.log, taus=10:85,
                              deltas=c(0,0.8,1), B=500, signif_level=0.05)


####Fitting a regression model by defining AO at times 10 and 23

my.interventions <- interv_covariate(n = length(Y), tau = c(10, 23), delta = c(0,0))
my.fit.int.log <- tsglm(Y, model=list(past_obs=1, past_mean=1), xreg=my.interventions, link="log")
summary(my.fit.int.log)

#####Corresponding predictions
my.int.fut <- matrix(0, nrow=7, ncol=2)
X=predict(my.fit.int.log, n.ahead =7, newxreg=my.int.fut,level = 0.95, global=TRUE, B=5000)
new.dates=c(dates[-length(Y)], seq(as.Date(dates[length(Y)]), length=8, by="day"))

data.plot <- c(Y,round(X$pred))

plot(as.Date(new.dates), data.plot, type = "o", ylab = "Number of Cases", xlab="Time")
lines(as.Date(dates), fitted(my.fit.int.log), col = "red", lty = "longdash", lwd = 2)
arrows(x0 = new.dates[time(X$interval)+length(Y)],
       y0 = X$interval[, "lower"],
       y1 = X$interval[, "upper"],
       angle = 90, code = 3, length = 0.04, col = "darkgrey", lwd = 2)
points(new.dates[(length(Y)+1):(length(Y)+7)], data.plot[(length(Y)+1):(length(Y)+7)], pch = 16, type = "o") 
lines(new.dates[(length(Y)):(length(Y)+7)], c( fitted(my.fit.int.log)[length(Y)],
                                               data.plot[(length(Y)+1):(length(Y)+7)]), col = "black",  lty = "solid", lwd = 2)

#####Method of Cori et al. 

library(EpiEstim)

########Start after some time to get more reliable results

Y.new=Y[10:length(Y)]
length(Y.new)

I=Y.new
In.df <-data.frame(I, dates[10:length(Y)])

res <- estimate_R(In.df, method = "parametric_si",
                  config = make_config(list(
                    mean_si = 6.48, std_si = 3.83)))

plot(res, what="R",options_R = list(ylab=expression(R[t])))

