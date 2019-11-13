#Part 1
#slide 8
1+1

#slide 10
year <- 2019	
catch_at_age <- c(10000, 30000)
year
catch_at_age
class(year)
class(catch_at_age)

#slide 11
gantan <- "2019-1-1"
class(gantan) #that's not the way
gantan <- as.Date("2019-1-1")
class(gantan)
gantan + 1 # date + numeric
gantan + "Christmas" #ERROR! date + character

#slide 13
a <- 7
b <- 2
a^b
a/b
a%/%b
a%%b
a <- c(1:10) # a is a vector
a+b #how about vector + numeric?
b <- c(10:1) # noe b is also a vector
a+b #how about vector + vector?
b <- c(1,2)
a+b #how about vector(10) + vector(2)?
b <- c(1,2,3)
a+b #how about vector(10) + vector(3)?

#slide 15
number <- 100
sum <- 0
for(i in 1:number) # change i from 1 to number (100)
{
  sum <- sum + i # self-assignment
}
sum
i

#slide 16
number <- 100
sum1 <- 0
sum2 <- 0
for(i in 1:number)
{
  if(i%%2 == 0) # if i is even number
  {
    sum1 <- sum1 + i
  }else
  {
    sum2 <- sum2 + i
  }
}
sum1
sum2

#slide 18
days <- c(1:100) # make a vector
weight <- -0.1*days + 65 + rnorm(length(days), 0, 0.4) #rnorm(): normal random number 
weight
plot(days, weight, type="o", xlab="Days", ylab="Weight") # plot 

#slide 19
?cor # how to use cor()?
?glm

#slide 20
list_factors <- function(x) # x is the argument
{
  factors <- c() # make an empty vector
  for(i in 1:x)
  {
    if(x%%i == 0) # if x is divisible by i
    {
      factors <- c(factors, i) # add i to the vector
    }
  }	
  factors # return factors
}

list_factors(36)
list_factors(1997)
factors



#Part 2
#slide 24, 25
x <- load("/Users/shin-ichironakayama/Dropbox/r_practice/beer.rda") # load data #rewrite the path!
x # confirm contents of the data
class(sales_data) # confirm class
head(sales_data) # first 6 row
tail(sales_data) # last 6 row
dim(sales_data) # confirm dimension
plot(sales_data) # overview

#slide 27
res <- lm(sales~day+temp-1, data=sales_data) #fit linear model
summary(res) #check the result

#slide 28
names(res)
res$coefficients
res$residuals
res$fitted.values

#slide 29
newdata <- data.frame("wed", 30) #create newdata
names(newdata) <- c("day", "temp") #give the same column name
predict(res, newdata=newdata) #predict the sales using the new data

#slide 30
pred <- res$fitted.values # extract the point estimates
xlim <- c(min(sales_data$temp), max(sales_data$temp)) #range of x
ylim <- c(min(c(sales_data$sales, pred)), max(c(sales_data$sales, pred))) #range of y
par(oma=c(3,3,3,3), mfrow=c(1,1))
plot(sales_data$temp, sales_data$sales, xlim=xlim, ylim=ylim, xlab="", ylab="")
par(new=T) 
plot(sales_data$temp, pred, pch=20, col="red", xlim=xlim, ylim=ylim, xlab="", ylab="")
mtext("Temperature", line=3, cex=2, side=1)
mtext("Sales", line=3, cex=2, side=2)

#slide 32
res2 <- lm(sales~day+temp-1+day*temp, data=sales_data) # a model with interactions
summary(res2)

#slide 33
pred2 <- res2$fitted.values
par(oma=c(3,3,3,3), mfrow=c(1,1))
plot(sales_data$temp, sales_data$sales, xlim=xlim, ylim=ylim, xlab="", ylab="")
par(new=T)
plot(sales_data$temp, pred, pch=20, col="red", xlim=xlim, ylim=ylim, xlab="", ylab="")
par(new=T)
plot(sales_data$temp, pred2, pch=20, col="blue", cex=0.5, xlim=xlim, ylim=ylim, xlab="", ylab="")
mtext("Temperature", line=3, cex=2, side=1)
mtext("Sales", line=3, cex=2, side=2)

#slide 34
logLik(res) # log-likelihood
logLik(res2)

#slide 38
AIC(res) # AIC
AIC(res2)



#Part 3
#slide 43, 44
data <- read.csv("/Users/shin-ichironakayama/Dropbox/r_practice/standardization.csv") #rewrite the path!
head(data)
tail(data)
dim(data)
plot(data)

#slide 45
# calculate spatial averaged CPUE
calc_nominal <- function(data)
{
  year <- unique(data$year) # unique(): remove overlapped elements from a vector
  nominal_cpue <- c() # make an empty vector
  for(i in year)
  {
    data_temp <- data[data$year==i,] # pick up the data in year i
    nominal_cpue <- c(nominal_cpue, mean(data_temp$cpue)) # calculate average
  }
  nominal_cpue/mean(nominal_cpue) # standardization
}
nominal_cpue <- calc_nominal(data)
nominal_cpue

#slide 46
year <- c(1990:2005)
par(mfrow=c(1,1), oma=c(3,3,3,3))
ylim <- c(0, 2.5)
plot(year, nominal_cpue, type="o", pch=20, xlab="", ylab="", ylim=ylim, col="black")
mtext("Year", side=1, line=3, cex=2)
mtext("CPUE (scaled)", side=2, line=3, cex=2)

#slide 52
#model1
model1 <- lm(log(cpue+0.1) ~ as.factor(year) + as.factor(lon%/%10) + as.factor(lat%/%10) - 1, data=data) # categorical lon and lat
summary(model1) 
#model2
model2 <- lm(log(cpue+0.1) ~ as.factor(year) + scale(lon) + scale(lat) + I(scale(lon)^2) + I(scale(lat)^2) - 1, data=data) # continuous lon and lat
summary(model2) 

#slide 53
lon <- c(132:194)
lat <- c(19:49)
dummy <- expand.grid(lon, lat, year) # make a data with all combinations of lon, lat, and year
names(dummy) <- c("lon", "lat", "year") # change the name according to the original data

#slide 54
#calculate standardized CPUE under a certain model
calc_std_cpue <- function(model, dummy)
{
  predict <- predict(model, newdata=dummy, se.fit=T) # predict using dummy data
  log_cpue <- predict$fit # extract the value of log(CPUE)
  std_error <- predict$se.fit # extract standard errors of the estimates 
  cpue <- exp(log_cpue + std_error^2/2) # convert to exp scale
  data <- cbind(dummy, cpue) # attach the new CPUE in the column direction
  calc_nominal(data) # calculate spatial average
}
cpue1 <- calc_std_cpue(model1, dummy)
cpue2 <- calc_std_cpue(model2, dummy)
cpue1
cpue2

#slide 56
par(mfrow=c(1,1), oma=c(3,3,3,3))
ylim <- c(0, 2.5)
plot(year, cpue1, type="o", pch=20, xlab="", ylab="", ylim=ylim, col="red")
par(new=T)
plot(year, cpue2, type="o", pch=20, xlab="", ylab="", ylim=ylim, col="blue")
par(new=T)
plot(year, nominal_cpue, type="o", pch=20, xlab="", ylab="", ylim=ylim, col="black")
mtext("Year", side=1, line=3, cex=2)
mtext("CPUE (scaled)", side=2, line=3, cex=2)

#slide 59
# generate a bootstrapped data from a set of fits and residuals
boot_sample <- function(data, fit, resid_list)
{
  resid <- sample(resid_list, length(fit), replace=T) # resampling of residuals
  log_cpue <- fit + resid # make log(CPUE)
  cpue <- exp(log_cpue) # convert to exp scale
  data$cpue <- cpue # replace 
  data
}

#slide 60
fit1 <- model1$fitted.values # extract fit 
fit2 <- model2$fitted.values
resid_list1 <- model1$residuals # extract residuals
resid_list2 <- model2$residuals
cpue_boot1 <- c()
cpue_boot2 <- c()
for(i in 1:1000) # 1000 bootstraps
{
  print(i) # print current number
  d1 <- boot_sample(data, fit1, resid_list1) # make bootstrap data for model 1
  d2 <- boot_sample(data, fit2, resid_list2) # make bootstrap data for model 2
  model_temp1 <- lm(log(cpue+0.1) ~ as.factor(year) + as.factor(lon%/%10) + as.factor(lat%/%10) - 1, data=d1) 
  model_temp2 <- lm(log(cpue+0.1) ~ as.factor(year) + scale(lon) + scale(lat) + I(scale(lon)^2) + I(scale(lat)^2) - 1, data=d2) 
  cpue_boot1 <- cbind(cpue_boot1, calc_std_cpue(model_temp1, dummy))
  cpue_boot2 <- cbind(cpue_boot2, calc_std_cpue(model_temp2, dummy))
}

#slide 61
upper1 <- c()
lower1 <- c()
upper2 <- c()
lower2 <- c()
for(i in 1:nrow(cpue_boot1))
{
  upper1 <- c(upper1, quantile(cpue_boot1[i,], 0.975)) # 97.5 percentile
  lower1 <- c(lower1, quantile(cpue_boot1[i,], 0.025)) # 0.25 percentile
  upper2 <- c(upper2, quantile(cpue_boot2[i,], 0.975))
  lower2 <- c(lower2, quantile(cpue_boot2[i,], 0.025))
}
x_list <- c(year, rev(year))
y1_list <- c(lower1, rev(upper1))
y2_list <- c(lower2, rev(upper2))
par(mfrow=c(1,1), oma=c(3,3,3,3))
ylim <- c(0, 2.5)
plot(year, cpue1, type="o", pch=20, xlab="", ylab="", ylim=ylim, col="red")
polygon(x_list, y1_list, col=rgb(1, 0, 0, alpha=0.2), border=F)
par(new=T)
plot(year, cpue1, type="o", pch=20, xlab="", ylab="", ylim=ylim, col="red")
par(new=T)
plot(year, cpue2, type="o", pch=20, xlab="", ylab="", ylim=ylim, col="blue")
polygon(x_list, y2_list, col=rgb(0, 0, 1, alpha=0.2), border=F)
par(new=T)
plot(year, cpue2, type="o", pch=20, xlab="", ylab="", ylim=ylim, col="blue")
par(new=T)
plot(year, nominal_cpue, type="o", pch=20, xlab="", ylab="", ylim=ylim, col="black")
mtext("Year", side=1, line=3, cex=2)
mtext("CPUE (scaled)", side=2, line=3, cex=2)

#slide 62
#model selection
AIC(model1)
AIC(model2)
