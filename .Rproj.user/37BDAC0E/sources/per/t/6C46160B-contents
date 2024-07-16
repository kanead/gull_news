# GAM analysis of newspaper articles over time 
# 4th March 2024

# clear everything
rm(list=ls())
graphics.off()

# load libraries
library(tidyverse)
library(mgcv)
library(gratia)
library(DHARMa)

# load the data
mydata <- read_csv("Combined_wDoubles_Nov15.csv")

# check it out
head(mydata)
glimpse(mydata)
head(mydata$Date)

# make sure country is set as a factor 
mydata$Country <- as.factor(mydata$Country)

# transform the date column into a date class
mydata$Date <- as.POSIXct(mydata$Date, format = "%d/%m/%Y")
head(mydata$Date)

# transform the dates to use in the model 
# mydata <- transform(
#   mydata,
#   ndate = as.numeric(Date),
#   nyear  = as.numeric(format(Date, '%Y')),
#   nmonth = as.numeric(format(Date, '%m')),
#   doy    = as.numeric(format(Date, '%j'))
# )

# check it worked 
# mydata$doy
# head(mydata)

# group the data by year, month and country, then count how many occurrences per group
# check this manually if it makes sense 
count_data <- mydata %>%
  group_by(Year, Month, Country) %>%
  summarise(count = n()) %>% arrange(Country, Year, Month)
count_data

# account for year of the gull
count_data <-
  count_data %>% mutate(special = if_else(Year == "2015", "yes", "no"))

# account for larger countries producing more articles by adding the population size
count_data <- count_data %>%
  mutate(
    population = case_when(
      Country == "England"  ~ 53,
      Country == "Ireland"  ~ 5.1,
      Country == "Northern Ireland"  ~ 1.8,
      Country == "Scotland"  ~ 5.3,
      Country == "Wales"  ~ 3
    )
  )

# fit the GAM
# two interactions, country with year and with month 
gam_mod1 <-
  gam(
    count ~ s(Year, k = 30, by = Country) +
      s(Month, k = 12, bs = "cc", by = Country) +
      Country + special,
    # potentially need to control for larger countries producing 
   #  more articles 
    offset = log(population), 
    family = nb,
    # other options for the family argument 
    # family = quasipoisson(),
    # family = poisson,
    method = "REML",
    data = count_data
  )

# check how the model went 
summary(gam_mod1)
plot(gam_mod1)
gam.check(gam_mod1)

# plot the model outputs on the original scale by taking the exponential 
# negative binomial family (nb) works on the log scale 
plot(gam_mod1, residuals = TRUE, pch=1, cex=1, seWithMean = TRUE,
     shift = coef(gam_mod1)[1],
     trans = exp)

# any autocorrelation?
acf(resid(gam_mod1), lag.max = 36, main = "ACF")
pacf(resid(gam_mod1), lag.max = 36, main = "pACF")


# fit the GAM
# one interaction with month but not with year
gam_mod2 <-
  gam(
    count ~ s(Year, k = 30) +
      s(Month, k = 12, bs = "cc", by = Country) +
      Country + special,
    # potentially need to control for larger countries producing 
    #  more articles 
    offset = log(population), 
    family = nb,
    # other options for the family argument 
    # family = quasipoisson(),
    # family = poisson,
    method = "REML",
    data = count_data
  )


# fit the GAM
# one interaction with year but not with month
gam_mod3 <-
  gam(
    count ~ s(Year, k = 30,  by = Country) +
      s(Month, k = 12, bs = "cc") +
      Country + special,
    # potentially need to control for larger countries producing 
    #  more articles 
    offset = log(population), 
    family = nb,
    # other options for the family argument 
    # family = quasipoisson(),
    # family = poisson,
    method = "REML",
    data = count_data
  )


# fit the GAM
# no interactions with country
gam_mod4 <-
  gam(
    count ~ s(Year, k = 30) +
      s(Month, k = 12, bs = "cc") +
      Country + special,
    # potentially need to control for larger countries producing 
    #  more articles 
    offset = log(population), 
    family = nb,
    # other options for the family argument 
    # family = quasipoisson(),
    # family = poisson,
    method = "REML",
    data = count_data
  )

# compare models
# gam_mod3 comes out as the best
AIC(gam_mod1, gam_mod2, gam_mod3, gam_mod4)

# take a look at it 
# check how the model went 
summary(gam_mod3)
plot(gam_mod3)
gam.check(gam_mod3)

# plot the model outputs on the original scale by taking the exponential 
# negative binomial family (nb) works on the log scale 
# this looks weird because of the huge uncertainty in earlier years when 
# we don't have much data to go on 
plot(gam_mod3, residuals = TRUE, pch=1, cex=1, seWithMean = TRUE,
     shift = coef(gam_mod3)[1],
     trans = exp, ylim = c(0,5))

# any autocorrelation?
acf(resid(gam_mod3), lag.max = 36, main = "ACF")
pacf(resid(gam_mod3), lag.max = 36, main = "pACF")

# another diagnostic test
res = simulateResiduals(gam_mod3)
plot(res)

# fit a more complex GAM
# interactions between year and month and country
gam_mod5 <-
  gam(
    count ~ s(Year, k = 30, by = Country) +
      s(Month, k = 12, bs = "cc" , by = Country) +
      ti(Month, Year, by = Country) +
      Country + special,
    # potentially need to control for larger countries producing 
    #  more articles 
    offset = log(population), 
    family = nb,
    # other options for the family argument 
    # family = quasipoisson(),
    # family = poisson,
    method = "REML",
    data = count_data
  )

# take a look at it 
# check how the model went 
summary(gam_mod5)
plot(gam_mod5)
gam.check(gam_mod5)

AIC(gam_mod1, gam_mod2, gam_mod3, gam_mod4, gam_mod5)


# subset to after some point 
count_data90s <- count_data %>% filter(Year > 1999)

# fit the GAM
# one interaction with year but not with month
gam_mod6 <-
  gam(
    count ~ s(Year,  by = Country, k = 20) +
      s(Month, bs = "cc") +
      Country + special,
    # potentially need to control for larger countries producing 
    #  more articles 
    offset = log(population), 
    family = nb,
    # other options for the family argument 
    # family = quasipoisson(),
    # family = poisson,
    method = "REML",
    data = count_data90s
  )

# take a look at it 
# check how the model went 
summary(gam_mod6)
plot(gam_mod6)
gam.check(gam_mod6)

# plot the model outputs on the original scale by taking the exponential 
# negative binomial family (nb) works on the log scale 
plot(gam_mod6, residuals = TRUE, pch=1, cex=1, seWithMean = TRUE,
     shift = coef(gam_mod6)[1],
     trans = exp)

# any autocorrelation?
acf(resid(gam_mod1), lag.max = 36, main = "ACF")
pacf(resid(gam_mod1), lag.max = 36, main = "pACF")


# playing around with gamm 
# gamm_mod <- gamm(
#   count ~ s(Year, k = 30) +
#     s(nmonth, k = 12, bs = "cc", by = Country) +
#     Country,
#   family = quasipoisson(),
#   method = 'REML',
#   correlation = corARMA(form = ~ 1 | Year, p = 1),
#  offset = log(population),
#   data = count_data
# )
# 
# summary(gamm_mod$gam)
# gam.check(gamm_mod$gam)
# plot(gamm_mod$gam)
# 
# layout(matrix(1:2, ncol = 2))
# acf(resid(gamm_mod$lme), lag.max = 36, main = "ACF")
# pacf(resid(gamm_mod$lme), lag.max = 36, main = "pACF")
# layout(1)
# 
# res = simulateResiduals(gamm_mod$gam)
# plot(res)

