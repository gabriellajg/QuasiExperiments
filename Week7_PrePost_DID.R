# Simple DID estimator 
# https://dss.princeton.edu/training/DID101R.pdf

# Getting sample data.

library(foreign)
mydata = read.dta("http://dss.princeton.edu/training/Panel101.dta")

head(mydata)

# scatterplot by country
library(car)
scatterplot(y~year|country, boxplots=FALSE, smooth=TRUE, data = mydata)

# Create a dummy variable to indicate the time when the treatment started. 
# Let's  assume that treatment started in 1994. In this case, years before 1994 will have a value of 0 and 1994+ a 1. If you already have this skip this step.

mydata$time = ifelse(mydata$year >= 1994, 1, 0)

# Create a dummy variable to identify the group exposed to the treatment. In this example let's assumed that countries with code 5,6, and 7 were treated (=1). 
# Countries 1-4 were not treated (=0). If you already have this skip this step.

mydata$treated = ifelse(mydata$country == "E" |
                          mydata$country == "F" |
                          mydata$country == "G", 1, 0)

# Create an interaction by multiplying time and treated. We will call this interaction ‘did’.

mydata$did = mydata$time * mydata$treated

# Estimating the DID estimator
didreg = lm(y ~ treated + time + did, data = mydata)
summary(didreg)

# The coefficient for ‘did’ is the differences-in-differences estimator. 
# The effect is significant at 10% with the treatment having a negative effect.



########################################################

# The "PoEdata"" package loads into R the data sets that accompany Principles of Econometrics 4e, by Carter Hill, William Griffiths, and Guay Lim.
# by Dr. Constantin Colonescu
# https://github.com/ccolonescu/PoEdata
# https://bookdown.org/ccolonescu/RPoE4/indvars.html#the-difference-in-differences-estimator

library(devtools)
#devtools::install_github('ccolonescu/PoEdata')
library(PoEdata)
library(stargazer)
data("njmin3", package="PoEdata")
?njmin3 # minimum wage example 

mod1 <- lm(fte~nj*d, data=njmin3)
mod2 <- lm(fte~nj*d+
             kfc+roys+wendys+co_owned, data=njmin3)
mod3 <- lm(fte~nj*d+
             kfc+roys+wendys+co_owned+
             southj+centralj+pa1, data=njmin3)

stargazer::stargazer(mod1, mod2, mod3, 
                     type = 'text', model.names = FALSE, 
                     header=FALSE, keep.stat="n",digits=2, 
                     column.labels = c('DID', 'DIDw/Cov', 'DIDw/All'))


# Write down the regression equation for model 2 and interpret the regression coefficients in the context. 
