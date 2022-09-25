library(weights)
library(survey)
library(twang)
library(CBPS)
library(cobalt)
library(jtools)
library(lmtest)
library(sandwich) #vcovCL
library(rbounds) #gamma
library(tidyr)
library(tidyverse)
library(janitor)


# load LaLonde from package Matching
data(LaLonde, package = 'CBPS')
?LaLonde

# outcome: re78

# Creating Employment Dummies: 1: employed; 0:unemployed 
LaLonde$un74 = ifelse(LaLonde$re74==0, 0, 1)
LaLonde$un75 = ifelse(LaLonde$re75==0, 0, 1)

# 1. Selection of covariates in X
fm1 = treat ~ age + educ + black + hisp + married + nodegr + un74 + un75 + re74 + re75
fm2 = treat ~ age + I(age^2) + I(age^3) + educ + black + hisp + married + nodegr + un74 + un75 + re74 + re75
fm3 = treat ~ age + I(age^2) + I(age^3) + educ + I(educ^2) + black + hisp + married + nodegr + un74 + un75 + re74 + re75

# 2. Calculation of propensity scores (p-scores)
# use glm() for logistic regression
pscore <- glm(fm1, data = LaLonde, family = 'binomial')
head(pscore$fitted.values)
hist(pscore$fitted.values[LaLonde$treat==0],xlim=c(0,1), main = 'Control')
hist(pscore$fitted.values[LaLonde$treat==1],xlim=c(0,1), main = 'Treated')
LaLonde$pscore = pscore$fitted.values
#View(LaLonde)
# try other fm?

# ------------------------------------------------
# weighting 
# ------------------------------------------------

# IPTW

# In the following code, the weight variable weightATE is created by using the ifelse function to obtain the inverse of the propensity score for treated units or the inverse of 1 minus the propensity score for control units

LaLonde$weightATE <- with(LaLonde, ifelse(treat==1, 1/pscore, 1/(1-pscore)))

#  A summary of the ATE weights for treated and untreated groups is obtained with the following: 
with(LaLonde, by(weightATE,treat,summary)) 

boxplot(LaLonde$weightATE ~ LaLonde$treat)

# the maximum weight for the ATE is 117.417
# Individuals with extreme weights for the ATE are those who are either very likely to participate in the treatment given their covariate values but did not or very unlikely to participate but did so.

#  perform weight truncation
# Truncation can be performed by assigning the weight at a cutoff percentile to observations with weights above the cutoff.

# The following code demonstrates weight truncation. More specifically, it uses the quantile function to calculate the weight at the 99th percentile and the ifelse function to assign this weight to any student whose weight exceeds the 99th percentile:

LaLonde$weightATETruncated <- with(LaLonde, ifelse(weightATE > quantile(weightATE, 0.99), quantile(weightATE, 0.99), weightATE))

# may truncate all above .05 quantiles if necessary

with(LaLonde, by(weightATETruncated,treat,summary)) 

boxplot(LaLonde$weightATETruncated ~ LaLonde$treat)


# 4. Balance check
# In this section, covariate balance evaluation is performed by comparing the standardized difference between the weighted means of the treated and untreated groups.

# The bal.stat function of the twang package (Ridgeway et al., 2013) is useful for covariate balance evaluation. 
#  If there are no sampling weights, then sampw=1.

covariateNames <-names(LaLonde)[3:14]
  
balanceTable <- bal.stat(LaLonde, vars= covariateNames,
                         treat.var = "treat",
                         w.all = LaLonde$weightATETruncated, 
                         get.ks=F,
                         sampw = 1,
                         estimand="ATE", multinom=F)
balanceTable <- balanceTable$results
round(balanceTable,3) 

# std.eff.sz quantifies effect size
# 0.2 small 0.5 medium 0.8 large
# most of them are medium to large - balance was not achieved 


# 5. Estimation of treatment effect
#  the final weights are divided by the mean of weights to make them sum to the sample size, which is a process known as normalization.

LaLonde$finalWeight <- LaLonde$weightATETruncated/mean(LaLonde$weightATETruncated)

#  Before the estimation can be performed, the surveyDesign object is created with the svydesign function to declare the names of the variables that contain cluster ids, strata ids, weights, and the data set: 

surveyDesign <- svydesign(ids=~1, weights=~finalWeight, data = LaLonde, nest=T) 

# Methods to obtain standard errors for propensity score weighted estimates include Taylor series linearization and resampling methods such as bootstrapping, jackknife, and balanced repeated replication (see review by Rodgers, 1999).

# To use bootstrap methods with the survey package, the surveyDesign object created above should be modified to include weights for each replication. The following code takes the surveyDesign object and adds weights for 1,000 bootstrapped samples:

set.seed(8)
surveyDesignBoot <- as.svrepdesign(surveyDesign, type=c("bootstrap"), replicates=1000) 

# First, the svyby function is used to apply the svymean function separately to treated and control units to obtain weighted outcome:

weightedMeans <- svyby(formula=~re78, by=~treat, design=surveyDesignBoot, FUN=svymean, covmat=TRUE)

weightedMeans


# The code to obtain the treatment effect with a regression model is shown as follows. 
# The model formula using R notation is re78~treat. This formula can be expanded to include any covariates and interaction effects of interest

outcomeModel <- svyglm(re78~treat,surveyDesignBoot)
summary(outcomeModel) # final model




# ------------------------------------------------
# Stratification 
# ------------------------------------------------

# The following code used the cut function to create five strata of approximately the same size based on the quintiles of the distribution of propensity scores for both treated and untreated groups. 
# The quintiles are obtained with the quantile function. The function levels is used to assign number labels from 1 to 5 to the strata, and then xtabs is used to display strata by treatment counts.

hist(LaLonde$pscore)
quantile(LaLonde$pscore, prob = seq(0, 1, 1/5))

LaLonde$subclass <- cut(x=LaLonde$pscore, breaks = quantile(LaLonde$pscore, prob = seq(0, 1, 1/5)), include.lowest=T)
levels(LaLonde$subclass) <- 1:length(levels(LaLonde$subclass))

ntable <- xtabs(~treat+subclass,LaLonde) 
ntable

surveyDesign <- svydesign(ids=~1, weights=~finalWeight, data = LaLonde, nest=T) 

# To use bootstrap methods with the survey package, the surveyDesign object created above should be modified to include weights for each replication. The following code takes the surveyDesign object and adds weights for 1,000 bootstrapped samples:

set.seed(8)
surveyDesignBoot <- as.svrepdesign(surveyDesign, type=c("bootstrap"), replicates=1000) 


# The following R code uses svyby to apply the svymean function:

head(LaLonde)

subclassMeans <- svyby(formula=~re78, by=~treat+subclass,
                       design=surveyDesignBoot, FUN=svymean, covmat=TRUE)

subclassMeans

# To obtain the ATE or ATT by pooling stratum-specific effects, svycontrast is used with the weights

# first, use ntable to obtain the weights
ntable
colSums(ntable)/sum(ntable)
ATTw = colSums(ntable)/sum(ntable)

# won't work 
# pooledEffects <- svycontrast(subclassMeans, contrasts = list(ATT=ATTw)) 

# ATTw needs to be specified for EVERY group within EVERY stratum
# in our case, 9 groups (no treatment in first stratum so 2*5-1 = 9)
# Control groups get negative weights and treatment group gets positive weights: 
ATTw2 <- c(-ATTw[1], # first stratum, just control
           -ATTw[2], ATTw[2],  # second stratum
           -ATTw[3], ATTw[3],  # third stratum
           -ATTw[4], ATTw[4],  # fourth stratum
           -ATTw[4], ATTw[5])  # fifth stratum
subclassMeans$ATTw2 <- ATTw2

# view ATTw2 here:
subclassMeans

pooledEffects <- svycontrast(subclassMeans, list(ATE=as.numeric(ATTw2))) 

pooledEffects

# ------------------------------------------------
# adjustment 
# ------------------------------------------------

adjustModel <- lm(re78~treat+pscore,LaLonde)
summ(adjustModel) 

# ------------------------------------------------
# Exercise: can you use fm2 to obtain propensity scores and see if that improves the balance and thus the estimate of treatment effect? 
# ------------------------------------------------

