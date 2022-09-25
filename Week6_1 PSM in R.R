library(Matching)
library(MatchIt)
library(optmatch)
library(weights)
library(cem)
library(tcltk2)
library(CBPS)
library(jtools)
library(Zelig)
library(cobalt)
library(lmtest)
library(sandwich) #vcovCL
library(rbounds) #gamma

# load LaLonde from package Matching
data(LaLonde, package = 'CBPS')
?LaLonde
lalonde = LaLonde

# data(lalonde, package = 'Matching')
#head(lalonde)

#View(lalonde)
dim(LaLonde)
names(LaLonde)
head(LaLonde)

# outcome: re78
# true causal effect: 1800

# Creating Unemployment Dummies
LaLonde$un74 = ifelse(LaLonde$re74==0, 0, 1)
LaLonde$un75 = ifelse(LaLonde$re75==0, 0, 1)

# 0. naive estimate 
reglm <- lm(re78 ~ treat, data = LaLonde)
summ(reglm)

reglm1 <- lm(re78 ~ treat + educ + age + black + hisp + married + nodegr + un74 + un75 + re74 + re75, data = LaLonde)
summ(reglm1)

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

# logistic regression model
summ(pscore)
  
# 3. Matching based on p-scores

# all in one function matchit()
# since matching process is random, set seed first
set.seed(42)
m.out <- matchit(data = LaLonde,
                 formula = fm1,
                 distance = "logit",
                 method = "nearest",
                 replace = TRUE,
                 caliper = 0.2, 
                 discard = 'both')
m.out
summary(m.out)
# m.out$match.matrix
# m.out$distance
# plot(m.out$distance, m.out$fitted.values) # same 
# method: exact, subclass, optimal, full, cem
# distance: pscore

# manual: https://imai.fas.harvard.edu/research/files/matchit.pdf
# starting from page 15
# Let's play with the arguments! 


# 4. Balance check
names(summary(m.out))
round(summary(m.out, standardize = T)$sum.all, 2)
round(summary(m.out, standardize = T)$sum.matched, 2)
plot(m.out, type = "hist", interactive = F)
plot(m.out, type = "jitter", interactive = F)
plot(summary(m.out,standardize = T), interactive = F)
plot(m.out, type = "QQ", interactive = F, which.xs = c("age", "educ", "re74", "re75"))

# find the problematic one
for (i in 1:length(covariatenames)){
  plot(m.out, type = "QQ", interactive = F, which.xs = covariatenames[i])
}

plot(m.out, type = "QQ", interactive = F, which.xs = covariatenames[-3])


love.plot(m.out, binary = "std")
bal.plot(m.out, var.name = "distance", which = "both",
         type = "histogram", mirror = TRUE)

# who matched to whom?
head(m.out$match.matrix, 10)

# 3b Optimal matching 
set.seed(101)
m.out.opt <- matchit(data = LaLonde,
                 formula = fm1,
                 distance = "logit",
                 method = "optimal",
                 ratio = 2)
m.out.opt
summary(m.out.opt)

# 5. Estimation of treatment effect
# extract matched data
m.data <- match.data(m.out)
m.data <- get_matches(m.out)

#Linear model without covariates
fit1 <- lm(re78 ~ treat, data = m.data)
summ(fit1)

#Cluster-robust standard errors
coeftest(fit1, vcov. = vcovCL, cluster = ~subclass)

#Linear model with covariates: double robust
fit2 <- lm(re78 ~ treat + age + educ + black + hisp + married + nodegr + un74 + un75 + re74 + re75, data = m.data)
summ(fit2)

#Cluster-robust standard errors
coeftest(fit2, vcov. = vcovCL, cluster = ~subclass)


# 6. Hidden bias analysis
# extract matched treatment and control units from m.data
psens(x = m.data$re78[m.data$treat==1], 
      y = m.data$re78[m.data$treat==0], 
      Gamma = 3, GammaInc=0.1)
# a gamma value with 1.4 or larger could lead to a change in ATT estimates

