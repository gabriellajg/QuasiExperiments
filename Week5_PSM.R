# read in data
# propensity score first - logistic
# add higher interactive terms
# other methods built in - boosting
# packages 
# examine overlap
# check the recal precision? bootstrap?
# explain matching strategy
# arguments in the r packages

# PSM designs with logistic regression 
# https://github.com/yishilin14/causal_playground/blob/master/causal_inference_intro2_propensity_score_matching.Rmd

library(Matching)
library(MatchIt)
library(optmatch)
library(weights)
# https://r.iq.harvard.edu/docs/matchit/2.4-20/matchit.pdf

data(lalonde)
View(lalonde)
dim(lalonde)
names(lalonde)

# outcome: re78
# true causal effect: 1800

# 1. Selection of covariates in X
fm = treat ~ age + educ + black + hispan + married + I(re74/1000) + I(re75/1000)
fm = treat ~ age + I(age^2) + I(age^3) + educ + black + hispan + married + I(re74/1000) + I(re75/1000)
fm = treat ~ age + I(age^2) + I(age^3) + educ + I(educ^2) + black + hispan + married + I(re74/1000) + I(re75/1000)

# 2. Calculation of propensity scores (p-scores)
pscore <- glm(fm, data = lalonde, family = 'binomial')
pscore$fitted.values
hist(pscore$fitted.values[lalonde$treat==0],xlim=c(0,1))
hist(pscore$fitted.values[lalonde$treat==1],xlim=c(0,1))
lalonde$pscore = pscore$fitted.values

# One-step
set.seed(42)
m.out <- matchit(data = lalonde,
                 formula = fm,
                 distance = "logit",
                 method = "nearest",
                 replace = FALSE,
                 caliper = 0.2, 
                 discard = 'both')
# m.out$match.matrix
# m.out$distance
# plot(m.out$distance, pscore$fitted.values) # same 
# method: exact, subclass, optimal, full, cem
# distance: pscore

# 3. Matching based on p-scores
# 4. Balance check
plot(summary(m.out,standardize = T))
plot(m.out, type = "hist", interactive = F)
plot(m.out, type = "QQ", interactive = F, which.xs = c("age", "I(re74/1000)", "I(re75/1000)"))
plot(m.out, type = "QQ", interactive = F)
summary(m.out, standardize = T)$sum.matched
round(summary(m.out, standardize = T)$sum.matched, 2)

# 5. Estimation of treatment effect
m.data <- match.data(m.out)
# Direct compare
res <- wtd.t.test(m.data$re78[m.data$treat == 1],
                  m.data$re78[m.data$treat == 0],
                  weight = m.data$weights[m.data$treat == 1],
                  weighty = m.data$weights[m.data$treat == 0])
print(res)
mu <- res$additional[1]
std <- res$additional[4]
cat("Confidence interval: ", sapply(qt(c(0.025, 0.975), coef(res)["df"]), function(x){return(mu+x*std)}), "\n")

res <- t.test(m.data$re78[m.data$treat == 1],
              m.data$re78[m.data$treat == 0])


# Fit - double robust - add covariates in the estimation stage
att.fml <- re78 ~ treat + age + educ + black + hispan + married + nodegree + re74 + re75
fit <- lm(att.fml, data = m.data, weights = m.data$weights)
summary(fit)
cat("Confidence interval: ", confint(fit, "treat", 0.95), "\n")


# 6. Post-hoc test for hidden bias

library(rbounds)
psm4 <- Match(Y = lalonde$re78, Tr = lalonde$treat, X = pscore$fitted.values, estimand = "ATT", M = 5, replace = TRUE)
summary(psm4)
psens(x = psm4, Gamma = 2, GammaInc = 0.1)
