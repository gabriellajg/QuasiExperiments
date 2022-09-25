# 1. Segmented regression for ITS designs
# Source: "Global Health Research: Design + Methods". By Dr. Eric Green

dat <- data.frame(month=seq(from=1, to=24), post=c(rep(0, 12), rep(1, 12)),
                  monthpost=c(rep(0, 12), seq(from=1, to=12)), 
                  outcome=c(50, 60, 45, 47, 50, 55, 53, 51, 48, 50, 47, 48,
                            55, 60, 63, 67, 65, 69, 72, 68, 74, 71, 76, 70))
lmFit1 <- lm(outcome ~ month + post + monthpost, data=dat) 
summary(lmFit1)

# regression diagnostics 
par(mfrow = c(2,2))
plot(lmFit1)
par(mfrow=c(1,1))
# nothing abnormal

########################################################
# Line plot
plot(dat$month, dat$outcome, type='b')
abline(v=12, lty = 2)

library(ggplot2)
p1 <- ggplot(data = dat, mapping = aes(x = month, y = outcome)) +
  geom_line(color = 'red') + theme_bw() + 
  geom_vline(xintercept = 12.5, color = 'blue', linetype = "longdash")
p1

# Overplot prediction
dat$predicted = predict(lmFit1, data = dat)
p1 + 
  geom_line(data = dat[dat$post==0,], aes(x=month, y=predicted), color="green") +
  geom_line(data = dat[dat$post==1,], aes(x=month, y=predicted), color="green")


########################################################
# Polynomial term
# Line plot
lmFit2 <- lm(outcome ~ month + I(month^2) + post + monthpost + I(monthpost^2), data=dat) 
summary(lmFit2)

# Overplot prediction
dat$predicted2 = predict(lmFit2, data = dat)
p1 + 
  geom_line(data = dat[dat$post==0,], aes(x=month, y=predicted2), color="green") +
  geom_line(data = dat[dat$post==1,], aes(x=month, y=predicted2), color="green")

########################################################
# ACF and PACF of residuals
par(mfrow=c(2,2), mar = c(2, 4, 4, 4))
acf(dat$outcome, lag.max = 23) # for Y
pacf(dat$outcome, lag.max = 23) # for Y
acf(resid(lmFit1)) # for e
pacf(resid(lmFit1)) # for e

library(forecast)
checkresiduals(lmFit1)

# the Breusch-Godfrey test for jointly testing up to 8th order autocorrelation.

# The residual plot shows some changing variation over time, but is not remarkable. 

# The histogram shows that the residuals seem to be slightly skewed, which may also affect the standard errors of the residuals.

# The autocorrelation plot shows no significant spike beyond the dashed blue line. Even up to lag 8, there is not quite enough evidence for the Breusch-Godfrey to be significant at the 5% level. The autocorrelations are not particularly large, and will be unlikely to have any noticeable impact on the forecasts or the prediction intervals.

# We are good with segmented regression...
# But for illustration purposes...

########################################################
# 2. Regression with HAC correction for standard errors

library(sandwich)
vcov(lmFit1) # original
vcovHAC(lmFit1) # corrected

# HAC
library(lmtest)
summary(lmFit1)
coeftest(lmFit1, vcov = vcovHAC(lmFit1))

