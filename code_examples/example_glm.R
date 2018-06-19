# https://www.zoology.ubc.ca/~schluter/R/fit-model/


library("MASS")
library("visreg")
library("emmeans")


## generalized linear model


# If the variance of the error distribution in the data is 
# greater than that expected under the binomial and poisson distributions 
# ("overdispersion"), try: error distributions model overdispersion 
# for binary and count data, 
# and additionally provide an estimate of the dispersion parameter 
# (a value greater than one indicates overdispersion)
family = quasibinomial(link = "logit") # or
family = quasipoisson(link = "log")

# model a binary response variable, specify a binomial error distribution and the logit link function
# logistic regression
z <- glm(response ~ explanatory, family = binomial(link="logit"), data = mydata)
# model count data, specify the Poisson error distribution and the log link function
# log-linear regression
z <- glm(response ~ explanatory, family = poisson(link="log"), data = mydata)




# no plot(z); visualize by adding fitted curves or means to scatterplot
summary(z)                    # parameter estimates and overall model fit
coef(z)                       # model coefficients
resid(z)                      # deviance residuals;  not the same as ordinary residuals
#                               a measure of the goodness of fit of the model to each data point
predict(z)                    # predicted values on the transformed scale
predict(z, se.fit = TRUE)     # Includes SE's of predicted values
fitted(z)                     # predicted values on the original scale
anova(z, test = "Chisq")      # Analysis of deviance - sequential
anova(z1, z2, test = "Chisq") # compare fits of 2 models, "reduced" vs "full"
anova(z, test = "F")          # Use F test for gaussian, quasibinomial or quasipoisson link

# calculating likelihood based confidence intervals for parameters
confint(z, level = 0.95) # approximate 95% confidence intervals
dose.p(z, p = 0.50)      # LD50 for a dose-response curve




# view transformed scale results
visreg(z, xvar = "x") # fit on the transformed scale
visreg(z, xvar = "x", ylim = range(y), scale = "response")  # fit on the original scale
points(y ~ x)                                               # add the data points

# generate predicted values
plot(jitter(y, amount = 0.02) ~ x, data = mydata)
yhat <- fitted(z)
lines(yhat[order(x)] ~ x[order(x)])

## Approximate confidence bands
# calculate the upper and lower limits
zhat <- predict(z, se.fit = TRUE)       # result on logit or log scale
zupper <- zhat$fit + 1.96 * zhat$se.fit # replace 1.96 with 1 to plot standard error bands
zlower <- zhat$fit - 1.96 * zhat$se.fit
# convert to original scale
yupper <- exp(zupper)/(1 + exp(zupper)) # for logit link
ylower <- exp(zlower)/(1 + exp(zlower))
# or
yupper <- exp(zupper)                   # for log link
ylower <- exp(zlower)
# plot
lines(yupper[order(x)] ~ x[order(x)], lty = 2)
lines(ylower[order(x)] ~ x[order(x)], lty = 2)

# slope, intercept, on logit/log scale, SE, 
# pvalues based on normal approx so not accurate for small sample size
# instead, used log likelihood ratio test
summary(z)

# analysis of deviance table
# gets log-likelihood ratio tests of most important par
# model terms tested sequentially; results depend on order var entered into formula
anova(z, test = "Chisq")





## categorical explanatory variable: same as single-factor ANOVA but response is 2+ categories

# order groups s.t. control level comes first
A = factor(A, levels=c("control","a","b","c"))

z <- glm(y ~ A, family = binomial(link="logit"), data = mydata)
z <- glm(y ~ A, family = poisson(link="log"), data = mydata)

# enerate model-based fitted means for treatments
emmeans(z, c("A"), data = mydata)  # means of "A" treatment on transformed scale

# you'll need to transform the means and the lower and upper limits 
# of the confidence intervals to the original scale 
# using the inverse of the link function




## modeling contingency tables
















