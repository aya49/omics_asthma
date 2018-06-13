# https://www.zoology.ubc.ca/~schluter/R/fit-model/

library("lm")
library("visreg") 
library("emmeans") #get model predicted means and tukey tests
library("car")

library("lmerTest") # for lmer() in lme4 package
library("nlme")     # for lme() in nlme package

options(a.action = na.exclude) # keep track of cases having missing values, in which case the residuals and predicted values will have NA's inserted for those cases. Otherwise R just drops missing cases.



## linear regression
z <- lm(response ~ explanatory, data = mydata)
z <- lm(y ~ x, data = mydata)

z <- lm(y ~ 1, data = mydata)               # slope=0
z <- lm(y ~ 1 + offset(1*x), data = mydata) # slope=1
z <- lm(y ~ 1 + offset(b*x), data = mydata) # slope=b (b must be a number)

z <- lm(y ~ x - 1, data = mydata)           # intercept=0
z <- lm(y ~ 0 + x, data = mydata)           # intercept=0

z <- lm(y ~ x + A + x:A, data = mydata)     # covariance models; different linear regression for each categorical variable group

z1 <- lm(y ~ x1 + x2, data = mydata)          # no interaction between x1 and x2
z2 <- lm(y ~ x1 + x2 + x1:x2, data = mydata)  # interaction term present
z2 <- lm(y ~ x1 * x2, data = mydata)          # interaction term present



coef(z)        # model coefficients (means, slopes, intercepts)
resid(z)       # residuals
predict(z)     # predicted valuess
predict(z, newdata = mynewdata) # used the model to predict values for new observations
fitted(z)      # predicted values
anova(z1, z2)  # compare fits of 2 models, "full" vs "reduced"
anova(z)       # ANOVA table (** terms tested sequentially **);  test the null hypothesis of zero slope with the ANOVA table


summary(z)            # parameter estimates and overall model fit; slope, intercept, SE's
confint(z, level=.95) # confidence intervals for parameters

# plots of residuals, normal quantiles, leverage
# "Leverage" calculates the influence that each data point has on the estimated parameters
# "Cook's distance" measures the effect of each data point on the predicted values for all the other data points
#   -> value > 1 is bad
plot(z)
hist(resid(z))        # histogram of residuals

# scatterplot + regression line
plot(y ~ x, data = mydata)
abline(z)
# or
ggplot(mydata, aes(y = y, x = x)) +
  geom_point(size = 2, col = "red") +
  geom_smooth(method = lm, se = FALSE) +
  theme(aspect.ratio = 0.80)

#scatterplot + 95% confidence bands
visreg(z, points.par = list(pch = 16, cex = 1.2, col = "red"))
# or
ggplot(x, aes(y = age, x = black)) +
  geom_point(size = 2, col = "red") +
  geom_smooth(method = lm, se = TRUE) +
  theme(aspect.ratio = 0.80)
# or
x.p <- predict(z, interval = "prediction")
x.p <- cbind.data.frame(x, x.p)
ggplot(x.p, aes(y = age, x = black)) +
  geom_point(size = 2, col = "red") +
  geom_smooth(method = lm, se = TRUE) +
  geom_line(aes(y = lwr), color = "red", linetype = "dashed") +
  geom_line(aes(y = upr), color = "red", linetype = "dashed") +
  theme(aspect.ratio = 0.80)




## single factor ANOVA

# order groups s.t. control level comes first
A = factor(A, levels=c("control","a","b","c"))

z <- lm(y ~ A, data = mydata)               # A=categorical factor levels; model=single factor ANOVA
head(model.matrix(z))

# stripchart
stripchart(y ~ A, vertical = TRUE, method = "jitter", pch = 16,
           col = "red", data = mydata)
# or
stripchart(fitted(z) ~ A, vertical = TRUE, add = TRUE, pch="------",
           method = "jitter", data = mydata)
# or + confint 95% for predicted means
visreg(z)                                           # basic plot
visreg(z, points.par = list(cex = 1.2, col = "red")) # with points options
visreg(z, whitespace = 0.4)
# or + predicted values
stripchart(y ~ A, vertical = TRUE, method = "jitter", pch = 16,
           col = "red", data = mydata)
yhat <- tapply(fitted(z), mydata$A, mean)
for(i in 1:length(yhat)){
  lines(rep(yhat[i], 2) ~ c(i-.2, i+.2))
}
# or + group means fitted by single factor ANOVA i.e. SE, confint, Tuckey btwn all mean pairs
# NOT SE's and confints calculated on each group separately;
# uses residual mean square from the model fitted to all the data. 
# generally result in smaller SE's and narrower confidence intervals. 
# valid if equal variances within different groups.
emmeans(z, "A", data = mydata)   # use the real name of your factor in place of "A"
# Tuckey test btwn all mean pairs
grpMeans <- emmeans(z, "A", data = mydata)
grpMeans
pairs(grpMeans) # table of Tukey pairwise comparisons between means
cld(grpMeans) # table depicting comparisons between means using compact-letter display (e.g, "1", "2")

# plot: check assumptions of single factor ANOVA
plot(z)        # residual plots, etc
hist(resid(z)) # histogram of residuals

# estimate parameters
# by default, intercept=mean of 1st group (set to control or reference group)
# other parameters estimate mean of each group - mean of 1st group
# ignore p values as they're invalid for hypothesis tests of differences (except in planned contrast)
summary(z)               # coefficients table with estimates
confint(z, level = 0.95) # conf. intervals for parameters

# test null hypothesis of equal group means with ANOVA tables
anova(z)




## multiple factors
z <- lm(y ~ A + B, data = mydata) # additive model, no interaction fitted
z <- lm(y ~ A * B, data = mydata) # main effects and interaction terms fitted
summary(z)                        # coefficients table
confint(z, level = 0.95)          # confidence intervals for parameters
anova(z)                          # A is tested before B; interaction tested last

# marginal fitting of terms (type III sum of squares); cars package
z <- lm(y ~ A * B, contrasts = c("contr.sum", "contr.poly"), data = mydata)
Anova(z, type = 3)

# plot
# When a subset of multiple variables is plotted, 
# the method visualizes the predicted values for the plotted variable(s) 
# while conditioning on the value of the other variable 
# (adjusted to the most common category for other factor)
# seeing only a subset of model only makes sense only when
# no interaction term is fitted.
# When all interaction terms are plotted (A and B) then predicted values are shown for all combination of levels
visreg(z, xvar = "B", whitespace = 0.4, 
       points.par = list(cex = 1.1, col = "red"))
visreg(z, xvar = "A", by = "B", whitespace = 0.4, 
       points.par = list(cex = 1.1, col = "red"))
visreg(z, xvar = "A", by = "B", whitespace = 0.5, overlay = TRUE, 
       band = FALSE, points.par = list(cex = 1.1))

# fitted group means
z <- lm(y ~ A + B, data = mydata) # additive model, no interaction
emmeans(z, "A", data = mydata)         # group means for A variable, averaged over levels of B
emmeans(z, c("A", "B"), data = mydata) # model fitted means* for all combinations of A and B groups

z <- lm(y ~ A * B, data = mydata) # main effects and interaction terms fitted
emmeans(z, c("A", "B"), data = mydata) # means for all combinations of A and B groups




## factor + continuous covariate; ancova or analysis of covariance
# note: order of terms in your formula affects the sums of squares and anova tests

# sequental fitting of terms
z <- lm(y ~ x + A, data = mydata) # no interaction term included, or
z <- lm(y ~ x * A, data = mydata) # interaction term present
summary(z)                        # coefficients table
confint(z, level = 0.95)          # confidence intervals for parameters
anova(z)                          # sequential fitting of terms

# check model assumptions
plot(z)        # residual plots, etc
hist(resid(z)) # histogram of residuals


# vs marginal fitting of terms ("type III sum of squares"); by overriding default contrasts
z <- lm(y ~ x * A, contrasts = c("contr.sum", "contr.poly"), data = mydata)
Anova(z, type = 3)

# plot
visreg(z, xvar = "A", whitespace = 0.4, 
       points.par = list(cex = 1.1, col = "red"))
visreg(z, xvar = "x", by = "A", 
       points.par = list(cex = 1.1, col = "red"))
visreg(z, xvar = "x", by = "A", whitespace = 0.4, overlay = TRUE, 
       band = FALSE, points.par = list(cex = 1.1))
# or add the separate regression lines for each group to the scatter plot yourself
plot(y ~ x, pch = as.numeric(A), data = mydata)
legend(locator(1), as.character(levels(A)),     
       pch = 1:length(levels(A)))              # click on plot to place legend
groups <- levels(A)                         # stores group names
for(i in 1:length(groups)){
  xi <- x[A==groups[i]]                   # grabs x-values for group i
  yhati <- fitted(z)[A==groups[i]]        # grabs yhat's for group i
  lines(yhati[order(xi)] ~ xi[order(xi)]) # connects the dots
}








## linear mixed effects model
## note: frequent error is to use the same labels to id different sampling units in different levels
##         e.g. same number codes 1-5 used to label 5 subplots of every plot;
##              instead use codes a1-a5 to label subplots from plot "a", b1-b5 to label subplots from plot "b" etc.

# lmer wants formula that includes fixed + random effects in ()
# e.g. single random effect B and only an intercept for a fixed effect
# "1" means fixed intercet is fitted; 
# "1|B" means a separate random intercept needs to be fitted for each level in B
z <- lmer(y ~ 1 + (1|B), data = mydata)

# lme wants 2 formulas, fixed effects + random effects
# e.g. single random effect B and only an intercept for a fixed effect
z <- lme(y ~ 1, random = ~ 1|B, data = mydata) # avoid mydata$y and mydata$B notation



summary(z)                # variances for random effects, fit metrics
plot(z)                   # plot of residuals against predicted values
VarCorr(z)                # variance components for random effects  
confint(z)                # lmer: conf. intervals for fixed effects and variances 
intervals(z)              # lme: conf. intervals for fixed effects and variances

resid(z)                  # residuals
fitted(z)                 # best linear unbiased predictors (BLUPs)

anova(z, type = 1)        # lmer: test fixed effects sequentially (Type I SS)
anova(z, type = 3)        # lmer: as above but using Type III Sums of Squares
anova(z)                  # lme: test fixed effects sequentially (Type I SS)





## e.g. 1 random factor; repeated measurements of individuals
z <- lmer(y ~ 1 + (1|B), data = mydata)         # lme4/lmerTest package
z <- lme (y ~ 1, random = ~ 1|B, data = mydata) # nlme package
# application: repeatability of a measurement is calculated as
# r = var_among / (var_among + var_within)
# var_among = var among the means of levels (individuals, in this case)
# var_within = var among repeat measurements within each level

# estimate parameters
summary(z)
VarCorr(z)   # estimated variance among groups is shown as the variance associated with the random intercepts
#              variance associated with the residual is the estimate of the within-group variance
confint(z)   # lmer: 95% confidence intervals fixed effects
intervals(z) # lme: 95% confidence intervals fixed effects



plot(z) # lot of the residual values against fitted values

# confirm what is happening, plot the original data 
# and superimpose the lme fitted (predicted) values
stripchart(y ~ B, vertical = TRUE, pch = 1)
stripchart(fitted(z) ~ B, vertical = TRUE, add = TRUE, pch = "---")
# lme is predicting that the "true" trait value for individuals 
# having extreme measurements is likely to be closer, on average, 
# to the grand mean than is the average of the repeat measurements. 
# This effect diminishes the more measurements you have taken of each individual. 





## e.g. 2 nested random factors; nested sampling design, no treatment fixed effects
## measure y of multiple offspring per mom (A) 
## and multiple moms are mated to each dad (B)
## no treatment fixed effects; so anova won't do anything
z <- lmer(y ~ 1 + (1|B/A), data = mydata)        # lme4/lmerTest
z <- lme(y ~ 1, random = ~ 1|B/A, data = mydata) # nlme
# B/A nesting: woodlot and transect within woodlot
# random intercepts are fitted to woodlots and to transects from each woodlot




## e.g. 1 fixed, 1 random factor; randomized block and subjects-by-treatment repeated measures
## RCB (randomized complete block) is a paired design for 2+ treatments
## i.e. each treatment (A) applied once in every block (B) (treatments are independant); so block is modeled as random effects
z <- lmer(y ~ A + (1|B), data = mydata)         # lme4/lmerTest package
z <- lme(y ~ A, random = ~ 1|B, data = mydata)  # nlme package
# similar to single random effect except fixed part has effect A explanetory variable

# how the mean y changes between treatment levels of A separately for each level of B
# parallel lines = lack of interaction btwn A, B
interaction.plot(A, B, y)

# if A is fixed (treatment) effect and B is random effect
visreg(z, xvar = "A")
visreg(z, xvar = "A", by = "B", scales=list(rot = 90))

summary(z)
VarCorr(z)
confint(z)   # lmer
intervals(z) # lme

# model-based estimates of A means
emmeans(z, "A", data = mydata)                            # uses kenward-roger degrees of freedom method
emmeans(z, "A", data = mydata, lmer.df = "satterthwaite") # uses satterthwaite degrees of freedom method

# test the fixed effects (the grand mean and the treatment A)
anova(z, type = 1)   # lmer: test fixed effects sequentially (Type I SS)
anova(z, type = 3)   # lmer: as above but using Type III Sums of Squares
anova(z)             # lme: test fixed effects sequentially (Type I SS)




## e.g. 1 fixed, 1 random factor; factorial design
## RBC except each combination of levels of fixed and random effect may be replicated
## 2 factor ANOVA is more complex when a factor is random; 
## because random sampling of groups 
## add extra sampling error = noise to measurements 
## of differences btwn group means for fixed factors 
## that interact with random factor

# intearction btwn fixed A and random B is a random effect
z <- lmer(y ~ A + (1|B/A), data = mydata)        # lme4/lmerTest package
z <- lme(y ~ A, random = ~ 1|B/A, data = mydata) # nlme package

# no interaction term
z <- lmer(y ~ A + (1|B), data = mydata)          # lme4/lmerTest package
z <- lme(y ~ A, random = ~ 1|B, data = mydata)   # nlme package

# plot residuals against fitted values
plot(z)

# test fixed effects
anova(z, type = 1)   # lmer: test fixed effects sequentially (Type I SS)
anova(z, type = 3)   # lmer: as above but using Type III Sums of Squares
anova(z)             # lme: test fixed effects sequentially (Type I SS)

# obtain model-based estimates of the fixed factor (A) means
emmeans(z, "A", data = mydata)


###### TBC examples #############
