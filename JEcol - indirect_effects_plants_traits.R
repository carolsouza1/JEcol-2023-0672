# IMPACT OF THE EFN TRAITS ON THE POTENTIAL INDIRECT EFFECTS BETWEEN PLANTS VIA DOMINANT AND SUBORDINATE ANT SPECIES #

# PACKAGES #
library(car)
library(bbmle)   ## To compare models in a table with AICs
library(nlme)     ## To run GLS models
library(piecewiseSEM)   ## To calculate partial residuals in models
library(report)
library(effectsize)
library(sjPlot)


data <- read.table("dataset_muller_planttraits.txt", header = TRUE, dec = ".")
str(data)

-------
# Tests to investigate the effect of two EFN traits 
# on potential indirect effects (promoter and receptor degree values) in plant species
-------

## ALL ANTS ##

# IMPACT OF EFNs ATTRIBUTES ON PROMOTER DEGREE


### Model 1 - lm with Gaussian distribution ###

mod1_promoter_normal <- glm(promoter_all_ants ~ scale(average_activenodes_perplot) + scale(Secretoryarea_perplantnode_mm2), data = data)
summary(mod1_promoter_normal)
par(mfrow = c(1, 1))
plot(mod1_promoter_normal)   ## observing the pattern of residuals and normality (see page 401 of Crawley's book!)

hist(data$promoter_all_ants)
residualPlot(mod1_promoter_normal)
shapiro.test(residuals(mod1_promoter_normal))
hist(residuals(mod1_promoter_normal))
vif(mod1_promoter_normal)  ## testing the assumption of absence of multicollinearity (see page 478 of Zuur's book!)

# Improving the model by excluding the species that most influences non-homogeneous variance pattern:
mod1_promoter_normal.updated <- update(mod1_promoter_normal, subset = (Species != "sp16"))
summary(mod1_promoter_normal.updated)
plot(mod1_promoter_normal.updated)   # not much improvement!

# Partial graphs:
par(mfrow = c(1, 2), pty = "s")
avPlots(mod1_promoter_normal, ylab = "Promoter degree (all ants)", pch = 16, col.lines = c("gray", "gray"), id = FALSE, main = "", grid = FALSE, layout = c(1, 1), las = 2, lwd = 2)  # With lines
avPlots(mod1_promoter_normal, ylab = "Promoter degree (all ants)", pch = 16, col.lines = c("gray", "gray"), id = FALSE, main = "", grid = FALSE, layout = c(1, 1), las = 2, lwd = 0)   # Without lines
par(mfrow = c(1, 1))

# Gamma and lognormal models did not work well when plotting the graphs
# Given that we don't have a normality problem but a variance homogeneity problem, 
# trying a strategy to model the heterogeneity using the gls() function as 
# implemented in Chapter 4 "Dealing with heterogeneity" - from Zuur et al 2009

# Traditional model
mod1_lm_promoter <- gls(promoter_all_ants ~ scale(average_activenodes_perplot) + scale(Secretoryarea_perplantnode_mm2), data = data)
summary(mod1_lm_promoter)

# Model 'modeling' the increase in variance (alternative 1)
vf1Fixed <- varFixed(~average_activenodes_perplot)
mod1_gls1_promoter <- gls(promoter_all_ants ~ scale(average_activenodes_perplot) + scale(Secretoryarea_perplantnode_mm2), weights = vf1Fixed, data = data)
summary(mod1_gls1_promoter)
plot(mod1_gls1_promoter$residuals ~ mod1_gls1_promoter$fitted)
plot(mod1_gls1_promoter$residuals ~ data$average_activenodes_perplot)
plot(mod1_gls1_promoter$residuals ~ data$Secretoryarea_perplantnode_mm2)


# Model 'modeling' the increase in variance (alternative 2)
vf1power <- varPower(form = ~average_activenodes_perplot)
mod1_gls2_promoter <- gls(promoter_all_ants ~ scale(average_activenodes_perplot) + scale(Secretoryarea_perplantnode_mm2), weights = vf1power, data = data)
summary(mod1_gls2_promoter)

# Model 'modeling' the increase in variance (alternative 3)
vf1Exp <- varExp(form = ~average_activenodes_perplot)
mod1_gls3_promoter <- gls(promoter_all_ants ~ scale(average_activenodes_perplot) + scale(Secretoryarea_perplantnode_mm2), weights = vf1Exp, data = data, na.action = na.omit)
summary(mod1_gls3_promoter)

# Model 'modeling' the increase in variance (alternative 4)
vf1ConstPower <- varConstPower(form = ~average_activenodes_perplot)
mod1_gls4_promoter <- gls(promoter_all_ants ~ scale(average_activenodes_perplot) + scale(Secretoryarea_perplantnode_mm2), weights = vf1ConstPower, data = data, na.action = na.omit)
summary(mod1_gls4_promoter)

# Checking which model fits better
anova (mod1_lm_promoter,mod1_gls1_promoter,mod1_gls2_promoter,mod1_gls3_promoter,mod1_gls4_promoter)
ICtab(mod1_lm_promoter,mod1_gls1_promoter,mod1_gls2_promoter,mod1_gls3_promoter,mod1_gls4_promoter, type=c("AICc"),weights=TRUE, delta=TRUE,sort=TRUE,nobs=21)

# Alternative 1 - best model

#### Defining the variance structure of the data:
## Modeling now to check if the two fixed factors should remain in the model

mod1_gls1_promoter

vf1Fixed <- varFixed(~average_activenodes_perplot)

mod1_gls1_promoter <- gls(promoter_all_ants ~ scale(average_activenodes_perplot) + scale(Secretoryarea_perplantnode_mm2), weights = vf1Fixed, data = data)
summary(mod1_gls1_promoter)
mod2_gls1_promoter <- gls(promoter_all_ants ~ scale(average_activenodes_perplot), weights = vf1Fixed, data = data)
mod3_gls1_promoter <- gls(promoter_all_ants ~ scale(Secretoryarea_perplantnode_mm2), weights = vf1Fixed, data = data)
mod0_gls1_promoter <- gls(promoter_all_ants ~ 1, weights = vf1Fixed, data = data)  ## Null model

ICtab(mod0_gls1_promoter, mod1_gls1_promoter, mod2_gls1_promoter, mod3_gls1_promoter, type = c("AICc"), weights = TRUE, delta = TRUE, sort = TRUE, nobs = 21)


# Constructing "partial residual graphs" for each factor:
# https://www.rdocumentation.org/packages/piecewiseSEM/versions/1.2.1/topics/partial.resid

# Best model above: mod1_gls1_promoter

mod1_gls1_promoter_final <- gls(promoter_all_ants ~ average_activenodes_perplot + Secretoryarea_perplantnode_mm2, weights = vf1Fixed, data = data)
summary(mod1_gls1_promoter_final)
tab_model(mod1_gls1_promoter_final)
str(mod1_gls1_promoter_final)
effectsize(mod1_gls1_promoter_final)




YX1 <- partialResid(promoter_all_ants ~ average_activenodes_perplot, mod1_gls1_promoter_final, data)
plot(YX1$yresid ~ YX1$xresid, las = 1, pch = 19, cex = 1, xlim = c(-0.3, 1.7), ylim = c(-0.6, 1.15), cex.lab = 0.9, ylab = "", xlab = "")
text(YX1$yresid ~ YX1$xresid, labels = rownames(data), data = data, cex = 0.8, font = 2, pos = 2)
abline(lm(YX1$yresid ~ YX1$xresid), lwd = 2)

YX2 <- partialResid(promoter_all_ants ~ Secretoryarea_perplantnode_mm2, mod1_gls1_promoter_final, data)
plot(YX2$yresid ~ YX2$xresid, las = 1, pch = 19, cex = 0.7, xlim = c(-10, 13), ylim = c(-0.6, 0.6), cex.lab = 0.9, ylab = "", xlab = "Partial [Secretory area of EFNs (mm²)]")
#abline(lm(YX2$yresid ~ YX2$xresid))

------
## ALL ANTS ##

# IMPACT OF EFNs ATTRIBUTES ON RECEPTOR DEGREE


# Model 1 - lm with Gaussian distribution
mod1_receptor_normal <- glm(receptor_all_ants + 0.01 ~ scale(average_activenodes_perplot) + scale(Secretoryarea_perplantnode_mm2), data = data)
summary(mod1_receptor_normal)
plot(mod1_receptor_normal)

# Model 2 - Gamma distribution
mod2_receptor_gamma <- glm(receptor_all_ants + 0.01 ~ scale(average_activenodes_perplot) + scale(Secretoryarea_perplantnode_mm2), data = data, family = Gamma)
summary(mod2_receptor_gamma)

# Model comparison
ICtab(mod1_receptor_normal, mod2_receptor_gamma, type = c("AICc"), weights = TRUE, delta = TRUE, sort = TRUE, nobs = 21)

# Visualizing data distribution
hist(data$receptor_all_ants)
residualPlot(mod1_receptor_normal)
shapiro.test(residuals(mod1_receptor_normal))
hist(residuals(mod1_receptor_normal))
vif(mod1_receptor_normal)

# Partial graphs (receptor)
par(mfrow = c(1, 1), pty = "s")
avPlots(mod1_receptor_normal, ylab = "Receptor degree (all ants)", pch = 16, col.lines = c("gray", "gray"), id = FALSE, main = "", grid = FALSE, layout = c(1, 1), las = 2, lwd = 0)   # Without lines

# Trying to model heterogeneity using gls() function as implemented in Zuur et al 2009

# Traditional model
mod1_receptor <- gls(receptor_all_ants ~ scale(average_activenodes_perplot) + scale(Secretoryarea_perplantnode_mm2), data = data)
summary(mod1_receptor)

# Model 'modeling' increased variance (alternative 1)
vf1Fixed <- varFixed(~average_activenodes_perplot)
mod2_receptor <- gls(receptor_all_ants ~ scale(average_activenodes_perplot) + scale(Secretoryarea_perplantnode_mm2), weights = vf1Fixed, data = data)
summary(mod2_receptor)

# Plotting residuals
plot(mod2_receptor$residuals ~ mod2_receptor$fitted)
plot(mod2_receptor$residuals ~ data$average_activenodes_perplot)
plot(mod2_receptor$residuals ~ data$Secretoryarea_perplantnode_mm2)

# Model 'modeling' increased variance (alternative 2)
vf1power <- varPower(form = ~average_activenodes_perplot)
mod3_receptor <- gls(receptor_all_ants ~ scale(average_activenodes_perplot) + scale(Secretoryarea_perplantnode_mm2), weights = vf1power, data = data)
summary(mod3_receptor)

# Model 'modeling' increased variance (alternative 3)
vf1Exp <- varExp(form = ~average_activenodes_perplot)
mod4_receptor <- gls(receptor_all_ants ~ scale(average_activenodes_perplot) + scale(Secretoryarea_perplantnode_mm2), weights = vf1Exp, data = data, na.action = na.omit)
summary(mod4_receptor)

# Model 'modeling' increased variance (alternative 4)
vf1ConstPower <- varConstPower(form = ~average_activenodes_perplot)
mod5_receptor <- gls(receptor_all_ants ~ scale(average_activenodes_perplot) + scale(Secretoryarea_perplantnode_mm2), weights = vf1ConstPower, data = data, na.action = na.omit)
summary(mod5_receptor)

# Checking which model fits better
anova(mod2_receptor_gamma, mod1_receptor, mod2_receptor, mod3_receptor, mod4_receptor, mod5_receptor)
ICtab(mod2_receptor_gamma, mod1_receptor, mod2_receptor, mod3_receptor, mod4_receptor, mod5_receptor, type = c("AICc"), weights = TRUE, delta = TRUE, sort = TRUE, nobs = 21)

# Constructing "partial residual graphs" by factor
# https://www.rdocumentation.org/packages/piecewiseSEM/versions/1.2.1/topics/partial.resid

# Best model: mod2_receptor_gamma
modified_data <- data
modified_data$Secretoryarea_perplantnode_mm2[13] <- 0.01

mod2_receptor_gamma_final <- glm(receptor_all_ants + 0.01 ~ average_activenodes_perplot + Secretoryarea_perplantnode_mm2, data = modified_data, family = Gamma)
summary(mod2_receptor_gamma_final)
effectsize(mod2_receptor_gamma_final)
tab_model(mod2_receptor_gamma_final)

data$Secretoryarea_perplantnode_mm2
YX1_r <- partialResid(receptor_all_ants ~ average_activenodes_perplot, mod2_receptor_gamma_final, data)
plot(YX1_r$yresid ~ YX1_r$xresid, las = 1, pch = 19, cex = 0.7, xlim = c(-2, 2), ylim = c(-2.5, 1.7), cex.lab = 0.9, ylab = "Partial [Receptor degree (all ants)]", xlab = "Partial [Number of active nodes per species]")
text(YX1_r$yresid ~ YX1_r$xresid, labels = rownames(data), data = data, cex = 0.8, font = 
1, pos = 1)
#abline(lm(YX1_r$yresid ~ YX1_r$xresid))

YX2_r <- partialResid(receptor_all_ants ~ Secretoryarea_perplantnode_mm2, mod2_receptor_gamma_final, modified_data)
plot(YX2_r$yresid ~ YX2_r$xresid, las = 1, pch = 19, cex = 0.7, xlim = c(-3.5, 2), ylim = c(-2.5, 2), cex.lab = 0.9, ylab = "", xlab = "Partial [Secretory area of EFNs (mm²)]")
#abline(lm(YX2_r$yresid ~ YX2_r$xresid))

------
### DOMINANTS ANTS ####

# IMPACT OF EFNs ATTRIBUTES ON PROMOTER DEGREE


### Model 1 - lm with Gaussian distribution ###
mod1_promoter_dominants <- glm(promoter_dominant ~ scale(average_activenodes_perplot) + scale(Secretoryarea_perplantnode_mm2), data = data)
summary(mod1_promoter_dominants)
par(mfrow = c(2, 2))
plot(mod1_promoter_dominants)  ## Observing the pattern of residuals and normality (see page 401 of Crawley's book!)

hist(data$promoter_dominant)
residualPlot(mod1_promoter_dominants)
shapiro.test(residuals(mod1_promoter_dominants))
hist(residuals(mod1_promoter_dominants))
vif(mod1_promoter_dominants)  ## Testing the assumption of no multicollinearity (see page 478 of Zuur's book!)

# Partial graphs:
par(mfrow = c(2, 1))
avPlots(mod1_promoter_dominants, ylab = "Promoter degree (for dominant ants)", pch = 16, col.lines = c("gray", "gray"), id = FALSE, main = "", grid = FALSE, layout = c(1, 1), las = 2, lwd = 2)  # With lines
avPlots(mod1_promoter_dominants, ylab = "Promoter degree (for dominant ants)", pch = 16, col.lines = c("gray", "gray"), id = FALSE, main = "", grid = FALSE, layout = c(1, 1), las = 2, lwd = 0)  # Without lines

## Traditional model
mod1_lm_promoter_dominants <- gls(promoter_dominant ~ scale(average_activenodes_perplot) + scale(Secretoryarea_perplantnode_mm2), data = data)
summary(mod1_lm_promoter_dominants)

## Model 'modeling' increased variance (alternative 1)
vf1Fixed <- varFixed(~average_activenodes_perplot)  # Variance fixed, the variance of residuals is linearly related to a variance covariate
mod1_gls1_promoter_dominants <- gls(promoter_dominant ~ scale(average_activenodes_perplot) + scale(Secretoryarea_perplantnode_mm2), weights = vf1Fixed, data = data)
summary(mod1_gls1_promoter_dominants)

# Intervals (mod1_gls1_promoter_dominants)
plot(mod1_gls1_promoter_dominants$residuals ~ mod1_gls1_promoter_dominants$fitted)
plot(mod1_gls1_promoter_dominants$residuals ~ data$average_activenodes_perplot)
plot(mod1_gls1_promoter_dominants$residuals ~ data$Secretoryarea_perplantnode_mm2)

## Model 'modeling' increased variance (alternative 2)
vf1power <- varPower(form = ~average_activenodes_perplot)
mod1_gls2_promoter_dominants <- gls(promoter_dominant ~ scale(average_activenodes_perplot) + scale(Secretoryarea_perplantnode_mm2), weights = vf1power, data = data)
summary(mod1_gls2_promoter_dominants)

## Model 'modeling' increased variance (alternative 3)
vf1Exp <- varExp(form = ~average_activenodes_perplot)
mod1_gls3_promoter_dominants <- gls(promoter_dominant ~ scale(average_activenodes_perplot) + scale(Secretoryarea_perplantnode_mm2), weights = vf1Exp, data = data, na.action = na.omit)
summary(mod1_gls3_promoter_dominants)

## Model 'modeling' increased variance (alternative 4)
vf1ConstPower <- varConstPower(form = ~average_activenodes_perplot)
mod1_gls4_promoter_dominants <- gls(promoter_dominant ~ scale(average_activenodes_perplot) + scale(Secretoryarea_perplantnode_mm2), weights = vf1ConstPower, data = data, na.action = na.omit)
summary(mod1_gls4_promoter_dominants)

## Checking which model fits better
anova(mod1_lm_promoter_dominants, mod1_gls1_promoter_dominants, mod1_gls2_promoter_dominants, mod1_gls3_promoter_dominants, mod1_gls4_promoter_dominants)
ICtab(mod1_lm_promoter_dominants, mod1_gls1_promoter_dominants, mod1_gls2_promoter_dominants, mod1_gls3_promoter_dominants, mod1_gls4_promoter_dominants, type = c("AICc"), weights = TRUE, delta = TRUE, sort = TRUE, nobs = 21)

# Constructing "partial residual graphs" by factor
# https://www.rdocumentation.org/packages/piecewiseSEM/versions/1.2.1/topics/partial.resid

# Best model above: mod1_gls1_promoter_dominants
mod1_gls1_promoter_final_dom <- gls(promoter_dominant ~ average_activenodes_perplot + Secretoryarea_perplantnode_mm2, weights = vf1Fixed, data = data)
summary(mod1_gls1_promoter_final_dom)
effectsize(mod1_gls1_promoter_final_dom)



YX1_dom <- partialResid(promoter_dominant ~ average_activenodes_perplot, mod1_gls1_promoter_final_dom, data)
plot(jitter(YX1_dom$yresid + 0.05) ~ YX1_dom$xresid, las = 1, pch = 19, cex = 1, xlim = c(-0.3, 1.7), ylim = c(-0.6, 1.15), cex.lab = 0.9, ylab = "Partial [Promoter degree (dominant ants)]", xlab = "Partial [Number of active nodes per species]")
text(YX1_dom$yresid ~ YX1_dom$xresid, labels = rownames(data), data = data, cex = 0.8, font = 2, pos = 2)
abline(lm(YX1_dom$yresid ~ YX1_dom$xresid), lwd = 2)

YX2_dom <- partialResid(promoter_dominant ~ Secretoryarea_perplantnode_mm2, mod1_gls1_promoter_final_dom, data)
plot(YX2_dom$yresid ~ YX2_dom$xresid, las = 1, pch = 19, cex = 0.7, xlim = c(-10, 13), ylim = c(-0.6, 0.6), cex.lab = 0.9, ylab = "", xlab = "Partial [Secretory area of EFNs (mm²)]")
#abline(lm(YX2_dom$yresid ~ YX2_dom$xresid))

------
## DOMINANT ANTS ##

# IMPACT OF EFNs ATTRIBUTES ON RECEPTOR DEGREE


### Model 1 - lm with Gaussian distribution ###
mod1_receptor_normal_dom = glm(receptor_dominant + 0.01 ~ scale(average_activenodes_perplot) + scale(Secretoryarea_perplantnode_mm2), data = data)  ## Corresponds to model_gls1 (normal, does not model variance)
summary(mod1_receptor_normal_dom)
plot(mod1_receptor_normal_dom)

mod2_receptor_gamma_dom = glm(receptor_dominant + 0.01 ~ scale(average_activenodes_perplot) + scale(Secretoryarea_perplantnode_mm2), data = data, family = Gamma)
summary(mod2_receptor_gamma_dom)
m1 = anova(mod2_receptor_gamma_dom)
m1
ICtab(mod1_receptor_normal_dom, mod2_receptor_gamma_dom, type = c("AICc"), weights = TRUE, delta = TRUE, sort = TRUE, nobs = 21)

hist(data$receptor_dominant)
residualPlot(mod1_receptor_normal_dom)
shapiro.test(residuals(mod1_receptor_normal_dom))
hist(residuals(mod1_receptor_normal_dom))
vif(mod1_receptor_normal_dom)

# Partial graphs (receptor):
par(mfrow = c(1, 1), pty = "s")
avPlots(mod1_receptor_normal_dom, ylab = "Receptor degree (for dominant ants)", pch = 16, col.lines = c("gray", "gray"), id = FALSE, main = "", grid = FALSE, layout = c(1, 1), las = 2, lwd = 0)   # Without lines

## Trying a strategy to model heterogeneity using the gls() function as implemented in ch. 4 "Dealing with heterogeneity" - from Zuur et al 2009 (read it!!!)

## Traditional model
mod1_receptor_dom <- gls(receptor_dominant + 0.01 ~ scale(average_activenodes_perplot) + scale(Secretoryarea_perplantnode_mm2), data = data)
summary(mod1_receptor_dom)

## Model 'modeling' increased variance (alternative 1)
vf1Fixed <- varFixed(~average_activenodes_perplot)
mod2_receptor_dom <- gls(receptor_dominant + 0.01 ~ scale(average_activenodes_perplot) + scale(Secretoryarea_perplantnode_mm2), weights = vf1Fixed, data = data)
summary(mod2_receptor_dom)

plot(mod2_receptor_dom$residuals ~ mod2_receptor_dom$fitted)
plot(mod2_receptor_dom$residuals ~ data$average_activenodes_perplot)
plot(mod2_receptor_dom$residuals ~ data$Secretoryarea_perplantnode_mm2)

## Model 'modeling' increased variance (alternative 2)
vf1power <- varPower(form = ~average_activenodes_perplot)
mod3_receptor_dom <- gls(receptor_dominant + 0.01 ~ scale(average_activenodes_perplot) + scale(Secretoryarea_perplantnode_mm2), weights = vf1power, data = data)
summary(mod3_receptor_dom)

## Model 'modeling' increased variance (alternative 3)
vf1Exp <- varExp(form = ~average_activenodes_perplot)
mod4_receptor_dom <- gls(receptor_dominant + 0.01 ~ scale(average_activenodes_perplot) + scale(Secretoryarea_perplantnode_mm2), weights = vf1Exp, data = data, na.action = na.omit)
summary(mod4_receptor_dom)

## Model 'modeling' increased variance (alternative 4)
vf1ConstPower <- varConstPower(form = ~average_activenodes_perplot)
mod5_receptor_dom <- gls(receptor_dominant + 0.01 ~ scale(average_activenodes_perplot) + scale(Secretoryarea_perplantnode_mm2), weights = vf1ConstPower, data = data, na.action = na.omit)
summary(mod5_receptor_dom)

## Checking which model fits better (which models variance differently)
anova(mod2_receptor_gamma_dom, mod1_receptor_dom, mod2_receptor_dom, mod3_receptor_dom, mod4_receptor_dom, mod5_receptor_dom)
ICtab(mod2_receptor_gamma_dom, mod1_receptor_dom, mod2_receptor_dom, mod3_receptor_dom, mod4_receptor_dom, mod5_receptor_dom, type = c("AICc"), weights = TRUE, delta = TRUE, sort = TRUE, nobs = 21)

# Constructing "partial residual graphs" by factor
# https://www.rdocumentation.org/packages/piecewiseSEM/versions/1.2.1/topics/partial.resid

# Best model above (considering the first part of the modeling): mod2_receptor_gamma_dom

data.modified_dom <- data 
data.modified_dom$receptor_dominant <- data$receptor_dominant + 0.01
data.modified_dom$Secretoryarea_perplantnode_mm2[13] <- 0.01

mod2_receptor_gamma_dom <- glm(receptor_dominant ~ average_activenodes_perplot + Secretoryarea_perplantnode_mm2, data = data.modified_dom, family = Gamma)  # Remove scale() before making the plots!
summary(mod2_receptor_gamma_dom)
effectsize(mod2_receptor_gamma_dom)


data.modified_dom$receptor_dominant
data.modified_dom$average_activenodes_perplot
data.modified_dom$Secretoryarea_perplantnode_mm2


YX1_dom <- partialResid(receptor_dominant ~ average_activenodes_perplot, mod2_receptor_gamma_dom, data.modified_dom)
plot(YX1_dom$yresid ~ YX1_dom$xresid, las = 1, pch = 19, cex = 
0.7, xlim = c(-2, 2), ylim = c(-2.5, 1.7), cex.lab = 0.9, ylab = "Partial [Receptor degree (dominant ants)]", xlab = "Partial [Number of active nodes per species]")
#abline(lm(YX1_dom$yresid ~ YX1_dom$xresid))

YX2_dom <- partialResid(receptor_dominant ~ Secretoryarea_perplantnode_mm2, mod2_receptor_gamma_dom, data.modified_dom)
plot(YX2_dom$yresid ~ YX2_dom$xresid, las = 1, pch = 19, cex = 0.7, xlim = c(-3.5, 2), ylim = c(-2.5, 2), cex.lab = 0.9, ylab = "", xlab = "Partial [Secretory area of NEFs (mm²)]")
#abline(lm(YX2_dom$yresid ~ YX2_dom$xresid))

------
### SUBORDINATE ANTS ###

# IMPACT OF EFNs ATTRIBUTES ON PROMOTER DEGREE


### Model 1 - lm with Gaussian distribution ###
mod1_promoter_sub = glm(promoter_subordinate ~ scale(average_activenodes_perplot) + scale(Secretoryarea_perplantnode_mm2), data = data)
summary(mod1_promoter_sub)
par(mfrow = c(2, 2))
plot(mod1_promoter_sub) ## observing the pattern of residuals and normality (see page 401 of Crawley's book!)

hist(data$promoter_subordinate)
residualPlot(mod1_promoter_sub)
shapiro.test(residuals(mod1_promoter_sub))
hist(residuals(mod1_promoter_sub))
vif(mod1_promoter_sub)  ## testing the assumption of absence of multicollinearity (see page 478 of Zuur's book!)

# Partial graphs:
par(mfrow = c(2, 1))
avPlots(mod1_promoter_sub, ylab = "Promoter degree (for subordinate ants)", pch = 16, col.lines = c("gray", "gray"), id = FALSE, main = "", grid = FALSE, layout = c(1, 1), las = 2, lwd = 2)  # With lines
avPlots(mod1_promoter_sub, ylab = "Promoter degree (for subordinate ants)", pch = 16, col.lines = c("gray", "gray"), id = FALSE, main = "", grid = FALSE, layout = c(1, 1), las = 2, lwd = 0)   # Without lines

## Traditional model
mod1_promoter_sub <- gls(promoter_subordinate ~ scale(average_activenodes_perplot) + scale(Secretoryarea_perplantnode_mm2), data = data)
summary(mod1_promoter_sub)  ## without modeling variance in relation to predictors

## Model 'modeling' increased variance (alternative 1)
vf1Fixed <- varFixed(~average_activenodes_perplot) # fixed variance, residuals variance is linearly related to a variance covariate #
mod2_promoter_sub <- gls(promoter_subordinate ~ scale(average_activenodes_perplot) + scale(Secretoryarea_perplantnode_mm2), weights = vf1Fixed, data = data)
summary(mod2_promoter_sub)

## Model 'modeling' increased variance (alternative 2)
vf1power <- varPower(form = ~average_activenodes_perplot)
mod3_promoter_sub <- gls(promoter_subordinate ~ scale(average_activenodes_perplot) + scale(Secretoryarea_perplantnode_mm2), weights = vf1power, data = data)
summary(mod3_promoter_sub)

## Model 'modeling' increased variance (alternative 3)
vf1Exp <- varExp(form = ~average_activenodes_perplot)
mod4_promoter_sub <- gls(promoter_subordinate ~ scale(average_activenodes_perplot) + scale(Secretoryarea_perplantnode_mm2), weights = vf1Exp, data = data, na.action = na.omit)
summary(mod4_promoter_sub)

## Model 'modeling' increased variance (alternative 4)
vf1ConstPower <- varConstPower(form = ~average_activenodes_perplot)
mod5_promoter_sub <- gls(promoter_subordinate ~ scale(average_activenodes_perplot) + scale(Secretoryarea_perplantnode_mm2), weights = vf1ConstPower, data = data, na.action = na.omit)
summary(mod5_promoter_sub)

## Checking which model fits better (which models variance differently)
anova(mod1_promoter_sub, mod2_promoter_sub, mod3_promoter_sub, mod4_promoter_sub, mod5_promoter_sub)
ICtab(mod1_promoter_sub, mod2_promoter_sub, mod3_promoter_sub, mod4_promoter_sub, mod5_promoter_sub, type = c("AICc"), weights = TRUE, delta = TRUE, sort = TRUE, nobs = 21)

# Constructing "partial residual graphs" by factor
# https://www.rdocumentation.org/packages/piecewiseSEM/versions/1.2.1/topics/partial.resid

# Best model above: mod3_promoter_sub
mod3_promoter_sub <- gls(promoter_subordinate ~ scale(average_activenodes_perplot) + scale(Secretoryarea_perplantnode_mm2), weights = vf1power, data = data)
summary(mod3_promoter_sub)
effectsize(mod3_promoter_sub)


par(mfrow = c(1, 1))


YX1 <- partialResid(promoter_subordinate ~ average_activenodes_perplot, mod3_promoter_sub, data)
plot(YX1$yresid ~ YX1$xresid, las = 1, pch = 19, cex.lab = 2, cex = 1, xlim = c(-0.5, 3.5), ylim = c(-0.1, 0.8), cex.lab = 0.9, ylab = "Partial [Promoter degree (subordinate ants)]", xlab = "Partial [Number of active nodes per species]")
text(YX1$yresid ~ jitter(YX1$xresid - 0.2), labels = row.names(data), data = data, cex = 0.8, font = 2, pos = 1)
abline(lm(YX1$yresid ~ YX1$xresid), lwd = 2)

YX2 <- partialResid(promoter_subordinate ~ Secretoryarea_perplantnode_mm2, mod2_promoter_sub, data)
plot(YX2$yresid ~ YX2$xresid, las = 1, pch = 19, cex = 1, xlim = c(-2, 2.6), ylim = c(-0.6, 0.75), cex.lab = 0.9, ylab = "", xlab = "Partial [Secretory area of NEFs (mm²)]")
#abline(lm(YX2$yresid~YX2$xresid))

------
## SUBORDINATE ANTS ##

# IMPACT OF EFNs ATTRIBUTES ON RECEPTOR DEGREE


### Model 1 - lm with Gaussian distribution ###
mod1_receptor_normal_sub = glm(receptor_subordinate + 0.01 ~ scale(average_activenodes_perplot) + scale(Secretoryarea_perplantnode_mm2), data = data)
summary(mod1_receptor_normal_sub)
plot(mod1_receptor_normal_sub)

mod2_receptor_gamma_sub = glm(receptor_subordinate + 0.01 ~ scale(average_activenodes_perplot) + scale(Secretoryarea_perplantnode_mm2), data = data, family = Gamma)
summary(mod2_receptor_gamma_sub)

ICtab(mod1_receptor_normal_sub, mod2_receptor_gamma_sub, type = c("AICc"), weights = TRUE, delta = TRUE, sort = TRUE, nobs = 21)

hist(data$receptor_subordinate)
residualPlot(mod1_receptor_normal_sub)
shapiro.test(residuals(mod1_receptor_normal_sub))
hist(residuals(mod1_receptor_normal_sub))
vif(mod1_receptor_normal_sub)

# partial graphs (receptor):
par(mfrow = c(1, 1), pty = "s")
avPlots(mod1_receptor_normal_sub, ylab = "Receptor degree (for subordinate ants)", pch = 16, col.lines = c("gray", "gray"), id = FALSE, main = "", grid = FALSE, layout = c(1, 1), las = 2, lwd = 0)   # Without lines

## Trying the strategy of modeling heterogeneity using the gls() function as 
## implemented in Chapter 4 "Dealing with heterogeneity" - Zuur et al 2009 (read it there!)

## Traditional model
mod1.receptor_sub <- gls(receptor_subordinate + 0.01 ~ scale(average_activenodes_perplot) + scale(Secretoryarea_perplantnode_mm2), data = data)
summary(mod1.receptor_sub)

## Model 'modeling' increased variance (alternative 1)
vf1Fixed <- varFixed(~average_activenodes_perplot)
mod2.receptor_sub <- gls(receptor_subordinate + 0.01 ~ scale(average_activenodes_perplot) + scale(Secretoryarea_perplantnode_mm2), weights = vf1Fixed, data = data)
summary(mod2.receptor_sub)

## Model 'modeling' increased variance (alternative 2)
vf1power <- varPower(form = ~average_activenodes_perplot)
mod3.receptor_sub <- gls(receptor_subordinate + 0.01 ~ scale(average_activenodes_perplot) + scale(Secretoryarea_perplantnode_mm2), weights = vf1power, data = data)
summary(mod3.receptor_sub)

## Model 'modeling' increased variance (alternative 3)
vf1Exp <- varExp(form = ~average_activenodes_perplot)
mod4.receptor_sub <- gls(receptor_subordinate + 0.01 ~ scale(average_activenodes_perplot) + scale(Secretoryarea_perplantnode_mm2), weights = vf1Exp, data = data, na.action = na.omit)
summary(mod4.receptor_sub)

## Model 'modeling' increased variance (alternative 4)
vf1ConstPower <- varConstPower(form = ~average_activenodes_perplot)
mod5.receptor_sub <- gls(receptor_subordinate + 0.01 ~ scale(average_activenodes_perplot) + scale(Secretoryarea_perplantnode_mm2), weights = vf1ConstPower, data = data, na.action = na.omit)
summary(mod5.receptor_sub)

## Checking which model fits better (which models variance differently)
anova(mod1.receptor_sub, mod2.receptor_sub, mod3.receptor_sub, mod4.receptor_sub, mod5.receptor_sub)
ICtab(mod1.receptor_sub, mod2.receptor_sub, mod3.receptor_sub, mod4.receptor_sub, mod5.receptor_sub, type = c("AICc"), weights = TRUE, delta = TRUE, sort = TRUE, nobs = 21)

# Constructing "partial residual graphs" by factor
# https://www.rdocumentation.org/packages/piecewiseSEM/versions/1.2.1/topics/partial.resid

# Best model above: mod2_receptor_sub

data.modified_sub <- data 
data.modified_sub$Secretoryarea_perplantnode_mm2[13] <- 0.01

mod2.receptor_sub <- gls(receptor_subordinate + 0.01 ~ scale(average_activenodes_perplot) + scale(Secretoryarea_perplantnode_mm2), weights = vf1Fixed, data = data)
summary(mod2.receptor_sub)
tab_model(mod2.receptor_sub)

ci_mod2.receptor_sub <- confint(mod2.receptor_sub)
print(ci_mod2.receptor_sub)

YX1.t <- partialResid(receptor_subordinate + 0.01 ~ average_activenodes_perplot, mod3.receptor_final_sub, data)
plot(YX1.t$yresid ~ YX1.t$xresid, las = 1, pch = 19, cex = 0.7, xlim = c(-0.2, 2.5), ylim = c(-0.2, 1), cex.lab = 0.9, ylab = "Partial [Receptor degree (subordinate ants)]", xlab = "Partial [Number of active nodes per species]")
text(YX1.t$yresid ~ YX1.t$xresid, labels = rownames(data), data = data, cex = 0.8, font = 2, pos = 2)

abline(lm(YX1.t$yresid ~ YX1.t$xresid))

YX2.t <- partialResid(receptor_subordinate + 0.01 ~ Secretoryarea_perplantnode_mm2, mod3.receptor_final_sub, data.modified_sub)
plot(YX2.t$yresid ~ YX2.t$xresid, las = 1, pch = 19, cex = 0.7, xlim = c(-10, 13), ylim = c(-0.4, 0.6), cex.lab = 0.9, ylab = "", xlab = "Partial [Secretory area of NEFs (mm²)]")

#abline(lm(YX2.t$yresid ~ YX2.t$xresid))
