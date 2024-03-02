### IMPACT OF INDIRECT EFFECTS ON LEAF HERBIVORY PATTERNS ###

# PACKAGES
library(car)
library(bbmle)   ## To compare models in a table with AICs
library(nlme)     ## To run GLS models
library(piecewiseSEM)   ## To calculate partial residuals in models
library(lmtest)
library(report)
library(effectsize)


dados <- read.table("herb_muller.txt", header = TRUE, dec = ".")
str(dados)

-------
# Test whether plant species that promote indirect effects by attracting more ants negatively (competition hypothesis)
# or positively (facilitation hypothesis) influence leaf herbivory of species receiving indirect effects.
-------
  
# PROMOTER DEGREE and RECEPTOR DEGREE - AVERAGE PROPORTION OF DAMAGED LEAFLETS #

### ALL ANTS ###


# Test correlation between promoter and receptor degrees
cor.test(dados$promoter_all_ants, dados$receptor_all_ants, method = "pearson")

# Fit GLM model for the average proportion of damaged leaflets
mod1_promoter_receptor = glm(mean_prop_leaflets_herb_sp ~ promoter_all_ants + receptor_all_ants, data = dados)
summary(mod1_promoter_receptor)
anova(mod1_promoter_receptor, test="Chisq")
report(mod1_promoter_receptor)
report_effectsize(mod1_promoter_receptor)
report_table(mod1_promoter_receptor)


# Check for multicollinearity using VIF
vif(mod1_promoter_receptor)

# Plot model diagnostics
plot(mod1_promoter_receptor)


# Histogram of the dependent variable
hist(dados$mean_prop_leaflets_herb_sp)

# Residual plot
residualPlot(mod1_promoter_receptor)

# Shapiro-Wilk test for normality of residuals
shapiro.test(residuals(mod1_promoter_receptor))

# Histogram of residuals
hist(residuals(mod1_promoter_receptor))

# Heteroscedasticity test using Breusch-Pagan
bptest(mod1_promoter_receptor)

# Partial residual plots
YX1 <- partialResid(mean_prop_leaflets_herb_sp ~ promoter_all_ants, mod1_promoter_receptor, dados)
plot(YX1$yresid ~ YX1$xresid, las = 1, pch = 19, cex = 1, cex.lab = 0.7, ylab = "Partial [Average proportion of damaged leaflets per species]", xlab = "Partial [Promoter degree (all ants)]")
text(YX1$yresid ~ YX1$xresid, labels = rownames(dados), data = dados, cex = 0.8, font = 2, pos = 2)

YX2 <- partialResid(mean_prop_leaflets_herb_sp ~ receptor_all_ants, mod1_promoter_receptor, dados)
plot(YX2$yresid ~ YX2$xresid, las = 1, pch = 19, cex = 0.7, cex.lab = 0.9, ylab = "Partial [Average proportion of damaged leaflets per species]", xlab = "Partial [Receptor degree (all ants)]")
-------

### DOMINANT ANTS ###

# PROMOTER DEGREE and RECEPTOR DEGREE - AVERAGE PROPORTION OF DAMAGED LEAFLETS #

# Test correlation between promoter and receptor degrees
cor.test(dados$promoter_dominant, dados$receptor_dominant, method = "pearson")

# Scatter plot of promoter and receptor degrees
plot(dados$promoter_dominant ~ dados$receptor_dominant, las = 1)

# Fit GLM model for the average proportion of damaged leaflets
mod1_promoter_receptor_dom = glm(mean_prop_leaflets_herb_sp ~ promoter_dominant + receptor_dominant, data = dados)
summary(mod1_promoter_receptor_dom)
anova(mod1_promoter_receptor_dom, test="Chisq")
report(mod1_promoter_receptor_dom)
report_effectsize(mod1_promoter_receptor_dom)


# Check for multicollinearity using VIF
vif(mod1_promoter_receptor_dom)

# Plot model diagnostics
plot(mod1_promoter_receptor_dom)


# Histogram of the dependent variable
hist(dados$media_prop_foliolos_herb_sp)

# Residual plot
residualPlot(mod1_promoter_receptor_dom)

# Shapiro-Wilk test for normality of residuals
shapiro.test(residuals(mod1_promoter_receptor_dom))
hist(residuals(mod1_promoter_receptor_dom))

# Heteroscedasticity test using Breusch-Pagan
bptest(mod1_promoter_receptor_dom)

# Partial residual plots
YX1 <- partialResid(mean_prop_leaflets_herb_sp ~ promoter_dominant, mod1_promoter_receptor_dom, dados)
plot(YX1$yresid ~ jitter(YX1$xresid + 0.02), las = 1, pch = 19, cex = 1, cex.lab = 0.7, ylab = "Partial [Average proportion of damaged leaflets per species]", xlab = "Partial [Promoter degree (dominant ants)]")
text(YX1$yresid ~ jitter(YX1$xresid), labels = rownames(dados), data = dados, cex = 0.8, font = 2, pos = 2)

YX2 <- partialResid(mean_prop_leaflets_herb_sp ~ receptor_dominant, mod1_promoter_receptor_dom, dados)
plot(YX2$yresid ~ YX2$xresid, las = 1, pch = 19, cex = 0.7, cex.lab = 0.9, ylab = "Partial [Average proportion of damaged leaflets per species]", xlab = "Partial [Receptor degree (dominant ants)]")
text(YX2$yresid ~ jitter(YX2$xresid), labels = rownames(dados), data = dados, cex = 0.8, font = 2, pos = 2)
-------

### SUBORDINATE ANTS ###

# PROMOTER DEGREE and RECEPTOR DEGREE - AVERAGE PROPORTION OF DAMAGED LEAFLETS #

# Test correlation between promoter and receptor degrees
cor.test(dados$promoter_subordinate, dados$receptor_subordinate, method = "pearson")

# Scatter plot of promoter and receptor degrees
plot(jitter(dados$promoter_subordinate, 30) ~ dados$receptor_subordinate, las = 1)

# Fit GLM model for the average proportion of damaged leaflets
mod1_promoter_receptor_sub = glm(mean_prop_leaflets_herb_sp ~ promoter_subordinate + receptor_subordinate, data = dados)
summary(mod1_promoter_receptor_sub)
report(mod1_promoter_receptor_sub)


# VIF for multicollinearity check
vif(mod1_promoter_receptor_sub)

# Plot model diagnostics
plot(mod1_promoter_receptor_sub)

# Confidence intervals for coefficients
confidance_intervals <- confint(mod1_promoter_receptor_sub)
confidance_intervals

# Histogram of the dependent variable
hist(dados$media_prop_foliolos_herb_sp)

# Shapiro-Wilk test for normality of residuals
shapiro.test(residuals(mod1_promoter_receptor_sub))
residualPlot(mod1_promoter_receptor_sub)
hist(residuals(mod1_promoter_receptor_sub))

# Heteroscedasticity test using Breusch-Pagan
bptest(mod1_promoter_receptor_sub)

# Partial residual plots
YX1 <- partialResid(mean_prop_leaflets_herb_sp ~ promoter_subordinate, mod1_promoter_receptor_sub, dados)
plot(YX1$yresid ~ YX1$xresid, las = 1, pch = 19, cex = 1, cex.lab = 0.9, ylim = c(-0.5, 0.5), ylab = "Partial [Average proportion of damaged leaflets per species]", xlab = "Partial [Promoter degree (subordinate ants)]")
text(YX1$yresid ~ YX1$xresid, labels = row.names(dados), data = dados, cex = 0.8, font = 2, pos = 1)

YX2 <- partialResid(mean_prop_leaflets_herb_sp ~ receptor_subordinate, mod1_promoter_receptor_sub, dados)
plot(YX2$yresid ~ YX2$xresid, las = 1, pch = 19, cex = 0.7, cex.lab = 0.9, ylab = "Partial [Average proportion of damaged leaflets per species]", xlab = "Partial [Receptor degree (subordinate ants)]")
text(YX2$yresid ~ YX2$xresid, labels = row.names(dados), data = dados, cex = 0.8, font = 2, pos = 1)
--------

# AVERAGE PERCENTAGE OF LEAF AREA CONSUMED BY SPECIES #

# ALL ANTS #

# MODEL 1 - NORMAL GLM #
  

m1_promoter_receptor = glm(mean_herb_sp ~ promoter_all_ants + receptor_all_ants, data = dados)  # Weak correlation between acting and target, and VIF <3
vif(m1_promoter_receptor)
summary(m1_promoter_receptor)

plot(m1_promoter_receptor)

hist(dados$mean_herb_sp)
residualPlot(m1_promoter_receptor)
shapiro.test(residuals(m1_promoter_receptor))
hist(residuals(m1_promoter_receptor))

bptest(m1_promoter_receptor)  # Accept null hypothesis of homoscedasticity (Breusch & Pagan 1979)

# As there are normality issues, try transforming the response variable #

# MODEL 2 - GLM WITH SQUARE ROOT TRANSFORMATION #

m2_promoter_receptor = glm(sqrt(mean_herb_sp) ~ promoter_all_ants + receptor_all_ants, data = dados)  # Weak correlation between acting and target, and VIF <3
summary(m2_promoter_receptor)
anova(m2_promoter_receptor, test="Chisq")
vif(m2_promoter_receptor)
plot(m2_promoter_receptor)
report(m2_promoter_receptor)


bptest(m2_promoter_receptor)

ICtab(m1_promoter_receptor, m2_promoter_receptor, type = c("AICc"), weights = TRUE, delta = TRUE, sort = TRUE, nobs = 21)
# m2_promoter_receptor is the better model #

# Partial residual plots #

YX1 <- partialResid(mean_herb_sp ~ promoter_all_ants, m2_promoter_receptor, dados)
plot(YX1$yresid ~ YX1$xresid, las = 1, pch = 19, cex = 0.7, cex.lab = 0.9, ylab = "Partial log[Mean damaged leaflets per species]", xlab = "Partial [Promoter degree (all ants)]")
text(YX1$yresid ~ jitter(YX1$xresid), labels = rownames(dados), data = dados, cex = 0.8, font = 2, pos = 2)
#abline(lm(YX1$yresid ~ YX1$xresid))

YX2 <- partialResid(mean_herb_sp ~ receptor_all_ants, m2_promoter_receptor, dados)
plot(jitter(YX2$yresid + 0.05) ~ YX2$xresid, las = 1, pch = 19, cex = 1, cex.lab = 0.9, ylab = "", xlab = "")
text(YX2$yresid ~ YX2$xresid, labels = rownames(dados), data = dados, cex = 0.8, font = 2, pos = 2)
#abline(lm(YX2$yresid ~ YX2$xresid))

-------
# AVERAGE PERCENTAGE OF LEAF AREA CONSUMED BY SPECIES #

### DOMINANT ANTS ###

m1_promoter_receptor_dom = glm(mean_herb_sp ~ promoter_dominant + receptor_dominant, data = dados)  # Weak correlation between acting and target, and VIF <3
summary(m1_promoter_receptor_dom)
plot(m1_promoter_receptor_dom)
vif(m1_promoter_receptor_dom)


hist(dados$mean_herb_sp)
residualPlot(m1_promoter_receptor_dom)
shapiro.test(residuals(m1_promoter_receptor_dom))
hist(residuals(m1_promoter_receptor_dom))

bptest(m1_promoter_receptor_dom)

# MODEL 2 GLM, SQUARE ROOT TRANSFORMATION #

m2_promoter_receptor_dom = glm(sqrt(mean_herb_sp) ~ promoter_dominant + receptor_dominant, data = dados)  # Weak correlation between acting and target, and VIF <3
summary(m2_promoter_receptor_dom)
anova(m2_promoter_receptor_dom, test="Chisq")
plot(m2_promoter_receptor_dom)
vif(m2_promoter_receptor_dom)
report(m2_promoter_receptor_dom)


residualPlot(m2_promoter_receptor_dom)
shapiro.test(residuals(m2_promoter_receptor_dom))
hist(residuals(m2_promoter_receptor_dom))

bptest(m2_promoter_receptor_dom)

bptest(m2_promoter_receptor_dom, mean_herb_sp ~ promoter_dominant * receptor_dominant + I(promoter_dominant^2) + I(receptor_dominant^2), data = dados) 

ICtab(m1_promoter_receptor_dom, m2_promoter_receptor_dom, type = c("AICc"), weights = TRUE, delta = TRUE, sort = TRUE, nobs = 21)

# Best model: m2_promoter_receptor_dom #

# Partial residual plots #

YX1 <- partialResid(mean_herb_sp ~ promoter_dominant, m2_promoter_receptor_dom, dados)
plot(YX1$yresid ~ YX1$xresid, las = 1, pch = 19, cex = 0.7, cex.lab = 0.9, ylab = "Partial log[Mean damaged leaflets per species]", xlab = "Partial [Promoter degree (dominant ants)]")
#abline(lm(YX1$yresid ~ YX1$xresid))

YX2 <- partialResid(mean_herb_sp ~ receptor_dominant, m2_promoter_receptor_dom, dados)
plot(YX2$yresid ~ YX2$xresid, las = 1, pch = 19, cex = 0.7, cex.lab = 0.9, ylab = "Partial log[Mean damaged leaflets per species]", xlab = "Partial [Receptor degree (dominant ants)]")
#abline(lm(YX2$yresid ~ YX2$xresid))

-------
# AVERAGE PERCENTAGE OF LEAF AREA CONSUMED BY SPECIES #

# SUBORDINANT ANTS #

m1_promoter_receptor_sub = glm(mean_herb_sp ~ promoter_subordinate + receptor_subordinate, data = dados)
summary(m1_promoter_receptor_sub)
vif(m1_promoter_receptor_sub)

plot(m1_promoter_receptor_sub)

hist(dados$mean_herb_sp)
residualPlot(m1_promoter_receptor_sub)
shapiro.test(residuals(m1_promoter_receptor_sub))
hist(residuals(m1_promoter_receptor_sub))

bptest(m1_promoter_receptor_sub)

# MODEL 2 GLM, SQUARE ROOT TRANSFORMATION #

m2_promoter_receptor_sub = glm(sqrt(mean_herb_sp) ~ promoter_subordinate + receptor_subordinate, data = dados) # Weak correlation between acting and target, and VIF <3
summary(m2_promoter_receptor_sub)
report(m2_promoter_receptor_sub)
vif(m2_promoter_receptor_sub)


plot(m2_promoter_receptor_sub)

mean_herb_sp_sqrt = sqrt(dados$mean_herb_sp)

hist(med_herb_sp_sqrt)
residualPlot(m2_promoter_receptor_sub)
shapiro.test(residuals(m2_promoter_receptor_sub))
hist(residuals(m2_promoter_receptor_sub))

bptest(m2_promoter_receptor_sub)

ICtab(m1_promoter_receptor_sub, m2_promoter_receptor_sub, type = c("AICc"), weights = TRUE, delta = TRUE, sort = TRUE, nobs = 21)

# Best model: m2_promoter_receptor_sub #

# Partial residual plots #

YX1 <- partialResid(mean_herb_sp ~ promoter_subordinate, m2_promoter_receptor_sub, dados)
plot(YX1$yresid ~ YX1$xresid, las = 1, pch = 19, cex = 0.7, cex.lab = 0.9, ylab = "Partial log[Mean damaged leaflets per species]", xlab = "Partial [Promoter degree (Subordinate ants)]")
#abline(lm(YX1$yresid ~ YX1$xresid))

YX2 <- partialResid(mean_herb_sp ~ receptor_subordinate, m2_promoter_receptor_sub, dados)
plot(YX2$yresid ~ YX2$xresid, las = 1, pch = 19, cex = 0.7, cex.lab = 0.9, ylab = "Partial log[Mean damaged leaflets per species]", xlab = "Partial [Receptor degree (Subordinate ants)]")
#abline(lm(YX2$yresid ~ YX2$xresid))




