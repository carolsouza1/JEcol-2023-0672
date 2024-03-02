#### Local impact of the promoter plant species on the ant 
# visitation pattern and plant herbivory in their vicinity ####


# Packages

library(tidyverse)
library(Rmisc)
library(DHARMa)
library(lme4)
library(lmtest)
library(ggeffects)
library(emmeans)
library(glmmTMB)
library(multcomp)
library(AER)
library(effectsize)
library(sjPlot)
library(report)


dat <- read.table("plots_promoter_ants_herb.txt", header = TRUE, stringsAsFactors = TRUE)
str(dat)

--------
# Tests to explore competition and facilitation hypotheses at a local scale, 
#explicitly considering the neighborhood of the promoting plant species.
--------

### PROBABILIT OF ANT OCCURRENCE ###

# DOMINANT ANTS #
mod2.dom.simple <- glm(occ_dom ~ plant, family = binomial, data = dat)
summary(mod2.dom.simple) # Probabilistic distribution seems ok (no high overdispersion)
anova(mod2.dom.simple, test="Chisq")

mod2.dom <- glmer(occ_dom ~ plant + (1|plot) + (1|specie), family = binomial, data = dat)
summary(mod2.dom) ## no relation found!
drop1(mod2.dom, test = "Chisq")
plot(mod2.dom)
report(mod2.dom)
effectsize(mod2.dom)

## Plot with corrected confidence intervals (unequal intervals):

ggemmeans(mod2.dom, ~ plant, type = "fixed", back.transform = TRUE, interval = "confidence") %>% 
  plot(rawdata = TRUE)

dat3 <- data.frame(ggemmeans(mod2.dom, ~ plant, type = "fixed", back.transform = TRUE, interval = "confidence"))
colnames(dat3)[1] <- "plant"
colnames(dat3)[2] <- "occ_dom"

ggplot(dat, aes(x = plant, y = occ_dom, color = "darkgray")) +
  expand_limits(y = c(0.2, 1)) + 
  geom_jitter(width = 0.08, height = 0.06, color = "grey45", pch = 19, size = 2, alpha = 0.3) + 
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high), size = 0.7, data = dat3) +
  scale_color_manual(values = c("black", "black")) +
  labs(x = "Plants", y = "Probability of occurrence of interaction (dominant ants)") +
  scale_x_discrete(labels = c("Promoter", "Neighbour", "Non-neighbour")) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(color = "black", size = 15),
        axis.text.y = element_text(size = 15, color = "black"),
        axis.title = element_text(size = 15),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(10, 0, 0, 0)))

# SUBORDINATE ANTS #

mod2.sub.simple <- glm(occ_sub ~ plant, family = binomial, data = dat)
summary(mod2.sub.simple) # Probabilistic distribution seems ok (no high overdispersion)

mod2.sub <- glmer(occ_sub ~ plant + (1|plot) + (1|specie), family = binomial, data = dat)
summary(mod2.sub) ## no relation found!
plot(mod2.sub)
drop1(mod2.sub, test = "Chisq")
report(mod2.sub)
effectsize(mod2.dom)


## Plot with corrected confidence intervals (unequal intervals):


ggemmeans(mod2.sub, ~ plant, type = "fixed", back.transform = TRUE, interval = "confidence") %>% 
  plot(rawdata = TRUE)

dat3 <- data.frame(ggemmeans(mod2.sub, ~ plant, type = "fixed", back.transform = TRUE, interval = "confidence"))
colnames(dat3)[1] <- "plant"
colnames(dat3)[2] <- "occ_sub"

ggplot(dat, aes(x = plant, y = occ_sub, color = "darkgray")) +
  expand_limits(y = c(0.2, 1)) + 
  geom_jitter(width = 0.08, height = 0.06, color = "grey45", pch = 19, size = 2, alpha = 0.3) + 
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high), size = 0.7, data = dat3) +
  scale_color_manual(values = c("black", "black")) +
  labs(x = "Plants", y = "Probability of occurrence of interaction (subordinate ants)") +
  scale_x_discrete(labels = c("Promoter", "Neighbour", "Non-Neighbour")) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(color = "black", size = 15),
        axis.text.y = element_text(size = 15, color = "black"),
        axis.title = element_text(size = 15),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(10, 0, 0, 0)))

----------------------------------------------------------------------------
## NUMBER (SUM) OF ANTS VISITING THE EFNs ##

## DOMINANT ANTS ##

mod.dom <- glm(sum_dom ~ plant, family = poisson(link = log), data = dat)
summary(mod.dom)


mod.2a.dom <- glmmTMB(sum_dom ~ plant + (1|plot) + (1|specie), family = nbinom2, data = dat)
summary(mod.2a.dom)  ## no relation found!
res = residuals(mod.2a.dom)
hist(res)
drop1(mod.2a.dom, test = "Chisq")
report(mod.2a.dom)
effectsize(mod.2a.dom)


## Plot with corrected confidence intervals (unequal intervals):

ggemmeans(mod.2a.dom, ~ plant, type = "fixed", back.transform = TRUE, interval = "confidence") %>% 
  plot(rawdata = TRUE)

dat.sub <- data.frame(ggemmeans(mod.2a.dom, ~ plant, type = "fixed", back.transform = TRUE, interval = "confidence"))
colnames(dat.sub)[1] <- "plant"
colnames(dat.sub)[2] <- "sum_dom"

ggplot(dat, aes(x = plant, y = sum_dom, color = "darkgray")) +
  expand_limits(y = c(0.2, 25)) + 
  geom_jitter(width = 0.08, height = 0.06, color = "grey38", pch = 19, size = 2, alpha = 0.3) + 
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high), size = 0.7, data = dat.sub) +
  scale_color_manual(values = c("black", "black")) +
  labs(x = "Plants", y = "Total abundance of dominant ants") +
  scale_x_discrete(labels = c("Promoter", "Neighbour", "Non-Neighbour")) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(color = "black", size = 15),
        axis.text.y = element_text(size = 15, color = "black"),
        axis.title = element_text(size = 15),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(10, 0, 0, 0)))


# SUBORDINATE ANTS #

mod.sub.poisson <- glm(sum_sub ~ plant, data = dat, family = poisson(link = log))
summary(mod.sub.poisson)  # Overdispersion!

mod.sub.bneg <- glm.nb(sum_sub ~ plant, data = dat)
summary(mod.sub.bneg)   # Improved... therefore continue with nb

mod.2.sub <- glmmTMB(sum_sub ~ plant + (1|plot) + (1|specie), family = nbinom1, data = dat)
summary(mod.2.sub) ## no relation found!
drop1(mod.2.sub, test = "Chisq")
effectsize(mod.2.sub)


## Plot with corrected confidence intervals (unequal intervals):

ggemmeans(mod.2.sub, ~ plant, type = "fixed", back.transform = TRUE, interval = "confidence") %>% 
  plot(rawdata = TRUE)

dat.sub <- data.frame(ggemmeans(mod.2.sub, ~ plant, type = "fixed", back.transform = TRUE, interval = "confidence"))
colnames(dat.sub)[1] <- "plant"
colnames(dat.sub)[2] <- "sum_sub"

ggplot(dat, aes(x = plant, y = sum_sub, color = "darkgray")) +
  expand_limits(y = c(0.2, 25)) + 
  geom_jitter(width = 0.08, height = 0.06, color = "grey38", pch = 19, size = 2, alpha = 0.3) + 
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high), size = 0.7, data = dat.sub) +
  scale_color_manual(values = c("black", "black")) +
  labs(x = "Plants", y = "Total abundance of subordinate ants") +
  scale_x_discrete(labels = c("Promoter", "Neighbour", "Non-Neighbour")) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(color = "black", size = 15),
        axis.text.y = element_text(size = 15, color = "black"),
        axis.title = element_text(size = 15),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(10, 0, 0, 0)))

----------------------------------------------------------------------------
########## HERBIVORY ##############

### Percentage of leaf herbivory per plant ###

mod0.herb2 <- glm(mean_herb ~ plant, family = gaussian, data = dat)
summary(mod0.herb2)

mod1.herb2 <- lmer(mean_herb ~ plant + (1|plot), data = dat)
summary(mod1.herb2)
drop1(mod1.herb2, test = "Chisq")

hist(dat$mean_herb)
dat$mean_herb


# Using arcsine transformation for herbivory data
mean_herb_01 <- asin(sqrt(dat$mean_herb/100)) # Transformation added as a new column in the data frame
dat <- cbind(dat, mean_herb_01)

mod2.herb2 <- lmer(mean_herb_01 ~ plant + (1|plot) + (1|specie), data = dat)
summary(mod2.herb2)
drop1(mod2.herb2, test = "Chisq")  ## no relation found 

mod.3.prop.herb2 <- glmmTMB(mean_herb_01 ~ plant + (1|plot) + (1|specie), family = gaussian, data = dat)
summary(mod.3.prop.herb2)
drop1(mod.3.prop.herb2, test = "Chisq")
report(mod.3.prop.herb2)
effectsize(mod.3.prop.herb2)
report_info(mod.3.prop.herb2)



mod.3.null <- glmmTMB(mean_herb_01 ~ 1 + (1|plot) + (1|specie), family = gaussian, data = dat)

AIC(mod.3.null, mod.3.prop.herb2)  # both are similar, so we stick with the simpler one - no effect either! 

bptest(mod.3.prop.herb2)
hist(residuals(mod.3.prop.herb2))
testDispersion(mod.3.prop.herb2)

## For the new plot with corrected confidence intervals (unequal intervals):

ggemmeans(mod.3.prop.herb2, ~ plant, type = "fixed", back.transform = FALSE, interval = "confidence") %>% 
  plot(rawdata = TRUE)

dat.prop.herb2 <- data.frame(ggemmeans(mod.3.prop.herb2, ~ plant, type = "fixed", back.transform = FALSE, interval = "confidence"))
colnames(dat.prop.herb2)[1] <- "plant"
colnames(dat.prop.herb2)[2] <- "mean_herb"


ggplot(dat, aes(x = plant, y = mean_herb)) +
  expand_limits(y = c(0, 1)) + 
  geom_jitter(width = 0.08, height = 0.06, color = "grey30", pch = 19, size = 2, alpha = 0.4) + 
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high), size = 0.7, data = dat.prop.herb2) +
  labs(x = "Plants", y = "Mean leaf area damaged (%)") +
  scale_x_discrete(labels = c("Promoter", "Neighbours", "Non-Neighbours")) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(color = "black", size = 15),
        axis.text.y = element_text(size = 16, color = "black"),
        axis.title = element_text(size = 15),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(10, 0, 0, 0)))
