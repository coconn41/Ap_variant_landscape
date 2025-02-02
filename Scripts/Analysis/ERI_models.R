#####
# Load libraries
#####
library(readr)
library(tidyverse)
library(mgcv) 
library(spdep)
library(glmmTMB)
library(performance)
library(spdep)
library(lme4)

#####
# Make ERI
#####

Regression_df_ERI = Regression_df_w_density %>%
  pivot_longer(cols = c(v1,ha)) %>%
  mutate(ERI = (value/tot_tested)*`Target density`)


#####
# Split between adult/nymph
#####
Adult_ERI_df = Regression_df_ERI %>%
  filter(Lifestage == "Adult",
         is.na(`Target density`) == F,
         is.na(UNIT) == F) %>%
  mutate(Year = as.numeric(substring(Date,1,4)))

Nymph_ERI_df = Regression_df_ERI %>%
  filter(Lifestage == "Nymph",
         is.na(`Target density`) == F) %>%
  mutate(Year = as.numeric(substring(Date,1,4)))  

#####
# Adult ha ERI model
#####

adult_ha_ERI_mod = glmmTMB(formula = ERI ~ metric + 
                                 Landscape_metric +
                                 (1|UNIT) + (1|UNIT:Site),
                              family = tweedie(link = 'log'),
                               data = Adult_ERI_df %>%
                             filter(name == "ha"))
adult_ha_ERI_mod2 = glmmTMB(formula = ERI ~ Landscape_metric + 
                             (1|UNIT) + (1|UNIT:Site),
                           family = tweedie(link = 'log'),
                           data = Adult_ERI_df %>%
                             filter(name == "ha"))
adult_ha_ERI_mod3 = glmmTMB(formula = ERI ~ metric + 
                              (1|UNIT) + (1|UNIT:Site),
                            family = tweedie(link = 'log'),
                            data = Adult_ERI_df %>%
                              filter(name == "ha"))
AIC1 = AIC(adult_ha_ERI_mod)
AIC2 = AIC(adult_ha_ERI_mod2)
AIC3 = AIC(adult_ha_ERI_mod3)
diff1 = AIC1-AIC2
diff2 = AIC1-AIC3
diff3 = AIC2-AIC3
adult_ha_ERI_best = adult_ha_ERI_mod2
summary(adult_ha_ERI_best)
#####
# Adult v1 ERI models
#####

adult_v1_ERI_mod = glmmTMB(formula = ERI ~ metric + 
                             Landscape_metric +
                             (1|UNIT) + (1|UNIT:Site),
                           family = tweedie(link = 'log'),
                           data = Adult_ERI_df %>%
                             filter(name == "v1"))
adult_v1_ERI_mod2 = glmmTMB(formula = ERI ~ Landscape_metric + 
                              (1|UNIT) + (1|UNIT:Site),
                            family = tweedie(link = 'log'),
                            data = Adult_ERI_df %>%
                              filter(name == "v1"))
adult_v1_ERI_mod3 = glmmTMB(formula = ERI ~ metric + 
                              (1|UNIT)+(1|UNIT:Site),
                            family = tweedie(link = 'log'),
                            data = Adult_ERI_df %>%
                              filter(name == "v1"))
AIC1 = AIC(adult_v1_ERI_mod)
AIC2 = AIC(adult_v1_ERI_mod2)
AIC3 = AIC(adult_v1_ERI_mod3)
diff1 = AIC1-AIC2
diff2 = AIC1-AIC3
diff3 = AIC2-AIC3
adult_v1_ERI_best = adult_v1_ERI_mod3 # Best picked here due to p not AIC
summary(adult_v1_ERI_best)
#####
# Nymph ha ERI models
#####

nymph_ha_ERI_mod = glmmTMB(formula = ERI ~ metric + 
                             Landscape_metric +
                             (1|UNIT) + (1|UNIT:Site),
                           family = tweedie(link = 'log'),
                           data = Nymph_ERI_df %>%
                             filter(name == "ha"))
nymph_ha_ERI_mod2 = glmmTMB(formula = ERI ~ Landscape_metric + 
                              (1|UNIT) + (1|UNIT:Site),
                            family = tweedie(link = 'log'),
                            data = Nymph_ERI_df %>%
                              filter(name == "ha"))
nymph_ha_ERI_mod3 = glmmTMB(formula = ERI ~ metric + 
                              (1|UNIT)+(1|UNIT:Site),
                            family = tweedie(link = 'log'),
                            data = Nymph_ERI_df %>%
                              filter(name == "ha"))
AIC1 = AIC(nymph_ha_ERI_mod)
AIC2 = AIC(nymph_ha_ERI_mod2)
AIC3 = AIC(nymph_ha_ERI_mod3)
diff1 = AIC1-AIC2
diff2 = AIC1-AIC3
diff3 = AIC2-AIC3
nymph_ha_ERI_best = nymph_ha_ERI_mod2
summary(nymph_ha_ERI_best)
#####
# Nymph v1 ERI models
#####

nymph_v1_ERI_mod = glmmTMB(formula = ERI ~ metric + 
                             Landscape_metric +
                             (1|UNIT) + (1|UNIT:Site),
                           family = tweedie(link = 'log'),
                           data = Nymph_ERI_df %>%
                             filter(name == "v1"))
nymph_v1_ERI_mod2 = glmmTMB(formula = ERI ~ Landscape_metric + 
                              (1|UNIT) + (1|UNIT:Site),
                            family = tweedie(link = 'log'),
                            data = Nymph_ERI_df %>%
                              filter(name == "v1"))
nymph_v1_ERI_mod3 = glmmTMB(formula = ERI ~ metric + 
                              (1|UNIT)+(1|UNIT:Site),
                            family = tweedie(link = 'log'),
                            data = Nymph_ERI_df %>%
                              filter(name == "v1"))
AIC1 = AIC(nymph_v1_ERI_mod)
AIC2 = AIC(nymph_v1_ERI_mod2)
AIC3 = AIC(nymph_v1_ERI_mod3)
diff1 = AIC1-AIC2
diff2 = AIC1-AIC3
diff3 = AIC2-AIC3

nymph_v1_ERI_best = nymph_v1_ERI_mod3
summary(nymph_v1_ERI_best)
final_ERI_models = data.frame(Lifestage = c("Adult","Nymph","Adult","Nymph"),
                              parameter = c("Landscape","Patch","Landscape","Patch"),
                              Genotype = c("ha","ha","v1","v1"),
                              Beta = c(summary(adult_ha_ERI_best)$coefficients$cond[2],
                                       summary(nymph_ha_ERI_best)$coefficients$cond[2],
                                       summary(adult_v1_ERI_best)$coefficients$cond[2],
                                       summary(nymph_v1_ERI_best)$coefficients$cond[2]),
                              # exp_beta = c(exp(summary(adult_ha_ERI_best)$coefficients$cond[2]),
                              #              exp(summary(nymph_ha_ERI_best)$coefficients$cond[2]),
                              #              exp(summary(adult_v1_ERI_best)$coefficients$cond[2]),
                              #              exp(summary(nymph_v1_ERI_best)$coefficients$cond[2])),
                              P_val = c(summary(adult_ha_ERI_best)$coefficients$cond[8],
                                        summary(nymph_ha_ERI_best)$coefficients$cond[8],
                                        summary(adult_v1_ERI_best)$coefficients$cond[8],
                                        summary(nymph_v1_ERI_best)$coefficients$cond[8]))





