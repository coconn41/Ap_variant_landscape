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

adult_ha_ERI_mod = lme4::lmer(formula = log(ERI) ~ metric + 
                                 Landscape_metric +
                                 (1|UNIT) + (1|UNIT:Site),
                               data = Adult_ERI_df %>%
                                        filter(name=="ha",
                                               ERI!=0))
summary(adult_ha_ERI_mod)
AIC(adult_ha_ERI_mod)

adult_ha_ERI_mod = lme4::lmer(formula = log(ERI) ~ metric + 
                                (1|UNIT) + (1|UNIT:Site),
                              data = Adult_ERI_df %>%
                                filter(name=="ha",
                                       ERI!=0))
summary(adult_ha_ERI_mod)
AIC(adult_ha_ERI_mod)
