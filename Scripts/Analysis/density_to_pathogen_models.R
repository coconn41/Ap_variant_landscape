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
# Divide for interpretation?
#####

divider = F

#####
# Rename density
#####
Regression_df_w_density$density = Regression_df_w_density$`Target density`
if(divider == T){
  Regression_df_w_density$Landscape_metric = Regression_df_w_density$Landscape_metric*10
  Regression_df_w_density$metric = Regression_df_w_density$metric*10}


#####
# Split between adult/nymph
#####
Adult_density_df = Regression_df_w_density %>%
  filter(Lifestage == "Adult",
         is.na(`Target density`) == F,
         is.na(UNIT) == F) %>%
  mutate(Year = as.numeric(substring(Date,1,4)))

Nymph_density_df = Regression_df_w_density %>%
  filter(Lifestage == "Nymph",
         is.na(`Target density`) == F) %>%
  mutate(Year = as.numeric(substring(Date,1,4)))

#####
# Build adult models
#####


adult_ha_density_mod = glmmTMB(formula = ha ~ log(density+0.0001) +
                                 (1|UNIT) + (1|UNIT:Site),
                               zi = ~1,
                            family = poisson,
                            offset = log(tot_tested),
                            data = Adult_density_df)
summary(adult_ha_density_mod)

adult_v1_density_mod = glmmTMB(formula = v1 ~ log(density+0.0001) +
                                 (1|UNIT) + (1|UNIT:Site),
                              zi = ~1,
                            family = poisson,
                            offset = log(tot_tested),
                            data = Adult_density_df)
summary(adult_v1_density_mod)

#####
# Build nymph models
#####

nymph_ha_density_mod = glmmTMB(formula = ha ~ log(density+0.0001) +
                                 (1|UNIT) + (1|UNIT:Site),
                              zi = ~1,
                            family = poisson,
                            offset = log(tot_tested),
                            data = Nymph_density_df)
summary(nymph_ha_density_mod)

nymph_v1_density_mod = glmmTMB(formula = v1 ~ log(density+0.0001) +
                                 (1|UNIT) + (1|UNIT:Site),
                               zi = ~1,
                            family = poisson,
                            offset = log(tot_tested),
                            data = Nymph_density_df)
summary(nymph_v1_density_mod)

final_d_to_p_models = data.frame(Lifestage = c("Adult","Nymph","Adult","Nymph"),
                                 Pathogen = c('ha','ha','v1','v1'),
                          Beta = c(summary(adult_ha_density_mod)$coefficients$cond[2],
                                   summary(nymph_ha_density_mod)$coefficients$cond[2],
                                   summary(adult_v1_density_mod)$coefficients$cond[2],
                                   summary(nymph_v1_density_mod)$coefficients$cond[2]),
                          exp_beta = c(exp(summary(adult_ha_density_mod)$coefficients$cond[2]),
                                       exp(summary(nymph_ha_density_mod)$coefficients$cond[2]),
                                       exp(summary(adult_v1_density_mod)$coefficients$cond[2]),
                                       exp(summary(nymph_v1_density_mod)$coefficients$cond[2])),
                          P_val = c(summary(adult_ha_density_mod)$coefficients$cond[8],
                                    summary(nymph_ha_density_mod)$coefficients$cond[8],
                                    summary(adult_v1_density_mod)$coefficients$cond[8],
                                    summary(nymph_v1_density_mod)$coefficients$cond[8]))

