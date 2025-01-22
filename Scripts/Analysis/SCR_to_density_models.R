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
# Adult density model
#####
adult_density_mod = glmmTMB(formula = density ~ metric + 
                             Landscape_metric +
                             (1|UNIT) + (1|UNIT:Site),
                           family = tweedie(link = 'log'),
                           data = Adult_density_df)
adult_density_mod2 = glmmTMB(formula = density ~ metric + 
                             # Landscape_metric +
                              (1|UNIT) + (1|UNIT:Site),
                            family = tweedie(link = 'log'),
                            data = Adult_density_df)
adult_density_mod3 = glmmTMB(formula = density ~# metric + 
                              Landscape_metric +
                              (1|UNIT) + (1|UNIT:Site),
                            family = tweedie(link = 'log'),
                            data = Adult_density_df)
AIC1=AIC(adult_density_mod)
AIC2=AIC(adult_density_mod2)
AIC3=AIC(adult_density_mod3)
diff1=AIC1-AIC2
diff2=AIC1-AIC3
diff3=AIC2-AIC3
adult_density_mod_best = adult_density_mod3
summary(adult_density_mod_best)
#####
# Nymph density model
#####

nymph_density_mod = glmmTMB(formula = density ~ metric + 
                              Landscape_metric +
                              (1|UNIT) + (1|UNIT:Site),
                            family = tweedie(link = 'log'),
                            data = Nymph_density_df)
nymph_density_mod2 = glmmTMB(formula = density ~ metric + 
                               # Landscape_metric +
                               (1|UNIT) + (1|UNIT:Site),
                             family = tweedie(link = 'log'),
                             data = Nymph_density_df)
nymph_density_mod3 = glmmTMB(formula = density ~# metric + 
                               Landscape_metric +
                               (1|UNIT) + (1|UNIT:Site),
                             family = tweedie(link = 'log'),
                             data = Nymph_density_df)
AIC1=AIC(nymph_density_mod)
AIC2=AIC(nymph_density_mod2)
AIC3=AIC(nymph_density_mod3)
diff1=AIC1-AIC2
diff2=AIC1-AIC3
diff3=AIC2-AIC3
nymph_density_mod_best = nymph_density_mod3
summary(nymph_density_mod_best)

final_scr_to_d_models = data.frame(Lifestage = c("Adult","Nymph"),
                          parameter = c("Landscape","Landscape"),
                          Beta = c(summary(adult_density_mod_best)$coefficients$cond[2],
                                   summary(nymph_density_mod_best)$coefficients$cond[2]),
                          P_val = c(summary(adult_density_mod_best)$coefficients$cond[8],
                                    summary(nymph_density_mod_best)$coefficients$cond[8]))






