# this is bunk


#####
# Load libraries
#####
library(readr)
library(tidyverse)
library(mgcv) 
library(spdep)
library(glmmTMB)
library(performance)
source(paste0(getwd(),"/Scripts/Analysis/Local_Moran_fn.R"))
######
# Competition models: Adult, a-spatial and test for spatial autocorrelation
#####
adult_competition = glmmTMB(formula = competition ~ metric + 
                                 Closest_UA +
                                 (1|UNIT) + (1|Site:UNIT),
                               offset = log(tot_tested+0.01),
                               data = Adult_df)
adult_competition2 = glmmTMB(formula = competition ~ metric + 
                                  (1|UNIT) + (1|Site:UNIT),
                                offset = log(tot_tested+0.01),
                               data = Adult_df)
if((AIC(adult_competition2) - AIC(adult_competition))<2){
  adult_competition = adult_competition2}
remove(adult_competition2)
summary(adult_competition)
Moran_maps(df = Nymph_df,
           model = adult_ha_competition,
           model_type = "competition",
           seed = T,
           seed_n = 1)


######
# Competition models: Adult, a-spatial v1 and test for spatial autocorrelation
#####

adult_v1_competition = glmmTMB(v1_competition ~ metric + 
                                 Closest_UA +
                                 (1|UNIT) + (1|Site:UNIT),
                               offset = log(tot_tested+0.01),
                               data = Adult_df)
adult_v1_competition2 = glmmTMB(v1_competition ~ metric + 
                                  (1|UNIT) + (1|Site:UNIT),
                                offset = log(tot_tested+0.01),
                               data = Adult_df)
if((AIC(adult_v1_competition2) - AIC(adult_v1_competition))<2){
  adult_v1_competition = adult_v1_competition2}
remove(adult_v1_competition2)
summary(adult_v1_competition)
Moran_maps(df = Adult_df,
           model = adult_v1_competition,
           model_type = "competition",
           seed = T,
           seed_n = 1)

#####
# Competition models: Nymph, a-spatial ha and test for spatial autocorrelation
#####
nymph_ha_competition = glmmTMB(ha_competition ~ metric + 
                                 Closest_UA +
                                 (1|UNIT) + (1|Site:UNIT),
                               offset = log(tot_tested+0.01),
                               data = Nymph_df)
nymph_ha_competition2 = glmmTMB(ha_competition ~ metric + 
                                  (1|UNIT) + (1|Site:UNIT),
                                offset = log(tot_tested+0.01),
                               data = Nymph_df)
if((AIC(nymph_ha_competition2) - AIC(nymph_ha_competition))<2){
  nymph_ha_competition = nymph_ha_competition2}
remove(nymph_ha_competition2)
summary(nymph_ha_competition)
Moran_maps(df = Nymph_df,
           model = nymph_ha_competition,
           model_type = "competition",
           seed = T,
           seed_n = 1)


#####
# Competition models: Nymph, a-spatial v1 and test for spatial autocorrelation
#####

nymph_v1_competition = glmmTMB(v1_competition ~ metric + 
                                 Closest_UA +
                                 (1|UNIT) + (1|Site:UNIT),
                               offset = log(tot_tested+0.01),
                               data = Nymph_df)
nymph_v1_competition2 = glmmTMB(v1_competition ~ metric + 
                                  (1|UNIT) + (1|Site:UNIT),
                                offset = log(tot_tested+0.01),
                               data = Nymph_df)
if((AIC(nymph_v1_competition2) - AIC(nymph_v1_competition))<2){
  nymph_v1_competition = nymph_v1_competition2}
remove(nymph_v1_competition2)
summary(nymph_v1_competition)
Moran_maps(df = Nymph_df,
           model = nymph_v1_competition,
           model_type = "competition",
           seed = T,
           seed_n = 1)

#####
# Combine a-spatial model results:
#####

a_spat_comp_mods = data.frame(Lifestage = c(),
                              Parameter = c(),
                              Estimate = c(),
                              P_val = c())

#####
# Spatial models:
#####
pos = numFactor(Adult_df$longitude, Adult_df$latitude)
group = factor(rep(1,nrow(Adult_df)))
adult_ha_competition_spatial = glmmTMB(ha_competition ~ metric + 
                                         Closest_UA + exp(pos + 0 | group)+
                 (1|Site),
               zi = ~1,
               family = poisson,
               offset = log(tot_tested),
               data = Adult_df)
adult_ha_competition_spatial2 = glmmTMB(ha_competition ~ metric + 
                                          exp(pos + 0 | group)+
                                         (1|Site),
                                       zi = ~1,
                                       family = poisson,
                                       offset = log(tot_tested),
                                       data = Adult_df)
summary(adult_ha_competition_spatial2)
if((AIC(adult_ha_competition_spatial2) - AIC(adult_ha_competition_spatial))<2){
  adult_ha_competition_spatial = adult_ha_competition_spatial2}
remove(adult_ha_competition_spatial2)
summary(adult_ha_competition_spatial)
Moran_maps(df = Adult_df,
           model = adult_ha_competition_spatial,
           model_type = "competition",
           seed = T,
           seed_n = 1)
