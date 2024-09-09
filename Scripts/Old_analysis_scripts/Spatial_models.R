#####
# Load libraries
#####
library(readr)
library(tidyverse)
library(mgcv) 
library(spdep)
library(glmmTMB)
library(performance)

#####
# 1.) Spatial adult ha model:
#####
adult_ha_spatial = mgcv::gam(ha ~ metric + Closest_UA + 
                               s(longitude,latitude) + 
                               s(Site,bs='re'),
                             family = ziP,
                             offset = log(tot_tested),
                             data = Adult_df)
adult_ha_spatial2 = mgcv::gam(ha ~ metric +  
                                s(longitude,latitude) + 
                                s(Site,bs='re'),
                              family = ziP,
                              offset = log(tot_tested),
                              data = Adult_df)
if((AIC(adult_ha_spatial2) - AIC(adult_ha_spatial))<2){
  adult_ha_spatial = adult_ha_spatial2}
remove(adult_ha_spatial2)
summary(adult_ha_spatial)
Moran_maps(df = Adult_df,
           model = adult_ha_spatial,
           seed = T,
           seed_n = 1)
# Comparing to aspatial:
if((AIC(adult_ha_aspatial) - AIC(adult_ha_spatial))<2){
  final_adult_ha_model = adult_ha_aspatial}
if((AIC(adult_ha_aspatial) - AIC(adult_ha_spatial))>=2){
  final_adult_ha_model = adult_ha_spatial}
summary(final_adult_ha_model)

#####
# 2.) Spatial adult v1 model:
#####
adult_v1_spatial = mgcv::gam(v1 ~ metric + Closest_UA + 
                               s(longitude,latitude) + 
                               s(Site,bs='re'),
                             family = ziP,
                             offset = log(tot_tested),
                             data = Adult_df)
adult_v1_spatial2 = mgcv::gam(v1 ~ metric + #Closest_UA + 
                                s(longitude,latitude) + 
                                s(Site,bs='re'),
                              family = ziP,
                              offset = log(tot_tested),
                              data = Adult_df)
if((AIC(adult_v1_spatial2) - AIC(adult_v1_spatial))<2){
  adult_v1_spatial = adult_v1_spatial2}
remove(adult_v1_spatial2)
summary(adult_v1_spatial)
Moran_maps(df = Adult_df,
           model = adult_v1_spatial,
           seed = T,
           seed_n = 1)
# Comparing to aspatial:
if((AIC(adult_v1_aspatial) - AIC(adult_v1_spatial))<2){
  final_adult_v1_model = adult_v1_aspatial}
if((AIC(adult_v1_aspatial) - AIC(adult_v1_spatial))>=2){
  final_adult_v1_model = adult_v1_spatial}
summary(final_adult_v1_model)

#####
# 3.) Spatial nymph ha model:
#####
nymph_ha_spatial = mgcv::gam(ha ~ metric + Closest_UA + 
                               s(longitude,latitude) + 
                               s(Site,bs='re'),
                             family = ziP,
                             offset = log(tot_tested),
                             data = Nymph_df)
nymph_ha_spatial2 = mgcv::gam(ha ~ metric + #Closest_UA + 
                                s(longitude,latitude) + 
                                s(Site,bs='re'),
                              family = ziP,
                              offset = log(tot_tested),
                              data = Nymph_df)
if((AIC(nymph_ha_spatial2) - AIC(nymph_ha_spatial))<2){
  nymph_ha_spatial = nymph_ha_spatial2}
remove(nymph_ha_spatial2)
summary(nymph_ha_spatial)
Moran_maps(df = Nymph_df,
           model = nymph_ha_spatial,
           seed = T,
           seed_n = 1)
# Comparing to aspatial:
if((AIC(nymph_ha_aspatial) - AIC(nymph_ha_spatial))<2){
  final_nymph_ha_model = nymph_ha_aspatial}
if((AIC(nymph_ha_aspatial) - AIC(nymph_ha_spatial))>=2){
  final_nymph_ha_model = nymph_ha_spatial}
summary(final_nymph_ha_model)

#####
# 4.) Spatial nymph v1 model:
#####
nymph_v1_spatial = mgcv::gam(v1 ~ metric + Closest_UA + 
                               s(longitude,latitude) + 
                               s(Site,bs='re'),
                             family = ziP,
                             offset = log(tot_tested),
                             data = Nymph_df)
nymph_v1_spatial2 = mgcv::gam(v1 ~ metric + #Closest_UA + 
                                s(longitude,latitude) + 
                                s(Site,bs='re'),
                              family = ziP,
                              offset = log(tot_tested),
                              data = Nymph_df)
if((AIC(nymph_v1_spatial2) - AIC(nymph_v1_spatial))<2){
  nymph_v1_spatial = nymph_v1_spatial2}
remove(nymph_v1_spatial2)
summary(nymph_v1_spatial)
Moran_maps(df = Nymph_df,
           model = nymph_v1_spatial,
           seed = T,
           seed_n = 1)
# Comparing to aspatial:
if((AIC(nymph_v1_aspatial) - AIC(nymph_v1_spatial))<2){
  final_nymph_v1_model = nymph_v1_aspatial}
if((AIC(nymph_v1_aspatial) - AIC(nymph_v1_spatial))>=2){
  final_nymph_v1_model = nymph_v1_spatial}
summary(final_nymph_v1_model)

#####
# Combine results from final models:
#####

summary(final_adult_ha_model)
summary(final_adult_v1_model)
summary(final_nymph_ha_model)
summary(final_nymph_v1_model)

final_patch_results = data.frame(Lifestage = c(rep("Adult",6),
                                            rep("Nymph",4)),
                              Genotype = c(rep('ha',4),
                                           "v1","v1",
                                           "ha","ha",
                                           "v1","v1"),
                              Variable = c("Metric","Metric","Closest_UA","Closest_UA",
                                           "Metric","Metric",
                                           "Metric","Metric",
                                           "Metric","Metric"),
                              Parameter = rep(c("Estimate","P-value"),5),
                              Value = c(summary(final_adult_ha_model)$p.coeff[2],
                                        summary(final_adult_ha_model)$p.pv[2],
                                        summary(final_adult_ha_model)$p.coeff[3],
                                        summary(final_adult_ha_model)$p.pv[3],
                                        summary(final_adult_v1_model)$p.coeff[2],
                                        summary(final_adult_v1_model)$p.pv[2],
                                        summary(final_nymph_ha_model)$p.coeff[2],
                                        summary(final_nymph_ha_model)$p.pv[2],
                                        summary(final_nymph_v1_model)$p.coeff[2],
                                        summary(final_nymph_v1_model)$p.pv[2])) %>%
  pivot_wider(names_from = Parameter,
              values_from = Value) %>%
  mutate(Significant = ifelse(`P-value`<0.05,1,0))
