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
#####
# Include local moran?
#####
local_moran = F

######
# A-spatial adult ha model:
#####
adult_ha_aspatial = glmmTMB(formula = ha ~ metric + 
                              Landscape_metric +
                              (1|UNIT) + (1|UNIT:Site),
                            zi = ~1,
                            family = poisson,
                            offset = log(tot_tested),
                            data = Adult_df)
adult_ha_aspatial2 = glmmTMB(formula = ha ~ metric + 
                               (1|Site),
                             zi = ~1,
                             offset = log(tot_tested),
                             family = poisson,
                             data = Adult_df)
adult_ha_aspatial3 = glmmTMB(formula = ha ~ Landscape_metric + 
                              (1|UNIT) + (1|UNIT:Site),
                            zi = ~1,
                            offset = log(tot_tested),
                            family = poisson,
                            data = Adult_df)
if(AIC(adult_ha_aspatial)==min(c(AIC(adult_ha_aspatial),
      AIC(adult_ha_aspatial2),
      AIC(adult_ha_aspatial3)))){
  adult_ha_best = adult_ha_aspatial}
if(AIC(adult_ha_aspatial2)==min(c(AIC(adult_ha_aspatial),
                                 AIC(adult_ha_aspatial2),
                                 AIC(adult_ha_aspatial3)))){
  adult_ha_best = adult_ha_aspatial2}
if(AIC(adult_ha_aspatial3)==min(c(AIC(adult_ha_aspatial),
                                 AIC(adult_ha_aspatial2),
                                 AIC(adult_ha_aspatial3)))){
  adult_ha_best = adult_ha_aspatial3}

summary(adult_ha_best)
Adult_df$aspatial_ha_residuals = residuals(adult_ha_best)
if(local_moran==T){Moran_maps(df = Adult_df,
           model = adult_ha_best,
           seed = T,
           seed_n = 1)}
res = data.frame(estimate = NA,
                 p.value = NA)
ind = 0
for(i in min(Adult_df$Year):max(Adult_df$Year)){
  ind = ind + 1
  df = Adult_df %>% 
    filter(Year == i)
try({res = moran.test(df$aspatial_ha_residuals,
                      nb2listw(knn2nb(knearneigh(x = df %>%
                                                   st_as_sf(.,
                                                            coords=c('longitude',
                                                                     'latitude'))))))})
if(ind==1){moran_df = data.frame(Lifestage = "Adult",
                                 Genotype = "ha",
                                 Year = i,
                                 Moran_I = as.numeric(res$estimate[1]),
                                 p_val = res$p.value)}
if(ind>1){moran_df2 = data.frame(Lifestage = "Adult",
                                 Genotype = "ha",
                                 Year = i,
                                 Moran_I = as.numeric(res$estimate[1]),
                                 p_val = res$p.value)
          moran_df = rbind(moran_df2,moran_df)}
}
######
# A-spatial adult v1 model:
#####

adult_v1_aspatial = glmmTMB(formula = v1 ~ metric + 
                              Landscape_metric +
                              (1|UNIT) + (1|UNIT:Site),
                            zi = ~1,
                            family = poisson,
                            offset = log(tot_tested),
                            data = Adult_df)
adult_v1_aspatial2 = glmmTMB(formula = v1 ~ metric + 
                               (1|Site),
                             zi = ~1,
                             offset = log(tot_tested),
                             family = poisson,
                             data = Adult_df)
adult_v1_aspatial3 = glmmTMB(formula = v1 ~ Landscape_metric + 
                               (1|UNIT) + (1|UNIT:Site),
                             zi = ~1,
                             offset = log(tot_tested),
                             family = poisson,
                             data = Adult_df)
if(AIC(adult_v1_aspatial)==min(c(AIC(adult_v1_aspatial),
                                 AIC(adult_v1_aspatial2),
                                 AIC(adult_v1_aspatial3)))){
  adult_v1_best = adult_v1_aspatial}
if(AIC(adult_v1_aspatial2)==min(c(AIC(adult_v1_aspatial),
                                  AIC(adult_v1_aspatial2),
                                  AIC(adult_v1_aspatial3)))){
  adult_v1_best = adult_v1_aspatial2}
if(AIC(adult_v1_aspatial3)==min(c(AIC(adult_v1_aspatial),
                                  AIC(adult_v1_aspatial2),
                                  AIC(adult_v1_aspatial3)))){
  adult_v1_best = adult_v1_aspatial3}
summary(adult_v1_best)
Adult_df$aspatial_v1_residuals = residuals(adult_v1_best)
if(local_moran==T){Moran_maps(df = Adult_df,
           model = adult_v1_best,
           seed = T,
           seed_n = 1)}
for(i in min(Adult_df$Year):max(Adult_df$Year)){
  ind = ind + 1
  df = Adult_df %>% 
    filter(Year == i)
  res = data.frame(estimate = NA,
                   p.value = NA)
try({res = moran.test(df$aspatial_ha_residuals,
                      nb2listw(knn2nb(knearneigh(x = df %>%
                                                   st_as_sf(.,
                                                            coords=c('longitude',
                                                                     'latitude'))))))})
  if(ind>1){moran_df2 = data.frame(Lifestage = "Adult",
                                   Genotype = "v1",
                                   Year = i,
                                   Moran_I = as.numeric(res$estimate[1]),
                                   p_val = res$p.value)
  moran_df = rbind(moran_df2,moran_df)}
}

######
# A-spatial nymph ha model:
#####

nymph_ha_aspatial = glmmTMB(formula = ha ~ metric + 
                              Landscape_metric +
                              (1|Site) + (1|UNIT:Site),
                            zi = ~1,
                            family = poisson,
                            offset = log(tot_tested),
                            data = Nymph_df)
nymph_ha_aspatial2 = glmmTMB(formula = ha ~ metric + 
                               Landscape_metric +
                               (1|Site),
                             zi = ~1,
                             offset = log(tot_tested),
                             family = poisson,
                             data = Nymph_df)
nymph_ha_aspatial3 = glmmTMB(formula = ha ~ Landscape_metric +
                               (1|UNIT) + (1|UNIT:Site),
                             zi = ~1,
                             offset = log(tot_tested),
                             family = poisson,
                             data = Nymph_df)
if(AIC(nymph_ha_aspatial)==min(c(AIC(nymph_ha_aspatial),
                                 AIC(nymph_ha_aspatial2),
                                 AIC(nymph_ha_aspatial3)))){
  nymph_ha_best = adult_v1_aspatial}
if(AIC(nymph_ha_aspatial2)==min(c(AIC(nymph_ha_aspatial),
                                  AIC(nymph_ha_aspatial2),
                                  AIC(nymph_ha_aspatial3)))){
  nymph_ha_best = nymph_ha_aspatial2}
if(AIC(nymph_ha_aspatial3)==min(c(AIC(nymph_ha_aspatial),
                                  AIC(nymph_ha_aspatial2),
                                  AIC(nymph_ha_aspatial3)))){
  nymph_ha_best = nymph_ha_aspatial3}
summary(nymph_ha_best)
Nymph_df$aspatial_ha_residuals = residuals(nymph_ha_best)
if(local_moran==T){Moran_maps(df = Nymph_df,
           model = nymph_ha_best,
           seed = T,
           seed_n = 1)}
for(i in min(Nymph_df$Year):max(Nymph_df$Year)){
  ind = ind + 1
  df = Nymph_df %>% 
    filter(Year == i)
  res = data.frame(estimate = NA,
                   p.value = NA)
  try({res = moran.test(df$aspatial_ha_residuals,
                        nb2listw(knn2nb(knearneigh(x = df %>%
                                                     st_as_sf(.,
                                                              coords=c('longitude',
                                                                       'latitude'))))))})
  if(ind>1){moran_df2 = data.frame(Lifestage = "Nymph",
                                   Genotype = "ha",
                                   Year = i,
                                   Moran_I = as.numeric(res$estimate[1]),
                                   p_val = res$p.value)
  moran_df = rbind(moran_df2,moran_df)}
}
######
# A-spatial nymph v1 model:
#####

nymph_v1_aspatial = glmmTMB(formula = v1 ~ metric + 
                              Landscape_metric +
                              (1|UNIT) + (1|UNIT:Site),
                            zi = ~1,
                            family = poisson,
                            offset = log(tot_tested),
                            data = Nymph_df)
nymph_v1_aspatial2 = glmmTMB(formula = v1 ~ metric + 
                               Landscape_metric +
                               (1|Site),
                             zi = ~1,
                             offset = log(tot_tested),
                             family = poisson,
                             data = Nymph_df)
nymph_v1_aspatial3 = glmmTMB(formula = v1 ~ Landscape_metric +
                               (1|UNIT) + (1|UNIT:Site),
                             zi = ~1,
                             offset = log(tot_tested),
                             family = poisson,
                             data = Nymph_df)
if(AIC(nymph_v1_aspatial)==min(c(AIC(nymph_v1_aspatial),
                                 AIC(nymph_v1_aspatial2),
                                 AIC(nymph_v1_aspatial3)))){
  nymph_v1_best = nymph_v1_aspatial}
if(AIC(nymph_v1_aspatial2)==min(c(AIC(nymph_v1_aspatial),
                                 AIC(nymph_v1_aspatial2),
                                 AIC(nymph_v1_aspatial3)))){
  nymph_v1_best = nymph_v1_aspatial2}
if(AIC(nymph_v1_aspatial3)==min(c(AIC(nymph_v1_aspatial),
                                 AIC(nymph_v1_aspatial2),
                                 AIC(nymph_v1_aspatial3)))){
  nymph_v1_best = nymph_v1_aspatial3}
summary(nymph_v1_best)
Nymph_df$aspatial_v1_residuals = residuals(nymph_v1_best)
if(local_moran==T){Moran_maps(df = Nymph_df,
           model = nymph_v1_best,
           seed = T,
           seed_n = 1)}
for(i in min(Nymph_df$Year):max(Nymph_df$Year)){
  ind = ind + 1
  df = Nymph_df %>% 
    filter(Year == i)
  res = data.frame(estimate = NA,
                   p.value = NA)
  try({res = moran.test(df$aspatial_ha_residuals,
                        nb2listw(knn2nb(knearneigh(x = df %>%
                                                     st_as_sf(.,
                                                              coords=c('longitude',
                                                                       'latitude'))))))})
  if(ind>1){moran_df2 = data.frame(Lifestage = "Nymph",
                                   Genotype = "v1",
                                   Year = i,
                                   Moran_I = as.numeric(res$estimate[1]),
                                   p_val = res$p.value)
  moran_df = rbind(moran_df2,moran_df)}
}
moran_df$p_adj = p.adjust(moran_df$p_val,
                          method = "BH")
moran_df$significant = ifelse(moran_df$p_adj<0.05,1,0)
#####
# A spatial results
#####
summary(adult_ha_best)$coefficients$cond
summary(adult_v1_best)$coefficients$cond
summary(nymph_ha_best)$coefficients$cond
summary(nymph_v1_best)$coefficients$cond


aspatial_results = data.frame(Lifestage = c(rep("Adult",6),
                                            rep("Nymph",6)),
                              Genotype = c("ha","ha",
                                           "v1","v1",
                                           "v1","v1",
                                           "ha","ha",
                                           "v1","v1",
                                           "v1","v1"),
                              Coeff_n = c(2,2,
                                          4,4,4,4,
                                          2,2,
                                          4,4,4,4),
                              Variable = c("Landscape metric","Landscape metric",
                                           "metric","metric",
                                           "Landscape metric","Landscape metric",
                                           "Landscape metric","Landscape metric",
                                           "metric","metric",
                                           "Landscape metric","Landscape metric"),
                              Parameter = rep(c("Estimate","P-value"),6),
                              Value = c(summary(adult_ha_best)$coefficients$cond[2],
                                        summary(adult_ha_best)$coefficients$cond[8],
                                        summary(adult_v1_best)$coefficients$cond[2],
                                        summary(adult_v1_best)$coefficients$cond[11],
                                        summary(adult_v1_best)$coefficients$cond[3],
                                        summary(adult_v1_best)$coefficients$cond[12],
                                        summary(nymph_ha_best)$coefficients$cond[2],
                                        summary(nymph_ha_best)$coefficients$cond[8],
                                        summary(nymph_v1_best)$coefficients$cond[2],
                                        summary(nymph_v1_best)$coefficients$cond[11],
                                        summary(nymph_v1_best)$coefficients$cond[3],
                                        summary(nymph_v1_best)$coefficients$cond[12])) %>%
  pivot_wider(names_from = Parameter,
              values_from = Value) %>%
  mutate(Significant = ifelse(`P-value`<0.05,1,0),
         IRR = exp(Estimate))

All_models = data.frame(Lifestage = c(rep("Adult",6),
                                      rep("Nymph",6)),
                        Genotype = c(rep("ha",3),
                                     rep("v1",3),
                                     rep("ha",3),
                                     rep("v1",3)),
                        Parameters = rep(c("metric + Landscape metric",
                                       "metric",
                                       "Landscape metric"),4),
                        AICs = c(AIC(adult_ha_aspatial),
                                 AIC(adult_ha_aspatial2),
                                 AIC(adult_ha_aspatial3),
                                 AIC(adult_v1_aspatial),
                                 AIC(adult_v1_aspatial2),
                                 AIC(adult_v1_aspatial3),
                                 AIC(nymph_ha_aspatial),
                                 AIC(nymph_ha_aspatial2),
                                 AIC(nymph_ha_aspatial3),
                                 AIC(nymph_v1_aspatial),
                                 AIC(nymph_v1_aspatial2),
                                 AIC(nymph_v1_aspatial3)))
                        
                        
                        
                        