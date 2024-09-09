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


adult_ha_density_mod = glmmTMB(formula = ha ~ density +
                                 (1|UNIT),
                               zi = ~1,
                            family = poisson,
                            offset = log(tot_tested),
                            data = Adult_density_df)
summary(adult_ha_density_mod)

Adult_density_df$residuals = abs(residuals(adult_ha_density_mod))

set.seed(1)
ind = 0
for(i in min(Adult_density_df$Year):max(Adult_density_df$Year)){
  ind = ind + 1
  df = Adult_density_df %>% 
    mutate(Year = as.numeric(substring(Date,1,4)),
           latitude = latitude+runif(n=n(),
                                     min = -0.000001,
                                     max = 0.000001),
           longitude = longitude+runif(n=n(),
                                       min = -0.000001,
                                       max = 0.000001)) %>%
    filter(Year == i)
  try({res = moran.test(df$residuals,
                        nb2listw(knn2nb(knearneigh(x = df %>%
                                                     st_as_sf(.,
                                                              coords=c('longitude',
                                                                       'latitude'))))))})
  if(ind==1){moran_df = data.frame(Lifestage = "Adult",
                                   Year = i,
                                   Genotype = "ha",
                                   Moran_I = as.numeric(res$estimate[1]),
                                   p_val = res$p.value)}
  if(ind>1){moran_df2 = data.frame(Lifestage = "Adult",
                                   Year = i,
                                   Genotype = "ha",
                                   Moran_I = as.numeric(res$estimate[1]),
                                   p_val = res$p.value)
  moran_df = rbind(moran_df2,moran_df)}
}



adult_v1_density_mod = glmmTMB(formula = v1 ~ density +
                              (1|UNIT),
                              zi = ~1,
                            family = poisson,
                            offset = log(tot_tested),
                            data = Adult_density_df)
summary(adult_v1_density_mod)

Adult_density_df$residuals = abs(residuals(adult_v1_density_mod))

set.seed(1)
ind = 0
for(i in min(Adult_density_df$Year):max(Adult_density_df$Year)){
  ind = ind + 1
  df = Adult_density_df %>% 
    mutate(Year = as.numeric(substring(Date,1,4)),
           latitude = latitude+runif(n=n(),
                                     min = -0.000001,
                                     max = 0.000001),
           longitude = longitude+runif(n=n(),
                                       min = -0.000001,
                                       max = 0.000001)) %>%
    filter(Year == i)
  try({res = moran.test(df$residuals,
                        nb2listw(knn2nb(knearneigh(x = df %>%
                                                     st_as_sf(.,
                                                              coords=c('longitude',
                                                                       'latitude'))))))})

  if(ind>1){moran_df2 = data.frame(Lifestage = "Adult",
                                   Year = i,
                                   Genotype = "v1",
                                   Moran_I = as.numeric(res$estimate[1]),
                                   p_val = res$p.value)
  moran_df = rbind(moran_df2,moran_df)}
}
#####
# Build nymph models
#####

nymph_ha_density_mod = glmmTMB(formula = ha ~ density +
                              (1|UNIT),
                              zi = ~1,
                            family = poisson,
                            offset = log(tot_tested),
                            data = Nymph_density_df)
summary(nymph_ha_density_mod)

Nymph_density_df$residuals = abs(residuals(nymph_ha_density_mod))

set.seed(1)
ind = 0
for(i in min(Nymph_density_df$Year):max(Nymph_density_df$Year)){
  ind = ind + 1
  df = Nymph_density_df %>% 
    mutate(Year = as.numeric(substring(Date,1,4)),
           latitude = latitude+runif(n=n(),
                                     min = -0.000001,
                                     max = 0.000001),
           longitude = longitude+runif(n=n(),
                                       min = -0.000001,
                                       max = 0.000001)) %>%
    filter(Year == i)
  try({res = moran.test(df$residuals,
                        nb2listw(knn2nb(knearneigh(x = df %>%
                                                     st_as_sf(.,
                                                              coords=c('longitude',
                                                                       'latitude'))))))})
  
  if(ind>1){moran_df2 = data.frame(Lifestage = "Nymph",
                                   Year = i,
                                   Genotype = "ha",
                                   Moran_I = as.numeric(res$estimate[1]),
                                   p_val = res$p.value)
  moran_df = rbind(moran_df2,moran_df)}
}

nymph_v1_density_mod = glmmTMB(formula = v1 ~ density +
                                 (1|UNIT),
                               zi = ~1,
                            family = poisson,
                            offset = log(tot_tested),
                            data = Nymph_density_df)
summary(nymph_v1_density_mod)

Nymph_density_df$residuals = abs(residuals(nymph_v1_density_mod))

set.seed(1)
ind = 0
for(i in min(Nymph_density_df$Year):max(Nymph_density_df$Year)){
  ind = ind + 1
  df = Nymph_density_df %>% 
    mutate(Year = as.numeric(substring(Date,1,4)),
           latitude = latitude+runif(n=n(),
                                     min = -0.000001,
                                     max = 0.000001),
           longitude = longitude+runif(n=n(),
                                       min = -0.000001,
                                       max = 0.000001)) %>%
    filter(Year == i)
  try({res = moran.test(df$residuals,
                        nb2listw(knn2nb(knearneigh(x = df %>%
                                                     st_as_sf(.,
                                                              coords=c('longitude',
                                                                       'latitude'))))))})
  
  if(ind>1){moran_df2 = data.frame(Lifestage = "Nymph",
                                   Year = i,
                                   Genotype = "v1",
                                   Moran_I = as.numeric(res$estimate[1]),
                                   p_val = res$p.value)
  moran_df = rbind(moran_df2,moran_df)}
}

moran_df$padjust = p.adjust(moran_df$p_val)

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

