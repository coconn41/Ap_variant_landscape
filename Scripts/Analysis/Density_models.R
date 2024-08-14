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

adult_density_mod = lme4::lmer(formula = log(density) ~ metric + 
                              Landscape_metric +
                              (1|UNIT) + (1|UNIT:Site),
                            data = Adult_density_df)
summary(adult_density_mod)
AIC(adult_density_mod)

adult_density_mod = lme4::lmer(formula = log(density) ~ metric + 
                                               (1|UNIT) + (1|UNIT:Site),
                               data = Adult_density_df)
summary(adult_density_mod)
AIC(adult_density_mod)

Adult_density_df$residuals = abs(residuals(adult_density_mod))

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
                                   Moran_I = as.numeric(res$estimate[1]),
                                   p_val = res$p.value)}
  if(ind>1){moran_df2 = data.frame(Lifestage = "Adult",
                                   Year = i,
                                   Moran_I = as.numeric(res$estimate[1]),
                                   p_val = res$p.value)
  moran_df = rbind(moran_df2,moran_df)}
}

#####
# Nymph density model
#####


nymph_density_mod = lme4::lmer(formula = log(density) ~ metric + 
                                 Landscape_metric +
                                 (1|UNIT) + (1|UNIT:Site),
                               data = Nymph_density_df)
summary(nymph_density_mod)
AIC(nymph_density_mod)

nymph_density_mod = lme4::lmer(formula = log(density) ~ metric + 
                                 (1|UNIT) + (1|UNIT:Site),
                               data = Nymph_density_df)
summary(nymph_density_mod)
AIC(nymph_density_mod)

Nymph_density_df$residuals = abs(residuals(nymph_density_mod))

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
  if(ind==1){moran_df3 = data.frame(Lifestage = "Nymph",
                                   Year = i,
                                   Moran_I = as.numeric(res$estimate[1]),
                                   p_val = res$p.value)}
  if(ind>1){moran_df4 = data.frame(Lifestage = "Nymph",
                                   Year = i,
                                   Moran_I = as.numeric(res$estimate[1]),
                                   p_val = res$p.value)
  moran_df3 = rbind(moran_df3,moran_df4)}
}

morandf_fin = rbind(moran_df,moran_df3)
morandf_fin$p_adjust = p.adjust(morandf_fin$p_val)

final_models = data.frame(Lifestage = c("Adult","Nymph"),
                          parameter = c("Patch","Patch"),
                          Beta = c(summary(adult_density_mod)$coefficients[2],
                                   summary(nymph_density_mod)$coefficients[2]),
                          exp_beta = c(exp(summary(adult_density_mod)$coefficients[2]),
                                       exp(summary(nymph_density_mod)$coefficients[2])),
                          T_val = c(summary(adult_density_mod)$coefficients[6],
                                    summary(nymph_density_mod)$coefficients[6]))






