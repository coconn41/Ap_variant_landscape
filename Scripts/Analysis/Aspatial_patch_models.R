#####
# Load libraries
#####
library(readr)
library(tidyverse)
library(mgcv) 
library(spdep)
library(glmmTMB)
library(performance)

######
# A-spatial adult ha model:
#####
adult_ha_aspatial = glmmTMB(formula = ha ~ metric + 
                              (1|Site),
                            zi = ~1,
                            family = poisson,
                            offset = log(tot_tested+0.01),
                            data = Adult_df)
summary(adult_ha_aspatial)

Adult_df$ha_patch_resid = abs(residuals(adult_ha_aspatial))
for(i in 2008:max(Adult_df$Year)){
  df = Adult_df %>%
    filter(Year==i)
  moran_test = moranfast(df$ha_landscape_resid,
                         df$longitude,
                         df$latitude)
  if(i==2008){morandf = data.frame(I = moran_test$observed,
                                   p = moran_test$p.value,
                                   Year = i,
                                   dependent = "ha",
                                   lifestage = "Adult")}
  if(i>2008){morandf2 = data.frame(I = moran_test$observed,
                                   p = moran_test$p.value,
                                   Year = i,
                                   dependent = "ha",
                                   lifestage = "Adult")
  morandf = rbind(morandf,morandf2)}}
morandf$p.adj = p.adjust(morandf$p)
morandf$sig = ifelse(morandf$p.adj<0.05,'sig','ns')
ggplot(data=morandf,
       aes(x=Year,y=I,alpha=sig))+
  geom_point()+
  scale_alpha_manual(values=c(0.3,1))

Moran_maps(df = Adult_df,
           model = adult_ha_aspatial,
           seed = T,
           seed_n = 1)

######
# A-spatial adult v1 model:
#####

adult_v1_aspatial = glmmTMB(formula = v1 ~ metric + 
                              (1|Site),
                            zi = ~1,
                            family = poisson,
                            offset = log(tot_tested+0.01),
                            data = Adult_df)
summary(adult_v1_aspatial)

Adult_df$v1_patch_resid = abs(residuals(adult_v1_aspatial))
for(i in 2008:max(Adult_df$Year)){
  df = Adult_df %>%
    filter(Year==i)
  moran_test = moranfast(df$v1_landscape_resid,
                         df$longitude,
                         df$latitude)
  if(i==2008){morandf = data.frame(I = moran_test$observed,
                                   p = moran_test$p.value,
                                   Year = i,
                                   dependent = "v1",
                                   lifestage = "Adult")}
  if(i>2008){morandf2 = data.frame(I = moran_test$observed,
                                   p = moran_test$p.value,
                                   Year = i,
                                   dependent = "v1",
                                   lifestage = "Adult")
  morandf = rbind(morandf,morandf2)}}
morandf$p.adj = p.adjust(morandf$p)
morandf$sig = ifelse(morandf$p.adj<0.05,'sig','ns')
ggplot(data=morandf,
       aes(x=Year,y=I,alpha=sig))+
  geom_point()+
  scale_alpha_manual(values=c(0.3,1))

Moran_maps(df = Adult_df,
           model = adult_v1_aspatial,
           seed = T,
           seed_n = 1)

######
# A-spatial nymph ha model:
#####

nymph_ha_aspatial = glmmTMB(formula = ha ~ metric + 
                              (1|Site),
                            zi = ~1,
                            family = poisson,
                            offset = log(tot_tested+0.01),
                            data = Nymph_df)

summary(nymph_ha_aspatial)

Nymph_df$ha_landscape_resid = abs(residuals(nymph_ha_aspatial))
starting_year = 2008
for(i in starting_year:max(Nymph_df$Year)){
  df = Nymph_df %>%
    filter(Year==i)
  moran_test = moranfast(df$ha_landscape_resid,
                         df$longitude,
                         df$latitude)
  if(i==starting_year){morandf = data.frame(I = moran_test$observed,
                                            p = moran_test$p.value,
                                            Year = i,
                                            dependent = "ha",
                                            lifestage = "Nymph")}
  if(i>starting_year){morandf2 = data.frame(I = moran_test$observed,
                                            p = moran_test$p.value,
                                            Year = i,
                                            dependent = "ha",
                                            lifestage = "Nymph")
  morandf = rbind(morandf,morandf2)}}
morandf$p.adj = p.adjust(morandf$p)
morandf$sig = ifelse(morandf$p.adj<0.05,'sig','ns')
ggplot(data=morandf,
       aes(x=Year,y=I,alpha=sig))+
  geom_point()+
  scale_alpha_manual(values=c(0.3,1))

Moran_maps(df = Nymph_df,
           model = nymph_ha_aspatial,
           seed = T,
           seed_n = 1)

######
# A-spatial nymph v1 model:
#####

nymph_v1_aspatial = glmmTMB(formula = v1 ~ metric + 
                              (1|Site),
                            zi = ~1,
                            family = poisson,
                            offset = log(tot_tested+0.01),
                            data = Nymph_df)
summary(nymph_v1_aspatial)

Nymph_df$v1_landscape_resid = abs(residuals(nymph_v1_aspatial))
starting_year = 2008
for(i in starting_year:max(Nymph_df$Year)){
  df = Nymph_df %>%
    filter(Year==i)
  moran_test = moranfast(df$v1_landscape_resid,
                         df$longitude,
                         df$latitude)
  if(i==starting_year){morandf = data.frame(I = moran_test$observed,
                                            p = moran_test$p.value,
                                            Year = i,
                                            dependent = "v1",
                                            lifestage = "Nymph")}
  if(i>starting_year){morandf2 = data.frame(I = moran_test$observed,
                                            p = moran_test$p.value,
                                            Year = i,
                                            dependent = "v1",
                                            lifestage = "Nymph")
  morandf = rbind(morandf,morandf2)}}
morandf$p.adj = p.adjust(morandf$p)
morandf$sig = ifelse(morandf$p.adj<0.05,'sig','ns')
ggplot(data=morandf,
       aes(x=Year,y=I,alpha=sig))+
  geom_point()+
  scale_alpha_manual(values=c(0.3,1))

Moran_maps(df = Nymph_df,
           model = nymph_v1_aspatial,
           seed = T,
           seed_n = 1)

#####
# A spatial results
#####
summary(adult_ha_aspatial)$coefficients$cond
summary(adult_v1_aspatial)$coefficients$cond
summary(nymph_ha_aspatial)$coefficients$cond
summary(nymph_v1_aspatial)$coefficients$cond


aspatial_results = data.frame(Lifestage = c(rep("Adult",4),
                                            rep("Nymph",4)),
                              Genotype = c("ha","ha",
                                           "v1","v1",
                                           "ha","ha",
                                           "v1","v1"),
                              Variable = rep("Metric",8),
                              Parameter = rep(c("Estimate","P-value"),4),
                              Value = c(summary(adult_ha_aspatial)$coefficients$cond[2],
                                        summary(adult_ha_aspatial)$coefficients$cond[8],
                                        summary(adult_v1_aspatial)$coefficients$cond[2],
                                        summary(adult_v1_aspatial)$coefficients$cond[8],
                                        summary(nymph_ha_aspatial)$coefficients$cond[2],
                                        summary(nymph_ha_aspatial)$coefficients$cond[8],
                                        summary(nymph_v1_aspatial)$coefficients$cond[2],
                                        summary(nymph_v1_aspatial)$coefficients$cond[8])) %>%
  pivot_wider(names_from = Parameter,
              values_from = Value) %>%
  mutate(Significant = ifelse(`P-value`<0.05,1,0))
