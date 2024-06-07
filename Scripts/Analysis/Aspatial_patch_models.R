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
